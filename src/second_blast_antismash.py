import multiprocessing
import re
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from Loggers import Logger
from first_blast_anti import Parttxt, SecondProcessData
from src.judgement import Judement
from src.antismash import Antismash, FindProteinInOrNotAntisamsh
from src.filter import FilterBaseAntiRes
from utils.config import *
from utils.util import *


class SecondPatchData(SecondProcessData):
    def bactie_except_blast_poteinprotein_stepblast(self, filter_df, bl):
        protein = list(set(filter_df["同源蛋白"].tolist()))
        sample_path = os.path.join(shell_path, bl)
        for faa in os.listdir(sample_path):
            if faa.endswith("_protein.faa"):
                faa_path = os.path.join(sample_path, faa)
                faa_dict = util_read_sequcence(faa_path)
                faa_dict = {k: v for k, v in faa_dict.items() if k in protein}
                with open(f"{os.path.join(sample_path, 'query_cluster.faa')}", mode="w") as f:
                    for k, v in faa_dict.items():
                        f.write(f">{k}" + "\n")
                        f.write(f"{v}" + "\n")
                with open(f"{os.path.join(sample_path, 'query_cluster.sh')}", mode="w") as f1:
                    f1.write(f"cd {sample_path}" + "\n")
                    f1.write(
                        f"blastp -db {os.path.join(sample_path, 'blast', 'ref.db')} -query {os.path.join(sample_path, 'query_cluster.faa')} "
                        f"-out {os.path.join(sample_path, 'blast_cluster.txt')} -outfmt 7  -num_threads 11")
                return os.path.join(sample_path, 'blast_cluster.txt'), os.path.join(sample_path,
                                                                                    'query_cluster.sh'), protein

    def sep_comman_second(self, sam, shell, data, logger):
        sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
        samples_parrt = Parttxt().readtxt(sample4w_path)
        if sam not in samples_parrt:
            shell_sam_path = os.path.join(shell, sam)
            data_sam_path = os.path.join(data, sam)
            shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
            shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
            if Judement().exist_gbff_faa(data_sam_path):
                if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                    logger.info(sam + "start")
                    try:
                        anti_res = Antismash().get_result(shell_sam_path_anti)
                    except Exception as e:
                        logger.error(shell_sam_path_anti + "not result is empty")
                        return []
                    else:
                        blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                        blast_df.columns = self.blast_columns
                        blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                        b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                       anti_df.copy())
                        b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                        if b_anti_res_genecluster.empty:
                            logger.info(f"{sam}  没有在基因簇的蛋白")
                            return []
                        b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                            b_anti_res_genecluster.copy(), sam)
                        try:
                            result = subprocess.run(f"sh {query_blast_sh}", shell=True, capture_output=True, text=True)
                            if result.stderr:
                                logger.error(f"{query_blast_sh}\n{result.stderr}")
                        except Exception as e:
                            logger.error(f"{query_blast_sh} blast error")
                            return []
                        else:
                            b_anti_blast_df_include_c = self.bactie_except_protein_incluster_notcluster(b_blast_df_path,
                                                                                                        b_anti_res_genecluster,
                                                                                                        b_anti_res)
                            if b_anti_blast_df_include_c.empty:
                                logger.info(f"{sam} 没有基因簇内的基因 比对到全部时 有符合基因簇外的基因 ")
                                return []
                            all_df = self.add_human_to_rea(b_anti_blast_df_include_c.copy(), blast_filter)
                            if all_df.empty:
                                logger.info(f"{sam} 关联人蛋白时，没有相关的符合的")
                                return []
                            all_df_t = self.slove_allcsv_addhumna(all_df.copy())
                            df1_clusternei_b_b1, t1_co = self.get_cluster_b_b1(b_blast_df_path,
                                                                               b_protein, b_anti_res_genecluster)
                            mer = pd.merge(all_df_t, t1_co, on=["subject acc.ver_a", "基因簇名称_b"], how="left")
                            mer1 = FilterBaseAntiRes().tmp2(mer.copy())
                            if mer1.empty:
                                logger.info(f"{sam} 全部合并时无符合的")
                                return []
                            mer1.to_csv(os.path.join(shell_sam_path, "merge.csv"), index=False)
                            group, id = self.slove_group_tmp(mer1.copy(), sam)
                            group.to_csv(os.path.join(shell_sam_path, "group.csv"), index=False)
                            group.to_csv(os.path.join(shell, "..", "data", "group", f"{id}.csv"), index=False)
                            return group

            elif Judement().exist_gbff_fna(data_sam_path):
                return []

    def get_res(self, shell, data):
        log = Logger(log_file=os.path.join(log_path, "all", "second.log"))
        logger = log.logger
        pool = multiprocessing.Pool(70)
        for sam in os.scandir(shell):
            if sam.is_dir():
                pool.apply_async(self.sep_comman_second, args=(sam.name, shell, data, logger),
                                 error_callback=self.error_callback)
        pool.close()
        pool.join()


if __name__ == '__main__':
    SecondPatchData().get_res(shell_path, data_path)
