import multiprocessing
import os
import sys
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from antismash import Antismash, ProkkaProteinInorNotAntismsh
from src.judgement import Judement
from src.first_blast_anti import SecondProcessData, FirstPatchData
from src.filter import FilterBaseAntiRes
from utils.config import *
from utils.util import util_read_sequcence, Command
from src.Loggers import Logger


class StringSample:
    def get_strings_with(self, string_list, name):
        name_list = [s for s in string_list if s.endswith(name)]
        return name_list

    def write_blast(self, path):
        prokka_faa_path = self.get_strings_with(os.listdir(os.path.join(path, "prokka")), ".faa")[0]
        with open(os.path.join(path, 'blast.sh'), mode="w") as f:
            f.write(f"cd {path}" + "\n")
            f.write(
                f"makeblastdb -in {os.path.join(path, 'prokka', prokka_faa_path)}  -dbtype prot -out ./blast/ref.db" + "\n")
            f.write(
                f"blastp -db ./blast/ref.db -query {targetgene_path} -out blastout.txt -outfmt 7  -num_threads 11" + "\n")
        return os.path.join(path, 'blast.sh')

    def write_antismatch(self, path, sam_data_path):
        gbff = [p for p in os.listdir(sam_data_path) if p.endswith("genomic.gbff.gz")][0]
        with open(os.path.join(path, "antismash.sh"), mode='w') as f:
            f.write(f"cd {path}" + "\n")
            f.write(f"source activate antismash7" + "\n")
            f.write(
                f'/public/Users/kongjind/anoconda/envs/antismash7/bin/antismash {os.path.join(sam_data_path, gbff)}  '
                f'--cc-mibig  '
                f'--cb-knownclusters  '
                f'--genefinding-tool prodigal --output-dir ./antismash' + '\n')
        return os.path.join(path, "antismash.sh")

    def write_prokka(self, path, data_sample):
        gene_path = self.get_strings_with(os.listdir(data_sample), "_genomic.fna")[0]
        with open(os.path.join(path, "prokka.sh"), mode='w') as f:
            if os.path.exists(os.path.join(path, "prokka")):
                f.write(f"rm -rf {os.path.join(path, 'prokka')}" + "\n")
            f.write(f'alias prokka="LC_ALL=en_US.UTF-8 prokka"' + "\n")
            f.write(f"cd {path}" + "\n")
            f.write(f"source activate prokka" + "\n")
            f.write(f"prokka {os.path.join(data_sample, gene_path)} --outdir ./prokka --kingdom Bacteria")
        return os.path.join(path, "prokka.sh")


class ProkkaSample(StringSample):

    def slove_prokka(self, datapath, shell_path):
        prokka_com = []
        blast_com = []
        ant_com = []
        for i in os.listdir(datapath):
            sample_path = os.path.join(datapath, i)
            shell_path_sample = os.path.join(shell_path, os.path.basename(sample_path))
            if Judement().exist_gbff_faa(sample_path):
                continue
            elif Judement().exist_gbff_fna(sample_path):
                blast_path = self.write_blast(shell_path_sample)
                antis_path = self.write_antismatch(shell_path_sample, sample_path)
                prokka_path = self.write_prokka(shell_path_sample, sample_path)
                prokka_com.append(f"sh {prokka_path}")
                ant_com.append(f"sh {antis_path}")
                blast_com.append(f"sh {blast_path}")
        return prokka_com, ant_com, blast_com

    def pool_sample(self, shell, data):
        prokka_com, ant_com, blast_com = self.slove_prokka(data, shell)
        pool = multiprocessing.Pool(processes=5)
        pool.map(Command().excute_prokka_anti_command, ant_com)
        pool.close()
        pool.join()


class SloveProkkaFirst(FirstPatchData):
    def bactie_except_blast_poteinprotein_stepblast(self, filter_df, bl):
        protein = list(set(filter_df["同源蛋白"].tolist()))
        sample_path = os.path.join(shell_path, bl)
        for faa in os.listdir(os.path.join(sample_path, 'prokka')):
            if faa.endswith("faa"):
                faa_path = os.path.join(sample_path, 'prokka', faa)
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

    def slove_main(self, datapath, shell):
        comms = []
        for sam in os.listdir(datapath):
            sample_path = os.path.join(datapath, sam)
            shell_path_sample = os.path.join(shell, os.path.basename(sample_path))
            if Judement().exist_gbff_faa(sample_path):
                continue
            elif Judement().exist_gbff_fna(sample_path):
                print(shell_path_sample)
                shell_sam_path_anti = os.path.join(shell_path_sample, "antismash")
                shell_sam_path_blast = os.path.join(shell_path_sample, "blastout.txt")
                if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                    try:
                        anti_res = Antismash().get_result(shell_sam_path_anti)
                    except Exception as e:
                        print(shell_sam_path_anti + "antismash没有结果")
                        continue
                    else:
                        blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                        blast_df.columns = self.blast_columns
                        blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                        b_anti_res = ProkkaProteinInorNotAntismsh().find_protein(shell_path_sample, pro, anti_df)
                        b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                        if b_anti_res_genecluster.empty:
                            print(f"{sam}  没有在基因簇的蛋白")
                            continue
                        b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                            b_anti_res_genecluster.copy(), sam)
                        comms.append(f"sh {query_blast_sh}")
        return comms

    def slove_pro_sample(self, datapath, shell, sam, que):
        comms = []
        sample_path = os.path.join(datapath, sam)
        shell_path_sample = os.path.join(shell, os.path.basename(sample_path))
        if Judement().exist_gbff_faa(sample_path):
            return []
        elif Judement().exist_gbff_fna(sample_path):
            print(shell_path_sample)
            shell_sam_path_anti = os.path.join(shell_path_sample, "antismash")
            shell_sam_path_blast = os.path.join(shell_path_sample, "blastout.txt")
            if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                try:
                    anti_res = Antismash().get_result(shell_sam_path_anti)
                except Exception as e:
                    print(shell_sam_path_anti + "antismash没有结果")
                    return []
                else:
                    blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                    blast_df.columns = self.blast_columns
                    blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                    b_anti_res = ProkkaProteinInorNotAntismsh().find_protein(shell_path_sample, pro, anti_df)
                    b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                    if b_anti_res_genecluster.empty:
                        print(f"{sam}  没有在基因簇的蛋白")
                        return []
                    b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                        b_anti_res_genecluster.copy(), sam)
                    comms.append(f"sh {query_blast_sh}")
        que.put(comms)

    def get_que_res(self, shell, data):
        queue = multiprocessing.Manager().Queue()
        pool = multiprocessing.Pool(20)
        for sam in os.listdir(shell):
            pool.apply_async(self.slove_pro_sample, args=(data, shell, sam, queue), error_callback=self.error_callback)
        pool.close()
        pool.join()
        commands = []
        while not queue.empty():
            commands.extend(queue.get())
        return commands

    def pool_sample(self, data, shell):
        comms = self.get_que_res(shell, data)
        pool = multiprocessing.Pool(50)
        pool.map(Command().excute_prokka_blastnei, comms)
        pool.close()
        pool.join()


class SloveProkkaBlastAntismah(SecondProcessData):

    def main_pro(self, datapath, shell, sam, logge):
        sample_path = os.path.join(datapath, sam)
        shell_path_sample = os.path.join(shell, os.path.basename(sample_path))
        if Judement().exist_gbff_faa(sample_path):
            return ""
        elif Judement().exist_gbff_fna(sample_path):
            logge.info(shell_path_sample)
            shell_sam_path_anti = os.path.join(shell_path_sample, "antismash")
            shell_sam_path_blast = os.path.join(shell_path_sample, "blastout.txt")
            if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                try:
                    anti_res = Antismash().get_result(shell_sam_path_anti)
                except Exception as e:
                    logge.error(shell_sam_path_anti + "antismash没有结果")
                    return ""
                else:
                    blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                    blast_df.columns = self.blast_columns
                    blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())

                    b_anti_res = ProkkaProteinInorNotAntismsh().find_protein(shell_path_sample, pro, anti_df)

                    b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())

                    if b_anti_res_genecluster.empty:
                        logge.info(f"{sam}  没有在基因簇的蛋白")
                        return ""
                    b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                        b_anti_res_genecluster.copy(), sam)

                    b_anti_blast_df_include_c = self.bactie_except_protein_incluster_notcluster(b_blast_df_path,
                                                                                                b_anti_res_genecluster,
                                                                                                b_anti_res)
                    if b_anti_blast_df_include_c.empty:
                        logge.info(f"{sam} 没有基因簇内的基因 比对到全部时 有符合基因簇外的基因 ")
                        return ""
                    all_df = self.add_human_to_rea(b_anti_blast_df_include_c.copy(), blast_filter)
                    if all_df.empty:
                        logge.info(f"{sam} 关联人蛋白时，没有相关的符合的")
                        return ""
                    all_df_t = self.slove_allcsv_addhumna(all_df.copy())
                    df1_clusternei_b_b1, t1_co = self.get_cluster_b_b1(b_blast_df_path,
                                                                       b_protein, b_anti_res_genecluster)
                    mer = pd.merge(all_df_t, t1_co, on=["subject acc.ver_a", "基因簇名称_b"], how="left")
                    mer1 = FilterBaseAntiRes().tmp2(mer.copy())
                    if mer1.empty:
                        logge.info(f"{sam} 全部合并时无符合的")
                        return ""
                    mer1.to_csv(os.path.join(shell_path_sample, "merge.csv"), index=False)
                    group, id = self.slove_group_tmp(mer1.copy(), sam)
                    group.to_csv(os.path.join(shell_path_sample, "group.csv"), index=False)
                    group.to_csv(os.path.join(shell, "..", "data", "group", f"{id}.csv"), index=False)
                    return group

    def main_al(self, shell, data):
        log = Logger(log_file=os.path.join(log_path, "prokka", "log8.log"))
        logger = log.logger
        pool = multiprocessing.Pool(40)
        for sam in os.scandir(shell):
            if sam.is_dir():
                pool.apply_async(self.main_pro, args=(data, shell, sam.name, logger),
                                 error_callback=self.error_callback)
        pool.close()
        pool.join()



if __name__ == '__main__':
    # ProkkaSample().pool_sample(shell_path, data_path)
    # SloveProkkaFirst().pool_sample(data_path, shell_path)
    SloveProkkaBlastAntismah().main_al(shell_path, data_path)
