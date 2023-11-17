"""
part sample (4w) first result
"""
import multiprocessing
import re
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import warnings

warnings.filterwarnings("ignore")
from Loggers import Logger
from antismash import Antismash, FindProteinInOrNotAntisamsh
from src.filter import FilterBlast, FilterBaseAntiRes, FilterBGeneBlast, FilterHumanGenePro

from src.judgement import Judement
from utils.config import *
from utils.util import *


class Parttxt:
    def readtxt(self, path):
        with open(path, 'r') as f:
            samples = [i.strip() for i in f.readlines()]
        return samples


class FirstPatchData:
    blast_columns = ["query acc.ver", "subject acc.ver", "identity",
                     "alignment length", "mismatches", "gap opens",
                     "q. start", "q. end", "s. start", "s. end",
                     "evalue", "bit score"]

    def blast_anti_params(self, df_blast, df_anti, evalue=1e-10, score=60, identity=20):
        df_blast_filter = FilterBlast().getScoreMax(df_blast)
        df_blast_scoremax = df_blast_filter.loc[(df_blast_filter["evalue"] <= evalue) \
                                                & (df_blast_filter["bit score"] >= score) & (
                                                        df_blast_filter["identity"] >= identity)]
        df_anti['From'] = df_anti['From'].astype(str).apply(lambda x: int(x.replace(',', '')))
        df_anti['To'] = df_anti['To'].astype(str).apply(lambda x: int(x.replace(',', '')))
        df_anti['gene_id'] = df_anti['gene_id'].fillna(method='ffill')
        pro_id = list(set(df_blast_scoremax['subject acc.ver'].tolist()))
        return df_blast_scoremax, df_anti, pro_id

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

    def writedata(self, shell, data):
        commands = []
        sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
        samples_parrt = Parttxt().readtxt(sample4w_path)
        for sam in os.listdir(shell):
            if sam in samples_parrt:
                print(sam)
                shell_sam_path = os.path.join(shell, sam)
                data_sam_path = os.path.join(data, sam)
                shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
                shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
                if Judement().exist_gbff_faa(data_sam_path):  # 61986(47257)
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
                            b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                           anti_df.copy())
                            b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                            if b_anti_res_genecluster.empty:
                                print(f"{sam}  没有在基因簇的蛋白")
                                continue
                            b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                                b_anti_res_genecluster.copy(), sam)
                            commands.append(f"sh {query_blast_sh}")
                elif Judement().exist_gbff_fna(data_sam_path):
                    continue
        return commands

    def sep_comman(self, sam, que, shell, data):
        commands = []
        sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
        samples_parrt = Parttxt().readtxt(sample4w_path)
        if sam in samples_parrt:
            print(sam)
            shell_sam_path = os.path.join(shell, sam)
            data_sam_path = os.path.join(data, sam)
            shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
            shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
            if Judement().exist_gbff_faa(data_sam_path):  # 61986(47257)
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
                        b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                       anti_df.copy())
                        b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                        if b_anti_res_genecluster.empty:
                            print(f"{sam}  没有在基因簇的蛋白")
                            return []
                        b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                            b_anti_res_genecluster.copy(), sam)
                        commands.append(f"sh {query_blast_sh}")
            elif Judement().exist_gbff_fna(data_sam_path):
                return []
        que.put(commands)

    def error_callback(self, error):
        print(f"Error info: {error}")

    def get_queue_res(self, shell, data):
        queue = multiprocessing.Manager().Queue()
        pool = multiprocessing.Pool(70)
        for sam in os.listdir(shell):
            pool.apply_async(self.sep_comman, args=(sam, queue, shell, data), error_callback=self.error_callback)
        pool.close()
        pool.join()
        commands = []
        while not queue.empty():
            commands.extend(queue.get())
        return commands

    def main_blast(self, shell, data):
        commands = self.get_queue_res(shell, data)
        pool = multiprocessing.Pool(processes=30)
        pool.map(Command().excute_query_command, commands)
        pool.close()
        pool.join()


class SecondProcessData(FirstPatchData):

    def bactie_except_blast_poteinprotein_stepblast(self, filter_df, bl):
        protein = list(set(filter_df["同源蛋白"].tolist()))
        sample_path = os.path.join(shell_path, bl)
        return os.path.join(sample_path, 'blast_cluster.txt'), os.path.join(sample_path,
                                                                            'query_cluster.sh'), protein

    def bactie_except_protein_incluster_notcluster(self, path, b_anti_res_genecluster, b_anti_res):
        df_blast_scoremax = FilterBGeneBlast().getBblastResult(path)
        b_anti_res1 = b_anti_res_genecluster.copy()
        b_anti_res1["query acc.ver"] = b_anti_res1["同源蛋白"]
        b_anti_blast_df = pd.merge(df_blast_scoremax, b_anti_res1, on="query acc.ver")
        b_anti_blast_df.drop("同源蛋白", axis=1, inplace=True)
        b_anti_blast_df = add_column_tail(b_anti_blast_df, "_b")
        b_anti_res2 = FilterBaseAntiRes().reverser_filterdf(b_anti_res).copy()
        b_anti_res2 = b_anti_res2.add_suffix("_c")
        b_anti_res2["subject acc.ver_b"] = b_anti_res2["同源蛋白_c"]
        b_anti_blast_df_include_c = pd.merge(b_anti_blast_df.copy(),
                                             b_anti_res2, on="subject acc.ver_b")
        b_anti_blast_df_include_c.drop("同源蛋白_c", axis=1, inplace=True)
        return b_anti_blast_df_include_c

    def add_human_to_rea(self, b_anti_blast_df_include_c, df_blast_filter):
        df_blast_filter_copy = df_blast_filter.copy()
        df_blast_filter = df_blast_filter.add_suffix("_a")
        df_blast_filter["query acc.ver_b"] = df_blast_filter["subject acc.ver_a"]
        df1 = pd.merge(df_blast_filter, b_anti_blast_df_include_c, on="query acc.ver_b")
        df_blast_filter_copy = df_blast_filter_copy.add_suffix("_ac")
        df_blast_filter_copy["subject acc.ver_b"] = df_blast_filter_copy["subject acc.ver_ac"]
        df_blast_filter_copy["query acc.ver_a"] = df_blast_filter_copy["query acc.ver_ac"]
        df2 = pd.merge(df1, df_blast_filter_copy, on=("subject acc.ver_b", "query acc.ver_a"))
        return df2

    def slove_allcsv_addhumna(self, allcsv):
        df, human_gene, humna_protein = read_human_extra_message()
        allcsv['Gene name'] = allcsv['query acc.ver_a'].apply(lambda x: human_gene.get(next(
            (key for key in human_gene.keys() if key in x), None)))
        allcsv['Protein name'] = allcsv['query acc.ver_a'].apply(lambda x: humna_protein.get(next(
            (key for key in humna_protein.keys() if key in x), None)))
        return allcsv

    def get_cluster_b_b1(self, b_blast_df_path, b_protein, b_anti_res_genecluster):
        df_blast_scoremax = FilterBGeneBlast().getBblastResult(b_blast_df_path)
        filter_df1_nei = df_blast_scoremax[df_blast_scoremax['subject acc.ver'].isin(b_protein)]
        filter_df1_nei = filter_df1_nei.add_suffix("_b1")
        filter_df1_nei_tmp = filter_df1_nei.copy()
        filter_df1_nei_tmp['subject acc.ver_a'] = filter_df1_nei_tmp['query acc.ver_b1']
        # 簇内的 同源,在同一基因簇内
        b_anti_res_genecluster_co = b_anti_res_genecluster.copy()
        b_anti_res_genecluster_co["subject acc.ver_b1"] = b_anti_res_genecluster_co["同源蛋白"]
        ttt = pd.merge(filter_df1_nei_tmp, b_anti_res_genecluster_co, on="subject acc.ver_b1")
        tttt1 = ttt.groupby(["query acc.ver_b1", "基因簇名称"])["subject acc.ver_b1"].agg(
            lambda x: ','.join(map(str, x))).reset_index()
        tttt1_co = tttt1.copy()
        try:
            tttt1_co["subject acc.ver_a"] = tttt1_co["query acc.ver_b1"]
        except Exception as e:
            tttt1_co["基因簇名称"] = ""
            tttt1_co["query acc.ver_b1"] = ""
            tttt1_co["subject acc.ver_a"] = tttt1_co["query acc.ver_b1"]
        tttt1_co["基因簇名称_b"] = tttt1_co["基因簇名称"]
        tttt1_co.rename(columns={"基因簇名称": "基因簇名称_b1"}, inplace=True)
        tttt1_co = tttt1_co.astype(str)
        return filter_df1_nei_tmp, tttt1_co

    def slove_group_tmp(self, mer, bl):
        mer["query acc.ver_a"] = mer["query acc.ver_a"].astype("str")
        mer["query acc.ver_b"] = mer["query acc.ver_b"].astype("str")
        mer["subject acc.ver_ac"] = mer["subject acc.ver_ac"].astype("str")
        mer["subject acc.ver_b1"] = mer["subject acc.ver_b1"].astype("str")
        group = mer.groupby(["query acc.ver_b", "基因簇名称_b", "基因簇起始位点_b", "基因簇终止位点_b",
                             "Type_b", "Most similar known cluster_b", "other_b", "Similarity_b"]).agg(
            a=("query acc.ver_a", lambda x: ",".join(x.unique())),
            c=("subject acc.ver_b", lambda x: ",".join(x.unique())),
            b1=("subject acc.ver_b1", lambda x: ",".join(x.unique()))
        ).reset_index().rename(columns={"基因簇名称_b": "基因簇名称", "基因簇起始位点_b": "基因簇起始位点",
                                        "基因簇终止位点_b": "基因簇终止位点", "Type_b": "Type",
                                        "Most similar known cluster_b": "Most similar known cluster",
                                        "other_b": "other", "query acc.ver_b": "b",
                                        "Similarity_b": "Similarity"})

        group['b1'] = group["b1"].apply(lambda x: ",".join([i for i in x.split(',') if i != "nan"]))
        bl = bl.replace("_blast.xlsx", "")
        id = re.search(r'GCF_\d+\.\d+', bl).group() if "GCF" in bl else re.search(r'GCA_\d+\.\d+', bl).group()
        junzhu_df, d11, d22 = read_junzhu_df()
        group["序列号"] = id
        group["orgin_name"] = group["序列号"].map(d11)
        group["infra_name"] = group["序列号"].map(d22)

        human_merge, dict1 = read_human_merge_slove()
        group1 = FilterHumanGenePro().filter1(group.copy(), dict1)

        return group1, id

    def slove_sample(self, sam, shell, data, logger):
        sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
        samples_parrt = Parttxt().readtxt(sample4w_path)
        if sam in samples_parrt:
            logger.info(sam)
            shell_sam_path = os.path.join(shell, sam)
            data_sam_path = os.path.join(data, sam)
            shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
            shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
            if Judement().exist_gbff_faa(data_sam_path):  # 61986(47257)
                if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                    try:
                        anti_res = Antismash().get_result(shell_sam_path_anti)
                    except Exception as e:
                        logger.info(sam + "antismash没有结果")
                        return ""
                    else:
                        blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                        blast_df.columns = self.blast_columns
                        blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                        b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                       anti_df.copy())
                        b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                        if b_anti_res_genecluster.empty:
                            logger.info(f"{sam}  没有在基因簇的蛋白")
                            return ""
                        b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                            b_anti_res_genecluster.copy(), sam)
                        b_anti_blast_df_include_c = self.bactie_except_protein_incluster_notcluster(b_blast_df_path,
                                                                                                    b_anti_res_genecluster,
                                                                                                    b_anti_res)
                        if b_anti_blast_df_include_c.empty:
                            logger.info(f"{sam} 没有基因簇内的基因 比对到全部时 有符合基因簇外的基因 ")
                            return ""
                        all_df = self.add_human_to_rea(b_anti_blast_df_include_c.copy(), blast_filter)
                        if all_df.empty:
                            logger.info(f"{sam} 关联人蛋白时，没有相关的符合的")
                            return ""
                        all_df_t = self.slove_allcsv_addhumna(all_df.copy())
                        df1_clusternei_b_b1, t1_co = self.get_cluster_b_b1(b_blast_df_path,
                                                                               b_protein, b_anti_res_genecluster)
                        mer = pd.merge(all_df_t, t1_co, on=["subject acc.ver_a", "基因簇名称_b"], how="left")
                        mer1 = FilterBaseAntiRes().tmp2(mer.copy())
                        if mer1.empty:
                            logger.info(f"{sam} 全部合并时无符合的")
                            return ""
                        mer1.to_csv(os.path.join(shell_sam_path, "merge.csv"), index=False)
                        group, id = self.slove_group_tmp(mer1.copy(), sam)
                        group.to_csv(os.path.join(shell_sam_path, "group.csv"), index=False)
                        group.to_csv(os.path.join(shell, "..", "data", "group", f"{id}.csv"), index=False)
                        return group
            elif Judement().exist_gbff_fna(data_sam_path):
                return ""

    def main(self, shell, data):
        log = Logger(log_file=os.path.join(log_path, "all", "log4w.log"))
        logger = log.logger
        pool = multiprocessing.Pool(80)
        for sam in os.scandir(shell):
            if sam.is_dir():
                pool.apply_async(self.slove_sample, args=(sam.name, shell, data, logger),
                                 error_callback=self.error_callback)
        pool.close()
        pool.join()


if __name__ == '__main__':
    # FirstPatchData().main_blast(shell_path, data_path)
    SecondProcessData().main(shell_path, data_path)
