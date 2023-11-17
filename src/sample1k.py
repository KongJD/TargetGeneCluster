import sys
import os

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append("/public/Users/kongjind/pipeline/TargetCluster")

from utils.config import *
from utils.util import *
import datetime
from src.antismash import *
from src.filter import *

target_date = datetime.datetime(2023, 8, 15)


class BlastAntismash1k(object):
    blast_columns = ["query acc.ver", "subject acc.ver", "identity",
                     "alignment length", "mismatches", "gap opens",
                     "q. start", "q. end", "s. start", "s. end",
                     "evalue", "bit score"]

    def hmmer2gbffres(self, data_path, dir, pid):
        for p in os.listdir(os.path.join(data_path, dir)):
            if p.endswith('_genomic.gbff.gz'):
                hmmerout = GbffReults2Blast.blast2locations(data_path, dir, p, pid)
                return hmmerout

    def get1k(self, shell, datapath):
        t = 0
        for sh in os.listdir(shell):
            sh_path = os.path.join(shell, sh)
            sh_path_antismah = os.path.join(sh_path, "antismash")
            blast_txt = os.path.join(sh_path, "blastout.txt")
            if os.path.exists(sh_path_antismah) and os.path.exists(blast_txt):
                blast_time = os.path.getctime(blast_txt)
                creation_date = datetime.datetime.fromtimestamp(blast_time)
                if creation_date > target_date and os.path.exists(sh_path_antismah):
                    if not is_empty_directory(sh_path_antismah):
                        print(sh_path_antismah)
                        try:
                            anti_res = Antismash().get_result(sh_path_antismah)
                        except Exception as e:
                            print(sh_path_antismah + "没有结果")
                            continue
                        blast_res = pd.read_csv(blast_txt, sep="\t", comment="#", header=None)
                        blast_res.columns = self.blast_columns
                        blast_res = blast_res.loc[blast_res["evalue"] < 1e-5]
                        anti_res.to_excel(os.path.join(shell, "..", "tmp", "sample1k", "antismash", f"{sh}_anti.xlsx"),
                                          index=False)
                        blast_res.to_excel(os.path.join(shell, "..", "tmp", "sample1k", "blast", f"{sh}_blast.xlsx"),
                                           index=False)
                        t += 1
                        if t == 1000:
                            break


class FindProteinInOrNotAntisamsh(BlastAntismash1k):
    def findProteinBaseAnti(self, datapath, bl, pro_id, df_anti, minbp=2000):
        blast_res = self.hmmer2gbffres(datapath, bl.replace("_blast.xlsx", ""), pro_id)
        blast_res_group = blast_res.groupby(['blast_id', 'id']). \
            agg({'start': 'min', 'end': 'max'}).reset_index()
        blast_res_group.rename(columns={"id": "gene_id"}, inplace=True)
        merged_df = pd.merge(blast_res_group,
                             df_anti[['Region', 'Type', 'From', 'To', 'Most similar known cluster', 'other',
                                      'Similarity', 'gene_id']], on='gene_id',
                             how='left')
        for index, row in merged_df.iterrows():
            if pd.isna(row['From']):
                merged_df.loc[index, "是否在基因簇中"] = "基因簇外"
                merged_df.loc[index, "基因簇名称"] = ""
                continue
            distance = []
            form = int(row['From'])
            to = int(row['To'])
            start = int(row['start'])
            end = int(row['end'])
            area = row['gene_id'] + " " + row['Region']
            if form < start < to or form < end < to:
                merged_df.loc[index, "是否在基因簇中"] = "基因簇中"
                merged_df.loc[index, "基因簇名称"] = area
            elif start <= form and end >= to:
                merged_df.loc[index, "是否在基因簇中"] = "基因簇中"
                merged_df.loc[index, "基因簇名称"] = area
            else:
                distance.append(min(abs(start - form), abs(start - to),
                                    abs(end - form), abs(end - to)))
                merged_df.loc[index, "基因簇名称"] = area
            if distance:
                if min(distance) <= minbp:
                    merged_df.loc[index, "是否在基因簇中"] = "邻域"
                else:
                    merged_df.loc[index, "是否在基因簇中"] = f"{min(distance)}bp"
        merged_df.rename(columns={"start": "蛋白起始位点", "end": "蛋白终止位点", "blast_id": "同源蛋白",
                                  "From": "基因簇起始位点", "To": "基因簇终止位点"}, inplace=True)

        merged_df = merged_df[['同源蛋白', '蛋白起始位点', '蛋白终止位点', '基因簇起始位点',
                               '基因簇终止位点', '基因簇名称', 'Type', 'Most similar known cluster', 'other',
                               'Similarity', '是否在基因簇中']].copy()
        return merged_df


class BlastAnti1kResult(BlastAntismash1k):
    tmp_blast = os.path.join(shell_path, "..", "tmp", "sample1k", "blast")
    tmp_anti = os.path.join(shell_path, "..", "tmp", "sample1k", "antismash")
    columns = ["同源蛋白", "蛋白起始位点", "蛋白终止位点", "基因簇起始位点", "基因簇终止位点",
               "是否在基因簇中", '基因簇名称']

    def filter_blast_readanti(self, bl, evalue=1e-10, score=60, identity=20):
        bl_sample = os.path.join(self.tmp_blast, bl)
        an_sample = os.path.join(self.tmp_anti, bl.replace("blast", "anti"))
        df_blast = pd.read_excel(bl_sample)
        df_blast_filter = df_blast.loc[(df_blast["evalue"] <= evalue) \
                                       & (df_blast["bit score"] >= score) & (
                                               df_blast["identity"] >= identity)]
        df_blast_scoremax = FilterBlast().getScoreMax(df_blast_filter.copy())
        df_anti = pd.read_excel(an_sample)
        df_anti['From'] = df_anti['From'].astype(str).apply(lambda x: int(x.replace(',', '')))
        df_anti['To'] = df_anti['To'].astype(str).apply(lambda x: int(x.replace(',', '')))
        df_anti['gene_id'] = df_anti['gene_id'].fillna(method='ffill')
        pro_id = list(set(df_blast_scoremax['subject acc.ver'].tolist()))
        return df_blast_scoremax, df_anti, pro_id

    def bactie_except_blast_poteinprotein_stepblast(self, filter_df, bl):
        protein = list(set(filter_df["同源蛋白"].tolist()))
        sample_path = os.path.join(shell_path, bl.replace("_blast.xlsx", ""))
        craete_path = os.path.join(shell_path, '..', "tmp", "sa1k", bl.replace("_blast.xlsx", ""))
        make_dir(craete_path)
        for faa in os.listdir(sample_path):
            if faa.endswith("_protein.faa"):
                faa_path = os.path.join(sample_path, faa)
                faa_dict = util_read_sequcence(faa_path)
                faa_dict = {k: v for k, v in faa_dict.items() if k in protein}
                with open(f"{os.path.join(craete_path, 'query_cluster.faa')}", mode="w") as f:
                    for k, v in faa_dict.items():
                        f.write(f">{k}" + "\n")
                        f.write(f"{v}" + "\n")
                with open(f"{os.path.join(craete_path, 'query_cluster.sh')}", mode="w") as f1:
                    f1.write(f"cd {craete_path}" + "\n")
                    f1.write(
                        f"blastp -db {os.path.join(sample_path, 'blast', 'ref.db')} -query {os.path.join(craete_path, 'query_cluster.faa')} "
                        f"-out {os.path.join(craete_path, 'blast_cluster.txt')} -outfmt 7  -num_threads 11")
                return os.path.join(craete_path, 'blast_cluster.txt'), protein, craete_path

    def bactie_except_blast_poteinprotein_stepgetouterluster(self, path, bl, b_anti_res):
        filter_blast_df = FilterBlast().getStep2(path)
        if filter_blast_df.empty:
            # print(bl + "没有满足要求的结果，基因簇内的基因 比对到 全部时")
            return filter_blast_df
        b_anti_res1 = b_anti_res.copy()
        b_anti_res1["query acc.ver"] = b_anti_res1["同源蛋白"]
        b_anti_blast_df = pd.merge(filter_blast_df, b_anti_res1, on="query acc.ver", how="left")
        b_anti_blast_df.drop("同源蛋白", axis=1, inplace=True)
        b_anti_blast_df = add_column_tail(b_anti_blast_df, "_b")
        b_anti_res2 = b_anti_res.copy()
        b_anti_res2 = b_anti_res2.add_suffix("_c")
        b_anti_res2["subject acc.ver_b"] = b_anti_res2["同源蛋白_c"]
        b_anti_blast_df_include_c = pd.merge(b_anti_blast_df.copy(),
                                             b_anti_res2, on="subject acc.ver_b", how="left")
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

    def main(self, datapath):
        for bl in os.listdir(self.tmp_blast):
            print(bl)
            df_blast_filter, df_anti, pro_id = self.filter_blast_readanti(bl)
            b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(datapath, bl, pro_id, df_anti)
            b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res)
            # b_anti_res_genecluster = add_column_tail(b_anti_res_genecluster, "_b")
            if b_anti_res_genecluster.empty:
                print(f"{bl} 无 在基因簇的蛋白, ")
                continue
            b_blast_df_path, b_protein, create_path = self.bactie_except_blast_poteinprotein_stepblast(
                b_anti_res_genecluster, bl)
            b_anti_blast_df_include_c = self.bactie_except_blast_poteinprotein_stepgetouterluster(
                b_blast_df_path, bl, b_anti_res)
            if b_anti_blast_df_include_c.empty:
                print(f"{bl} 没有满足要求的结果，基因簇内的基因 比对到 全部时")
                continue
            all_df = self.add_human_to_rea(b_anti_blast_df_include_c.copy(), df_blast_filter.copy())
            all_df.to_csv(os.path.join(create_path, "all.csv"))


if __name__ == '__main__':
    BlastAnti1kResult().main(data_path)
