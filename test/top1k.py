import multiprocessing
import sys
import os

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append("/public/Users/kongjind/pipeline/TargetCluster")

from utils.config import *
from utils.util import *
import datetime
from src.antismash import *

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
                        # pro_id = list(set(blast_res['subject acc.ver'].tolist()))
                        # print(pro_id)
                        # print(len(pro_id))
                        # hmmerout = self.hmmer2gbffres(datapath, sh, pro_id)
                        # print(hmmerout)


class BlastAnti1kResult(BlastAntismash1k):
    tmp_blast = os.path.join(shell_path, "..", "tmp", "sample1k", "blast")
    tmp_anti = os.path.join(shell_path, "..", "tmp", "sample1k", "antismash")
    columns = ["同源蛋白", "蛋白起始位点", "蛋白终止位点", "基因簇起始位点", "基因簇终止位点",
               "是否在基因簇中", '基因簇名称']
    columns1 = ["同源蛋白", "基因簇起始位点", "基因簇终止位点",
                "是否在基因簇中", '基因簇名称']

    def findProtein_BaseAntiRes(self, datapath, bl, pro_id, df_anti, df_blast, minbp=2000):
        blast_res = self.hmmer2gbffres(datapath, bl.replace("_blast.xlsx", ""), pro_id)
        print(blast_res)
        print(df_blast)

        re1 = {}
        re2 = {}
        re3 = {}
        re4 = {}

        df = pd.DataFrame(columns=self.columns)
        blast_res_group = blast_res.groupby(['blast_id', 'id']).agg({'start': 'min', 'end': 'max'}).reset_index()
        print(blast_res_group)

        df_blast["blast_id"] = df_blast['subject acc.ver']
        res1 = pd.merge(blast_res, df_blast, on='blast_id', how='left')
        print(res1)
        df.iloc[:, 0] = pro_id
        df_anti['id'] = df_anti["gene_id"]
        print(df_anti.columns)
        res2 = pd.merge(blast_res_group, df_anti, on="id", how="left")
        print(res2)

        for pid in pro_id:
            distance = []
            names_left = []
            names_right = []
            area_all = []
            ge_id = res1.loc[res1["blast_id"] == pid]['id'].tolist()
            ge_id_full = res1.loc[res1["blast_id"] == pid]
            if len(ge_id_full) > 1:
                ge_id = list(set(ge_id))
                if len(ge_id) > 1:
                    print(f"{pid} have problem")
            anti_res = df_anti.loc[df_anti['gene_id'] == ge_id[0]]
            start = ge_id_full['start'].tolist()[0]
            end = ge_id_full['end'].tolist()[0]
            for index, row in anti_res.iterrows():
                form = int(row['From'])
                to = int(row['To'])
                area = row['gene_id'] + " " + row['Region']
                if form < start < to or form < end < to:
                    re1[pid] = "基因簇里"
                    re2 = Filutil().get_res(pid, re2, form)
                    re3 = Filutil().get_res(pid, re3, to)
                    re4 = Filutil().get_res(pid, re4, area)
                elif start <= form and end >= to:
                    re1[pid] = '基因簇里'
                    re2 = Filutil().get_res(pid, re2, form)
                    re3 = Filutil().get_res(pid, re3, to)
                    re4 = Filutil().get_res(pid, re4, area)
                else:
                    distance.append(min(abs(start - form), abs(start - to),
                                        abs(end - form), abs(end - to)))
                    names_left.append(form)
                    names_right.append(to)
                    area_all.append(area)

            if pid not in re1:
                if distance:
                    if min(distance) <= minbp:
                        re1[pid] = "邻域"
                        re2 = Filutil().get_res(pid, re2, names_left[distance.index(min(distance))])
                        re3 = Filutil().get_res(pid, re3, names_right[distance.index(min(distance))])
                        re4 = Filutil().get_res(pid, re4, area_all[distance.index(min(distance))])
                    else:
                        re1[pid] = f"{min(distance)}bp"
                else:
                    re1[pid] = '基因簇外'

        df.loc[:, "是否在基因簇中"] = df.iloc[:, 0].map(re1)
        df.loc[:, "基因簇起始位点"] = df.iloc[:, 0].map(re2).fillna(""). \
            apply(lambda x: ",".join(map(str, x)) if x else '')

        df.loc[:, "基因簇终止位点"] = df.iloc[:, 0].map(re3).fillna(""). \
            apply(lambda x: ",".join(map(str, x)) if x else '')

        df.loc[:, "基因簇名称"] = df.iloc[:, 0].map(re4).fillna(""). \
            apply(lambda x: ",".join(map(str, x)) if x else "")

        print(df)
        df["blast_id"] = df["同源蛋白"]
        gg = pd.merge(res1, df, on="blast_id")[['同源蛋白', 'start', 'end',
                                                '基因簇起始位点', '基因簇终止位点', '是否在基因簇中', '基因簇名称']]
        gg.drop_duplicates(inplace=True)
        gg.rename(columns={"start": "蛋白起始位点", "end": "蛋白终止位点"}, inplace=True)
        print(gg)
        gg.to_excel("1.xlsx", index=False)

    def filter_blast_readanti(self, bl, evalue=1e-10, score=60, identity=20):
        bl_sample = os.path.join(self.tmp_blast, bl)
        an_sample = os.path.join(self.tmp_anti, bl.replace("blast", "anti"))
        df_blast = pd.read_excel(bl_sample)
        df_blast_filter = df_blast.loc[(df_blast["evalue"] <= evalue) \
                                       & (df_blast["bit score"] >= score) & (
                                               df_blast["identity"] >= identity)]
        df_anti = pd.read_excel(an_sample)
        df_anti['From'] = df_anti['From'].astype(str).apply(lambda x: int(x.replace(',', '')))
        df_anti['To'] = df_anti['To'].astype(str).apply(lambda x: int(x.replace(',', '')))
        df_anti['gene_id'] = df_anti['gene_id'].fillna(method='ffill')
        pro_id = list(set(df_blast_filter['subject acc.ver'].tolist()))
        return df_blast_filter, df_anti, pro_id

    def blaste_filter2(self, df, evalue=1e-5, score=60, identity=40):
        df1 = df.loc[(df["evalue"] <= evalue) \
                     & (df["bit score"] >= score) & (
                             df["identity"] >= identity)]
        filtered_df = df1[df1['query acc.ver'] != df1['subject acc.ver']]
        return filtered_df

    def findProtein_BaseAntiRes1(self, datapath, bl, pro_id, df_anti, minbp=2000):
        blast_res = self.hmmer2gbffres(datapath, bl.replace("_blast.xlsx", ""), pro_id)
        blast_res_group = blast_res.groupby(['blast_id', 'id']). \
            agg({'start': 'min', 'end': 'max'}).reset_index()
        blast_res_group.rename(columns={"id": "gene_id"}, inplace=True)
        merged_df = pd.merge(blast_res_group, df_anti[['Region', 'Type', 'From', 'To', 'gene_id']], on='gene_id',
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
                               '基因簇终止位点', '基因簇名称', '是否在基因簇中']].copy()
        return merged_df

    def filter_mergedf(self, df):
        def check_condition(row):
            value = row['是否在基因簇中']
            if value.endswith('bp'):
                value = int(value[:-2])
                return value <= 2000
            else:
                return value == '基因簇中' or value == "邻域"

        merge_df = df[df.apply(check_condition, axis=1)]
        return merge_df

    def reverser_filterdf(self, df):
        def check(row):
            value = row['是否在基因簇中']
            if value.endswith('bp'):
                value = int(value[:-2])
                return value > 2000
            else:
                return value == '基因簇外'

        merge_df = df[df.apply(check, axis=1)]
        return merge_df

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

    def bactie_except_blast_poteinprotein_stepgetouterluster(self, path, datapath, bl, df_anti, protein):
        df = pd.read_csv(path, comment="#", sep="\t", header=None)
        df.columns = self.blast_columns
        filter_df = self.blaste_filter2(df)
        filter_df1 = filter_df[~filter_df['subject acc.ver'].isin(protein)]
        if filter_df1.empty:
            print(bl + "没有基因簇外的基因")
            return filter_df1, filter_df1
        pro_list = list(set(filter_df1["subject acc.ver"].tolist()))
        tmp_merge_df = self.findProtein_BaseAntiRes1(datapath, bl, pro_list, df_anti)
        tmp_merge_df1 = self.reverser_filterdf(tmp_merge_df)
        return tmp_merge_df1, filter_df1

    def merge_df_(self, filter_df1, tmp_merge_df1):
        filter_df1["同源蛋白"] = filter_df1["subject acc.ver"]
        mer_df = pd.merge(filter_df1, tmp_merge_df1, on="同源蛋白", how="left")
        return mer_df

    def find_proteinB12(self, path, protein, datapath, bl, df_anti):
        df = pd.read_csv(path, comment="#", sep="\t", header=None)
        df.columns = self.blast_columns
        filter_df = self.blaste_filter2(df)
        filter_df1 = filter_df[filter_df['subject acc.ver'].isin(protein)]
        pro_list2 = list(set(filter_df1["subject acc.ver"].tolist()))
        b_merge = self.findProtein_BaseAntiRes1(datapath, bl, pro_list2, df_anti)
        mer = self.merge_df_(filter_df1.copy(), b_merge)
        return mer

    def cu_combline(self, merge_df_filter, cuwai_blastdf):
        pro_list = cuwai_blastdf["query acc.ver"].to_list()
        cu_tongyuandf = merge_df_filter[merge_df_filter["同源蛋白"].isin(pro_list)]
        return cu_tongyuandf

    def human_potein_protein(self, df_blaster_filter, cuwai_blastdf, merger_df_filter, cunei_tongyuan,
                             cuwai_tongyuandf):
        protein_c = cuwai_blastdf["同源蛋白"].tolist()
        protein_b = cuwai_blastdf["query acc.ver"].tolist()
        blastc_df = df_blaster_filter[df_blaster_filter['subject acc.ver'].isin(protein_c)]
        blastb_df = df_blaster_filter[df_blaster_filter['subject acc.ver'].isin(protein_b)]
        b_df = self.merge_df_(blastb_df.copy(), cunei_tongyuan)
        b_df.drop("同源蛋白", axis=1, inplace=True)
        # b_df.to_excel("u.xlsx")
        c_df = self.merge_df_(blastc_df.copy(), cuwai_tongyuandf)
        c_df.drop("同源蛋白", axis=1, inplace=True)
        # merge_bc = pd.merge(blastb_df, blastc_df, on="query acc.ver", suffixes=("_b", "_c"))
        # merge_bc.drop_duplicates(inplace=True)
        merge_bc = pd.merge(b_df, c_df, on="query acc.ver", suffixes=("_b", "_c"))
        return merge_bc

    def main(self, datapath):
        gg = 0
        for bl in os.listdir(self.tmp_blast):
            df_blast_filter, df_anti, pro_id = self.filter_blast_readanti(bl)
            print(bl)
            df_blast_filter.drop_duplicates(inplace=True)
            merge_df = self.findProtein_BaseAntiRes1(datapath, bl, pro_id, df_anti)
            merge_df_filter = self.filter_mergedf(merge_df)
            blast_gene_in_gencu, protein, create_path = self.bactie_except_blast_poteinprotein_stepblast(
                merge_df_filter, bl)
            # genecluster_tongyuna_blast_df = self.find_proteinB12(blast_gene_in_gencu, protein, datapath, bl, df_anti)

            cuwai_tongyuandf, cuwai_blastdf = self.bactie_except_blast_poteinprotein_stepgetouterluster(
                blast_gene_in_gencu, datapath, bl, df_anti, protein)

            if cuwai_tongyuandf.empty:
                gg += 1
                continue

            cunei_tongyuan_df = self.cu_combline(merge_df_filter, cuwai_blastdf)
            cuwai_tongyuan_blast_merge = self.merge_df_(cuwai_blastdf, cuwai_tongyuandf)
            all_bc = self.human_potein_protein(df_blast_filter, cuwai_tongyuan_blast_merge, merge_df_filter,
                                               cunei_tongyuan_df, cuwai_tongyuandf)
            all_bc.to_csv(os.path.join(create_path, "all.csv"), index=False)
        print(gg)

    def run_sh(self, datapath):
        all_sh = []
        for bl in os.listdir(self.tmp_blast):
            df_blast_filter, df_anti, pro_id = self.filter_blast_readanti(bl)
            df_blast_filter.drop_duplicates(inplace=True)
            merge_df = self.findProtein_BaseAntiRes1(datapath, bl, pro_id, df_anti)
            merge_df_filter = self.filter_mergedf(merge_df)
            blast_gene_in_gencu, protein, create_parh = self.bactie_except_blast_poteinprotein_stepblast(
                merge_df_filter, bl)
            all_sh.append(f"sh {blast_gene_in_gencu.replace('blast_cluster.txt', 'query_cluster.sh')}")
        return all_sh

    def run_run_sh(self, datapath):
        commands = self.run_sh(datapath)
        pool = multiprocessing.Pool(processes=4)
        pool.map(Command().execute_command_process, commands)
        pool.close()
        pool.join()


if __name__ == '__main__':
    # BlastAntismash1k().get1k(shell_path, data_path)
    BlastAnti1kResult().main(data_path)
    # BlastAnti1kResult().run_run_sh(data_path)
