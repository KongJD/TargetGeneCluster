import pandas as pd
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))


class FilterBlast(object):
    blast_columns = ["query acc.ver", "subject acc.ver", "identity",
                     "alignment length", "mismatches", "gap opens",
                     "q. start", "q. end", "s. start", "s. end",
                     "evalue", "bit score"]

    def getScoreMax(self, df):
        df1 = df.loc[df.groupby(["query acc.ver", 'subject acc.ver'])['bit score'].idxmax()]
        return df1

    def blaste_filter2(self, df, evalue=1e-5, score=60, identity=40):
        df1 = df.loc[(df["evalue"] <= evalue) \
                     & (df["bit score"] >= score) & (
                             df["identity"] >= identity)]
        filtered_df = df1[df1['query acc.ver'] != df1['subject acc.ver']]
        return filtered_df

    def getStep2(self, path):
        df = pd.read_csv(path, comment="#", sep="\t", header=None)
        df.columns = self.blast_columns
        filter_df = self.blaste_filter2(df)
        return filter_df


class FilterBaseAntiRes:
    columns = ["同源蛋白", "蛋白起始位点", "蛋白终止位点", "基因簇起始位点", "基因簇终止位点",
               "是否在基因簇中", '基因簇名称']

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

    def tmp(self, df):
        def check1(row):
            value = row['是否在基因簇中_c']
            if value.endswith('bp'):
                value = int(value[:-2])
                return value > 2000
            else:
                return value == '基因簇外'

        merge_df = df[df.apply(check1, axis=1)]
        return merge_df

    def tmp2(self, df):

        def check2(row):
            c = row["subject acc.ver_b"]
            b1 = row["subject acc.ver_b1"]
            if c in b1:
                return False
            else:
                return True

        df["query acc.ver_b"] = df["query acc.ver_b"].astype("str")
        df["subject acc.ver_b1"] = df["subject acc.ver_b1"].astype("str")
        mer_df = df[df.apply(check2, axis=1)]
        return mer_df


class FilterHumanGenePro(object):

    def filter1(self, df, my_dict):

        def map_dict_value(row):
            key = row['a']
            key_list = key.split(",")
            m_all = []
            for k in key_list:
                for k2 in my_dict:
                    if k in k2:
                        m_all.append(my_dict[k2])
            return ";".join(list(set(m_all)))

        df['gene_protein'] = df.apply(map_dict_value, axis=1)

        return df


class FilterBGeneBlast(FilterBlast):
    def filter(self, path):
        df = pd.read_csv(path, comment="#", sep="\t", header=None)
        df.columns = self.blast_columns
        df1 = self.getScoreMax(df)
        return df1

    def getBblastResult(self, path):
        data = self.filter(path)
        data_b = self.blaste_filter2(data)
        return data_b




