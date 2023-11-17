import pandas as pd
import os


def filter1(df, dict1):
    def _map_dict_values(row):
        key = row["a"]
        key_list = key.split(",")
        res = []
        for k in key_list:
            for k2 in dict1:
                if k in k2:
                    res.append(k)
        return ";".join(list(set(res)))

    df["靶标上市蛋白"] = df.apply(_map_dict_values, axis=1)
    return df


def main():
    df = pd.read_excel("../data/all1.xlsx")
    bai_biao = pd.read_excel("../data/整理上市药物靶标信息.xlsx")
    bai_biao['Gene name'] = bai_biao['Gene name'].apply(lambda x: x.replace('\n', '').replace('\t', ' '))
    bai_biao['人蛋白靶标'] = bai_biao['人蛋白靶标'].apply(
        lambda x: ",".join([i.split(' ')[0].replace('>', '').strip() for i in x.split(';')]))
    bai_biao["gene_protein"] = bai_biao.apply(lambda row: f"{row['Gene name']} ({row['Protein name']})", axis=1)
    dict1 = dict(zip(bai_biao["人蛋白靶标"], bai_biao["gene_protein"]))
    group = filter1(df.copy(), dict1)
    group.to_excel("../data/all_2.xlsx", index=False)


def add():
    df = pd.read_csv("/public/Users/kongjind/pipeline/TargetCluster/data/summary/assembly_summary.txt",
                     sep="\t", header=[1])[["assembly_accession", "organism_name",
                                            "infraspecific_name", "gbrs_paired_asm", ]]

    print(df.columns)
    print(df)
    df_add = pd.read_csv("/public/Users/kongjind/pipeline/TargetCluster/data/summary/assembly_summary.txt.genbank",
                         sep="\t", header=[1])[["assembly_accession", "organism_name",
                                                "infraspecific_name", "gbrs_paired_asm", ]]

    print(df_add)

    dict1 = dict(zip(df["assembly_accession"], df["organism_name"]))
    dict2 = dict(zip(df["assembly_accession"], df["infraspecific_name"]))
    dict3 = dict(zip(df["gbrs_paired_asm"], df["organism_name"]))
    dict4 = dict(zip(df["gbrs_paired_asm"], df["infraspecific_name"]))
    dict1.update(dict3)
    dict2.update(dict4)

    dict1_add = dict(zip(df_add["assembly_accession"], df_add["organism_name"]))
    dict3_add = dict(zip(df_add["gbrs_paired_asm"], df_add["organism_name"]))
    dict1.update(dict3_add)
    dict1.update(dict1_add)

    dict2_add = dict(zip(df_add["assembly_accession"], df_add["infraspecific_name"]))
    dict4_add = dict(zip(df_add["gbrs_paired_asm"], df_add["infraspecific_name"]))
    dict2.update(dict2_add)
    dict2.update(dict4_add)

    group = pd.read_excel("../data/all_2.xlsx")
    group["orgin_name"] = group["序列号"].map(dict1)
    group["infra_name"] = group["序列号"].map(dict2)
    group.to_excel("../data/all_4.xlsx", index=False)


def slove():
    path = "/public/Users/kongjind/pipeline/GeneP/data/"
    df = pd.read_excel("../data/all_5.xlsx")
    df_copy = df[df[["infra_name"]].isnull().T.any()]
    print(df_copy)

    df_copy1 = df[df[["infra_name"]].notnull().T.any()]
    id_list = df_copy["序列号"].to_list()
    # all = pd.concat([df_copy,df_copy1])
    # print(all)
    dict1 = {}
    dict2 = {}
    for sam in os.listdir(path):
        for id in id_list:
            if sam.startswith(id):
                sample_path = os.path.join(path, sam)
                report_txt = [i for i in os.listdir(sample_path) if i.endswith("_report.txt")][0]
                f = open(os.path.join(sample_path, report_txt), mode="r")
                lines = f.readlines()
                for line in lines:
                    if line.startswith("# Organism name:"):
                        wuzhong = line.split(":")[1]
                        dict1[id] = wuzhong.strip()
                    elif line.startswith("# Infraspecific name:"):
                        starin = line.split(":")[1]
                        dict2[id] = starin.strip()
    df_copy["orgin_name"] = df_copy["序列号"].map(dict1)
    df_copy["infra_name"] = df_copy["序列号"].map(dict2)
    all = pd.concat([df_copy1, df_copy])
    all.to_excel("../data/all_6.xlsx", index=False)


if __name__ == '__main__':
    # main()
    # add()
    # slove()
    df = pd.read_excel("../data/all_6.xlsx")
    # df["orgin_name"] = df["orgin_name"].str.replace(r'\(.*\)', '')
    df_copy = df[df[["infra_name"]].isnull().T.any()]
    print(df_copy)
    # df.to_excel("../data/all_7.xlsx")


