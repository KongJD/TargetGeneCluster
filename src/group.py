import pandas as pd
import os

root = "/public/Users/kongjind/pipeline/TargetCluster/data"
group = os.path.join(root, "group")
group1 = os.path.join(root, "group1")


def main():
    for sam in os.scandir(group):
        sample = sam.name
        sample_path = os.path.join(group, sample)
        df = pd.read_csv(sample_path)
        df_filter = df.dropna(subset=["gene_protein"])
        if df_filter.empty:
            continue
        df_filter.to_csv(os.path.join(group1, sample), index=False)


def slove():
    # all = pd.DataFrame()
    # all.columns = ['b', '基因簇名称', '基因簇起始位点', '基因簇终止位点', 'Type',
    #                'Most similar known cluster', 'other', 'Similarity', 'a', 'c', 'b1',
    #                '序列号', 'orgin_name', 'infra_name', 'gene_protein']
    cluster = 0
    for sam in os.listdir(group1):
        sam_path = os.path.join(group1, sam)
        df = pd.read_csv(sam_path)
        cluster += len(list(set(df["基因簇名称"].tolist())))
        # all = pd.concat([all, df])

    print(cluster)
    # print(all)
    # try:
    #     all.to_excel("/public/Users/kongjind/pipeline/TargetCluster/data/all.xlsx", index=False)
    # except Exception as e:
    #     all.to_csv("/public/Users/kongjind/pipeline/TargetCluster/data/all.csv", index=False)


if __name__ == '__main__':
    # main()
    slove()
    # df = pd.read_excel("/public/Users/kongjind/pipeline/TargetCluster/data/all.xlsx")
    # df = df[['序列号'] + [col for col in df.columns if col != '序列号']]
    # df.to_excel("/public/Users/kongjind/pipeline/TargetCluster/data/all1.xlsx", index=False)
