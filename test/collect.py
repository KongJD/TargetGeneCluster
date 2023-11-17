import pandas as pd
import os

path = "/public/Users/kongjind/pipeline/TargetCluster/tmp/sa1k"

for sam in os.listdir(path):
    sam_path = os.path.join(path, sam)
    csv_path = os.path.join(sam_path, "all.csv")
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        df1 = df.drop(['Unnamed: 0'], axis=1)
        df1.to_csv(os.path.join("/public/Users/kongjind/pipeline/TargetCluster/tmp/1kall",
                                f"{sam}.csv"), index=False)
