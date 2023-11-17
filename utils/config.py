import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

root_path = "/public/Users/kongjind/pipeline/TargetCluster"

data_path = "/public/Users/kongjind/pipeline/GeneP/data"

targetgene_path = os.path.join(root_path, 'target_gene', 'data_tmp', 'GCF_009914755.1_T2T-CHM13v2.0_protein.faa')

log_path = os.path.join(root_path, "log")

shell_path = os.path.join(root_path, 'shell')

babiao_message = "/public/Users/kongjind/pipeline/TargetCluster/data/整理上市药物靶标信息.xlsx"

junzhu_path = os.path.join(root_path, "data", "summary", "assembly_summary.refseq.txt")
