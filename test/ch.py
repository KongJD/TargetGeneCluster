import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append("/public/Users/kongjind/pipeline/TargetCluster")

import multiprocessing
import subprocess

from utils.config import shell_path, data_path, log_path
from src.judgement import Judement


#
# tmp_log_path = "/public/Users/kongjind/pipeline/TargetCluster/test"
#
# comms = []
# m = 0
# for sam in os.listdir(shell_path):
#     sam_path = os.path.join(shell_path, sam)
#     data_sam_path = os.path.join(data_path, sam)
#     if Judement().exist_gbff_faa(data_sam_path):
#         blast_txt = os.path.join(sam_path, "blastout.txt")
#         if os.path.getsize(blast_txt) < 10000000:
#             m += 1
#             print(sam_path)
#             comms.append(f"sh {os.path.join(sam_path, 'blast.sh')}")
#     elif Judement().exist_gbff_fna(data_sam_path):
#         continue
#
# print(m)
#
#
# def excute_command(command):
#     try:
#         result = subprocess.run(command, shell=True, capture_output=True, text=True)
#         if result.stderr:
#             with open(os.path.join(tmp_log_path, "0blast.e.txt"), 'a') as f1:
#                 f1.write(f"{command}\n{result.stderr}\n")
#     except Exception as e:
#         with open(os.path.join(tmp_log_path, "0blast.problem.txt"), 'a') as f2:
#             f2.write(f"{command}\nError executing command\n")

# pool = multiprocessing.Pool(processes=30)
# pool.map(excute_command, comms)
# pool.close()
# pool.join()

class Parttxt:
    def readtxt(self, path):
        with open(path, 'r') as f:
            samples = [i.strip() for i in f.readlines()]
        return samples


sample4w_path = os.path.join(shell_path, "..", "data", "partsample.txt")

samples_parrt = Parttxt().readtxt(sample4w_path)
print(len(samples_parrt))
m = 0
all = 0
t = 0


# 48916
# 13070
def get_strings_with_gbff_gz(string_list):
    gbff_gz_strings = [s for s in string_list if s.endswith("genomic.gbff.gz")]
    return gbff_gz_strings


commandas = []

for sam in os.listdir(shell_path):
    data_sam_path = os.path.join(data_path, sam)
    if Judement().exist_gbff_faa(data_sam_path):
        if sam not in samples_parrt:
            all += 1
            if os.path.exists(os.path.join(shell_path, sam, "antismash", "index.html")):
                m += 1
            else:
                gbff_gz = get_strings_with_gbff_gz(os.listdir(data_sam_path))[0]
                # print(os.path.join(shell_path, sam))
                with open(os.path.join(shell_path, sam, "antismash.sh"), mode="w") as f:
                    if os.path.exists(os.path.join(shell_path, sam, "antismash")):
                        f.write(f"rm -rf {os.path.join(shell_path, sam, 'antismash')}" + "\n")
                    f.write(f"source activate antismash7" + "\n")
                    f.write(f"cd {os.path.join(shell_path, sam)}" + "\n")
                    f.write(
                        f'antismash {os.path.join(data_sam_path, gbff_gz)}  '
                        f'--cc-mibig  '
                        f'--cb-knownclusters  '
                        f'--genefinding-tool none  --output-dir ./antismash' + '\n')
                commandas.append(f"sh {os.path.join(shell_path, sam, 'antismash.sh')}")
    elif Judement().exist_gbff_fna(data_sam_path):
        t += 1

print(m)
print(all)
print(all - m)


def excute_command(command):
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.stderr and result.stderr.startswith('ERROR'):
            with open(os.path.join(log_path, "anti", "anti3.e.txt"), 'a') as f1:
                f1.write(f"{command}\n{result.stderr}\n")
    except Exception as e:
        with open(os.path.join(log_path, "anti", "anti3.problem.txt"), 'a') as f2:
            f2.write(f"{e}\n{command}\n")


# pool = multiprocessing.Pool(processes=5)
# pool.map(excute_command, commandas)
# pool.close()
# pool.join()
