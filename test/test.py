# from utils.config import *
#
# f = open('all.txt', 'w')
# for i in os.listdir(data_path):
#     f.write(f"{i}\n")
# f.close()
import multiprocessing
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))


# from utils.config import *
# import os
# import subprocess
#
# file_path = "/public/Users/kongjind/pipeline/TargetCluster/log/blast_prokka.e.txt"
# target_prefix = "/public/Users/kongjind/pipeline/TargetCluster/shell"
# count = 0
# commands = []
#
# with open(file_path, "r") as file:
#     for line in file:
#         if line.startswith(target_prefix):
#             commands.append(f"sh {line.strip().split(':')[0]}")
#             count += 1
#
# commands = list(set(commands))
# print(commands[:5])
# print("Number of lines starting with '{}' is: {}".format(target_prefix, count))

def is_empty_directory(directory_path):
    try:
        contents = os.listdir(directory_path)
        contents = [item for item in contents if item not in ['.', '..']]
        if len(contents) == 0:
            return True
        else:
            return False
    except OSError:
        return False


root_path = "/public/Users/kongjind/pipeline/TargetCluster/shell"
m = 0
for sample in os.listdir(root_path):
    sample_path = os.path.join(root_path, sample)
    sample_path_anti = os.path.join(sample_path, "antismash")
    if os.path.exists(sample_path_anti):
        if not is_empty_directory(sample_path_anti):
            # print(sample_path_anti)
            m += 1
print(m)
