import os
file_path = "/public/Users/kongjind/pipeline/TargetCluster/log/blast_extra_excption.o.txt"
target_prefix = "New DB name:   /public/Users/kongjind/pipeline"

commands = []

with open(file_path, "r") as file:
    for line in file:
        if line.startswith(target_prefix):
            commands.append(f"{os.path.dirname(os.path.dirname(line.strip().split(':')[1].strip()))}/blastout.txt")



awk = """grep -v "#" blastout.txt | awk '$11<0.00001 {print $0}' > blastout_filter.txt"""
print(commands[0])

with open("/public/Users/kongjind/test/blast_tmp/awk.sh", "w") as file:
    for i in range(20):
        file.write(f"""grep -v "#" {commands[i]} | awk '$11<0.00001 {{print $0}}' > ./data_query/blastout_filter{i}.txt""" + "\n")





