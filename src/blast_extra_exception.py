import multiprocessing
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))


from utils.config import *
import os
import subprocess

file_path = "/public/Users/kongjind/pipeline/TargetCluster/log/blast_prokka.e.txt"
target_prefix = "/public/Users/kongjind/pipeline/TargetCluster/shell"


def execute_command_process(command):
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        with open(os.path.join(log_path, 'blast_extra_excption.o.txt'), 'a') as f:
            f.write(result.stdout)
        with open(os.path.join(log_path, 'blast_extra_excption.e.txt'), 'a') as f:
            f.write(result.stderr)
    except Exception as e:
        with open(os.path.join(log_path, 'blast_extra_excption.problem.txt'), 'a') as f:
            f.write(f"Error executing command: {command}\n{e}\n")


count = 0
commands = []

with open(file_path, "r") as file:
    for line in file:
        if line.startswith(target_prefix):
            commands.append(f"sh {line.strip().split(':')[0]}")
            count += 1


# print(commands[:5])
# print("Number of lines starting with '{}' is: {}".format(target_prefix, count))


def main(commands_):
    pool = multiprocessing.Pool(processes=20)
    pool.map(execute_command_process, commands_)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main(commands)
