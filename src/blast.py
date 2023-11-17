import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from utils.config import *
from utils.util import *
import pandas as pd
import multiprocessing
import subprocess
from multiprocessing import cpu_count
import shutil
from Bio import SeqIO


@Singleton
class Readdir:
    def readconfig(self, dir):
        df = pd.read_csv(dir, sep='\t', header=None)
        return df


@Singleton
class Judement:
    def exist_gbff_faa(self, file_path):
        dirs = [dir for dir in os.listdir(file_path)]
        gbff_gz = any(j.endswith('genomic.gbff.gz') for j in dirs)
        faa_gz = any(j.endswith('protein.faa.gz') for j in dirs)
        if gbff_gz and faa_gz:
            return True
        return False

    def exist_gbff_fna(self, file_path):
        dirs = [dir for dir in os.listdir(file_path)]
        gbff_gz = any(j.endswith('genomic.gbff.gz') for j in dirs)
        fna = any(j.endswith('genomic.fna') for j in dirs)
        if fna and gbff_gz:
            return True
        return False

    def get_latest_folder(self, folder_path):
        folders = [f for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]
        latest_folder = max(folders, key=lambda x: os.path.getmtime(os.path.join(folder_path, x)))
        return os.path.join(folder_path, latest_folder)


class BlastWay:
    def blast_write(self, dir_name, target):
        path1 = os.path.join(data_path, os.path.basename(dir_name))
        make_dir(path1)

    def write_blast_antisamsh(self, dir_name, shellpath, target):
        path1 = os.path.join(shellpath, os.path.basename(dir_name))
        make_dir(path1)
        pp = self.get_strings_with_faa_gz_extension(os.listdir(dir_name))[0]
        gbff_gz = self.get_strings_with_gbff_gz(os.listdir(dir_name))[0]

        with open(os.path.join(path1, "blast.sh"), mode='w') as f:
            f.write(f"cd {path1}" + "\n")
            f.write(f"source activate base" + "\n")
            f.write(
                f"gzip -dc {os.path.join(dir_name, pp)}>{pp.replace('.gz', '')}" + "\n")
            f.write(
                f"makeblastdb -in {pp.replace('.gz', '')} -dbtype prot -out ./blast/ref.db" + "\n")
            f.write(
                f"blastp -db ./blast/ref.db -query {target} -out blastout.txt -outfmt 7  -num_threads 11" + "\n")

        with open(os.path.join(path1, "antismash.sh"), mode='w') as anti:
            anti.write(f"source activate antismash7" + "\n")
            anti.write(f"cd {path1}" + "\n")
            anti.write(f'antismash {os.path.join(dir_name, gbff_gz)}  '
                       f'--cc-mibig --asf '
                       f'--cb-knownclusters  '
                       f'--genefinding-tool none  --output-dir ./antismash' + '\n')
        return path1

    def write_prokka_blast_antismash(self, filepath, target):
        path1 = os.path.join(data_path, os.path.basename(filepath))
        make_dir(path1)
        gbff_gz = self.get_strings_with_gbff_gz(os.listdir(filepath))[0]
        ge_path = self.get_strings_with_gene(os.listdir(filepath))[0]
        with open(os.path.join(path1, "prokka.sh"), mode='w') as pro:
            pro.write(f"cd {path1}" + "\n")
            pro.write(f"source activate prokka" + "\n")
            pro.write(f"prokka {os.path.join(filepath, ge_path)} --outdir ./prokka --kingdom Bacteria" + "\n")

        prokka_gff_path = self.get_strings_with(os.listdir(os.path.join(path1, 'prokka')), ".gff")[0]
        with open(os.path.join(path1, "antismash.sh"), mode='w') as anti:
            anti.write(f"source activate antismash7" + "\n")
            anti.write(f"cd {path1}" + "\n")
            anti.write(f'antismash {os.path.join(filepath, ge_path)}  '
                       f'--cc-mibig --asf '
                       f'--cb-knownclusters  '
                       f'--genefinding-tool none  --output-dir ./antismash --genefinding-gff3 {os.path.join(path1, "prokka", prokka_gff_path)}' + '\n')

        # make_dir(os.path.join(path1, "prokka"))
        prokka_faa_path = self.get_strings_with(os.listdir(os.path.join(path1, "prokka")), ".faa")[0]
        with open(os.path.join(path1, "blast.sh"), mode='w') as bla:
            bla.write(f"source activate base" + "\n")
            bla.write(f"cd {path1}" + "\n")
            bla.write(
                f"makeblastdb -in {os.path.join(path1, 'prokka', prokka_faa_path)} -dbtype prot -out ./blast/ref.db" + "\n")
            bla.write(f"blastp -db ./blast/ref.db -query {target} -out blastout.txt -outfmt 7  -num_threads 11" + "\n")
        return path1

    def write_prokka(self, filepath, shellpath):
        path1 = os.path.join(shellpath, os.path.basename(filepath))
        make_dir(path1)
        ge_path = self.get_strings_with_gene(os.listdir(filepath))[0]
        with open(os.path.join(path1, "prokka.sh"), mode='w') as pro:
            pro.write(f"cd {path1}" + "\n")
            pro.write(f"source activate prokka" + "\n")
            pro.write(f"prokka {os.path.join(filepath, ge_path)} --outdir ./prokka --kingdom Bacteria" + "\n")
        return path1

    def get_strings_with_faa_gz_extension(self, string_list):
        faa_gz_strings = [string for string in string_list if string.endswith("protein.faa.gz")]
        return faa_gz_strings

    def get_strings_with_gbff_gz(self, string_list):
        gbff_gz_strings = [s for s in string_list if s.endswith("genomic.gbff.gz")]
        return gbff_gz_strings

    def get_strings_with_gene(self, string_list):
        gene_strings = [s for s in string_list if s.endswith("genomic.fna")]
        return gene_strings

    def get_strings_with(self, string_list, name):
        name_list = [s for s in string_list if s.endswith(name)]
        return name_list

    def main(self, datapath, shellpath, target):
        for i in os.listdir(datapath):
            sample_path = os.path.join(datapath, i)
            if Judement().exist_gbff_faa(sample_path):
                self.write_blast_antisamsh(sample_path, shellpath, target)
            elif Judement().exist_gbff_fna(sample_path):
                self.write_prokka(sample_path, shellpath)

    def step1_commands_prokka_blast(self, shellpath):
        commands = []
        for sample in os.listdir(shellpath):
            sample_path = os.path.join(shellpath, sample)
            if any(dir.endswith('prokka.sh') for dir in os.listdir(sample_path)):
                commands.append(f'sh {os.path.join(sample_path, "prokka.sh")}')
            else:
                commands.append(f"sh {os.path.join(sample_path, 'blast.sh')}")
        return commands

    def blast_runsh(self, shellpath):
        commands = self.step1_commands_prokka_blast(shellpath)
        pool = multiprocessing.Pool(processes=cpu_count())
        pool.map(Command().execute_command_process, commands)
        pool.close()
        pool.join()


class Command:
    def execute_command(self, command):
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        return output.decode(), error.decode()

    def process_command(self, commands):
        # commands = ["ls", "echo 'Hello, World!'", "pwd"]
        processes = []
        for command in commands:
            process = multiprocessing.Process(target=self.execute_command, args=(command,))
            process.start()
            processes.append(process)
        for process in processes:
            process.join()

    def pro_command(self, command, logpath):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        output, error = process.communicate()

        with open(logpath, 'a') as file:
            file.write(output.decode())
            file.write(error.decode())

    def execute_command_process(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            with open(os.path.join(log_path, 'blast_prokka.o.txt'), 'a') as f:
                f.write(result.stdout)
            with open(os.path.join(log_path, 'blast_prokka.e.txt'), 'a') as f:
                f.write(result.stderr)
        except Exception as e:
            with open(os.path.join(log_path, 'blast_prokka.problem.txt'), 'a') as f:
                f.write(f"Error executing command: {command}\n{e}\n")

    def execute_command_process1(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            with open(os.path.join(log_path, 'blast_prokka1.o.txt'), 'a') as f:
                f.write(result.stdout)
            with open(os.path.join(log_path, 'blast_prokka1.e.txt'), 'a') as f:
                f.write(result.stderr)
        except Exception as e:
            with open(os.path.join(log_path, 'blast_prokka1.problem.txt'), 'a') as f:
                f.write(f"Error executing command: {command}\n{e}\n")

    def execute_command_process2(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            with open(os.path.join(log_path, 'blast_prokka3.o.txt'), 'a') as f:
                f.write(result.stdout)
            with open(os.path.join(log_path, 'blast_prokka3.e.txt'), 'a') as f:
                f.write(result.stderr)
        except Exception as e:
            with open(os.path.join(log_path, 'blast_prokka3.problem.txt'), 'a') as f:
                f.write(f"Error executing command: {command}\n{e}\n")

    def execute_prokka(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            with open(os.path.join(log_path, 'prokka.o.txt'), 'a') as f:
                f.write(result.stdout)
            with open(os.path.join(log_path, 'prokka.e.txt'), 'a') as f:
                f.write(result.stderr)
        except Exception as e:
            with open(os.path.join(log_path, 'prokka.problem.txt'), 'a') as f:
                f.write(f"Error executing command: {command}\n{e}\n")

    def execute_anti(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            with open(os.path.join(log_path, 'anti.o.txt'), 'a') as f:
                f.write(result.stdout)
            with open(os.path.join(log_path, 'anti.e.txt'), 'a') as f:
                f.write(result.stderr)
        except Exception as e:
            with open(os.path.join(log_path, 'anti.problem.txt'), 'a') as f:
                f.write(f"Error executing command: {command}\n{e}\n")


class AntismashWay(BlastWay):
    def write_prokka_blast_antismash(self, filepath, target):
        path1 = os.path.join(data_path, os.path.basename(filepath))
        # make_dir(path1)
        ge_path = self.get_strings_with_gene(os.listdir(filepath))[0]
        prokka_gff_path = self.get_strings_with(os.listdir(os.path.join(path1, 'prokka')), ".gff")[0]
        with open(os.path.join(path1, "antismash.sh"), mode='w') as anti:
            anti.write(f"source activate antismash7" + "\n")
            anti.write(f"cd {path1}" + "\n")
            anti.write(f'antismash {os.path.join(filepath, ge_path)}  '
                       f'--cc-mibig --asf '
                       f'--cb-knownclusters  '
                       f'--genefinding-tool none  --output-dir ./antismash --genefinding-gff3 {os.path.join(path1, "prokka", prokka_gff_path)}' + '\n')

        # make_dir(os.path.join(path1, "prokka"))
        prokka_faa_path = self.get_strings_with(os.listdir(os.path.join(path1, "prokka")), ".faa")[0]
        with open(os.path.join(path1, "blast.sh"), mode='w') as bla:
            bla.write(f"source activate base" + "\n")
            bla.write(f"cd {path1}" + "\n")
            bla.write(
                f"makeblastdb -in {os.path.join(path1, 'prokka', prokka_faa_path)} -dbtype prot -out ./blast/ref.db" + "\n")
            bla.write(f"blastp -db ./blast/ref.db -query {target} -out blastout.txt -outfmt 7  -num_threads 11" + "\n")
        return path1

    def write_prokka(self, filepath, shellpath):
        path1 = os.path.join(shellpath, os.path.basename(filepath))
        ge_path = self.get_strings_with_gene(os.listdir(filepath))[0]
        with open(os.path.join(path1, "prokka.sh"), mode='w') as pro:
            pro.write(f"cd {path1}" + "\n")
            pro.write(f"conda activate prokka" + "\n")
            pro.write(f"prokka {os.path.join(filepath, ge_path)} --outdir ./prokka1 --kingdom Bacteria" + "\n")
        return path1

    def writer_antisamsh(self, shellpath, dir_name):
        path1 = os.path.join(shellpath, os.path.basename(dir_name))
        gbff_gz = self.get_strings_with_gbff_gz(os.listdir(dir_name))[0]
        if os.path.exists(os.path.join(path1, 'antismash')):
            if not self.is_empty_directory(os.path.join(path1, 'antismash')):
                shutil.rmtree(os.path.join(path1, 'antismash'))
        with open(os.path.join(path1, "antismash.sh"), mode='w') as anti:
            # anti.write(f"source /public/Users/kongjind/anoconda/bin/activate" + "\n")
            # anti.write(f"conda activate antismash7" + "\n")
            anti.write(f"cd {path1}" + "\n")
            anti.write(
                f'/public/Users/kongjind/anoconda/envs/antismash7/bin/antismash {os.path.join(dir_name, gbff_gz)}  '
                f'--cc-mibig --asf '
                f'--cb-knownclusters  '
                f'--genefinding-tool none  --output-dir ./antismash' + '\n')

    def is_empty_directory(self, directory_path):
        try:
            contents = os.listdir(directory_path)
            contents = [item for item in contents if item not in ['.', '..']]
            if len(contents) == 0:
                return True
            else:
                return False
        except OSError:
            return False

    def main_anti(self, datapath):
        anti_com = []
        prokka_com = []
        for i in os.listdir(datapath):
            sample_path = os.path.join(datapath, i)
            shell_path_sample = os.path.join(shell_path, os.path.basename(sample_path))
            if Judement().exist_gbff_faa(sample_path):
                if any(dir.endswith('antismash.sh') for dir in os.listdir(shell_path_sample)):
                    # print(shell_path_sample)
                    self.writer_antisamsh(shell_path, sample_path)
                    anti_com.append(f'sh {os.path.join(shell_path_sample, "antismash.sh")}')
            elif Judement().exist_gbff_fna(sample_path):
                self.write_prokka(sample_path, shell_path)
                prokka_com.append(f"sh {os.path.join(shell_path_sample, 'prokka.sh')}")
        # pool = multiprocessing.Pool(processes=100)
        # pool.map(Command().execute_prokka, prokka_com)
        # pool.close()
        # pool.join()

        pool = multiprocessing.Pool(processes=10)
        pool.map(Command().execute_anti, anti_com)
        pool.close()
        pool.join()
        return anti_com, prokka_com

    def finsh_problem_sample(self):
        with open("/public/Users/kongjind/pipeline/TargetCluster/log/blast_prokka2.problem.txt", mode="r") as f:
            lines = [i for i in f.readlines() if i != "[Errno 11] Resource temporarily unavailable\n"]
            commands = [i.replace("Error executing command:", '').strip() for i in lines]
            pool = multiprocessing.Pool(processes=128)
            pool.map(Command().execute_command_process2, commands)
            pool.close()
            pool.join()


class ProkkaSample(AntismashWay):
    def write_gff_faa_dict(self, shellpath):
        prokka_faa_path = self.get_strings_with(os.listdir(os.path.join(shellpath, "prokka1")), ".faa")[0]
        prokka_gff_path = self.get_strings_with(os.listdir(os.path.join(shellpath, "prokka1")), ".gff")[0]
        print(prokka_faa_path)
        print(prokka_gff_path)

        seqid_list = []
        attrs_list = []
        with open(os.path.join(shellpath, 'prokka1', prokka_gff_path), mode='r') as ggf:
            for line in ggf:
                if line.startswith('#'):
                    continue
                if line.startswith('>'):
                    break
                fields = line.strip().split('\t')
                seqid = fields[0]
                attributes = fields[8].split(';')[0].split('=')[1]
                seqid_list.append(seqid)
                attrs_list.append(attributes)
        map_dict = dict(zip(attrs_list, seqid_list))
        # print(map_dict)

    def slove_prokka(self, datapath):
        m = 0
        for i in os.listdir(datapath):
            sample_path = os.path.join(datapath, i)
            shell_path_sample = os.path.join(shell_path, os.path.basename(sample_path))
            if Judement().exist_gbff_faa(sample_path):
                if any(dir.endswith('antismash.sh') for dir in os.listdir(shell_path_sample)):
                    m += 1
                    # continue
            elif Judement().exist_gbff_fna(sample_path):
                # print(shell_path_sample)
                # print(sample_path)
                # self.write_gff_faa_dict(shell_path_sample)
                pass
        print(m)


if __name__ == '__main__':
    # BlastWay().main(data_path, shell_path, targetgene_path)
    # BlastWay().blast_runsh(shell_path)
    # AntismashWay().main_anti(data_path)
    ProkkaSample().slove_prokka(data_path)
