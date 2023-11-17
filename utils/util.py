import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import os
import subprocess
from threading import Thread
from multiprocessing import Process
from utils.config import junzhu_path, babiao_message, log_path
import pandas as pd


class Singleton(object):
    def __init__(self, cls):
        self._cls = cls
        self._instance = {}

    def __call__(self):
        if self._cls not in self._instance:
            self._instance[self._cls] = self._cls()
        return self._instance[self._cls]


def make_dir(inpath):
    if not os.path.exists(inpath):
        os.makedirs(inpath)


class MyThread(Thread):
    def __init__(self, func, args):
        Thread.__init__(self)
        self.func = func
        self.args = args

    def run(self):
        self.func(*self.args)


class MyProcess(Process):
    def __init__(self, func, args):
        Process.__init__(self)
        self.func = func
        self.args = args

    def run(self):
        self.func(*self.args)


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


@Singleton
class Filutil:
    def get_res(self, pid, re, form):
        if pid not in re:
            re[pid] = []
        re[pid].append(form)
        return re


def util_read_sequcence(path):
    bla_dict = {}
    seq = ""
    f = open(path, mode="r")
    lines = f.readlines()
    for line in lines:
        if line.startswith(">"):
            if seq:
                bla_dict[header] = seq
                seq = ""
            header = line.strip()[1:].split(" ")[0]
        else:
            seq += line.strip()
    bla_dict[header] = seq
    return bla_dict


tmp_log_path = "/public/Users/kongjind/pipeline/TargetCluster/test"


class Command:
    def execute_command_process(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            with open(os.path.join(tmp_log_path, 'blast.o.txt'), 'a') as f:
                f.write(result.stdout)
            with open(os.path.join(tmp_log_path, 'blast.e.txt'), 'a') as f:
                f.write(result.stderr)
        except Exception as e:
            with open(os.path.join(tmp_log_path, 'blast.problem.txt'), 'a') as f:
                f.write(f"Error executing command: {command}\n")

    def excute_query_command(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            # with open(os.path.join(log_path, 'query', "query.o.txt"), 'a') as f:
            #     f.write(result.stdout)
            if result.stderr:
                with open(os.path.join(log_path, 'query', "query.e.txt"), 'a') as f1:
                    f1.write(f"{command}\n{result.stderr}\n")
        except Exception as e:
            with open(os.path.join(log_path, 'query', "query.problem.txt"), 'a') as f2:
                f2.write(f"{command}\nError executing command\n")

    def excute_prokka_command(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            if result.stderr:
                with open(os.path.join(log_path, 'prokka', "prokka.e.txt"), 'a') as f1:
                    f1.write(f"{command}\n{result.stderr}\n")
        except Exception as e:
            with open(os.path.join(log_path, 'prokka', "prokka.problem.txt"), 'a') as f2:
                f2.write(f"{command} problem\n")

    def excute_prokka_blast_command(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            if result.stderr:
                with open(os.path.join(log_path, 'prokka', "blast.e.txt"), 'a') as f1:
                    f1.write(f"{command}\n{result.stderr}\n")
        except Exception as e:
            with open(os.path.join(log_path, 'prokka', "blast.problem.txt"), 'a') as f2:
                f2.write(f"{command} problem\n")

    def excute_prokka_blastnei(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            if result.stderr:
                with open(os.path.join(log_path, 'prokka', "blastnei.e.txt"), 'a') as f1:
                    f1.write(f"{command}\n{result.stderr}\n")
        except Exception as e:
            with open(os.path.join(log_path, 'prokka', "blastnei.problem.txt"), 'a') as f2:
                f2.write(f"{command} problem\n")

    def excute_prokka_anti_command(self, command):
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            if result.stderr:
                with open(os.path.join(log_path, 'prokka', "anti.e.txt"), 'a') as f1:
                    f1.write(f"{command}\n{result.stderr}\n")
        except Exception as e:
            with open(os.path.join(log_path, 'prokka', "anti.problem.txt"), 'a') as f2:
                f2.write(f"{command} problem\n{e}\n")


def add_column_tail(df, tail=""):
    df = df.rename(columns=lambda x: x + tail)
    return df


def read_junzhu_df():
    df = pd.read_csv(junzhu_path,
                     sep="\t", comment="#")[["assembly_accession", "organism_name",
                                             "infraspecific_name", "gbrs_paired_asm"]]
    dict1 = dict(zip(df["assembly_accession"], df["organism_name"]))
    dict2 = dict(zip(df["assembly_accession"], df["infraspecific_name"]))
    dict3 = dict(zip(df["gbrs_paired_asm"], df["organism_name"]))
    dict1.update(dict3)
    return df, dict1, dict2


def read_human_merge_slove():
    df = pd.read_excel(babiao_message)
    df['Gene name'] = df['Gene name'].apply(lambda x: x.replace('\n', '').replace('\t', ' '))
    df['人蛋白靶标'] = df['人蛋白靶标'].apply(
        lambda x: ",".join([i.split(' ')[0].replace('>', '').strip() for i in x.split(';')]))
    df["gene_protein"] = df.apply(lambda row: f"{row['Gene name']} ({row['Protein name']})", axis=1)
    dict1 = dict(zip(df["人蛋白靶标"], df["gene_protein"]))
    return df, dict1


def read_human_extra_message():
    df = pd.read_excel(babiao_message)
    df['Gene name'] = df['Gene name'].apply(lambda x: x.replace('\n', '').replace('\t', ' '))
    df['人蛋白靶标'] = df['人蛋白靶标'].apply(
        lambda x: ",".join([i.split(' ')[0].replace('>', '').strip() for i in x.split(';')]))
    human_gene = dict(zip(df['人蛋白靶标'], df['Gene name']))
    humna_protein = dict(zip(df['人蛋白靶标'], df['Protein name']))
    return df, human_gene, humna_protein
