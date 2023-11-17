import multiprocessing


def process_sample(shell, data, sam):
    commands = []
    sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
    samples_parrt = Parttxt().readtxt(sample4w_path)
    if sam in samples_parrt:
        print(sam)
        shell_sam_path = os.path.join(shell, sam)
        data_sam_path = os.path.join(data, sam)
        shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
        shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
        if Judement().exist_gbff_faa(data_sam_path):  # 61986(47257)
            if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                try:
                    anti_res = Antismash().get_result(shell_sam_path_anti)
                except Exception as e:
                    print(shell_sam_path_anti + "antismash没有结果")
                    return commands
                else:
                    blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                    blast_df.columns = self.blast_columns
                    blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                    b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                   anti_df.copy())
                    b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                    if b_anti_res_genecluster.empty:
                        print(f"{sam}  没有在基因簇的蛋白")
                        return commands
                    b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                        b_anti_res_genecluster.copy(), sam)
                    commands.append(f"sh {query_blast_sh}")
        elif Judement().exist_gbff_fna(data_sam_path):
            return commands

    return commands


def writedata(self, shell, data):
    commands = []
    pool = multiprocessing.Pool()
    for sam in os.listdir(shell):
        result = pool.apply_async(process_sample, (shell, data, sam))
        commands.extend(result.get())
    pool.close()
    pool.join()
    return commands


###
import multiprocessing
from queue import Queue


def process_sample(shell, data, sam, queue):
    sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
    samples_parrt = Parttxt().readtxt(sample4w_path)
    if sam in samples_parrt:
        print(sam)
        shell_sam_path = os.path.join(shell, sam)
        data_sam_path = os.path.join(data, sam)
        shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
        shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
        if Judement().exist_gbff_faa(data_sam_path):  # 61986(47257)
            if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                try:
                    anti_res = Antismash().get_result(shell_sam_path_anti)
                except Exception as e:
                    print(shell_sam_path_anti + "antismash没有结果")
                    return
                else:
                    blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                    blast_df.columns = self.blast_columns
                    blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                    b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                   anti_df.copy())
                    b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                    if b_anti_res_genecluster.empty:
                        print(f"{sam}  没有在基因簇的蛋白")
                        return
                    b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                        b_anti_res_genecluster.copy(), sam)
                    queue.put(f"sh {query_blast_sh}")
        elif Judement().exist_gbff_fna(data_sam_path):
            return


def writedata(self, shell, data):
    queue = multiprocessing.Queue()
    processes = []

    for sam in os.listdir(shell):
        process = multiprocessing.Process(target=process_sample, args=(shell, data, sam, queue))
        process.start()
        processes.append(process)

    for process in processes:
        process.join()

    commands = []
    while not queue.empty():
        commands.append(queue.get())

    return commands

####



def process_sample(shell, data, sam, queue):
    commands = []
    sample4w_path = os.path.join(shell, "..", "data", "partsample.txt")
    samples_parrt = Parttxt().readtxt(sample4w_path)
    if sam in samples_parrt:
        print(sam)
        shell_sam_path = os.path.join(shell, sam)
        data_sam_path = os.path.join(data, sam)
        shell_sam_path_anti = os.path.join(shell_sam_path, "antismash")
        shell_sam_path_blast = os.path.join(shell_sam_path, "blastout.txt")
        if Judement().exist_gbff_faa(data_sam_path):  # 61986(47257)
            if os.path.exists(os.path.join(shell_sam_path_anti, 'index.html')):
                try:
                    anti_res = Antismash().get_result(shell_sam_path_anti)
                except Exception as e:
                    print(shell_sam_path_anti + "antismash没有结果")
                    return []
                else:
                    blast_df = pd.read_csv(shell_sam_path_blast, sep="\t", comment="#", header=None)
                    blast_df.columns = self.blast_columns
                    blast_filter, anti_df, pro = self.blast_anti_params(blast_df.copy(), anti_res.copy())
                    b_anti_res = FindProteinInOrNotAntisamsh().findProteinBaseAnti(data, sam, pro,
                                                                                   anti_df.copy())
                    b_anti_res_genecluster = FilterBaseAntiRes().filter_mergedf(b_anti_res.copy())
                    if b_anti_res_genecluster.empty:
                        print(f"{sam}  没有在基因簇的蛋白")
                        return []
                    b_blast_df_path, query_blast_sh, b_protein = self.bactie_except_blast_poteinprotein_stepblast(
                        b_anti_res_genecluster.copy(), sam)
                    commands.append(f"sh {query_blast_sh}")
        elif Judement().exist_gbff_fna(data_sam_path):
            return []

    queue.put(commands)


def writedata(self, shell, data):
    queue = multiprocessing.Queue()
    pool = multiprocessing.Pool()

    for sam in os.listdir(shell):
        pool.apply_async(process_sample, args=(shell, data, sam, queue))

    pool.close()
    pool.join()

    commands = []
    while not queue.empty():
        commands.extend(queue.get())

    return commands

