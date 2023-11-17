import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils.util import Singleton


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
