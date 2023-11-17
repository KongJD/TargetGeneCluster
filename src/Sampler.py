import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from utils.util import Singleton
from abc import *
import gzip
from Bio import SeqIO
from collections import OrderedDict
import re
import os


class GZ(metaclass=ABCMeta):

    def __init__(self, file_path):
        super(GZ, self).__init__()
        self.file_path = file_path

    @abstractmethod
    def read_gz_file(self):
        pass


class Fna_GZ(GZ):
    def read_gz_file(self):
        with gzip.open(self.file_path, 'rt') as f:
            seq = ''
            seq_dict = {}
            for line in f:
                if line.startswith('>'):
                    if seq:
                        seq_dict[header] = seq
                        seq = ''
                    header = line.strip()[1:]
                else:
                    seq += line.strip()
            seq_dict[header] = seq
        return seq_dict


class Gff_GZ(GZ):
    def read_gz_file(self):
        with gzip.open(self.file_path, 'rt') as f:
            gff_dict = {}
            for line in f:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    seqid = line[0]
                    source = line[1]
                    feature = line[2]
                    start = int(line[3])
                    end = int(line[4])
                    score = line[5]
                    strand = line[6]
                    phase = line[7]
                    attributes = line[8]
                    if seqid not in gff_dict:
                        gff_dict[seqid] = []
                    gff_dict[seqid].append({
                        'source': source,
                        'feature': feature,
                        'start': start,
                        'end': end,
                        'score': score,
                        'strand': strand,
                        'phase': phase,
                        'attributes': attributes
                    })
        return gff_dict


class GB_gz(GZ):
    def read_gz_file(self):
        gb_records = SeqIO.parse(self.file_path, "genbank")

        genbank_dict = {}

        for record in gb_records:

            dm_list = []

            for feature in record.features:
                dm_dict = OrderedDict()
                if feature.type == "source":
                    dm_dict["source_location"] = str(feature.location)
                    dm_dict["source_db_xref"] = feature.qualifiers.get("db_xref", [""])[0]
                    dm_dict["source_mol_type"] = feature.qualifiers.get("mol_type", [""])[0]
                    dm_dict["source_organism"] = feature.qualifiers.get("organism", [""])[0]
                    dm_dict["source_strain"] = feature.qualifiers.get("strain", [""])[0]
                    dm_dict["source_sub_species"] = feature.qualifiers.get("sub_species", [""])[0]
                    dm_list.append(dm_dict)
                    continue

                if feature.type == "CDS":
                    dm_dict["cds_location"] = str(feature.location)
                    dm_dict["cds_codon_start"] = feature.qualifiers.get("codon_start", [""])[0]
                    dm_dict["cds_gene"] = feature.qualifiers.get("gene", [""])[0]
                    dm_dict["cds_note"] = feature.qualifiers.get("note", [""])[0]
                    dm_dict["cds_product"] = feature.qualifiers.get("product", [""])[0]
                    dm_dict["cds_protein_id"] = feature.qualifiers.get("protein_id", [""])[0]
                    dm_dict["cds_transl_table"] = feature.qualifiers.get("transl_table", [""])[0]
                    dm_dict["cds_translation"] = feature.qualifiers.get("translation", [""])[0]
                    dm_list.append(dm_dict)
                    continue

                if feature.type == "gene":
                    dm_dict["gene_location"] = str(feature.location)
                    dm_dict["gene_gene"] = feature.qualifiers.get("gene", [""])[0]
                    dm_list.append(dm_dict)
                    continue

            id = record.id
            if id not in genbank_dict:
                genbank_dict[id] = []
            genbank_dict[id].append({
                "name": record.name,
                "annotations": {k: v for k, v in record.annotations.items() if k != "references"},
                "features": dm_list,
                "seq": str(record.seq)
            })

        return genbank_dict


class Gbff_GZ(GZ):
    def read_gz_file(self):
        genbank_dict = {}
        with gzip.open(self.file_path, "rt") as f:
            for record in SeqIO.parse(f, "genbank"):
                dm_list = []
                for feature in record.features:
                    dm_dict = OrderedDict()

                    if feature.type == "CDS":
                        dm_dict["cds_location"] = str(feature.location)
                        dm_dict["cds_codon_start"] = feature.qualifiers.get("codon_start", [""])[0]
                        dm_dict["cds_gene"] = feature.qualifiers.get("gene", [""])[0]
                        dm_dict["cds_note"] = feature.qualifiers.get("note", [""])[0]
                        dm_dict["cds_product"] = feature.qualifiers.get("product", [""])[0]
                        dm_dict["cds_protein_id"] = feature.qualifiers.get("protein_id", [""])[0]
                        dm_dict["cds_transl_table"] = feature.qualifiers.get("transl_table", [""])[0]
                        dm_dict["cds_translation"] = feature.qualifiers.get("translation", [""])[0]
                        dm_dict["cds_locus_tag"] = feature.qualifiers.get("locus_tag", [""])[0]
                        dm_dict["cds_old_locus_tag"] = feature.qualifiers.get("old_locus_tag", [""])[0]
                        dm_list.append(dm_dict)
                        continue

                id = record.id
                if id not in genbank_dict:
                    genbank_dict[id] = []
                genbank_dict[id].append({
                    "name": record.name,
                    "annotations": {k: v for k, v in record.annotations.items() if k != "references"},
                    "features": dm_list,
                    "seq": str(record.seq)
                })

            return genbank_dict


class Faa_GZ(GZ):
    def read_gz_file(self):
        with gzip.open(self.file_path, 'rt') as f:
            seq = ''
            seq_dict = {}
            for line in f:
                if line.startswith('>'):
                    if seq:
                        seq_dict[header] = seq
                        seq = ''
                    header = line.strip()[1:]
                else:
                    seq += line.strip()
            seq_dict[header] = seq
        return seq_dict


class ProkkaFaa(GZ):
    def read_gz_file(self):
        with open(self.file_path, mode="r") as f:
            lines = f.readlines()
            seq = ''
            seq_dict = {}
            for line in lines:
                if line.startswith('>'):
                    if seq:
                        seq_dict[header] = seq
                        seq = ''
                    header = line.strip()[1:]
                else:
                    seq += line.strip()
            seq_dict[header] = seq
        return seq_dict


class ExtraGb6(GZ):
    def read_gz_file(self):
        genbank_dict = {}
        with open(self.file_path, mode="r") as f:
            for record in SeqIO.parse(f, "genbank"):
                dm_list = []
                for feature in record.features:
                    dm_dict = OrderedDict()
                    if feature.type == "CDS":
                        dm_dict["cds_location"] = str(feature.location)
                        dm_dict["cds_codon_start"] = feature.qualifiers.get("codon_start", [""])[0]
                        dm_dict["cds_gene"] = feature.qualifiers.get("gene", [""])[0]
                        dm_dict["cds_note"] = feature.qualifiers.get("note", [""])[0]
                        dm_dict["cds_product"] = feature.qualifiers.get("product", [""])[0]
                        dm_dict["cds_protein_id"] = feature.qualifiers.get("protein_id", [""])[0]
                        dm_dict["cds_transl_table"] = feature.qualifiers.get("transl_table", [""])[0]
                        dm_dict["cds_translation"] = feature.qualifiers.get("translation", [""])[0]
                        dm_dict["cds_locus_tag"] = feature.qualifiers.get("locus_tag", [""])[0]
                        dm_dict["cds_old_locus_tag"] = feature.qualifiers.get("old_locus_tag", [""])[0]
                        dm_list.append(dm_dict)
                        continue
                id = record.id
                if id not in genbank_dict:
                    genbank_dict[id] = []
                genbank_dict[id].append({
                    "name": record.name,
                    "annotations": {k: v for k, v in record.annotations.items() if k != "references"},
                    "features": dm_list,
                    "seq": str(record.seq)
                })

            return genbank_dict


class ProkkaTmpGZ(GZ):
    def read_gz_file(self):
        genbank_dict = {}
        for record in SeqIO.parse(self.file_path, "genbank"):
            dm_list = []
            for feature in record.features:
                dm_dict = OrderedDict()

                if feature.type == "CDS":
                    dm_dict["cds_location"] = str(feature.location)
                    dm_dict["cds_codon_start"] = feature.qualifiers.get("codon_start", [""])[0]
                    dm_dict["cds_gene"] = feature.qualifiers.get("gene", [""])[0]
                    dm_dict["cds_note"] = feature.qualifiers.get("note", [""])[0]
                    dm_dict["cds_product"] = feature.qualifiers.get("product", [""])[0]
                    dm_dict["cds_protein_id"] = feature.qualifiers.get("protein_id", [""])[0]
                    dm_dict["cds_translation"] = feature.qualifiers.get("translation", [""])[0]
                    dm_dict["cds_locus_tag"] = feature.qualifiers.get("locus_tag", [""])[0]
                    dm_list.append(dm_dict)
                    continue

            id = record.id
            if id not in genbank_dict:
                genbank_dict[id] = []
            genbank_dict[id].append({
                "name": record.name,
                "features": dm_list,
                "seq": str(record.seq)
            })

        return genbank_dict


@Singleton
class GBFF_collect:
    def __init__(self):
        pass

    def tmp(self, genbank_dict, sre, output):
        location, translation, seq_all = self.get_location_translation_list(genbank_dict, sre)

        # print()
        # if sre == "50S ribosomal protein L11":
        #     # print(cds_translation)
        #     print(cds_location)
        # if cds_location != None:
        #     cs_locat = cds_location.split(",")
        #     for loc in cs_locat:
        #         loc_left = int("".join(re.findall("\d+", loc.strip().split(":")[0]))) + 1
        #         loc_right = int("".join(re.findall("\d+", loc.strip().split(":")[1])))
        #         seq1 = seq[int(loc_left) - 1:int(loc_right)]
        #         if sre not in output:
        #             output[sre] = []
        #
        #         output[sre].append({
        #             "start": loc_left,
        #             "end": loc_right,
        #             "seq": seq1,
        #             "translation": cds_translation
        #         })
        pass

    def seqperate_seqloc(self, v3, v1, location, translation, seq_all):
        cds_location = v3["cds_location"]
        cds_translation = v3["cds_translation"]
        location.append(cds_location)
        translation.append(cds_translation)
        seq_all.append(v1["seq"])
        return location, translation, seq_all

    def get_location_translation_list(self, genbank_dict, sre):

        location = []
        translation = []
        seq_all = []

        for k, v in genbank_dict.items():
            for v1 in v:
                for k2, v2 in v1.items():
                    if k2 == "features":
                        for v3 in v2:
                            if "cds_location" in v3.keys():
                                try:
                                    if sre == v3["cds_product"]:
                                        location, translation, seq_all = self.seqperate_seqloc(v3, v1, location,
                                                                                               translation, seq_all)

                                    elif sre == v3["cds_protein_id"]:
                                        location, translation, seq_all = self.seqperate_seqloc(v3, v1, location,
                                                                                               translation, seq_all)

                                    elif sre == v3["cds_locus_tag"]:
                                        location, translation, seq_all = self.seqperate_seqloc(v3, v1, location,
                                                                                               translation, seq_all)

                                    elif sre in v3["cds_old_locus_tag"]:
                                        location, translation, seq_all = self.seqperate_seqloc(v3, v1, location,
                                                                                               translation, seq_all)

                                    elif sre == v3["cds_gene"] or sre == v3["cds_note"]:
                                        location, translation, seq_all = self.seqperate_seqloc(v3, v1, location,
                                                                                               translation, seq_all)

                                    elif sre in v3["cds_product"]:
                                        location, translation, seq_all = self.seqperate_seqloc(v3, v1, location,
                                                                                               translation, seq_all)

                                except Exception as e:
                                    raise Exception(f"{sre} not in cds_gene and cds_product")

        return location, translation, seq_all

    def collect_gbffdict(self, genbank_dict, sre, output):

        location, translation, seq_all = self.get_location_translation_list(genbank_dict, sre)
        for lo in location:
            cs_locat = lo.split(",")
            for csl in cs_locat:
                loc_left = int("".join(re.findall("\d+", csl.strip().split(":")[0]))) + 1
                loc_right = int("".join(re.findall("\d+", csl.strip().split(":")[1])))
                seq1 = seq_all[location.index(lo)][int(loc_left) - 1:int(loc_right)]
                if sre not in output:
                    output[sre] = []

                output[sre].append({
                    "start": loc_left,
                    "end": loc_right,
                    "seq": seq1,
                    "translation": translation[location.index(lo)]
                })
        return output


@Singleton
class ExtraUtils(object):
    def __init__(self):
        pass

    def getGenebankDict(self, rawdir, path):
        for i in os.listdir(os.path.join(rawdir, os.path.basename(path))):
            if i.endswith(".gbff.gz"):
                genebank_dict = Gbff_GZ(os.path.join(rawdir,
                                                     os.path.basename(path), i)).read_gz_file()
                return genebank_dict
