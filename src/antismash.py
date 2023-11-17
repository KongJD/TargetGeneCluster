import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from collections import OrderedDict
from bs4 import BeautifulSoup
import pandas as pd
import os
import re

from src.Sampler import Gbff_GZ, ProkkaTmpGZ


class Antismash(object):

    def get_result(self, sample_path):
        content = []
        index_path = os.path.join(sample_path, 'index.html')
        try:
            with open(index_path, 'r') as file:
                html_content = file.read()
        except:
            print(index_path + "不存在")
        soup = BeautifulSoup(html_content, 'html.parser')
        records = list(
            map(lambda x: x.text.strip().split()[0], soup.find_all('div', class_='record-overview-header')))
        column = ['Region', 'Type', 'From', 'To', 'Most similar known cluster', 'other', 'Similarity']
        df = pd.DataFrame(columns=column)
        table = soup.find_all('table', class_='region-table')
        table = list(OrderedDict.fromkeys(table))
        for row in table:
            cells = row.find_all(['td', 'th'])
            for cell in cells:
                text1 = cell.get_text().strip()
                if text1 in column: continue
                if text1.startswith('Region'):
                    if content:
                        if len(content) != len(column):
                            for i in range(len(column) - len(content)):
                                content.append('')
                        df.loc[len(df.index)] = content
                    content = [text1]
                else:
                    if cell.a and 'https' in cell.a['href'] and len(content) > 3:
                        href = cell.a['href'].split('/')[-2]
                        href = '(' + href + ')'
                        text1 += href
                        content.append(text1)
                    else:
                        content.append(text1)
        if content:
            if len(content) != len(column):
                for i in range(len(column) - len(content)):
                    content.append('')
        df.loc[len(df.index)] = content
        df['Region'] = df['Region'].apply(lambda x: x.replace("&nbsp", ''))
        df['number'] = df['Region'].apply(lambda x: x.split('.')[0].replace('Region', ''))
        number = df['number'].tolist()
        num = list(OrderedDict.fromkeys(number))
        num_record = dict(zip(num, records))
        df['gene_id'] = df['number'].map(num_record)
        df = df.drop_duplicates()
        if (df['To'] == '').any():
            for ind, row in df.iterrows():
                if row['To'] == "":
                    row['To'] = row['From']
                    row['From'] = row['Type']
                    row['Type'] = 'other'

        for ind, row in df.iterrows():
            if row['Similarity'] == "" and row['other'] != "":
                row['Similarity'] = row['other']
                row['other'] = row['Most similar known cluster']
                row['Most similar known cluster'] = row['To']
                row['To'] = row['From']
                row['From'] = row['Type']
                row['Type'] = 'other'
        df.drop(['number'], inplace=True, axis=1)
        return df


class GbffReults2Blast:

    def collect_results(self, ndr, dir, p, p_id):
        locate = []
        nc = []
        gbff_dict = Gbff_GZ(os.path.join(ndr, dir, p)).read_gz_file()
        for k, v in gbff_dict.items():
            for v1 in v:
                for k2, v2 in v1.items():
                    if k2 == "features":
                        for v3 in v2:
                            if "cds_location" in v3.keys():
                                for sre in p_id:
                                    if sre == v3["cds_protein_id"]:
                                        cds_location = v3["cds_location"]
                                        nc.append(k)
                                        locate.append({
                                            sre: cds_location
                                        })

        return locate, nc

    @classmethod
    def blast2locations(cls, ndr, dir, p, p_id):
        locate, nc = cls().collect_results(ndr, dir, p, p_id)
        df = pd.DataFrame(columns=["id", "blast_id", "start", "end"])
        output = {}
        for loc in locate:
            for k, v in loc.items():
                if "join" in v:
                    v1 = v.split(",")
                    for v2 in v1:
                        loc_left = int("".join(re.findall("\d+", v2.strip().split(":")[0]))) + 1
                        loc_right = int("".join(re.findall("\d+", v2.strip().split(":")[1])))
                        if nc[locate.index(loc)] not in output:
                            output[nc[locate.index(loc)]] = []
                        output[nc[locate.index(loc)]].append({
                            k: [loc_left, loc_right]
                        }
                        )
                        df.loc[(len(df.index))] = [nc[locate.index(loc)], k, loc_left, loc_right]
                else:
                    loc_left = int("".join(re.findall("\d+", v.strip().split(":")[0]))) + 1
                    loc_right = int("".join(re.findall("\d+", v.strip().split(":")[1])))
                    if nc[locate.index(loc)] not in output:
                        output[nc[locate.index(loc)]] = []
                    output[nc[locate.index(loc)]].append({
                        k: [loc_left, loc_right]
                    }
                    )
                    df.loc[(len(df.index))] = [nc[locate.index(loc)], k, loc_left, loc_right]
        return df


class FindProteinInOrNotAntisamsh:
    def hmmer2gbffres(self, data_path, dir, pid):
        for p in os.listdir(os.path.join(data_path, dir)):
            if p.endswith('_genomic.gbff.gz'):
                hmmerout = GbffReults2Blast.blast2locations(data_path, dir, p, pid)
                return hmmerout

    def findProteinBaseAnti(self, datapath, bl, pro_id, df_anti, minbp=2000):
        blast_res = self.hmmer2gbffres(datapath, bl, pro_id)
        blast_res_group = blast_res.groupby(['blast_id', 'id']). \
            agg({'start': 'min', 'end': 'max'}).reset_index()
        blast_res_group.rename(columns={"id": "gene_id"}, inplace=True)
        merged_df = pd.merge(blast_res_group,
                             df_anti[['Region', 'Type', 'From', 'To', 'Most similar known cluster', 'other',
                                      'Similarity', 'gene_id']], on='gene_id',
                             how='left')
        for index, row in merged_df.iterrows():
            if pd.isna(row['From']):
                merged_df.loc[index, "是否在基因簇中"] = "基因簇外"
                merged_df.loc[index, "基因簇名称"] = ""
                continue
            distance = []
            form = int(row['From'])
            to = int(row['To'])
            start = int(row['start'])
            end = int(row['end'])
            area = row['gene_id'] + " " + row['Region']
            if form < start < to or form < end < to:
                merged_df.loc[index, "是否在基因簇中"] = "基因簇中"
                merged_df.loc[index, "基因簇名称"] = area
            elif start <= form and end >= to:
                merged_df.loc[index, "是否在基因簇中"] = "基因簇中"
                merged_df.loc[index, "基因簇名称"] = area
            else:
                distance.append(min(abs(start - form), abs(start - to),
                                    abs(end - form), abs(end - to)))
                merged_df.loc[index, "基因簇名称"] = area
            if distance:
                if min(distance) <= minbp:
                    merged_df.loc[index, "是否在基因簇中"] = "邻域"
                else:
                    merged_df.loc[index, "是否在基因簇中"] = f"{min(distance)}bp"
        merged_df.rename(columns={"start": "蛋白起始位点", "end": "蛋白终止位点", "blast_id": "同源蛋白",
                                  "From": "基因簇起始位点", "To": "基因簇终止位点"}, inplace=True)

        merged_df = merged_df[['同源蛋白', '蛋白起始位点', '蛋白终止位点', '基因簇起始位点',
                               '基因簇终止位点', '基因簇名称', 'Type', 'Most similar known cluster', 'other',
                               'Similarity', '是否在基因簇中']].copy()
        return merged_df


class ProkkaGbffResult:
    def collect_results(self, path, p_id):
        locate = []
        nc = []
        gbff_dict = ProkkaTmpGZ(path).read_gz_file()
        for k, v in gbff_dict.items():
            for v1 in v:
                for k2, v2 in v1.items():
                    if k2 == "features":
                        for v3 in v2:
                            if "cds_location" in v3.keys():
                                for sre in p_id:
                                    if sre == v3["cds_locus_tag"]:
                                        cds_location = v3["cds_location"]
                                        nc.append(k)
                                        locate.append({
                                            sre: cds_location
                                        })

        return locate, nc

    @classmethod
    def blast2locations(cls, path, p_id):
        locate, nc = cls().collect_results(path, p_id)
        df = pd.DataFrame(columns=["id", "blast_id", "start", "end"])
        output = {}
        for loc in locate:
            for k, v in loc.items():
                if "join" in v:
                    v1 = v.split(",")
                    for v2 in v1:
                        loc_left = int("".join(re.findall("\d+", v2.strip().split(":")[0]))) + 1
                        loc_right = int("".join(re.findall("\d+", v2.strip().split(":")[1])))
                        if nc[locate.index(loc)] not in output:
                            output[nc[locate.index(loc)]] = []
                        output[nc[locate.index(loc)]].append({
                            k: [loc_left, loc_right]
                        }
                        )
                        df.loc[(len(df.index))] = [nc[locate.index(loc)], k, loc_left, loc_right]
                else:
                    loc_left = int("".join(re.findall("\d+", v.strip().split(":")[0]))) + 1
                    loc_right = int("".join(re.findall("\d+", v.strip().split(":")[1])))
                    if nc[locate.index(loc)] not in output:
                        output[nc[locate.index(loc)]] = []
                    output[nc[locate.index(loc)]].append({
                        k: [loc_left, loc_right]
                    }
                    )
                    df.loc[(len(df.index))] = [nc[locate.index(loc)], k, loc_left, loc_right]
        return df


class ProkkaProteinInorNotAntismsh:

    def hmmer2gbffres(self, sample_path, pid):
        for p in os.listdir(os.path.join(sample_path, 'prokka')):
            if p.endswith('gbk'):
                gbff = os.path.join(sample_path, 'prokka', p)
                hmmerout = ProkkaGbffResult.blast2locations(gbff, pid)
                return hmmerout

    def find_protein(self, sample_path, pro_id, df_anti, minbp=2000):
        blast_res = self.hmmer2gbffres(sample_path, pro_id)
        blast_res_group = blast_res.groupby(['blast_id', 'id']). \
            agg({'start': 'min', 'end': 'max'}).reset_index()
        blast_res_group.rename(columns={"id": "gene_id"}, inplace=True)
        merged_df = pd.merge(blast_res_group,
                             df_anti[['Region', 'Type', 'From', 'To', 'Most similar known cluster', 'other',
                                      'Similarity', 'gene_id']], on='gene_id',
                             how='left')
        for index, row in merged_df.iterrows():
            if pd.isna(row['From']):
                merged_df.loc[index, "是否在基因簇中"] = "基因簇外"
                merged_df.loc[index, "基因簇名称"] = ""
                continue
            distance = []
            form = int(row['From'])
            to = int(row['To'])
            start = int(row['start'])
            end = int(row['end'])
            area = row['gene_id'] + " " + row['Region']
            if form < start < to or form < end < to:
                merged_df.loc[index, "是否在基因簇中"] = "基因簇中"
                merged_df.loc[index, "基因簇名称"] = area
            elif start <= form and end >= to:
                merged_df.loc[index, "是否在基因簇中"] = "基因簇中"
                merged_df.loc[index, "基因簇名称"] = area
            else:
                distance.append(min(abs(start - form), abs(start - to),
                                    abs(end - form), abs(end - to)))
                merged_df.loc[index, "基因簇名称"] = area
            if distance:
                if min(distance) <= minbp:
                    merged_df.loc[index, "是否在基因簇中"] = "邻域"
                else:
                    merged_df.loc[index, "是否在基因簇中"] = f"{min(distance)}bp"
        merged_df.rename(columns={"start": "蛋白起始位点", "end": "蛋白终止位点", "blast_id": "同源蛋白",
                                  "From": "基因簇起始位点", "To": "基因簇终止位点"}, inplace=True)

        merged_df = merged_df[['同源蛋白', '蛋白起始位点', '蛋白终止位点', '基因簇起始位点',
                               '基因簇终止位点', '基因簇名称', 'Type', 'Most similar known cluster', 'other',
                               'Similarity', '是否在基因簇中']].copy()
        return merged_df
