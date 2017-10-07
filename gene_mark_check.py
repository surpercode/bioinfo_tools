"""Mark the rows by patient and its id for patient gene mutation localtion"""

import openpyxl as opx
import os
import re
import sys

def get_info_of_varified_gene(summary_xls_file):
    """Extract the patient id and gene from summary_xls_file.
    Argv:
        summary_xls_file -- the summary of patient infos (1 col: id; 3 col: gene; 4 col: name)
    Return:
        info -- a dict {'varified': [(id, name, gene), ...],
                        'uncertain': [(id, name, gene), ...]}
    """
    wb = opx.load_workbook(summary_xls_file)
    sheet_names = wb.get_sheet_names()
    table = wb.get_sheet_by_name(sheet_names[0]) 
    nrows, ncols = table.max_row, table.max_column
    info = {'varified':[], 'uncertain':[]}
    for row in range(3, nrows): # start from 3 row
        row = str(row)
        pid, gene, name = table['A' + row].value, table['C' + row].value, table['D' + row].value
        #print(pid, gene, name)
        if gene and gene.strip() != '无':
            if gene.strip()[-1] == '?':
                info['uncertain'].append((str(pid), gene, name))
            else:
                gene = gene[:-1]
                info['varified'].append((str(pid), gene, name))
    return info

def read_patient_gene_loc(gene_loc_xls_file):
    wb = opx.load_workbook(gene_loc_xls_file)
    sheet_names = wb.get_sheet_names()
    table = wb.get_sheet_by_name(sheet_names[1]) # sheet2
    nrows, ncols = table.max_row, table.max_column
    patient_gene = {}
    for row in range(4, nrows):
        row = str(row)
        pid, gene, geno, alt, ref = table['B' + row].value, table['E' + row].value, table['J' + row].value, table['K' + row].value, table['L' + row].value
        try:
            patient_gene[(str(pid), gene)].append((row, str(geno).strip(), alt, ref))
        except KeyError as e:
            patient_gene[(str(pid), gene)] = []
            patient_gene[(str(pid), gene)].append((row, str(geno).strip(), alt, ref))
    return table, patient_gene
  
def mark_patient(patient_gene, info):
    varified_pids = set([x[0] for x in info['varified']])
    uncertain_pids = set([x[0] for x in info['uncertain']])
    patient_mark = {}
    for k, v in patient_gene.items():
        pid, gene = k
        mark = None
        try:
            pid_elements = re.sub(r'[-_.]', '-', str(pid)).split('-')
            if len(pid_elements) == 3:
                pid_elements.pop(1)
            #print(pid_elements)
            if set(pid_elements) & varified_pids:
                mark = 'varified'
            elif set(pid_elements) & uncertain_pids:
                mark = 'uncertain'
        except TypeError as e:
            pass
        patient_mark[pid] = mark
    return patient_mark

def patient_gene_select(patient_gene):
    patient_gene_match = {}
    for k, v in patient_gene.items():
        geno_match = {'0/1':[], '1/1':[], 'sp':[]}
        for row, geno, alt, ref in v:
            #print([type(x) for x in (row, geno, alt, ref)])
            try:
                int(alt)
            except ValueError as e:
                geno_match['sp'].append(row)
                continue
            if int(alt) >=5:
                if geno == '1/1':
                    geno_match[geno].append(row)
                elif geno == '0/1':
                    if int(ref) / int(alt) <=5:
                        geno_match[geno].append(row)
        print((row, geno, alt, ref),geno_match)
        if geno_match['1/1'] or len(geno_match['0/1']) >= 2 or geno_match['sp']:
            patient_gene_match[k] = 1
        else:
            patient_gene_match[k] = 0
    return patient_gene_match

def read_gene_check(gene_check_file):
    gene_list = {}
    with open(gene_check_file, 'r') as f:
        
def filter_gene(patient_gene):
    pass
        
if __name__ == '__main__':
    info = get_info_of_varified_gene('Summary_of_immature_medical_records.xlsx')
    table, patient_gene = read_patient_gene_loc('待验证复杂合与纯合位点.xlsx')
    #print(patient_gene)
    patient_mark = mark_patient(patient_gene, info)
    #print([(k, v) for k, v in patient_mark.items() if v != None])
    patient_gene_match = patient_gene_select(patient_gene)
    print([(k, v) for k, v in patient_gene_match.items() if v == 1])
