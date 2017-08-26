#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Finished by Douyaguang at 2017/8/26
# It is a annovar extraction script by some column filters.
# It can also output the result with only the selected columns.
# Details see the info in usage. You can see it by run this
# script by the -h or -help argv.
# This script is created for my love wenjing.

"""ANNOVAR file extract by filterings.
    Argv:
        annovar_file:
            A file with the snv annotation with many collumns
            which represent different features.
            More details, see:
            https://brb.nci.nih.gov/seqtools/colexpanno.html
            and
            http://annovar.openbioinformatics.org/en/latest/
            especially:
            http://annovar.openbioinformatics.org/en/latest/user-guide/filter/
            and
            http://varianttools.sourceforge.net/Annotation/DbNSFP
        filters:
            Filters to extract the subset of the annovar_file.
        output_file:
            The file to store the annovar_file after filtering.
    Format:
        filters:
            See the details in the usage block.
        output_file:
            Similar to the annovar_file
"""

usage = """
    Usage: [-in input_annovar_file] [-out output_extract_file] [-h help] [-e empty] [filters]
           [-col columns to output] [-s sort the columns by -col para]
        -col arg                   Columns to output.
                                       1  first column
                                       :10  first 10 columns
                                       5:10   5 to 10 columns
                                       exf  Exonic function column
                                       exonFunc  Exonic function column
                                       mutiple choices: link by comma
        -s arg                     Sort the columns by -col argv parameters.
                                       1: sort the columns by -col argv parameters
                                       2: by the original column order

        -e arg                     Whether include the empty data or blank data. Usually 
                                   represented by . or -
                                       1: yes (default)
                                       2: no
    filters_format:
        -chr arg (Dis)             Chromosome.
                                       1~22: corresponding to the same autosome.
                                       X: X chromosome
                                       Y: Y chromosome
                                       S: sex chromosomes
                                       A: autosomes
                                       mutiple choices: link by comma

        -zyg arg (Dis)             Zygosity type.
                                       1 or hom: homozygous
                                       2 or het: heterozygous
                                       multiple choices: link by comma

        -reg arg (Con)             FuncRegion type.
                                       default: 1
                                       1 or exon: exonic
                                       2 or intro: intronic
                                       3 or utr3: UTR3
                                       4 or utr5: UTR5
                                       5 or inter: intergenic
                                       6 or splicing: splicing
                                       7 or nc: ncRNA_exonic, ncRNA_intronic, ncRNA_splicing
                                       8 or up: upstream
                                       9 or down: downstream
                                       multiple choices: link by comma

        -exf arg (Dis)             Exonic function of SNV.
                                       1 or syn: synonymous SNV
                                       2 or nonsyn: nonsynonymous SNV
                                       3 or fsdel: frameshift deletion
                                       4 or fsins: frameshift insertion
                                       5 or nfsdel: nonframeshift deletion
                                       6 or nfsins: nonframeshift insertion
                                       7 or stopg: stopgain
                                       8 or stopl: stoploss
                                       9 or uk: unkown
                                       multiple choices: link by comma

        -snp arg (Dis)             Snp type.
                                       1: in the dbSNP
                                       2: not in the dbSNP

        Mutation rate options:

        -exac_eas arg (Con)        Mutation rate of EAS(east Asian) of The Exome Aggregation Consortium.

        -exac_all arg (Con)        Mutation rate of all populations of The Exome Aggregation Consortium.

        -1kg_eas arg (Con)         Mutation rate of EAS(east Asian) of the 1000 Genomes Project.

        -1kg_all arg (Con)         Mutation rate of all populations of the 1000 Genomes Project.

        -cg (Con)                  Mutation rate of CG (complete genomics) frequency annotations.

        -esp_all arg (Con)         Mutation rate of all populations of NHLBI GO Exome Sequencing Project.

                                   arg:
                                       threshold: the threshold of mutation rate
                                           e.g:
                                               0.01  equal to or less than 0.01
                                               :0.01  equal to or less than 0.01
                                               0.01:  larger than 0.01
                                               0.01:0.02  between 0.01 and 0.02
                                               0.05, 0.09:0.2  less than 0.05 or between 0.09 and 0.2
                                       mutiple choices: link by comma

        Functional Predictions and Annotations for Human Nonsynonymous and Splice-Site SNVs.

        -sift_s arg (Con)          Sift score of SNV. The smaller the more deleterious.

        -sift_p arg (Dis)          Sift prediction of SNV.
                                       1 or D: damaging
                                       2 or T:tolerated

        -pp_hdiv_s arg (Con)       Polyphen score of SNV based on HumDiv. The larger the deleterious.

        -pp_hdiv_p arg (Dis)       Polyphen predictions of SNV based on HumDiv.
                                       1 or D: probably damaging
                                       2 or P: possibly damaging
                                       3 or B: benign

        -pp_hvar_s arg (Con)       Polyphen score of SNV based on HumVar. The larger the deleterious.

        -pp_hvar_p arg (Dis)       Polyphen predictions of SNV based on HumVar.
                                       1 or D: probably damaging
                                       2 or P: possibly damaging
                                       3 or B: benign

        -lrt_s arg (Con)           The original LRT two-sided p-value. The smaller the more coffidence.

        -lrt_p arg (Dis)           LRT predictions.
                                       1 or D: deleterious
                                       2 or N: neutral
                                       3 or U: unkown

        -tas_s arg (Con)           MutationTaster score represent the probability.

        -tas_p arg (Dis)           MutationTaster predictions.
                                       1 or A: disease_causing_automatic - i.e. known to be deleterious
                                       2 or D: disease causing - i.e. probably deleterious
                                       3 or P: polymorphism automatic - i.e. known to be harmless
                                       4 or N: polymorphism - i.e. probably harmless

        -ass_s arg (Con)           MutationAssessor score represent the probability. The larger the more
                                   coffidence it has potential functional impact.

        -ass_p arg (Dis)           MutationAssessor predictions.
                                       1 or H: high
                                       2 or M: medium
                                       3 or L: low
                                       4 or N: neutral

        -fm_s arg (Con)            FATHMM (Functional Analysis through Hidden Markov Models) score.
                                   The smaller the more deleterious.

        -fm_p arg (Dis)            FATHMM predictions.
                                       1 or D: damaging
                                       2 or T: tolerated

        -pv_s arg (Con)            PROVEAN (Protein Variation Effect Analyzer) score. The smaller the 
                                   more deleterious is.

        -pv_p arg (Dis)            PROVEAN predictions.
                                       1 or D: damaging
                                       2 or N: neutral

        -vest arg (Con)            VEST (Variant Effect Scoring Tool) version3 score. The larger the 
                                   larger posibility its pathogenic functional significance is.
                                   It is a machine learning method that predicts the functional 
                                   significance of missense mutations based on the probability that 
                                   they are pathogenic.

        -cadd_r arg (Con)          CADD raw C-scores. The larger the more deleterious is.
                                   CADD can quantitatively prioritize functional, deleterious, and 
                                   disease causal variants across a wide range of functional categories,
                                   effect sizes and genetic architectures and can be used prioritize 
                                   causal variation in both research and clinical settings.

        -cadd_p arg (Con)          CADD phred score: normalized score that can be comparable. 
                                   The larger the more deleterious is.
                                       range: 1-99

        -dann arg (Con)            DANN score. The larger the more deleterious is.
                                       range: 0-1

        -fmm_s arg (Con)           FATHMM score integrates with functional annotations from ENCODE with 
                                   nucleotide-based HMMS. The larger the more deleterious is.

        -fmm_p arg (Dis)           FATHMM predictions integrates with functional annotations from ENCODE 
                                   with nucleotide-based HMMS.
                                       1 or D: damaging
                                       2 or N: neutral

        -metasvm_s (Con)           METASVM score: Ensemble scores of almost all the different popular 
                                   scores methods for deleterious missense mutations. The larger the 
                                   more deleterious is.

        -metasvm_p (Dis)           METASVM predictions.
                                       1 or D: deleterious
                                       2 or T: tolerated

        -metalr_s (Con)            METALR score. The larger the more deleterious is.

        -metalr_p (Dis)            METALR predictions.
                                       1 or D: deleterious
                                       2 or T: tolerated

        -intf (Con)                Integrated fitCons score: fitness consequences of functional 
                                   annotation. The larger the more deleterious is.

        -intc (Con)                Integrated confidence value: fitness cosequences of functional 
                                   annotation. The larger the more deleterious is.

        -gr (Con)                  GERP++ RS score, the larger the score, the more conserved or 
                                   deleterious the site.

        -ppv (Con)                 PhyloP7way vertebrate: phylogentic p-values calculated from a LRT, 
                                   score-based test, GERP test Use 7 species. The larger the more 
                                   deleterious is.

        -ppm (Con)                 PhyloP20way mammalian: a phylogenetic hidden Markov model (phylo-HMM) 
                                   use 20 species. The larger the more deleterious is.

        -pcv (Con)                 PhastCons7way vertebrate: A phylogenetic hidden Markov model 
                                   (phylo-HMM) Use 7 species. The larger the more deleterious is.

        -pcm (Con)                 PhastCons20way mammalian: A phylogenetic hidden Markov model 
                                   (phylo-HMM) Use 20 species. The larger the more deleterious is.

        -sp (Con)                  SiPhy 29way logOdds: Probablistic framework, HMM Use 29 species. 
                                   The larger the more deleterious is.

        -re (Con)                  REVEL score: An Ensemble Method for Predicting the Pathogenicity of 
                                   Rare Missense Variants. The larger the more deleterious is.

        -mcap_s (Con)              Mendelian Clinically Applicable Pathogenicity (M-CAP) Score for rare 
                                   missense variants in the human genome. The larger the deleterious is.

        -mcap_p (Con)              M-CAP predictions.
                                       1 or P: (Likely)Pathogenic
                                       2 or B: (Likely)Benign

        -car (Dis)                 Clinvar: it archives and aggregates information about relationships 
                                   among variation and human health.
                                       1 or P: pathogenic
                                       2 or LP: likely pathogenic
                                       3 or B: benign
                                       4 or LB: likely benign
                                       5 or DR: drug response
                                       6 or U: uncertain significance
                                       7 or N: not provided
                                       8 or O: other

        -omim (Dis)                The OMIM (Online Mendelian Inheritance in Man) dataabase annotation.
                                       1: is correlated with a known mendelian disease
                                       2: not provided
                                       key word: find the OMIM which contain this key word

        -dbsc_ada (Con)            dbscSNV: prediction score of splice-altering single nucleotide 
                                   variants in the human genome. The adaboost model. The larger the 
                                   more possibility it can influence the splicing.

        -dbsc_rf (Con)             dbscSNV with the random forest model. The larget the more
                                   possibility it can influence the splicing.


"""

# map the header of annovar file to the argvs
header_para_dict = {
    'chr':'chr',
    'funcRegion':'reg',
    'exonFunc':'exf',
    'zygosity':'zyg',
    'dbSNP':'snp',
    'ExAC_EAS':'exac_eas',
    'ExAX_all':'exac_all',
    '1000g_EAS':'1kg_eas',
    '1000g_all':'1kg_all',
    'cg':'cg',
    'esp_all':'esp_all',
    'dbscSNV':('dbsc_ada', 'dbsc_rf'),
    'REVEL_score':'re',
    'M_CAP_score':'mcap_s',
    'ClinVar':'car',
    'M_CAP_predict':'mcap_p',
    'OMIM':'omim',
    'SIFT_score':'sift_s',
    'SIFT_pred':'sift_p',
    'Polyphen2_HDIV_score':'pp_hdiv_s',
    'Polyphen2_HDIV_pred':'pp_hdiv_p',
    'Polyphen2_HVAR_score':'pp_hvar_s',
    'Polyphen2_HVAR_pred':'pp_hvar_p',
    'LRT_score':'lrt_s',
    'LRT_pred':'lrt_p',
    'MutationTaster_score':'tas_s',
    'MutationTaster_pred':'tas_p',
    'MutationAssessor_score':'ass_s',
    'MutationAssessor_pred':'ass_p',
    'FATHMM_score':'fm_s',
    'FATHMM_pred':'fm_p',
    'PROVEAN_score':'pv_s',
    'PROVEAN_pred':'pv_p',
    'VEST3_score':'vest',
    'CADD_raw':'cadd_r',
    'CADD_phred':'cadd_p',
    'DANN_score':'dann',
    'fathmm-MKL_coding_score':'fmm_s',
    'fathmm-MKL_coding_pred':'fmm_p',
    'MetaSVM_score':'metasvm_s',
    'MetaSVM_pred':'metasvm_p',
    'MetaLR_score':'metalr_s',
    'MetaLR_pred':'metalr_p',
    'integrated_fitCons_score':'intf',
    'integrated_confidence_value':'intc',
    'GERP++_RS':'gr',
    'phyloP7way_vertebrate':'ppv',
    'phyloP20way_mammalian':'ppm',
    'phastCons7way_vertebrate':'pcv',
    'phastCons20way_mammalian':'pcm',
    'SiPhy_29way_logOdds':'sp',
}

# construct the discrete para trans dictionary
discrete_paras = ['chr', 'zyg', 'reg', 'exf', 'snp', 'sift_p', 'pp_hdiv_p', 'pp_hvar_p', 'lrt_p',
                  'tas_p', 'ass_p', 'fm_p', 'pv_p', 'fmm_p', 'metasvm_p', 'metalr_p', 'mcap_p',
                  'car', 'omim']
discrete_para_options = {x:{} for x in discrete_paras}

discrete_para_options['chr'] = {str(x):'chr' + str(x) for x in range(1, 23)}
discrete_para_options['chr']['X'] = discrete_para_options['chr']['x'] = 'chrX'
discrete_para_options['chr']['Y'] = discrete_para_options['chr']['y'] = 'chrY'
discrete_para_options['chr']['S'] = discrete_para_options['chr']['s'] = ('chrX', 'chrY')
discrete_para_options['chr']['A'] = discrete_para_options['chr']['a'] = tuple(['chr' + str(x) for x in range(1, 23)])

discrete_para_options['zyg']['1'] = discrete_para_options['zyg']['hom'] = 'hom'
discrete_para_options['zyg']['2'] = discrete_para_options['zyg']['het'] = 'het'

discrete_para_options['reg']['1'] = discrete_para_options['reg']['exon'] = 'exonic'
discrete_para_options['reg']['2'] = discrete_para_options['reg']['intro'] = 'intronic'
discrete_para_options['reg']['3'] = discrete_para_options['reg']['utr3'] = 'UTR3'
discrete_para_options['reg']['4'] = discrete_para_options['reg']['utr5'] = 'UTR5'
discrete_para_options['reg']['5'] = discrete_para_options['reg']['inter'] = 'intergenic'
discrete_para_options['reg']['6'] = discrete_para_options['reg']['splicing'] = 'splicing'
discrete_para_options['reg']['7'] = discrete_para_options['reg']['nc'] = ('ncRNA_exonic', 'ncRNA_intronic', 'ncRNA_splicing')
discrete_para_options['reg']['8'] = discrete_para_options['reg']['up'] = 'upstream'
discrete_para_options['reg']['9'] = discrete_para_options['reg']['down'] = 'downstream'

discrete_para_options['exf']['1'] = discrete_para_options['exf']['syn'] = 'synonymous SNV'
discrete_para_options['exf']['2'] = discrete_para_options['exf']['nonsyn'] = 'nonsynonymous SNV'
discrete_para_options['exf']['3'] = discrete_para_options['exf']['fsdel'] = 'frameshift deletion'
discrete_para_options['exf']['4'] = discrete_para_options['exf']['fsins'] = 'frameshift insertion'
discrete_para_options['exf']['5'] = discrete_para_options['exf']['nfsdel'] = 'nonframeshift deletion'
discrete_para_options['exf']['6'] = discrete_para_options['exf']['nfsins'] = 'nonframeshift insertion'
discrete_para_options['exf']['7'] = discrete_para_options['exf']['stopg'] = 'stopgain'
discrete_para_options['exf']['8'] = discrete_para_options['exf']['stopl'] = 'stoploss'
discrete_para_options['exf']['9'] = discrete_para_options['exf']['uk'] = 'unkown'

discrete_para_options['sift_p']['1'] = discrete_para_options['sift_p']['D'] = 'D'
discrete_para_options['sift_p']['2'] = discrete_para_options['sift_p']['T'] = 'T'

discrete_para_options['pp_hdiv_p']['1'] = discrete_para_options['pp_hdiv_p']['D'] = 'D'
discrete_para_options['pp_hdiv_p']['2'] = discrete_para_options['pp_hdiv_p']['P'] = 'P'
discrete_para_options['pp_hdiv_p']['3'] = discrete_para_options['pp_hdiv_p']['B'] = 'B'

discrete_para_options['pp_hvar_p']['1'] = discrete_para_options['pp_hvar_p']['D'] = 'D'
discrete_para_options['pp_hvar_p']['2'] = discrete_para_options['pp_hvar_p']['P'] = 'P'
discrete_para_options['pp_hvar_p']['3'] = discrete_para_options['pp_hvar_p']['B'] = 'B'

discrete_para_options['lrt_p']['1'] = discrete_para_options['lrt_p']['D'] = 'D'
discrete_para_options['lrt_p']['2'] = discrete_para_options['lrt_p']['N'] = 'N'
discrete_para_options['lrt_p']['3'] = discrete_para_options['lrt_p']['U'] = 'U'

discrete_para_options['tas_p']['1'] = discrete_para_options['tas_p']['A'] = 'A'
discrete_para_options['tas_p']['2'] = discrete_para_options['tas_p']['D'] = 'D'
discrete_para_options['tas_p']['3'] = discrete_para_options['tas_p']['P'] = 'P'
discrete_para_options['tas_p']['4'] = discrete_para_options['tas_p']['N'] = 'N'

discrete_para_options['ass_p']['1'] = discrete_para_options['ass_p']['H'] = 'H'
discrete_para_options['ass_p']['2'] = discrete_para_options['ass_p']['M'] = 'M'
discrete_para_options['ass_p']['3'] = discrete_para_options['ass_p']['L'] = 'L'
discrete_para_options['ass_p']['4'] = discrete_para_options['ass_p']['N'] = 'N'

discrete_para_options['fm_p']['1'] = discrete_para_options['fm_p']['D'] = 'D'
discrete_para_options['fm_p']['2'] = discrete_para_options['fm_p']['T'] = 'T'

discrete_para_options['pv_p']['1'] = discrete_para_options['pv_p']['D'] = 'D'
discrete_para_options['pv_p']['2'] = discrete_para_options['pv_p']['N'] = 'N'

discrete_para_options['fmm_p']['1'] = discrete_para_options['fmm_p']['D'] = 'D'
discrete_para_options['fmm_p']['2'] = discrete_para_options['fmm_p']['N'] = 'N'

discrete_para_options['metasvm_p']['1'] = discrete_para_options['metasvm_p']['D'] = 'D'
discrete_para_options['metasvm_p']['2'] = discrete_para_options['metasvm_p']['T'] = 'T'

discrete_para_options['metalr_p']['1'] = discrete_para_options['metalr_p']['D'] = 'D'
discrete_para_options['metalr_p']['2'] = discrete_para_options['metalr_p']['T'] = 'T'

discrete_para_options['mcap_p']['1'] = discrete_para_options['mcap_p']['P'] = '(Likely)Pathogenic'
discrete_para_options['mcap_p']['2'] = discrete_para_options['mcap_p']['B'] = '(Likely)Benign'

discrete_para_options['car']['1'] = discrete_para_options['car']['P'] = 'Pathogenic'
discrete_para_options['car']['2'] = discrete_para_options['car']['LP'] = 'Likely pathogenic'
discrete_para_options['car']['3'] = discrete_para_options['car']['B'] = 'Benign'
discrete_para_options['car']['4'] = discrete_para_options['car']['LB'] = 'Likely benign'
discrete_para_options['car']['5'] = discrete_para_options['car']['DR'] = 'drug response'
discrete_para_options['car']['6'] = discrete_para_options['car']['U'] = 'Uncertain significance'
discrete_para_options['car']['7'] = discrete_para_options['car']['N'] = 'not provided'
discrete_para_options['car']['8'] = discrete_para_options['car']['O'] = 'other'


para_names = header_para_dict.keys()
valid_paras = []
for k, v in header_para_dict.items():
    if type(v) == type(''):
        valid_paras.append(v)
    else:
        valid_paras.extend(list(v))
continuous_paras = [x for x in valid_paras if x not in discrete_paras]
valid_paras += ['in', 'out', 'e', 'col', 's']

import re
import sys

def argv_process(argv):
    """parse the argv into the argv:para dict
        Format:
            argv type: str
            return type: dict
    """
    para_help = re.compile(r'-h|-H|-help|--help')
    res_help = para_help.search(argv)
    if res_help:
        print(usage)
        sys.exit()
    all_para_names = re.findall(r'-(\w+)', argv)
    argv += '-'
    para_dict = {}
    para = re.compile(r'-(\w+)\s+([\w\s,:.]+)(?=-)')
    while 1:
        res = para.search(argv)
        if not res:
            break
        para_key = res.group(1)
        para_value = res.group(2)
        para_value = re.sub(r'[\s,]+','\t', para_value).strip().split('\t')
        para_dict[para_key] = para_value
        argv = argv[:res.start()] + argv[res.end():]
    if 'in' not in para_dict or 'out' not in para_dict:
        print('Lack of Input file or Output file')
        sys.exit()
    unmatch_para_names = [x for x in all_para_names if x not in valid_paras]
    if unmatch_para_names:
        print('Unvalid Argv: ' + ' '.join(unmatch_para_names))
        sys.exit()
    empty_para_names = [x for x in all_para_names if x not in para_dict]
    if empty_para_names:
        print('Empty Argv: ' + ' '.join(empty_para_names))
        sys.exit()
    return para_dict


def para_judge(para_dict, para_name, content):
    """judge whether the content corresponding to the para_name meet the constriants
        Format:
            para_dict: dict
            para_name, content: str
    """
    empty_include = True
    if para_dict.get('e', None) == 2:
        empty_include = False
    paras = header_para_dict[para_name]
    if type(paras) == type(''):
        paras = (paras,) # trans the str to tuple format to unify paras format
    contents = content.split(',')
    if 1 == len(contents) < len(paras): # empty content, span it to the same length of paras
        contents *= len(paras)
    match = True
    for i in range(len(paras)):
        p, c = paras[i], contents[i]
        cs = c.replace('\\x2c', '|').replace(',', '|').split('|')
        if p not in para_dict:
            continue
        if p in discrete_paras:
            op_match = False
            for x in para_dict[p]:
                if p == 'snp':
                    if x == '1':
                        op_match |= c not in ['.', '-', '*', '']
                    elif x == '2':
                        op_match |= c in ['.', '-', '*', '']
                elif p == 'omim':
                    if x == '1':
                        op_match |= c not in ['.', '-', '*', '']
                    elif x == '2':
                        op_match |= c in ['.', '-', '*', '']
                    else:
                        op_match |= x in c # not use cs because the OMIM content and the key word is difficult to split.
                else:
                    option = discrete_para_options[p][x]
                    if type(option) == type(''):
                        option = (option, )
                    for op in option:
                        op_match |= op in cs # use cs because the content has been split and op is well to match.
                    if empty_include:
                        op_match |= c in ['.', '-', '*', '']
            match &= op_match
            if not match:
                return False
        else:
            op_match = False
            for x in para_dict[p]:
                interval = [- float('inf'), float('inf')]
                if ':' in x:
                    interval = [float(y) if y else interval[j] for j, y in enumerate(x.strip().split(':'))]
                else:
                    interval[1] = float(x.strip())
                try:
                    op_match |= (interval[0] <= float(c) <= interval[1])
                except ValueError as e:
                    pass
                if empty_include:
                    op_match |= c in ['.', '-', '*', '']
            match &= op_match
            if not match:
                return False
    return match


def unique(lst1):
   """Remove the duplicates in lst
       Format:
           lst: list
           return: list
   """
   dig_uni = []
   for x in lst1:
        if x in dig_uni:
            continue
        else:
            dig_uni.append(x)
   return dig_uni


def union(lst1, lst2):
    """Union of two list.
        Format:
            lst1, lst2: list
            return: list
    """
    lst_uni = lst1 + lst2
    return unique(lst_uni)


def col_choose(para_dict, header_ind, headers):
    """Choose the select columns indexs to output
        Format:
            para_dict: dict
            header_ind: list of integer
            headers: list of string
    """
    if 'col' not in para_dict:
        return list(range(len(headers)))
    else:
        select = []
        for p in para_dict['col']:
            if ':' in p:
                ps = p.strip().split(':')
                start = int(ps[0]) if ps[0] else header_ind[0] + 1
                end = int(ps[1]) if ps[1] else header_ind[-1] + 1
                select = union(select, list(range(start -1, end)))
            elif p in headers:
                p_ind = headers.index(p)
                if p_ind in header_ind and p_ind not in select:
                    select.append(p_ind)
            elif p in valid_paras:
                for k, v in header_para_dict.items():
                    if (type(v) == type('') and p == v) or (type(v) == type((1,)) and p in v):
                        p_ind = headers.index(k)
                        if p_ind in header_ind and p_ind not in select:
                            select.append(p_ind)
            else:
                try:
                    p_ind = int(p) - 1
                    if p_ind not in select:
                        select.append(p_ind)
                except ValueError as e:
                    pass
        return select
            

def main():
    para_dict = argv_process(' '.join(sys.argv))
    annovar_file = para_dict['in'][0]
    output_file = para_dict['out'][0]
    fi_annovar = open(annovar_file)
    fo = open(output_file, 'w')
    header = fi_annovar.readline()
    headers = header.strip('\n').split('\t')
    header_all_ind = [i for i, x in enumerate(headers) if x in para_names]
    output_col_ind = col_choose(para_dict, header_all_ind, headers)
    if 's' not in para_dict or para_dict['s'] == '2':
        output_col_ind = [x for x in list(range(len(headers))) if x in output_col_ind]
    print(output_col_ind)
    para_use_ind = []
    for i in header_all_ind:
        para = header_para_dict[headers[i]]
        if type(para) == type('') and para in para_dict:
            para_use_ind.append(i)
        else:
            if set(para_dict) & set(para):
                para_use_ind.append(i)
    fo.write('\t'.join([headers[x] for x in output_col_ind]) + '\n')
    for line in fi_annovar:
        items = line.strip('\n').split('\t')
        match = True
        for i in para_use_ind:
             match &= para_judge(para_dict, headers[i], items[i])
             if not match:
                 break
        if match:
             fo.write('\t'.join([items[x] for x in output_col_ind]) + '\n')
    fi_annovar.close()
    fo.close()


if __name__ == '__main__':
    main()
