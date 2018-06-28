import argparse
import csv
from datetime import datetime
from itertools import chain
from collections import defaultdict

#----- MyLaptop
# /home/matteo/git/VCF_UKB_Merger/VCF_UKB_Merger.py \
# -I /home/matteo/Scrivania/UKB-Results/UKB-22193-cleaned.tsv \
# -VCF /home/matteo/Scrivania/UKB-Results/GWAS_Results/rs78378222_fetched.vcf \
# -O /home/matteo/Scrivania/UKB-Results/GWAS_Results/UKB-22193-merged-genetic.tsv \
# -HomPhen /home/matteo/Scrivania/UKB-Results/GWAS_Results/rs78378222_Homozyg_UKB.tsv
# -Skip

if __name__ == '__main__':

    parser = argparse.ArgumentParser('This tool merge a vcf file with the UKB dataset. Developed by Matteo Di Giovannantonio for LICR and Department of Oncology of Oxford University.\nFor a detailed description abput the usage, please look at the README file or at the github page.\n\n')
    
    parser.add_argument('-I','--input',help="Path to the UK Biobank input file in tsv format.")
    parser.add_argument('-O','--out',help="Path to output file in tsv format.")
    parser.add_argument('-VCF','--VCFfile',help="Path to the vcf file")
    parser.add_argument('-HomPhen','--HomozygousPhentype',help="Path to the tsv file containing all the homozygous samples for the SNP rs78378222")
    parser.add_argument('-Skip','--SkipID',help="Report in a file in the same path of the out fle the number of skipping ID (no genetic data in .sample file)")

    global opts

    print('\n\nBioinformatic to rule them all, Python to find them, Python to bring them all, and in the darkness bind them\n\n')

    opts = parser.parse_args()

    out_file = open(opts.out,'w')

    if opts.HomozygousPhentype:
        phen_file = open(opts.HomozygousPhentype,'w')
        count_header = 0

    # We save genotypes in a dict
    Sample_list = []
    SNP_ID = {}
    # Skip column containing one fo these fields
    Skip = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    Count_skipping = 0

    with open(opts.VCFfile) as VCF:
        for var in VCF:
            if var.startswith('##'):
                continue
            elif var.startswith('#CHROM'):
                var = var.rstrip().split('\t')
                for elem in var:
                   if elem in Skip:
                       continue
                   else:
                       #Insert the sample ID as key. Then in the next for loop, I will assign to each key all the variants extracted.
                       Sample_list += [elem]

            else:
                var = var.rstrip().split('\t')
                SNP_ID[str(var[2])]={}
                count = 0

                # For loop over the line containing the genotypes (skip the other elements of the line)
                for genot in var[len(Skip):]:

                    genot = genot.split(':')[0]

                    GT = '.'
                    
                    if genot == '0/0':
                        GT = '0'
                    elif genot == '0/1':
                        GT = '1'
                    elif genot == '1/1':
                        GT = '2'

                    # Here I generate a dict in a dict. The first key is a SNP ID and the second the sample ID --> {'rs78378222,17:7571752_T_G': {'3416542': '0/0:1,0,0', '2088321': '0/0:1,0,0', '1828079': '0/0:1,0,0'}}
                    SNP_ID[str(var[2])][Sample_list[count]] = GT
                    count+=1

    with open(opts.input) as UKB:
        for line in UKB:
            line = line.rstrip().split('\t')
            # Here I add the SNP ID to the UKB header
            if line[0].startswith('eid'):
                header = line
                for keys in SNP_ID:
                    header += [keys]
                out_file.write('\t'.join(header)+'\n')

            else:
                for keys in SNP_ID:
                    # Here I extract the genotype of the sample for each snp in the SNP_ID dict and then I merge them in a new line
                    try:
                        line += [SNP_ID[keys][line[0]]]
                    except:
                        Count_skipping += 1
                        line += ['']
                out_file.write('\t'.join(line)+'\n')
                if opts.HomozygousPhentype:
                    if count_header == 0:
                        phen_header = ['eid','Genotype']
                        phen_file.write('\t'.join(phen_header+line[header.index('ICD_Code_0'):header.index('Cancer_Phenotype_16')+1])+'\n')
                        count_header = 1
                    elif line[header.index('rs78378222,17:7571752_T_G')] == '2':
                        new_line = [line[0] + line[header.index('rs78378222,17:7571752_T_G')]] + line[header.index('ICD_Code_0'):header.index('Cancer_Phenotype_16')+1]
                        phen_file.write('\t'.join(new_line)+'\n')



                    #else:
                    #    # Insert the sample ID as key. Then in the next for loop, I will assign to each key all the variants extracted.
                    #    Genotypes[elem] = []
