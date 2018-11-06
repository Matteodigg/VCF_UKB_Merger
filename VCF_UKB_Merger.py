import argparse
import csv
from datetime import datetime
from itertools import chain
from collections import defaultdict

#----- MyLaptop
# python3 /home/matteo/git/VCF_UKB_Merger/VCF_UKB_Merger.py \
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
    #parser.add_argument('-Skip','--SkipID',help="Report in a file in the same path of the out file the number of skipping ID (no genetic data in .sample file)")

    global opts

    print('\n\nBioinformatic to rule them all, Python to find them, Python to bring them all, and in the darkness bind them\n\n')

    opts = parser.parse_args()

    Index_collection_41204 = []
    Index_collection_41203 = []
    Index_collection_41202 = []
    start_41202 = []
    stop_41202 = []
    start_41203 = []
    stop_41203 = []
    start_41204 = []
    stop_41204 = []

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

                    GT = '-1'
                    
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

            new_line = list(line)

            # Here I add the SNP ID to the UKB header
            if line[0].startswith('eid'):
                header = line
                for keys in SNP_ID:
                    header += [keys]
                out_file.write('\t'.join(header)+'\n')

                for elem in line:
                    if '41202' in elem:
                        Index_collection_41202 += [line.index(elem)] 
                        if '41202-0.0' in elem:
                            #Here I can identify the FIRST element having ID of 40006
                            start_41202 = header.index(elem)
                        stop_41202 = header.index(elem)+1
                    

                    if '41203' in elem:
                        Index_collection_41203 += [line.index(elem)] 
                        if '41203-0.0' in elem:
                            #Here I can identify the FIRST element having ID of 40006
                            start_41203 = header.index(elem)
                        stop_41203 = header.index(elem)+1

                    if '41204' in elem:
                        Index_collection_41204 += [line.index(elem)] 
                        if '41204-0.0' in elem:
                            #Here I can identify the FIRST element having ID of 40006
                            start_41204 = header.index(elem)
                        stop_41204 = header.index(elem)+1

                #print(start_41202,start_41203)

            else:
                for keys in SNP_ID:
                    # Here I extract the genotype of the sample for each snp in the SNP_ID dict and then I merge them in a new line
                    if line[0] in SNP_ID[keys]:
                        line += [SNP_ID[keys][line[0]]]
                        #print([SNP_ID[keys]])
                        #print(len(line))
                    else:
                        # If the genotype is missing, I set the value as -2 (I want all the values in the column to be an integer)
                        Count_skipping += 1
                        line += ['-2']
                        #print(line[0],len(line))

                out_file.write('\t'.join(line)+'\n')

                if opts.HomozygousPhentype:
                    if count_header == 0:
                        phen_header = ['eid','Genotype']
                        #print(header.index('ICD_Code_0'))
                        phen_file.write('\t'.join(phen_header+header[header.index('ICD_Code_0'):header.index('Cancer_Phenotype_16')+1] + header[start_41202:stop_41202+1] + header[start_41203:stop_41203+1]+header[start_41204:stop_41204+1])+'\n')
                        count_header = 1
                    elif line[header.index('rs78378222,17:7571752_T_G')] == '2':
                        new_line = [line[0] + '\t' + line[header.index('rs78378222,17:7571752_T_G')]] + line[header.index('ICD_Code_0'):header.index('Cancer_Phenotype_16')+1] + line[start_41202:stop_41202+1] + line[start_41203:stop_41203+1] + line[start_41204:stop_41204+1]
                        phen_file.write('\t'.join(new_line)+'\n')



                    #else:
                    #    # Insert the sample ID as key. Then in the next for loop, I will assign to each key all the variants extracted.
                    #    Genotypes[elem] = []
