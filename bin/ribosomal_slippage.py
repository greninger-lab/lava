import subprocess 
import argparse
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio import Entrez
import shutil
import sys
import time
from datetime import datetime
import re
import os.path
import pandas as pd
import sys
import math

correction_number = 0
residue_correction_number = 0

# -1 for coronavirus, HIV
slippage_number = -1

# Grabs the start of the CDS for the original protein before ribosomal slip
slippage_cd_start = open("ribosomal_start.txt").read()

for line in open("proteins.csv"):
    if '_ribosomal_slippage' in line:
        slip_site = line.split(',')[1]
        residue_correction_number = int(slip_site) - int(slippage_cd_start) - slippage_number

        # For some reason, nucleotide counting is off but residue number is correct.
        # For now, adding back protein start is vaguely correct.
        correction_number = int(slip_site) - 1


# Ribosomal_corrected has corrected ribosomal slippage annotations.
# Visualization is what will be fed into Bokeh, which will include old and new annotations.
ribosomal_corrected = open('final_corrected_slippage.csv', 'w')
visualization = open('visualization.csv', 'w')

correction_number = int(correction_number)
residue_correction_number = int(residue_correction_number)

#visualization.write("STUDYID,USUBJID,NGSPL,Visit,SAMPLEID,clade,gene,AAPOS,AAREF,AASUB,TCOV,VCOV,AAFREQ,snpid,nsp,NSPPOS,NSPREF,NSPSUB\n")

# Looping through all mutations we found
with open("filtered_variants.txt") as f:
    # Adds headers as appropriate to each file
    #header = f.readline()
    #ribosomal_corrected.write(header.rstrip() + ",MatPeptide" + '\n')
    #visualization.write(header.rstrip() + ",NucCorrect,AminoCorrect,MatPeptide" + '\n')

    for row in f:
        line = row.rstrip()
        # Finds the nucleotide and amino acid numbers that need to be changed.
        # Formatting is different for deletions because of extra 'del.'
        type = line.split(',')[-2]
        amino_ref = line.split(',')[4]
        amino_alt = line.split(',')[5]
        nuc_num = line.split(',')[-1]
        position = int(line.split(',')[2])
        nuc = line.split(',')[6]

        # position = line.split(',')[2]
            # if position < int(start_num_list[index]):
            #     index = index + 1

        mat_peptide = "-"
        mat_peptide2 = ""
        # Going through and checking which mature peptide it falls under
        for mature_peptide in open("mat_peptides_additions.txt"):
            mat_name = mature_peptide.split(',')[0]
            mat_start = int(mature_peptide.split(',')[1])
            mat_end = int(mature_peptide.split(',')[2])
            mat_correction = int(mature_peptide.split(',')[3])

            # Check to see if mutation falls within this mature peptide
            if position >= mat_start and position <= mat_end:
                # Subtracts the difference between mature peptide and start of protein
                mat_nuc_num = position - mat_start + 1
                if (mat_name == "RNA-dependent_RNA_polymerase_rib_26"):
                    mat_nuc_num = mat_nuc_num + 27
                mat_aa_num = math.ceil(mat_nuc_num/3)
                
                # Grabs correct amino acid mutation.
                mat_aa = amino_ref + str(mat_aa_num) + amino_alt

                # Writes full mature peptide annotation.
                mat_peptide = mat_name + "," + str(mat_aa_num) + "," + amino_ref + "," + amino_alt
                mat_peptide2 = mat_name
                print(mat_peptide)

        if (line.split(",")[1]=="ORF1ab_polyprotein_ribosomal_slippage"):
            gene_name = "ORF1ab_polyprotein"
        else:
            gene_name = line.split(",")[1]

        
        # Corrects for ribosomal slippage by adding correction_number to 
        # original nucleotide/residue number.
        if '_ribosomal_slippage' in line:
            # Round up if decimal, which should only be if ribosomal slippage happens?
            amino_replacement = int(nuc_num) + residue_correction_number
            amino_replacement = math.ceil(amino_replacement / 3)
            #amino_replacement = amino[0] + str(amino_replacement) + amino[-1]
            #new_line = new_line.replace(amino, amino_replacement)

            #visualization_line = visualization_line.replace(amino, amino_replacement)
            #visualization_line = visualization_line + "," + amino_replacement

            #ribosomal_corrected.write(new_line + "," + mat_peptide + '\n')
            #STUDYID,USUBJID,NGSPL,Visit,SAMPLEID,clade,gene,AAPOS,AAREF,AASUB,TCOV,VCOV,AAFREQ,snpid,nsp,NSPPOS,NSPREF,NSPSUB
            amino_pos = str(amino_replacement)
        else:
            amino_pos = line.split(',')[3]
        
        #SAMPLEID,gene,AAPOS,AAREF,AASUB,TCOV,VCOV,AAFREQ,snpid,NTPOS,nsp,NSPPOS,NSPREF,NSPSUB    
        line = line.split(",")[0] + "," + line.split(",")[1] + "," + amino_pos + "," + amino_ref + "," + amino_alt + "," + line.split(",")[8] + "," + line.split(",")[10] + "," + line.split(",")[7] + "," + str(position) + "," + nuc
        visualization.write(line + "," + mat_peptide + "\n")


### Previous code prior to Novavax variant changes

# correction_number = 0
# residue_correction_number = 0

# # -1 for coronavirus, HIV
# slippage_number = -1

# # Grabs the start of the CDS for the original protein before ribosomal slip
# slippage_cd_start = open("ribosomal_start.txt").read()

# for line in open("proteins.csv"):
#     if '_ribosomal_slippage' in line:
#         slip_site = line.split(',')[1]
#         residue_correction_number = int(slip_site) - int(slippage_cd_start) - slippage_number

#         # For some reason, nucleotide counting is off but residue number is correct.
#         # For now, adding back protein start is vaguely correct.
#         correction_number = int(slip_site) - 1


# # Ribosomal_corrected has corrected ribosomal slippage annotations.
# # Visualization is what will be fed into Bokeh, which will include old and new annotations.
# ribosomal_corrected = open('final_corrected_slippage.csv', 'w')
# visualization = open('visualization.csv', 'w')

# correction_number = int(correction_number)
# residue_correction_number = int(residue_correction_number)

# visualization.write("Sample,Position,Protein,AAChange,NucleotideChange,AlleleFreq,Depth,Type,MatPeptide,MatPeptideAAChange,MatPeptideNucChange\n")

# # Looping through all mutations we found
# with open("filtered_variants.txt") as f:
#     # Adds headers as appropriate to each file
#     #header = f.readline()
#     #ribosomal_corrected.write(header.rstrip() + ",MatPeptide" + '\n')
#     #visualization.write(header.rstrip() + ",NucCorrect,AminoCorrect,MatPeptide" + '\n')

#     for row in f:
#         print(line)
#         line = row.rstrip()
#         # Finds the nucleotide and amino acid numbers that need to be changed.
#         # Formatting is different for deletions because of extra 'del.'
#         type = line.split(',')[-1]
#         nuc = line.split(',')[4]
#         if 'del' in nuc or 'dup' in nuc:
#             nuc_num = int(nuc[0:-4])
#         elif type == "frameshift insertion" or "ins" in nuc:
#             nuc_num = int(nuc.split("_")[0])
#         else:
#             nuc_num = int(nuc[1:-1])

#         nuc_ref = line.split(',')[5]
#         nuc_alt = line.split(',')[6]
#         amino = line.split(',')[3]
#         # fs
#         if type == "frameshift deletion" or type == "frameshift insertion":
#             amino_num = int(amino[1:-2])
#         elif "delins" in amino:
#             amino_num = int(amino.split("delins")[0][1:])
#         elif type == "nonframeshift deletion":
#             amino_num = int(amino.split("_")[0][1:])
#         else:
#             amino_num = int(amino[1:-1])
        
#         position = int(line.split(',')[2])

#         snpid = nuc_ref + str(position) + nuc_alt
#         print(snpid)

#         # position = line.split(',')[2]
#             # if position < int(start_num_list[index]):
#             #     index = index + 1

#         mat_peptide = "-"
#         # Going through and checking which mature peptide it falls under
#         for mature_peptide in open("mat_peptides_additions.txt"):
#             mat_name = mature_peptide.split(',')[0]
#             mat_start = int(mature_peptide.split(',')[1])
#             mat_end = int(mature_peptide.split(',')[2])
#             mat_correction = int(mature_peptide.split(',')[3])

#             # Check to see if mutation falls within this mature peptide
#             if position >= mat_start and position <= mat_end:
#                 # Subtracts the difference between mature peptide and start of protein
#                 mat_nuc_num = position - mat_start + 1
#                 if (mat_name == "RNA-dependent_RNA_polymerase_rib_26"):
#                     mat_nuc_num = mat_nuc_num + 27
#                 mat_aa_num = math.ceil(mat_nuc_num/3)

#                 # Grabs correct nucleotide annotation.
#                 if 'del' in nuc:
#                     mat_nuc = str(mat_nuc_num) + nuc[-4:]
#                 else:
#                     mat_nuc = nuc[0] + str(mat_nuc_num) + nuc[-1]
                
#                 # Grabs correct amino acid mutation.
#                 mat_aa = amino[0] + str(mat_aa_num) + amino[-1]

#                 # Writes full mature peptide annotation.
#                 mat_peptide = mat_name + ": " + mat_aa + "; " + mat_nuc

#         if (line.split(",")[1]=="ORF1ab_polyprotein_ribosomal_slippage"):
#             gene_name = "ORF1ab_polyprotein"
#         else:
#             gene_name = line.split(",")[1]

#         # Corrects for ribosomal slippage by adding correction_number to 
#         # original nucleotide/residue number.
#         if '_ribosomal_slippage' in line:
#             if 'del' in nuc:
#                 del_nuc_replacement = str(nuc_num + correction_number) + nuc[-4:]
#                 new_line = line.replace(nuc, del_nuc_replacement)
#                 visualization_line = new_line + "," + del_nuc_replacement
#             else:
#                 nuc_replacement = nuc[0] + str(nuc_num + correction_number) + nuc[-1]
#                 new_line = line.replace(nuc, nuc_replacement)
#                 visualization_line = new_line + "," + nuc_replacement

#             # Round up if decimal, which should only be if ribosomal slippage happens?
#             amino_replacement = int(nuc_num) + residue_correction_number
#             amino_replacement = math.ceil(amino_replacement / 3)
#             amino_replacement = amino[0] + str(amino_replacement) + amino[-1]
#             new_line = new_line.replace(amino, amino_replacement)

#             visualization_line = visualization_line.replace(amino, amino_replacement)
#             visualization_line = visualization_line + "," + amino_replacement

#             ribosomal_corrected.write(new_line + "," + mat_peptide + '\n')
#             visualization.write(visualization_line + "," + mat_peptide + '\n')
#         else:
#             ribosomal_corrected.write(line + "," + mat_peptide + '\n')
#             visualization.write(line + "," + nuc + "," + amino + "," + mat_peptide + '\n')