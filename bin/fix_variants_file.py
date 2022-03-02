# This file takes the bcftools_variants.csv and fixes the complex mutations.

import subprocess 
import argparse
import shutil
import sys
import time
from datetime import datetime
import re
import os.path
import pandas as pd
import sys
import math
import numpy as np

def correct_deletions(df, group_to_correct, correct_depth):
	for row_num, row in (group_to_correct[group_to_correct['Depth'] <= 0.4* correct_depth]).iterrows():
		#print(row)
		df.at[row_num,'Depth'] = correct_depth
		# print(df.at[row_num,'Depth'])
		#df.at[row_num,'AF'] = row['VCOV'] / correct_depth
	return df 

def translate(seq):
	table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X', 
		'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W', 
	} 
	protein = "" 
	if "-" in seq:
		return "fs"
	if len(seq)%3 == 0: 
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3] 
			protein+= table[codon] 
	return protein 

def process_fasta(fasta):
	gene_info = pd.DataFrame(columns=['gene_name','gene_start','nt_seq','aa_seq'])
	gene_name = ""
	gene_start = 0
	for line_num,line in enumerate(open(fasta)):	
		if line_num % 2 == 0:
			gene_name = line.split("transcript:")[1].split("Comment")[0].rstrip()
			gene_start = int(line.split(")")[0].split(":")[-1])
		else:
			gene_info.loc[line_num-1] = [gene_name] + [gene_start] + [line.rstrip()] + [translate(line.rstrip())]
	return gene_info

def find_new_nts(aa_nts, list_of_snpids, first_nt_pos):
	new_nts = list(aa_nts)
	for snpid in list_of_snpids:
		nt_pos = int(snpid[1:-1])
		new_nts[nt_pos - first_nt_pos] = snpid[-1]
	return(''.join(new_nts))

def find_aa_ref(protein,residue):
	rel_gene = gene_info.loc[gene_info['gene_name'] == protein]

	if(residue > len(rel_gene['aa_seq'].values[0])):
		return "?"
	return rel_gene['aa_seq'].values[0][residue-1]

def find_passage(metadata, file_name):
	for line in open(metadata):
		if file_name in line:
			return line.split(",")[1].rstrip()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-name', help="Provide sample name.")
	parser.add_argument('-fasta', help="Provide AT_refGeneMrna.fa file with gene sequences.")
	parser.add_argument('-metadata', help="Provide metadata path")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	sample_name = args.name

	# sample_name = args.name.split(".fastq")[0]
	# if "_R" in sample_name:
	# 	sample_name = sample_name.split("_R")[0]
	passage = find_passage(args.metadata, args.name)

	gene_info = process_fasta(args.fasta)

	######
	## First, we have to correct the AFs and pull out the relevant information from the variants.txt file
	######

	fixed_file = open("filtered_variants.txt", "w+")

	for line in open("variants.txt"):
		dp4=line.split("DP4=")[1].split(";")[0]
		ad = line.split("AD=")[1].split(";")[0]
		#allele_ref = int(ad.split(",")[0]) + int(ad.split(",")[1])
		#allele_alt = int(ad.split(",")[2]) + int(ad.split(",")[3])
		allele_ref = int(ad.split(",")[0])
		if (int(dp4.split(",")[2]) + int(dp4.split(",")[3])) == 0:
			allele_alt = 0
		else:
			allele_alt = int(ad.split(",")[1])
		
		fixed_depth = int(dp4.split(",")[0]) + int(dp4.split(",")[1]) + int(dp4.split(",")[2]) + int(dp4.split(",")[3])

		if (allele_ref + allele_alt) > 0 and allele_alt > 0:
			if("IMF" in line):
				af = float(line.split("IMF=")[1].split(";")[0])
			else:
				af = allele_alt / fixed_depth
			
			type = line.split("\t")[1]

			if(af >= 0.01 and len(line.split("\t")[6])<400 and "wholegene" not in line.split("\t")[2]):
				line_parts = line.split("\t")
				nuc_ref = (line_parts[6])
				nuc_alt = (line_parts[7])
				fixed_aa_change = line_parts[2].split(":p.")[1].split(",")[0]
				fixed_protein = line_parts[2].split(":")[1] 
				fixed_nuc_change = line_parts[2].split(":c.")[1].split(":")[0]
				nuc = fixed_nuc_change

				if '_' in nuc:
					nuc_num = nuc.split("_")[0]
				elif 'del' in nuc or 'dup' in nuc:
					nuc_num = int(nuc[0:-4])
				elif type == "frameshift insertion" or "ins" in nuc:
					nuc_num = int(nuc.split("_")[0])
				else:
					nuc_num = int(nuc[1:-1])
				
				nuc_num=str(nuc_num)
									
				# Time to break up some indels into multiple lines
				if(type=="nonframeshift deletion"):
					aa_start_pos = int(fixed_aa_change.split("_")[0])
					nuc_pos = int(line_parts[4])
					
					for codon in range(3,len(nuc_ref) + 1, 3):
						split_nuc_ref = (nuc_ref[codon-3:codon])
						split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)

						nuc_change = split_nuc_ref + str(nuc_pos) + "del"
						#                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF              AASUB           NUCCHANGE             AF               DEPTH                                                                   TYPE
						fixed_file.write(sample_name + "," + str(fixed_protein) + "," + str(nuc_pos) + "," + str(aa_start_pos) + "," + split_amino_ref + "," + "-" + "," + nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1] + "," + str(nuc_num) + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
						aa_start_pos += 1
						nuc_pos +=3
						nuc_num = int(nuc_num) + 3 
				elif(type=="nonframeshift insertion"):
					aa_start_pos = int(fixed_aa_change.split("delins")[0][1:])
					split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
					split_amino_alt = fixed_aa_change.split("delins")[1]
					nuc_change = line_parts[4] + "ins" + fixed_nuc_change.split("ins")[1]
					#                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB             NUCCHANGE             AF               DEPTH                                                                   TYPE
					fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," +  nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1]  + "," + nuc_num + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
				
				elif(type == "frameshift insertion"):
					#split_amino_ref = fixed_aa_change[0]
					aa_start_pos = int(fixed_aa_change.split("fs")[0][1:])
					split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
					split_amino_alt = "fs"
					if("dup" in nuc):
						nuc_change = line_parts[4] + "dup" + fixed_nuc_change.split("dup")[1]
					else:
						nuc_change = line_parts[4] + "ins" + fixed_nuc_change.split("ins")[1]
					#                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB               NUCCHANGE             AF               DEPTH                                                                   TYPE
					fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," + nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1] + "," + nuc_num  + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
				
				elif(type == "frameshift deletion"):
					#split_amino_ref = fixed_aa_change[0]
					aa_start_pos = int(fixed_aa_change.split("fs")[0][1:])
					split_amino_ref = find_aa_ref(fixed_protein,aa_start_pos)
					split_amino_alt = "fs"
					nuc_change = nuc_ref + line_parts[4] + "del"

					#                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB               NUCCHANGE             AF               DEPTH                                                                   TYPE
					fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," + nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + line_parts[1] + "," + nuc_num  + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 
				
				else:
					split_amino_ref = fixed_aa_change[0]
					split_amino_alt = fixed_aa_change[-1]
					# Some stoplosses do not have an aasub, so check for this first...
					if fixed_aa_change[-1].isnumeric():
						aa_start_pos = fixed_aa_change[1:]
					else:
						aa_start_pos = fixed_aa_change[1:-1]
					nuc_change = nuc_ref + line_parts[4] + nuc_alt

					if(split_amino_ref=="X"):
						type = "stoploss"
						split_amino_alt = "fs"
					else:
						type = line_parts[1]

					#                SAMPLE_ID              GENE                  GENPOS                   AAPOS                      AAREF                         AASUB           NUCCHANGE             AF               DEPTH                                                             TYPE
					fixed_file.write(sample_name + "," + str(fixed_protein) + "," + line_parts[4] + "," + str(aa_start_pos) + "," + split_amino_ref + "," + split_amino_alt + "," +  nuc_change + "," + str(af) + "," + str(fixed_depth) + "," + str(allele_ref) + "," + str(allele_alt) + "," + type + "," + nuc_num  + "\n")#+ "," + mat_peptide + "," + mat_peptide_nuc_change + "," + mat_peptide_aa_change + "\n") 

	fixed_file.close() 

	######
	## Now, correct for ribosomal slippage 
	######

	correction_number = 0
	residue_correction_number = 0

	# -1 for coronavirus, HIV
	slippage_number = -1

	# Grabs the start of the CDS for the original protein before ribosomal slip
	slippage_cd_start = open("ribosomal_start.txt").read()

	for line in open("proteins.csv"):
		if '_ribosomal_slippage' in line:
			slip_site = line.split(',')[1]
			if "ORF1ab_polyprotein_ribosomal_slippage" in line:
				slippage_cd_start = 266
			residue_correction_number = int(slip_site) - int(slippage_cd_start) - slippage_number

			# For some reason, nucleotide counting is off but residue number is correct.
			# For now, adding back protein start is vaguely correct.
			correction_number = int(slip_site) - 1

	# Ribosomal_corrected has corrected ribosomal slippage annotations.
	# Visualization is what will be fed into Bokeh, which will include old and new annotations.
	ribosomal_corrected = open('final_corrected_slippage.csv', 'w')
	previsualization = open('previsualization.csv', 'w')

	previsualization.write("Sample,AminoAcidChange,Position,AF,Change,Protein,AAPOS,AAREF,AASUB,NucleotideChange,LetterChange,Syn,Depth,VCOV,Passage,MatPeptide,nsp,NSPPOS,NSPREF,NSPSUB\n")

	correction_number = int(correction_number)
	residue_correction_number = int(residue_correction_number)


	# Looping through all mutations we found
	with open("filtered_variants.txt") as f:
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
				amino_pos = str(amino_replacement)
			else:
				amino_pos = line.split(',')[3]
			
			#SAMPLEID,gene,AAPOS,AAREF,AASUB,Depth,VCOV,AF,snpid,NTPOS,nsp,NSPPOS,NSPREF,NSPSUB    
			line = line.split(",")[0] + "," + line.split(",")[1] + "," + amino_pos + "," + amino_ref + "," + amino_alt + "," + line.split(",")[8] + "," + line.split(",")[10] + "," + line.split(",")[7] + "," + str(position) + "," + nuc
			
			# Sample,AminoAcidChange,Position,AF,Change,Protein,AAREF,AASUB,AAPOS,NucleotideChange,LetterChange,Syn,Depth,VCOV,nsp,NSPPOS,NSPREF,NSPSUB
			aa_change = line.split(",")[1] + " " + amino_ref + amino_pos + amino_alt
			change = amino_ref + amino_pos + amino_alt
			letterchange = nuc.split(str(position))[0] + " to " + nuc.split(str(position))[1]

			line = ",".join([line.split(",")[0],aa_change,str(position),line.split(",")[7],change,line.split(",")[1],amino_pos,amino_ref,amino_alt,nuc,letterchange,"",line.split(",")[5],line.split(",")[6],passage,""])
			previsualization.write(line + "," + mat_peptide + "\n")

	previsualization.close()

	######
	## Finally, fix some incorrect things
	######

	output_variants_file = args.name + ".csv"

	df = pd.read_csv("previsualization.csv").fillna("-")

	# First, let's correct some variants that are obviously wrong
	# Grabbing variants at the same aa position
	for name,group in df.groupby(["Protein","AAPOS"]):
		if group.shape[0] >1:
			for nuc,nuc_group in group.groupby(["Position"]):
				# Here we are protecting against real deletions, which mess up depths/afs of surrounding variants
				correct_depth = max(nuc_group["Depth"])
				df = correct_deletions(df, nuc_group, correct_depth)
				# print(df.iloc[row_num])
			
			# Here we protect against whole amino acid deletions, such as the spike 143-144 deletion
			if "-" in group['AASUB'].unique():
				group = group.sort_values(by=['VCOV'], ascending=False)
				max_depth = max(group['Depth'])
				real_depth = group[group['AASUB']=="-"]['Depth'].values[0]
				if real_depth == max_depth:
					df = correct_deletions(df, group, real_depth)

	# Now we recalculate AF and refilter by 0.01
	df['AF'] = df['VCOV'] / df['Depth']
	df = df.loc[df['AF'] >= 0.01]

	# Now for the actual complex mutations part!
	for name,group in df.groupby(["Protein","AAPOS"]):
		# We don't want to bother with super low-frequency variants
		filtered_group = group.loc[group['AF'] >= 0.05]
		filtered_group = filtered_group.loc[group['AASUB']!= "-"]
		filtered_group = filtered_group.loc[group['AASUB']!= "fs"]
		filtered_group = filtered_group.loc[~(filtered_group.AASUB.str.len()==1)]
		
		if filtered_group.shape[0] > 1:
			filtered_group = filtered_group.sort_values(['AF'], ascending=False)
			for i in range(len(filtered_group.index) - 1):
				j = i+1
				list_to_combine = pd.DataFrame()
				list_to_combine = list_to_combine.append(filtered_group.iloc[i])
				given_nuc = filtered_group.iloc[i]["Position"]

				# This list is sorted from high to low AF. So go through and add mutations that could combine with
				# the given mutation at i that are within 5% frequency.
				while(j < len(filtered_group.index) and filtered_group.iloc[i]['AF'] - filtered_group.iloc[j]['AF'] <=0.05):
					compare_nuc = filtered_group.iloc[j]["Position"]
					if given_nuc != compare_nuc:
						list_to_combine = list_to_combine.append(filtered_group.iloc[j])

					j+=1

				# If there are multiple mutations at the same nucleotide, we have no way
				# of telling which one is the complex one. So we drop these and report separately
				list_to_combine = list_to_combine.drop_duplicates(subset=["Position"],keep=False).sort_values(["Position"])

				# Ok we're done filtering... now the actual complex mutation part
				if list_to_combine.shape[0] > 1:
					# print(list_to_combine)
					# Grab relevant gene info
					rel_gene = gene_info.loc[gene_info['gene_name'] == name[0]]

					# Grab the nt position of first nt in relevant aa
					nt_pos_in_gene = int(name[1]) * 3
					first_nt_pos = nt_pos_in_gene + (rel_gene['gene_start'] - 2).values[0]
					if name[0] == "ORF1ab_polyprotein_ribosomal_slippage":
						first_nt_pos = first_nt_pos - 13203
						nt_pos_in_gene = nt_pos_in_gene - 13203

					# Grab nts
					aa_nts = rel_gene['nt_seq'].values[0][nt_pos_in_gene-3:nt_pos_in_gene]
					
					# Grab new changes
					new_nts = find_new_nts(aa_nts,(list(list_to_combine["NucleotideChange"])), first_nt_pos)
					new_aa = translate(new_nts)
					new_snpid = aa_nts + str(first_nt_pos) + new_nts 
					
					# Implement new changes in original dataframe
					row_to_change = list_to_combine.iloc[0].name
					df.at[row_to_change,'NSPSUB']=new_aa
					df.at[row_to_change,'AASUB']=new_aa
					df.at[row_to_change,"NucleotideChange"]=new_snpid
					df.at[row_to_change,"Syn"]="complex"
					df.at[row_to_change,"LetterChange"]=aa_nts + " to " + new_nts

					# Get rid of other rows old dataframe
					for i in range(1,list_to_combine.shape[0]):
						df = df.drop([list_to_combine.iloc[i].name], errors = "ignore")
					
					# print(df.loc[row_to_change])
					#print(df.loc[df['AAPOS']==name[1]])

	df['AF'] = df['AF'] * 100
	conditions = [
		(df['Syn']=="complex"),
		(df['AAREF']==df['AASUB']),
		(df['AAREF']=="X"),
		(df['AASUB']=="X"),
		(df['AASUB']=="fs"),
		(df['AASUB']=="-"),
		(df['AASUB'].str.len()>1),
		(df['AAREF']!=df['AASUB'])
		]

	values = ["complex","synonymous SNV","stoploss","stopgain","frameshift","deletion","insertion","nonsynonymous SNV"]

	df['Syn'] = np.select(conditions,values)
	df['MatPeptide'] = df['nsp'] + ": " + df['NSPREF'] + df['NSPPOS'].astype(str) + df['NSPSUB']
	df.loc[df.nsp=="-", "MatPeptide"] = "-"
	df['MatPeptide'] = df['MatPeptide'].str.replace(r'.0','')

	df = df.drop(columns=['AAPOS','AAREF','AASUB','VCOV','NSPSUB','NSPPOS','NSPREF','nsp'])

	# Drop rows with empty values - happens with 3-nucleotide complex changes that have already been covered by the 2-nt case
	df = df.dropna(subset=['Sample'])

	df.to_csv(output_variants_file, index=False)
