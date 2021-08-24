#for RHP: investigate R->K muts

import pandas as pd
import sys
import_path_base = '/storage1/hezscha/src/'
#import_path_base = '/home/henrike/Documents/PD_AS/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
from PrismData import PrismParser, VariantData#VariantParser, VariantData
from Bio import AlignIO
import numpy as np
import glob
import os
import argparse

def read_from_prism(primsfile):
	parser = PrismParser()
	dataframe = parser.read(primsfile).dataframe
	meta_data = parser.read_header(primsfile)
	return meta_data, dataframe

def write_prism(metadata, dataframe, prism_file, comment=''):
	variant_dataset = VariantData(metadata, dataframe)
	parser = PrismParser()
	parser.write(prism_file, variant_dataset, comment_lines=comment)

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-outdir', dest="outdir", help="Output directory for 20x20 csv files.")
parser.add_argument('-indir', dest="indir", help="Input directory for merge files.")
parser.add_argument('-swiss', dest="swiss", action='store_true', help="Does the input folder have a substructure by uniprot ID?")
args = parser.parse_args()
################################################################################

order_AA_WT = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'X']
order_AA_sub = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*']
#df that lists the counts of i.e. how many gnomad R->G vars there are 
df_all = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
#df that lists i.e. how many pathogenic gnomad R->G vars there are 
df_patho = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
df_ben = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
#freq(R->G_patho) = count(df_patho_gnomad[R,G])/count(df_all_gnomad[R,G]) . So out of all R->G muts, which fraction is known as pathogenic

#and then later, we will make a matrix where each field is divided by 1/freq(R->K_patho). So we have a field freq(R->G_patho)/freq(R->K_patho) and so on
#the difference to freq(S_to_Y_patho_inside_domain) = counts_S_to_Y_patho_inside_domain/sum(all_patho_vars_inside_domain) is that in that case I just divided each field by the sum of the matrix. Now, I want to divide each field by something different. So I'll use two matrices for that
#freq(R->K_patho) is just a special field inside df_patho_gnomad/df_all_gnomad 

#https://trello.com/c/vlUHag3A/62-interpretation-of-the-question
#get:
#count(R->K_patho), count(R->K_benign)
#freq(R->K_patho) = count(R->K_patho)/count(all_R->K_vars), freq(R->K_benign)
#count(R->K_patho) is all vars that are R->K and have a pathogenic rating, even if they do not have a gnomad freq! 
#count(all_R->K_vars) is all vars that are R->K and either have a gnomad freq or a clinvar rating.  

c_R_K_all = 0
c_R_K_patho = 0
c_R_K_ben = 0


#indir = "/storage1/shared/data/prism_merge"
if args.swiss:
	for d1 in os.scandir(args.indir): 
		try:
			for d2 in os.scandir(d1.path):
				try:
					for d3 in os.scandir(d2.path):
						for prism_file in os.scandir(d3.path):
							#skip over temporary files that exist due to merging
							if prism_file.name.endswith('_filled.txt') or prism_file.name.endswith('_SNPinv.txt') or prism_file.name.endswith('log'):
								continue
							uniprot = prism_file.name.split('_')[-1].split('.')[0]
							print(uniprot)
							
							try:
								metadata,df = read_from_prism(prism_file.path)
							except:
								print('Could not read prism file:', prism_file.path)
								continue	
								
							#get the substutions
							#identify which files have been merged in
							uni_file_nr = ''
							cv_file_nr = ''
							gnom_file_nr = ''
							sai_file_nr = ''
							nsp_file_nr = ''
							
							for File in metadata['merged']:
								if 'prism_uniprot_002' in metadata['merged'][File]:
									uni_file_nr = '_'+File.split('_')[-1]
								if 'prism_clinvar' in metadata['merged'][File]:
									cv_file_nr = '_'+File.split('_')[-1]
								if 'prism_gnomad_001' in metadata['merged'][File]:
									gnom_file_nr = '_'+File.split('_')[-1]
								if 'prism_spliceai' in metadata['merged'][File]:
									sai_file_nr = '_'+File.split('_')[-1]
									#print(File, metadata['merged'][File])
								elif 'prism_netsurfp' in metadata['merged'][File]:
									nsp_file_nr = '_'+File.split('_')[-1]
									
							if gnom_file_nr and cv_file_nr:
								#subset to vars that have either a gnomad freq or a clinvar rating:
								all_vars = df.loc[(df['clinvar;Clinvar_signifiance'+cv_file_nr].notnull()) | (df['gnomad;AF_tot'+gnom_file_nr].notnull())]
								
								for index, row in all_vars.iterrows():
									if row['aa_var'][0] == 'K' and row['aa_ref'][0] == 'R':
										c_R_K_all += 1
									if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
										df_all.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
								
								#subset to only pathogenic vars
								p = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'pathogenic']
								for index, row in p.iterrows():
									if row['aa_var'][0] == 'K' and row['aa_ref'][0] == 'R':
										c_R_K_patho += 1
									if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
										df_patho.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
								
								#subset to only benign vars
								b = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'benign']
								for index, row in b.iterrows():
									if row['aa_var'][0] == 'K' and row['aa_ref'][0] == 'R':
										c_R_K_ben += 1
									if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
										df_ben.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
									
										
				except NotADirectoryError:
					continue				
		except NotADirectoryError:
			continue

else:
	for prism_file in os.scandir(args.indir):
		#skip over temporary files that exist due to merging
		if prism_file.name.endswith('_filled.txt') or prism_file.name.endswith('_SNPinv.txt') or prism_file.name.endswith('log'):
			continue
		uniprot = prism_file.name.split('_')[-1].split('.')[0]
		print(uniprot)
		
		try:
			metadata,df = read_from_prism(prism_file.path)
		except:
			print('Could not read prism file:', prism_file.path)
			continue	
			
		#get the substutions
		#identify which files have been merged in
		uni_file_nr = ''
		cv_file_nr = ''
		gnom_file_nr = ''
		sai_file_nr = ''
		nsp_file_nr = ''
		
		for File in metadata['merged']:
			if 'prism_uniprot_002' in metadata['merged'][File]:
				uni_file_nr = '_'+File.split('_')[-1]
			if 'prism_clinvar' in metadata['merged'][File]:
				cv_file_nr = '_'+File.split('_')[-1]
			if 'prism_gnomad_001' in metadata['merged'][File]:
				gnom_file_nr = '_'+File.split('_')[-1]
			if 'prism_spliceai' in metadata['merged'][File]:
				sai_file_nr = '_'+File.split('_')[-1]
				#print(File, metadata['merged'][File])
			elif 'prism_netsurfp' in metadata['merged'][File]:
				nsp_file_nr = '_'+File.split('_')[-1]
				
		if gnom_file_nr and cv_file_nr:
			#subset to vars that have either a gnomad freq or a clinvar rating:
			all_vars = df.loc[(df['clinvar;Clinvar_signifiance'+cv_file_nr].notnull()) | (df['gnomad;AF_tot'+gnom_file_nr].notnull())]
			
			for index, row in all_vars.iterrows():
				if row['aa_var'][0] == 'K' and row['aa_ref'][0] == 'R':
					c_R_K_all += 1
				if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
					df_all.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
			
			#subset to only pathogenic vars
			p = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'pathogenic']
			for index, row in p.iterrows():
				if row['aa_var'][0] == 'K' and row['aa_ref'][0] == 'R':
					c_R_K_patho += 1
				if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
					df_patho.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
			
			#subset to only benign vars
			b = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'benign']
			for index, row in b.iterrows():
				if row['aa_var'][0] == 'K' and row['aa_ref'][0] == 'R':
					c_R_K_ben += 1
				if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
					df_ben.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
										
f_R_K_patho = c_R_K_patho/c_R_K_all
f_R_K_ben = c_R_K_ben/c_R_K_all

print('c_R_K_all',c_R_K_all)
print('c_R_K_patho',c_R_K_patho)
print('f_R_K_patho',f_R_K_patho)
print('c_R_K_ben',c_R_K_ben)
print('f_R_K_ben',f_R_K_ben)

outfile = os.path.join(args.outdir, "RHP_patho_20x20subs.csv")
df_patho.to_csv(outfile)
outfile = os.path.join(args.outdir, "RHP_ben_20x20subs.csv")
df_ben.to_csv(outfile)
outfile = os.path.join(args.outdir, "RHP_all_20x20subs.csv")
df_all.to_csv(outfile)

