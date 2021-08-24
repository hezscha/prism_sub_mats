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
#parser.add_argument('-indir', dest="indir", help="Input directory for merge files.")
parser.add_argument('-swiss', dest="swiss", action='store_true', help="Does the input folder have a substructure by uniprot ID?")
args = parser.parse_args()
################################################################################

order_AA_WT = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'X']
order_AA_sub = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*']
#df that lists the counts of i.e. how many gnomad R->G vars there are 
df_all_gnomad = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
#df that lists i.e. how many pathogenic gnomad R->G vars there are 
df_patho_gnomad = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
#freq(R->G_patho) = count(df_patho_gnomad[R,G])/count(df_all_gnomad[R,G]) . So out of all R->G muts, which fraction is known as pathogenic

#and then later, we will make a matrix where each field is divided by 1/freq(R->K_patho). So we have a field freq(R->G_patho)/freq(R->K_patho) and so on
#the difference to freq(S_to_Y_patho_inside_domain) = counts_S_to_Y_patho_inside_domain/sum(all_patho_vars_inside_domain) is that in that case I just divided each field by the sum of the matrix. Now, I want to divide each field by something different. So I'll use two matrices for that
#freq(R->K_patho) is just a special field inside df_patho_gnomad/df_all_gnomad 

indir = "/storage1/shared/data/prism_merge"
for d1 in os.scandir(indir): 
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
							for index, row in df.iterrows():
								if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
									df_all_gnomad.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
							
							#subset to only pathogenic vars
							p = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'pathogenic']
							for index, row in p.iterrows():
								if not row['aa_var'][0] == '*' and not row['aa_var'][0] == '=' and not (row['resi'][0] == 1):
									df_patho_gnomad.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
									
			except NotADirectoryError:
				continue				
	except NotADirectoryError:
		continue

outfile = os.path.join(args.outdir, "RHP_gnomad_patho_20x20subs.csv")
df_patho_gnomad.to_csv(outfile)
outfile = os.path.join(args.outdir, "RHP_gnomad_all_20x20subs.csv")
df_all_gnomad.to_csv(outfile)

