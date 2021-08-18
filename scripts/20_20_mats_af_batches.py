
#20x20 mats for the different gnomad allel freq batches

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
df_sub_b1 = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
df_sub_b2 = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
df_sub_b3 = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
#df_sub_folded = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)

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
									
							if gnom_file_nr:
								b1 = df.loc[df['gnomad;AF_tot'+gnom_file_nr] < 1e-04]
								for index, row in b1.iterrows():
									if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
										df_sub_b1.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
											
								b2 = df.loc[(df['gnomad;AF_tot'+gnom_file_nr] >= 1e-04) & (df['gnomad;AF_tot'+gnom_file_nr] < 0.01)]
								for index, row in b2.iterrows():
									if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
										df_sub_b2.loc[row['aa_var'][0],row['aa_ref'][0]] += 1

								b3 = df.loc[df['gnomad;AF_tot'+gnom_file_nr] >= 0.01]
								for index, row in b3.iterrows():
									if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
										df_sub_b3.loc[row['aa_var'][0],row['aa_ref'][0]] += 1

				except NotADirectoryError:
					continue				
		except NotADirectoryError:
			continue

outfile = os.path.join(args.outdir, "gnomad_b1_nss_20x20subs.csv")
df_sub_b1.to_csv(outfile)
outfile = os.path.join(args.outdir, "gnomad_b2_nss_20x20subs.csv")
df_sub_b2.to_csv(outfile)
outfile = os.path.join(args.outdir, "gnomad_b3_nss_20x20subs.csv")
df_sub_b3.to_csv(outfile)
