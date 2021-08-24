
#go through all marks genes merge files and count subs in all gnomad vars, across all allel freqs

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
args = parser.parse_args()
################################################################################

prism_dir = args.indir
#prism_dir = "/storage1/hezscha/marks_disease_genes/data/merge_files"
#prism_dir = "/storage1/hezscha/marks_disease_genes/test/"
#find all prism file in the dir
listing = glob.glob(os.path.join(prism_dir, '*.all.txt'))
outfile = os.path.join(args.outdir, "allAFs_20x20subs.csv")

#df for counting subs
order_AA_WT = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
order_AA_sub = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*']
df_sub = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
#df_sub_patho = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)

for prism_file in listing:
	#skip over temporary files that exist due to merging
	if prism_file.endswith('_filled.txt') or prism_file.endswith('_SNPinv.txt'):
		continue
	p_name = prism_file.split('/')[-1]
	uniprot = p_name.split('_')[-1].split('.')[0]
	print(p_name)

	metadata,df = read_from_prism(prism_file)
	
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
			
	#if there is a clinvar component file		
	if gnom_file_nr:
		#subset to vars on which we have gnomad data, i.e. tot_AF not NA
		df.loc[df['AF_tot'].notnull()]
		
		#################################
		b = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'benign']
		p = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'pathogenic']

		#iterate through vars in subset df and count which subst we observe how many times. Each var counts for 1
		for index, row in b.iterrows():
			df_sub_benign.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
		for index, row in p.iterrows():
			df_sub_patho.loc[row['aa_var'][0],row['aa_ref'][0]] += 1

#write csv
outfile = os.path.join(args.outdir, "benign_20x20subs.csv")
df_sub_benign.to_csv(outfile)#,header = False, index = False)
outfile = os.path.join(args.outdir, "patho_20x20subs.csv")
df_sub_patho.to_csv(outfile)#,header = False, index = False)
