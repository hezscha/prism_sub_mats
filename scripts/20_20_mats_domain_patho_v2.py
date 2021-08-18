
#v2: add option to omit start_lost and early termination vars since these are different from substitutions
#count substitutions for patho vars in domains and outside domains to see if there are difference in the types of substitutions
# logodds: log(count_MtoV_domain/count_MtoV_outside)

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
parser.add_argument('-omit_start_stop', dest="omit_start_stop", action='store_true', help="Omit counting start lost variants (M1 vars) and early termination variants (anything to *)?")
# ~ parser.add_argument('-omit_start_lost', dest="omit_start_lost", action='store_true', help="Omit counting start lost variants (M1 vars)?")
# ~ parser.add_argument('-omit_early_stop', dest="omit_early_stop", action='store_true', help="Omit early termination vars(anything to *).")

args = parser.parse_args()
################################################################################

#prism_dir = args.indir
#prism_dir = "/storage1/hezscha/marks_disease_genes/data/merge_files"
#prism_dir = "/storage1/hezscha/marks_disease_genes/test/"
#find all prism file in the dir
#if not args.swiss:
#	listing = glob.glob(os.path.join(prism_dir, '*.all.txt'))

#df for counting subs
order_AA_WT = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'X']
order_AA_sub = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*']
df_sub_dom = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
df_sub_outside = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)

#def update_subst_mats():
#if args.omit_start_lost:
	
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
								if 'prism_netsurfp' in metadata['merged'][File]:
									nsp_file_nr = '_'+File.split('_')[-1]
									
							if uni_file_nr and cv_file_nr:
								#is there a domain column?
								if 'uniprot;pfam_name'+uni_file_nr in df.columns:
									#subset to only pathogenic vars
									p = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'pathogenic']
									
									#subset to resis inside domains
									d_dom = p.loc[p['uniprot;pfam_name'+uni_file_nr].notnull()]
									for index, row in d_dom.iterrows():
										if args.omit_start_stop:
											if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
												df_sub_dom.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
										else:
											df_sub_dom.loc[row['aa_var'][0],row['aa_ref'][0]] += 1		
									
									#subset to resis outside domains
									d_out = p.loc[p['uniprot;pfam_name'+uni_file_nr].isnull()]
									for index, row in d_out.iterrows():
										if args.omit_start_stop:
											if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
												df_sub_outside.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
										else:
											df_sub_outside.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
							
				except NotADirectoryError:
					continue				
		except NotADirectoryError:
			continue						

#else all merge files are in the same top dir
else:
	for prism_file in os.scandir(args.indir):
		#skip over temporary files that exist due to merging
		if prism_file.name.endswith('_filled.txt') or prism_file.name.endswith('_SNPinv.txt'):
			continue
		uniprot = prism_file.name.split('/')[-1].split('_')[-1].split('.')[0]

		try:
			metadata,df = read_from_prism(prism_file.name)
		except:
			print('Could not read prism file:', prism_file.name)	
			
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
			if 'prism_netsurfp' in metadata['merged'][File]:
				nsp_file_nr = '_'+File.split('_')[-1]
			
		if uni_file_nr and cv_file_nr:
		#is there a domain column?
			if 'uniprot;pfam_name'+uni_file_nr in df.columns:
				#subset to only pathogenic vars
				p = df.loc[df['clinvar;Clinvar_signifiance'+cv_file_nr] == 'pathogenic']
				
				#subset to resis inside domains
				d_dom = p.loc[p['uniprot;pfam_name'+uni_file_nr].notnull()]
				for index, row in d_dom.iterrows():
					if args.omit_start_stop:
						if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
							df_sub_dom.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
					else:		
						df_sub_dom.loc[row['aa_var'][0],row['aa_ref'][0]] += 1
						
				#subset to resis outside domains
				d_out = p.loc[p['uniprot;pfam_name'+uni_file_nr].isnull()]
				for index, row in d_out.iterrows():
					if args.omit_start_stop:
						if not row['aa_var'][0] == '*' and not (row['resi'][0] == 1):
							df_sub_outside.loc[row['aa_var'][0],row['aa_ref'][0]] += 1	
					else:
						df_sub_outside.loc[row['aa_var'][0],row['aa_ref'][0]] += 1	
		


if args.omit_start_stop:
	outfile = os.path.join(args.outdir, "patho_dom_nss_20x20subs.csv")
else:
	outfile = os.path.join(args.outdir, "patho_dom_20x20subs.csv")
df_sub_dom.to_csv(outfile)

if args.omit_start_stop:
	outfile = os.path.join(args.outdir, "patho_outside_dom_nss_20x20subs.csv")
else:
	outfile = os.path.join(args.outdir, "patho_outside_dom_20x20subs.csv")

df_sub_outside.to_csv(outfile)


