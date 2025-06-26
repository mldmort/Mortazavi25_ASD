# you should be in an conda/micromamba environment with pandas
import pandas as pd
import argparse
import numpy as np

def get_args():
	# Create ArgumentParser object
	parser = argparse.ArgumentParser(description='Make a combined table from selected variants and scores.')
	
	# Add arguments
	parser.add_argument('--file_out', help='Output file name.', default='combined_table.tsv', type=str)
	parser.add_argument('--file_scores', help='Input scores file name.', default='scores.tsv', type=str)
	parser.add_argument('--file_vep', help='Input vep file name.', default='vep.tsv', type=str)
	
	# Parse the command-line arguments
	args = parser.parse_args()

	return args

args = get_args()
output_file = args.file_out
score_file = args.file_scores
vep_file = args.file_vep
print(f'output_file: {output_file}')
print(f'score_file: {score_file}')
print(f'vep_file: {vep_file}')

meta_file = '../../data/REACH_sample_info.tsv'
df_meta = pd.read_table(meta_file, sep='\t', header=0)
print(df_meta)

aff_dict = {}
for sample, aff in zip(df_meta['Sample_ID'].tolist(), df_meta['Affected'].tolist()):
	aff_dict[sample] = aff
#print(aff_dict)

def get_case(Ser):
	Ser_split = Ser.str.split(',')
	Ser_new = Ser_split.apply(lambda x: [aff_dict[xx.strip('_IL').strip('_PB').strip('_ONT')] if xx.strip('_IL').strip('_PB').strip('_ONT') in aff_dict != "." else "." for xx in x])
	return Ser_new.str.join(sep=',')

def get_case_2(Ser):
	Ser_split = Ser.str.split(',')
	Ser_new = Ser_split.apply(lambda x: [aff_dict[xx.split(':')[0].strip('_IL').strip('_PB').strip('_ONT')] if xx.split(':')[0].strip('_IL').strip('_PB').strip('_ONT') in aff_dict else "." for xx in x])
	return Ser_new.str.join(sep=',')

df_v = pd.read_table(vep_file, sep='\t', header=0)
df_v['key'] = df_v['CHROM'] + "_" + df_v['POS'].astype(str) + "_" + df_v['END'].astype(str) + "_" + df_v['ID']
print('vep table:')
print(df_v)

df_s = pd.read_table(score_file, sep='\t', header=0)
df_s['key'] = df_s['CHROM'] + "_" + df_s['POS'].astype(str) + "_" + df_s['END'].astype(str) + "_" + df_s['ID']
print('score table:')
print(df_s)

df = pd.merge(df_v, df_s, on='key', how='inner')
print('merged table:')
print(df)

print('reduce redundant variant annotations...')
# for variant columns, with redundant ','-separated repeated values, reduce the column to one.

# this must be present even when you have only one ALT, these bed-based annotations are repeated with "," for each gene.
var_cols = ['GENCODE', 'GNOCCHI', 'FB_PR_ENH_M', 'FB_PR_ENH_F', 'FANTOM_ENH', 'EV_CONS_EL']

for var_col in var_cols:
	df[var_col] = df[var_col].str.split(',',expand=False).str.get(0)

print('prioritize GENCODE/ENCODE annotations...')
# for GENCODE column, with functional annotations output the highest priority annotation.
# the rank of annotations is according to the following list, with decreasing priority.
#func_el_rank = ['start_codon', 'stop_codon', 'stop_codon_redefined_as_selenocysteine', 'CDS', 'TSS', 'exon', 'five_prime_UTR', 'three_prime_UTR', 'DNase-H3K4me3', 'PLS', 'pELS', 'dELS', 'CTCF-only']
func_el_rank = ['start_codon', 'stop_codon', 'stop_codon_redefined_as_selenocysteine', 'CDS', 'TSS', 'five_prime_UTR', 'three_prime_UTR', 'exon', 'DNase-H3K4me3', 'PLS', 'pELS', 'dELS', 'CTCF-only']
func2rank = {el: i for i, el in enumerate(func_el_rank)}

def get_priority(func_str):
	if func_str == '.':
		return '.'
	func_list = func_str.split('&')
	#func_rank_list = [(f, func2rank[f]) for f in func_list]
	func_rank_list = [(f, func2rank[f]) for f in func_list if f != 'gene']
	func_rank_list_sorted = sorted(func_rank_list, key=lambda x: x[1], reverse=False) # sorts from small to large
	# this is for cases with only "gene" annotation
	if len(func_rank_list_sorted) == 0:
		return '.'
	return func_rank_list_sorted[0][0]

df['GENCODE'] = df['GENCODE'].apply(get_priority)

print('reduce redundant FB and FANTOM annotations...')
# collapse repeated annotation for the columns below
df['FB_PR_ENH_M'] = df['FB_PR_ENH_M'].apply(lambda x: '&'.join(set(x.split('&'))))
df['FB_PR_ENH_F'] = df['FB_PR_ENH_F'].apply(lambda x: '&'.join(set(x.split('&'))))
df['FANTOM_ENH'] = df['FANTOM_ENH'].apply(lambda x: '&'.join(set(x.split('&'))))

print('rename and drop extra columns...')
rename_mapper = {'CHROM_x': 'CHROM', 'POS_x': 'POS', 'END_x': 'END', 'ID_x': 'ID'}
df.rename(columns=rename_mapper, inplace=True)

drop_cols = ['key' , 'CHROM_y', 'POS_y', 'END_y', 'ID_y']
df.drop(drop_cols, inplace=True, axis=1)

print('choose maximum/mean Gnocchi value...')
# choose maximum Gnocchi value for the variant
def get_gnocchi_max(gnocchi_str):
	if gnocchi_str == '.':
		return '.'
	return max([float(x.split('_')[-1]) for x in gnocchi_str.split('&')])

# choose mean Gnocchi value for the variant
def get_gnocchi_mean(row):
	if row.GNOCCHI == '.':
		return '.'
	pos0 = int(row.POS) - 1
	end = int(row.END)
	# this is the total interval which may not be annotated by GNOCCHI in some places
	#interval = end - pos0
	gnocchi = [float(x.split('_')[-1]) for x in row.GNOCCHI.split('&')]
	poss = [max([int(x.split('_')[1]), pos0]) for x in row.GNOCCHI.split('&')]
	ends = [min([int(x.split('_')[2]), end]) for x in row.GNOCCHI.split('&')]
	# this is the tatal lenght which is annotated by GNOCCHI
	interval = np.sum([e - p for p, e in zip(poss, ends)])
	mean = np.sum([g * (e-p) / interval for g, p, e in zip(gnocchi, poss, ends)])
	mean = f'{mean:.2f}'
	return mean

df['GNOCCHI_MAX'] = df['GNOCCHI'].apply(lambda gno: get_gnocchi_max(gno))
df['GNOCCHI_MEAN'] = df[['POS', 'END', 'GNOCCHI']].apply(lambda row: get_gnocchi_mean(row), axis=1)

print('process evolotionary conserved elements...')
def get_ev_cons_el(row):
	if row.EV_CONS_EL == '.':
		return '.', '.'
	pos0 = int(row.POS) - 1
	end = int(row.END)
	interval = end - pos0
	poss = [max([int(x.split('_')[1]), pos0]) for x in row.EV_CONS_EL.split('&')]
	ends = [min([int(x.split('_')[2]), end]) for x in row.EV_CONS_EL.split('&')]
	sum_bp = np.sum([e - p for p, e in zip(poss, ends)])
	frac = sum_bp / interval
	return sum_bp, f'{frac:.2f}'
	
df[['EV_CONS_EL_BP', 'EV_CONS_EL_FRAC']] = df[['POS', 'END', 'EV_CONS_EL']].apply(lambda row: get_ev_cons_el(row), axis=1, result_type='expand')

df[['case_DENOVO', 'case_DENOVO_MAT', 'case_DENOVO_PAT']] = df[['DENOVO', 'DENOVO_MAT', 'DENOVO_PAT']].apply(get_case, axis=0)

df[['case_ZS_SAMPLES', 'case_ZS_SAMPLES_PB', 'case_ZS_SAMPLES_ONT', 'case_LARGE25_DEV_GB', 'case_LARGE50_DEV_GB', 'case_LARGE100_DEV_GB', 'case_LARGE500_DEV_GB']] = df[['ZS_SAMPLES', 'ZS_SAMPLES_PB', 'ZS_SAMPLES_ONT', 'LARGE25_DEV_GB', 'LARGE50_DEV_GB', 'LARGE100_DEV_GB', 'LARGE500_DEV_GB']].apply(get_case_2, axis=0)

print('write the final table...')
df.to_csv(output_file, sep='\t', header=True, index=False)

