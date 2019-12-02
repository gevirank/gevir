import os
import sys
import progressbar
import pymongo
import csv
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from collections import OrderedDict
from decimal import Decimal
from scipy import stats
from scipy.stats import fisher_exact
from matplotlib_venn import venn3, venn3_circles
from sklearn.metrics import auc as sklearn_auc
from sklearn.utils import resample
from scipy.stats import spearmanr, pearsonr
from csv_reader import CsvReader
from terminaltables import AsciiTable
from common import MongoDB, is_int, sort_dict_by_values, write_table_to_csv, float_to_sci_str, calculate_confidence_interval
from common import get_gene_canonical_transcript_by_name, calculate_percentiles, proportion_to_percents_str
from common import report_2_x_2_test_table_as_obs_exp_str, is_float
from gnomad_utils import get_xpos, xpos_to_pos, worst_csq_from_csq
from export_data import export_gene_scores
from gene_sets import ScoreSets, OmimSets, EssentialSets, replace_nan_or_string_values_with_max_in_dictionary

# Comment these lines if you don't have Helvetica font
import matplotlib as mpl
mpl.rc('font',family='Helvetica')

#################
### CONSTANTS ###
#################

INCLUDE_GNOMAD_OUTLIERS = True

SOURCE_DATA_FOLDER = './source_data/'
OMIM_TSV = SOURCE_DATA_FOLDER + 'genemap2.txt'

REGIONS_COLLECTION = 'variant_regions'

'''
Variants in segmental duplication (SEGDUP) and low copy repeat (LCR) regions
are marked as non-pass quality in gnomAD 2.0, but in new gnomAD releases (e.g. 2.0.1, 2.1)
these variants are pass quality. Moreover, if whole gene is located inside 
segmental duplication, then all variants in it are marked as "non-pass" quality.
Therefore we decided to include these variants into the analysis
'''
VALID_FILTERS = set(['PASS', 'SEGDUP', 'LCR'])

VALID_CSQS = set(['stop_gained',
				  'frameshift_variant',
				  'stop_lost',
				  'start_lost',
				  'inframe_insertion',
				  'inframe_deletion',
				  'missense_variant',])

FIGURES_FOLDER = './figures/'
OUTPUT_FOLDER = './tables/'

# Pallete was taken from here: 
# http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
'''
Color blind friendly
Orange = #e69d00
Blue = #0072b2
Green = #009e73
'''

# Color blind palette: http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png
# Color mixer: https://trycolors.com/
# Color checker: https://www.color-blindness.com/coblis-color-blindness-simulator/ 

COLOR_PALETTE = OrderedDict()
COLOR_PALETTE['B'] = '#0072b2'  # Blue; Alternative: '#1f78b4'
COLOR_PALETTE['O'] = '#e69d00'  # Orange;  Alternative: '#ff7f00' 
COLOR_PALETTE['G'] = '#009e73'  # Green; Alternative: '#33a02c'
COLOR_PALETTE['R'] = '#e31a1c'  # Red
COLOR_PALETTE['P'] = '#6a3d9a'  # Purple
COLOR_PALETTE['Y'] = '#ffff99'  # Yellow
COLOR_PALETTE['PI'] = '#ff66cc' # Pink
COLOR_PALETTE['BR'] = '#b15928' # Brown;  Alternative: a05d2c
COLOR_PALETTE['DP'] = '#581845'

COLOR_PALETTE['BL'] = '#a6cee3' # Blue Light
COLOR_PALETTE['OL'] = '#fdbf6f' # Orange Light
COLOR_PALETTE['GL'] = '#b2df8a' # Green Light
COLOR_PALETTE['RL'] = '#ff7374' # Red Light
COLOR_PALETTE['PL'] = '#cab2d6' # Purple Light

COLOR_PALETTE['RD'] = '#a02c2c' # Red Dark

# Color blind friendly palette
COLOR_PALETTE['SB'] = '#56b3e9' # Sky Blue
COLOR_PALETTE['Y'] = '#f0e442' # Yellow
COLOR_PALETTE['V'] = '#d55e00' # Vermillion
COLOR_PALETTE['RP'] = '#cc79a7' # Reddish purple

# Additional colours
C_BLACK = '#000000'
C_GRAY = '#999999'
C_LIGHT_GRAY = '#f2f2f2'

# Fold Enrichment colours
C_DARK_GREEN = '#5c7943'
C_LIGHT_GREEN = '#abcb42'
C_YELLOW = '#fee71b'
C_ORANGE = '#feae17'
C_RED = '#f35001'

SCORE_COLORS = OrderedDict()

# Gene Score names used in figure legends.
MY_NAME = 'GeVIR'
MY_NO_GERP_NAME = 'GeVIR\n(without GERP++)'
PLI_NAME = 'pLI'
MISS_Z_NAME = 'Missense z-score'
MISS_OE_NAME = 'MOEUF'
LOF_OE_NAME = 'LOEUF'
COMBINED_NAME = 'VIRLoF'

SHET_NAME = 'sHet'
RVIS_NAME = 'RVIS'
DECIPHER_NAME = 'p(HI)'
EPISCORE_NAME = 'EPISCORE'
DOMINO_NAME = 'DOMINO'

CCRS_NAME = 'CCRS'
UNEECON_G_NAME = 'UNEECON-G'

# Gene Score colours used in figures.
SCORE_COLORS[MY_NAME] = COLOR_PALETTE['O']
SCORE_COLORS[MY_NO_GERP_NAME] = COLOR_PALETTE['V']
SCORE_COLORS[MISS_Z_NAME] = COLOR_PALETTE['RP']
SCORE_COLORS[MISS_OE_NAME] = COLOR_PALETTE['SB']
SCORE_COLORS[LOF_OE_NAME] = COLOR_PALETTE['B']
SCORE_COLORS[COMBINED_NAME] = COLOR_PALETTE['G']

SCORE_COLORS[PLI_NAME] = COLOR_PALETTE['BR']
SCORE_COLORS[SHET_NAME] = COLOR_PALETTE['PI']
SCORE_COLORS[RVIS_NAME] = COLOR_PALETTE['BR']
SCORE_COLORS[DECIPHER_NAME] = COLOR_PALETTE['BL']
SCORE_COLORS[EPISCORE_NAME] = COLOR_PALETTE['OL']
SCORE_COLORS[DOMINO_NAME] = COLOR_PALETTE['GL']

SCORE_COLORS[CCRS_NAME] = COLOR_PALETTE['B']

SCORE_COLORS[UNEECON_G_NAME] = COLOR_PALETTE['BR']

# Subplot number title letter. 
#SUBPLOT_LETTERS = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E', 6: 'F', 7: 'G', 8: 'H', 9: 'I'}
SUBPLOT_LETTERS = {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g', 8: 'h', 9: 'i'}

# We assume that in "unaffected" regions with coverage below these thresholds 
# variants absence might be at least partially explained by low coverage. 
MIN_AUTOSOMAL_COVERAGE = 50
MIN_XY_COVERAGE = 45 # 40

# required to analyse fold enrichment of known AD and AR genes by "export_gene_scores" method
# 5% offset = 18,352 / 20 ~= 918
# 5% offset = 19,361 / 20 ~= 968
AD_AR_ENRICHMENT_OFFSET = 968
AD_AR_ENRICHMENT_OFFSET_NO_OUTLIERS = 918

###############
### FIGURES ###
###############

def bar_chart_get_bar_width_and_x_paddings(n_bars):
	bar_width = 0.8 / n_bars

	if n_bars % 2 == 0:
		x_paddings = []
		left_bars_n = n_bars / 2
		for x in range(0, n_bars / 2):
			x_paddings.append(-1 * bar_width * left_bars_n + 0.5 * bar_width)
			left_bars_n -= 1

		right_bars_n = 1
		for x in range(0, n_bars / 2):
			x_paddings.append(bar_width * right_bars_n - 0.5 * bar_width)
			right_bars_n += 1
	else:
		x_paddings = []
		left_bars_n = (n_bars - 1) / 2
		for x in range(0, (n_bars - 1) / 2):
			x_paddings.append(-1 * bar_width * left_bars_n)
			left_bars_n -= 1

		# middle bar
		x_paddings.append(0)
		right_bars_n = 1
		for x in range(0, n_bars / 2):
			x_paddings.append(bar_width * right_bars_n)
			right_bars_n += 1
	return bar_width, x_paddings


def get_precision_recall_f1(ad, ar, all_ad):
	if (ad + ar) == 0:
		return [0,0,0]
	precision = float(ad) / (ad + ar) * 100
	recall = float(ad) / all_ad * 100
	if (precision + recall) > 0:
		f1 = 2 * (precision * recall) / (precision + recall)
	else:
		f1 = 0
	return  [precision, recall, f1]

#######################################################################################
### Figure 1a, Figure 2: Unaffected Regions, Pathogenic Variants (ClinVar) and GERP ###
#######################################################################################

class RegionGroup():
	def __init__(self):
		self.transcripts = []
		self.length_range = ''
		self.region_count = 0
		self.sum_length = 0

		self.miss_region_pathogenic_count = 0
		self.miss_region_pathogenic_percent = 0.0
		self.miss_pathogenic_count = 0
		self.miss_pathogenic_percent = 0.0
		self.miss_per_position = 0.0

		self.lof_region_pathogenic_count = 0
		self.lof_region_pathogenic_percent = 0.0
		self.lof_pathogenic_count = 0
		self.lof_pathogenic_percent = 0.0
		self.lof_per_position = 0.0

	def get_dictionary(self):
		dictionary = OrderedDict()
		if len(self.transcripts) > 0:
			dictionary['transcripts'] = self.transcripts
		dictionary['length_range'] = self.length_range
		dictionary['region_count'] = self.region_count
		dictionary['sum_length'] = self.sum_length

		dictionary['miss_region_pathogenic_count'] = self.miss_region_pathogenic_count
		dictionary['miss_region_pathogenic_percent'] = self.miss_region_pathogenic_percent
		dictionary['miss_pathogenic_count'] = self.miss_pathogenic_count
		dictionary['miss_pathogenic_percent'] = self.miss_pathogenic_percent
		dictionary['miss_per_position'] = self.miss_per_position

		dictionary['lof_region_pathogenic_count'] = self.lof_region_pathogenic_count
		dictionary['lof_region_pathogenic_percent'] = self.lof_region_pathogenic_percent
		dictionary['lof_pathogenic_count'] = self.lof_pathogenic_count
		dictionary['lof_pathogenic_percent'] = self.lof_pathogenic_percent
		dictionary['lof_per_position'] = self.lof_per_position

		return dictionary


def calculate_clin_var_figure_stats(db, region_clin_vars, bootstrap=False):
	overall_pathogenic_miss = 0
	overall_pathogenic_lof = 0

	regions_1_5_length = 0
	regions_6_10_length = 0
	regions_11_15_length = 0
	regions_16_20_length = 0
	regions_21_25_length = 0

	regions_1_5_count = 0
	regions_6_10_count = 0
	regions_11_15_count = 0
	regions_16_20_count = 0
	regions_21_25_count = 0

	regions_1_5_miss = 0
	regions_6_10_miss = 0
	regions_11_15_miss = 0
	regions_16_20_miss = 0
	regions_21_25_miss = 0

	regions_1_5_pathogenic_miss_count = 0
	regions_6_10_pathogenic_miss_count = 0
	regions_11_15_pathogenic_miss_count = 0
	regions_16_20_pathogenic_miss_count = 0
	regions_21_25_pathogenic_miss_count = 0

	regions_1_5_lof = 0
	regions_6_10_lof = 0
	regions_11_15_lof = 0
	regions_16_20_lof = 0
	regions_21_25_lof = 0

	regions_1_5_pathogenic_lof_count = 0
	regions_6_10_pathogenic_lof_count = 0
	regions_11_15_pathogenic_lof_count = 0
	regions_16_20_pathogenic_lof_count = 0
	regions_21_25_pathogenic_lof_count = 0

	genes_with_pathogenic_vars = 0
	for region_clin_var in region_clin_vars:
		# Skip genes with no pathogenic variants in VIRs
		if len(region_clin_var['regions']) == 0 and len(region_clin_var['regions_lof']) == 0:
			continue
		genes_with_pathogenic_vars += 1

		transcript_id = region_clin_var['_id']
		for region_length in region_clin_var['all_regions']:
			if region_length < 6:
				regions_1_5_length += region_length
				regions_1_5_count += 1
			elif region_length < 11:
				regions_6_10_length += region_length
				regions_6_10_count += 1
			elif region_length < 16:
				regions_11_15_length += region_length
				regions_11_15_count += 1
			elif region_length < 21:
				regions_16_20_length += region_length
				regions_16_20_count += 1
			else:
				regions_21_25_length += region_length
				regions_21_25_count += 1

		for region_id, variants_count in region_clin_var['regions'].iteritems():
			start, end = region_id.split('-')
			start = int(start)
			end = int(end)
			region_length = end - start - 1 # subtract 1 so that regions like 15-16 had zero length
			overall_pathogenic_miss += variants_count
			if region_length < 6:
				regions_1_5_miss += variants_count
				regions_1_5_pathogenic_miss_count += 1
			elif region_length < 11:
				regions_6_10_miss += variants_count
				regions_6_10_pathogenic_miss_count += 1
			elif region_length < 16:
				regions_11_15_miss += variants_count
				regions_11_15_pathogenic_miss_count += 1
			elif region_length < 21:
				regions_16_20_miss += variants_count
				regions_16_20_pathogenic_miss_count += 1
			else:
				regions_21_25_miss += variants_count
				regions_21_25_pathogenic_miss_count += 1

		for region_id, variants_count in region_clin_var['regions_lof'].iteritems():
			start, end = region_id.split('-')
			start = int(start)
			end = int(end)
			region_length = end - start - 1 # subtract 1 so that regions like 15-16 had zero length
			overall_pathogenic_lof += variants_count
			if region_length < 6:
				regions_1_5_lof += variants_count
				regions_1_5_pathogenic_lof_count += 1
			elif region_length < 11:
				regions_6_10_lof += variants_count
				regions_6_10_pathogenic_lof_count += 1
			elif region_length < 16:
				regions_11_15_lof += variants_count
				regions_11_15_pathogenic_lof_count += 1
			elif region_length < 21:
				regions_16_20_lof += variants_count
				regions_16_20_pathogenic_lof_count += 1
			else:
				regions_21_25_lof += variants_count
				regions_21_25_pathogenic_lof_count += 1

	groups_region_length = [regions_1_5_length, regions_6_10_length, regions_11_15_length, regions_16_20_length, regions_21_25_length]
	groups_region_miss = [regions_1_5_miss, regions_6_10_miss, regions_11_15_miss, regions_16_20_miss, regions_21_25_miss]
	groups_region_count = [regions_1_5_count, regions_6_10_count, regions_11_15_count, regions_16_20_count, regions_21_25_count]
	groups_region_pathogenic_miss_count = [regions_1_5_pathogenic_miss_count, regions_6_10_pathogenic_miss_count, regions_11_15_pathogenic_miss_count, regions_16_20_pathogenic_miss_count, regions_21_25_pathogenic_miss_count]
	groups_region_lof = [regions_1_5_lof, regions_6_10_lof, regions_11_15_lof, regions_16_20_lof, regions_21_25_lof,]
	groups_region_pathogenic_lof_count = [regions_1_5_pathogenic_lof_count, regions_6_10_pathogenic_lof_count, regions_11_15_pathogenic_lof_count, regions_16_20_pathogenic_lof_count, regions_21_25_pathogenic_lof_count,]

	if not bootstrap:
		print 'Genes with ClinVar variants:', genes_with_pathogenic_vars
		print 'Pathogenic Miss in VIRs:', sum(groups_region_miss)
		print 'Pathogenic LoF in VIRs:', sum(groups_region_lof)
		print 'All Pathogenic in VIRs:', sum(groups_region_miss) + sum(groups_region_lof)
		print 'Short regions proportion %.2f%%' % (float(regions_1_5_count) * 100 / sum(groups_region_count))
		print 'Miss in short regions proportion  %.2f%%' % (float(regions_1_5_miss) * 100 / sum(groups_region_miss))

	length_ranges = []
	miss_percents = []
	miss_per_positions = []
	miss_pathogenic_regions_percents = []

	lof_percents = []
	lof_per_positions = []
	lof_pathogenic_regions_percents = []

	ex_region_group = RegionGroup()
	ex_region_group = ex_region_group.get_dictionary()
	table = [ex_region_group.keys()]

	bootstrap_results = OrderedDict()
	for x in range(0, 5):
		region_group = RegionGroup()
		if x < 4:
			length_range = str(x*5+1) + '-' + str(x*5+5)
		else:
			length_range = '21+'
		length_ranges.append(length_range)

		region_group.length_range = length_range
		region_group.region_count = groups_region_count[x]
		region_group.sum_length = groups_region_length[x]

		region_group.miss_region_pathogenic_count = groups_region_pathogenic_miss_count[x]
		miss_pathogenic_regions_percents.append(float(groups_region_pathogenic_miss_count[x]) * 100 / groups_region_count[x])
		region_group.miss_region_pathogenic_percent = miss_pathogenic_regions_percents[x]
		region_group.miss_pathogenic_count = groups_region_miss[x]
		miss_percents.append(float(groups_region_miss[x]) * 100 / overall_pathogenic_miss)
		region_group.miss_pathogenic_percent = miss_percents[x]
		miss_per_positions.append(float(groups_region_miss[x]) / groups_region_length[x])
		region_group.miss_per_position = miss_per_positions[x]

		lof_pathogenic_regions_percents.append(float(groups_region_pathogenic_lof_count[x]) * 100 / groups_region_count[x])
		region_group.lof_region_pathogenic_count = lof_pathogenic_regions_percents[x]
		region_group.lof_region_pathogenic_percent = lof_pathogenic_regions_percents[x]
		region_group.lof_pathogenic_count = groups_region_lof[x]
		lof_percents.append(float(groups_region_lof[x]) * 100 / overall_pathogenic_lof)
		region_group.lof_pathogenic_percent = lof_percents[x]
		lof_per_positions.append(float(groups_region_lof[x]) / groups_region_length[x])
		region_group.lof_per_position = lof_per_positions[x]

		region_group = region_group.get_dictionary()
		table.append(region_group.values())
		bootstrap_results[length_range] = region_group

	if bootstrap:
		clin_var_transcripts = []
		for region_clin_var in region_clin_vars:
			clin_var_transcripts.append(region_clin_var['_id'])
		bootstrap_results['transcripts'] = clin_var_transcripts
		db.gevir.temp_regions_clin_var_bootsrap.insert(bootstrap_results)
	else:
		# Report table stats
		output_csv = OUTPUT_FOLDER + 'clin_var_regions.csv'
		write_table_to_csv(table, output_csv)	

		# Report 1-5 vs 21+ region bin comparison stats
		miss_short = [groups_region_miss[0], groups_region_length[0]]
		miss_long = [groups_region_miss[4], groups_region_length[4]]

		lof_short = [groups_region_lof[0], groups_region_length[0]]
		lof_long = [groups_region_lof[4], groups_region_length[4]]

		miss_fe, miss_p_value = fisher_exact([miss_long, miss_short])
		lof_fe, lof_p_value = fisher_exact([lof_long, lof_short])
		
		# p-valie = 0.0 means p-value < 2.23E-308 (to check try "print sys.float_info.min")

		headers  = ['Group', 'VIRs length group 1', 'Variants 1', 'Total length 1', 'VIRs length group 2', 'Variants 2', 'Total length 2', 'Fold-enrichment', 'p-value']
		miss_row = ['Missense', '1-5', groups_region_miss[0], groups_region_length[0], '21+', groups_region_miss[4], groups_region_length[4], miss_fe, miss_p_value]
		lof_row = ['Loss-of-Function', '1-5', groups_region_lof[0], groups_region_length[0], '21+', groups_region_lof[4], groups_region_length[4], lof_fe, lof_p_value]
		table = [headers, miss_row, lof_row]
		output_csv = OUTPUT_FOLDER + 'clin_var_short_and_long_virs_comparison.csv'
		write_table_to_csv(table, output_csv)

	return length_ranges, groups_region_count, groups_region_miss, miss_percents, miss_per_positions, groups_region_lof, lof_percents, lof_per_positions


def calculate_clin_var_bootstrap_stats_confidence_interval(db, stat_name):
	bin_names = ['1-5', '6-10', '11-15', '16-20', '21+']
	bin_stats = OrderedDict()
	for bin_name in bin_names:
		bin_stats[bin_name] = []

	bootstrap_clin_var_runs = db.gevir.temp_regions_clin_var_bootsrap.find({})
	for bootstrap_clin_var_run in bootstrap_clin_var_runs:
		for bin_name in bin_names:
			bin_stats[bin_name].append(bootstrap_clin_var_run[bin_name][stat_name])

	lower_cis = []
	upper_cis = []
	for bin_name, stats in bin_stats.iteritems():
		lower_ci, upper_ci = calculate_confidence_interval(stats)
		lower_cis.append(lower_ci)
		upper_cis.append(upper_ci)
	return [lower_cis, upper_cis], bin_stats


def analyse_clin_var(db, clean_temp_data=False, include_example_genes_subplot=False):
	if clean_temp_data:
		db.gevir.temp_regions_clin_var.drop()

	temp_clin_var_table = db.gevir.temp_regions_clin_var.find_one({})

	gnomad_transcript_ids = []

	if INCLUDE_GNOMAD_OUTLIERS:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True }) # "no_issues": True, 
	else:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	for gnomad_gene in gnomad_genes:
		gnomad_transcript_ids.append(gnomad_gene['_id'])

	clin_variants_num = 0

	if not temp_clin_var_table:
		transcript_variants = {}
		transcript_variants_lof = {}
		clin_vars = db.gevir.clin_vars.find({})

		miss_num = 0
		stop_gained_num = 0
		frameshift_num = 0

		total_lines = clin_vars.count()
		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for clin_var in clin_vars:
			# Skip variants for which no assertion criteria was provided (i.e. should have at least 1 ClinVar star)
			# Not used in the original analysis
			
			clin_sigs = clin_var['clin_sig'].split(', ')
			if 'Likely pathogenic' in clin_sigs or 'Pathogenic/Likely pathogenic' in clin_sigs or 'Pathogenic' in clin_sigs:
				if 'canonical_transcripts' in clin_var:
					for transcript_id in clin_var['canonical_transcripts']:
						if transcript_id not in gnomad_transcript_ids:
							continue

						if transcript_id in clin_var['csqs']:
							xpos = get_xpos(clin_var['chrom'], clin_var['start'])
							if clin_var['csqs'][transcript_id] == 'missense_variant': # in set(['stop_gained', 'frameshift_variant']):
								if transcript_id not in transcript_variants:
									transcript_variants[transcript_id] = []
								transcript_variants[transcript_id].append(xpos)
								miss_num += 1
							elif clin_var['csqs'][transcript_id] in set(['stop_gained', 'frameshift_variant']):
								if transcript_id not in transcript_variants_lof:
									transcript_variants_lof[transcript_id] = []
								transcript_variants_lof[transcript_id].append(xpos)
								if clin_var['csqs'][transcript_id] == 'stop_gained':
									stop_gained_num += 1
								elif clin_var['csqs'][transcript_id] == 'frameshift_variant':
									frameshift_num += 1

			line_number += 1
			bar.update((line_number + 0.0) / total_lines)
		bar.finish()

		db.gevir.temp_regions_clin_var.drop()
		print 'Missense Transcripts', len(transcript_variants)
		print 'LoF Transcripts', len(transcript_variants_lof)
		print 'Missense or LoF Transcripts', len(set(transcript_variants.keys() + transcript_variants_lof.keys()))

		print 'Missense', miss_num
		print 'Stop Gained', stop_gained_num
		print 'Frameshift', frameshift_num

		print 'All pathogenic', miss_num + stop_gained_num + frameshift_num

		# Pathogenic Missense
		transcript_var_regions = {}
		transcript_all_regions = {}
		transcript_all_regions_gerp = {}

		total_lines = len(transcript_variants)
		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for transcript_id, xposes in transcript_variants.iteritems():
			transcript_info = db.exac.transcripts.find_one({'transcript_id': transcript_id})

			coverage_threshold = MIN_AUTOSOMAL_COVERAGE
			if transcript_info['chrom'] == 'X' or transcript_info['chrom'] == 'Y':
				coverage_threshold = MIN_XY_COVERAGE

			transcript_var_regions[transcript_id] = {}

			for xpos in xposes:
				region = db.gevir[REGIONS_COLLECTION].find_one({ "transcript_id": transcript_id, "exome_coverage": { "$gte": coverage_threshold }, "xstart": { "$lt": xpos }, "xstop": { "$gt": xpos }, "lenght": { "$gte": 1 } })
				if region:
					region_id = str(region['start_variant']['protein_pos']) + '-' + str(region['stop_variant']['protein_pos'])
					if region_id not in transcript_var_regions[transcript_id]:
						transcript_var_regions[transcript_id][region_id] = 1
					else:
						transcript_var_regions[transcript_id][region_id] += 1

			transcript_all_regions[transcript_id] = []
			transcript_all_regions_gerp[transcript_id] = []
			regions = db.gevir[REGIONS_COLLECTION].find({ "transcript_id": transcript_id, "exome_coverage": { "$gte": coverage_threshold }, "lenght": { "$gte": 1 }})
			for region in regions:
				transcript_all_regions[transcript_id].append(region['lenght'])
				transcript_all_regions_gerp[transcript_id].append(region['gerp_mean'])

			db.gevir.temp_regions_clin_var.insert({'_id': transcript_id, 'xposes': xposes, 'regions': transcript_var_regions[transcript_id], 'regions_num': len(transcript_var_regions[transcript_id]),
												'all_regions': transcript_all_regions[transcript_id], 'all_regions_gerp': transcript_all_regions_gerp[transcript_id],
												'xposes_lof': [], 'regions_lof': {}, 'regions_lof_num': 0})

			line_number += 1
			bar.update((line_number + 0.0) / total_lines)
		bar.finish()

		# LoF
		transcript_var_regions_lof = {}
		total_lines = len(transcript_variants_lof)
		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for transcript_id, xposes in transcript_variants_lof.iteritems():
			transcript_info = db.exac.transcripts.find_one({'transcript_id': transcript_id})

			coverage_threshold = MIN_AUTOSOMAL_COVERAGE
			if transcript_info['chrom'] == 'X' or transcript_info['chrom'] == 'Y':
				coverage_threshold = MIN_XY_COVERAGE

			transcript_var_regions_lof[transcript_id] = {}

			for xpos in xposes:
				region = db.gevir[REGIONS_COLLECTION].find_one({ "transcript_id": transcript_id, "exome_coverage": { "$gte": coverage_threshold }, "xstart": { "$lt": xpos }, "xstop": { "$gt": xpos }, "lenght": { "$gte": 1 }, "not_in_cds": False })
				if region:
					region_id = str(region['start_variant']['protein_pos']) + '-' + str(region['stop_variant']['protein_pos'])
					if region_id not in transcript_var_regions_lof[transcript_id]:
						transcript_var_regions_lof[transcript_id][region_id] = 1
					else:
						transcript_var_regions_lof[transcript_id][region_id] += 1

			# If transcript also contain pathogenic missense variants, then update existing record,
			# Else calculate GERP data and add new record
			if transcript_id in transcript_variants:
				db.gevir.temp_regions_clin_var.update_one({'_id': transcript_id}, {'$set': {'xposes_lof': xposes, 'regions_lof': transcript_var_regions_lof[transcript_id]}})
			else:
				transcript_all_regions[transcript_id] = []
				transcript_all_regions_gerp[transcript_id] = []
				regions = db.gevir[REGIONS_COLLECTION].find({ "transcript_id": transcript_id, "exome_coverage": { "$gte": coverage_threshold }, "lenght": { "$gte": 1 }, "not_in_cds": False})
				for region in regions:
					transcript_all_regions[transcript_id].append(region['lenght'])
					transcript_all_regions_gerp[transcript_id].append(region['gerp_mean'])

				db.gevir.temp_regions_clin_var.insert({'_id': transcript_id, 'xposes': [], 'regions': {}, 
													'all_regions': transcript_all_regions[transcript_id], 'all_regions_gerp': transcript_all_regions_gerp[transcript_id],
													'xposes_lof': xposes, 'regions_lof': transcript_var_regions_lof[transcript_id], 'regions_lof_num': len(transcript_var_regions_lof[transcript_id]), 'regions_num': 0})

			line_number += 1
			bar.update((line_number + 0.0) / total_lines)
		bar.finish()

		# Remove transcripts with no pathogenic variants in VIR
		region_clin_vars = db.gevir.temp_regions_clin_var.find({})
		for region_clin_var in region_clin_vars:
			if len(region_clin_var['regions']) == 0 and len(region_clin_var['regions_lof']) == 0:
				db.gevir.temp_regions_clin_var.remove({'_id': region_clin_var['_id']})


	if clean_temp_data:
		db.gevir.temp_regions_clin_var_bootsrap.drop()

	# Run analysis 10,000 times on 50% of the data
	bootstrap_clin_var_table = db.gevir.temp_regions_clin_var_bootsrap.find_one({})
	if not bootstrap_clin_var_table:
		region_clin_vars = db.gevir.temp_regions_clin_var.find({})
		region_clin_vars_list = []
		for region_clin_var in region_clin_vars:
			region_clin_vars_list.append(region_clin_var)

		n_size = int(len(region_clin_vars_list) * 0.50)
		runs = 10000
		run = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for x in range(runs):
			region_clin_vars_subset = resample(region_clin_vars_list, n_samples=n_size)
			calculate_clin_var_figure_stats(db, region_clin_vars_subset, bootstrap=True)
			run += 1
			bar.update((run + 0.0) / runs)
		bar.finish()

	region_clin_vars = db.gevir.temp_regions_clin_var.find({})
	length_ranges, groups_region_count, groups_region_miss, miss_percents, miss_per_positions, groups_region_lof, lof_percents, lof_per_positions = calculate_clin_var_figure_stats(db, region_clin_vars)

	# Calculate 95% Confidence Intervals from bootsrap data
	miss_per_position_cis, miss_per_position_raw = calculate_clin_var_bootstrap_stats_confidence_interval(db, 'miss_per_position')
	lof_per_position_cis, lof_per_position_raw = calculate_clin_var_bootstrap_stats_confidence_interval(db, 'lof_per_position')

	draw_clin_var_figure(db, length_ranges, groups_region_count, groups_region_miss, miss_percents, miss_per_positions, miss_per_position_cis, miss_per_position_raw,
	                     groups_region_lof, lof_percents, lof_per_positions, lof_per_position_cis, lof_per_position_raw, include_example_genes_subplot=include_example_genes_subplot)


def get_x_padding_region_count(num):
	if float(num) < 1:
		return 0.15
	n = len(str(int(num)))

	if n < 3:	
		return 0.15
	elif n < 4:
		return 0.12
	elif n < 5:
		return 0.15
	elif n < 6:
		return 0.17
	else:
		return 0.2


def get_unique_rounded_numbers(numbers, ndigits):
	unique_numbers = set()
	for number in numbers:
		unique_numbers.add(round(number, 4))
	return list(unique_numbers)


def draw_region_count_subplot(subplot_num, x_labels, ys, annotations, cis=[], raw=[], ylabel='', title='', multi_bar_labels=[], legend_loc='upper right', subplot_letter='', row_num=3):
	ax = plt.subplot(row_num,2,subplot_num)

	if multi_bar_labels:
		xs = np.arange(0, len(x_labels))
		n_bars = len(ys)
		bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(n_bars)
		for x in range(0, n_bars):
			max_ys = max([z for y in ys for z in y])
			y_padding = max_ys / 30.0 # get maximum y out of all bar lists
			bar_ys = ys[x]
			bar_annotations = annotations[x]
			if cis:
				print 'Confidence Intervals'
				print cis[x]
				lower_errs = []
				upper_errs = []
				for i in range(0, len(ys[x])):
					y = ys[x][i]
					lower_errs.append(ys[x][i] - cis[x][0][i])
					upper_errs.append(cis[x][1][i] - ys[x][i])

				plt.bar(xs + x_paddings[x], ys[x], width=bar_width, label=multi_bar_labels[x], yerr=[lower_errs, upper_errs], capsize=3, error_kw={'elinewidth': 1}) # cis[x]
			else:
				plt.bar(xs + x_paddings[x], ys[x], width=bar_width, label=multi_bar_labels[x])


			plt.ylim(top=max_ys + y_padding * 10)
			for z in range(0, len(bar_annotations)):
				if cis:
					if x == 0:
						color = 'white'
					else:
						color = 'black'

					text_size = 6
					if len(raw) > 0:
						text_size = 5
					plt.text(x=xs[z] + x_paddings[x], y=0.0065, s=bar_annotations[z], size=text_size, color=color, ha='center',  rotation=45)
				else:
					num = ys[x][z]
					print num
					if num > 60:
						y=ys[x][z] + y_padding + 11
					elif num > 10:
						y=ys[x][z] + y_padding + 9
					else:
						y=ys[x][z] + y_padding + 5

					plt.text(x=xs[z] + x_paddings[x], y=y, s="{:,}".format(bar_annotations[z]), size=6, ha='center',  rotation=45)

		
	else:
		xs = np.arange(0, len(x_labels))

		# convert ys to %
		ys_percent = []
		all_ys = float(sum(ys))
		for y in ys:
			y_percent = y * 100 / all_ys
			ys_percent.append(y_percent)

		plt.bar(xs, ys_percent)
		y_padding = 3
		plt.ylim(0, 100)
		for x in range(0, len(annotations)):
			plt.text(x=xs[x], y=ys_percent[x]+y_padding, s="{:,}".format(annotations[x]), size=6, ha='center') # fontweight='bold'
	
		ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

	plt.xticks(xs, x_labels, fontsize=7)
	plt.yticks(fontsize=7)
	if subplot_letter == 'a' or subplot_letter == 'b':
		plt.yticks([0, 25, 50, 75, 100])


	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	if subplot_letter:
		subplot_letter_x_padding = 0
		if row_num == 2:
			if subplot_letter == 'a':
				subplot_letter_x_padding = -0.18 #-0.20
			elif subplot_letter == 'b':	
				subplot_letter_x_padding = -0.14
			elif subplot_letter == 'c':	
				subplot_letter_x_padding = -0.19 #-0.125 #; To be consistent with 'b' use -0.20
		elif row_num == 3:
			if subplot_letter == 'b':
				subplot_letter_x_padding = -0.16 #-0.20
			elif subplot_letter == 'c':	
				subplot_letter_x_padding = -0.13
			elif subplot_letter == 'd':	
				subplot_letter_x_padding = -0.16 #-0.125 #; To be consistent with 'b' use -0.20
		ax.text(subplot_letter_x_padding, 1.11, subplot_letter, transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')
	else:
		ax.text(-0.08, 1.13, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')

	plt.xlabel('VIR length bins', fontsize=7)
	if ylabel:
		plt.ylabel(ylabel, fontsize=7)

	if len(raw) > 0:
		s_x = 0
		bin_names = ['1-5', '6-10', '11-15', '16-20', '21+']
		all_s_ys = []
		for bin_name in bin_names:
			# miss
			s_ys = raw[0][bin_name]
			# To reduce figure size, we do not plot dots closely overlapping each other
			s_ys = get_unique_rounded_numbers(s_ys, 4)

			s_xs = [s_x] * len(s_ys)
			plt.scatter(s_xs, s_ys, zorder=4, s=1, c=COLOR_PALETTE['G'])
			all_s_ys += s_ys

			# LoF
			s_ys = raw[1][bin_name]
			# To reduce figure size, we do not plot dots closely overlapping each other
			s_ys = get_unique_rounded_numbers(s_ys, 4)

			s_xs = [s_x + 0.4] * len(s_ys)
			plt.scatter(s_xs, s_ys, zorder=4, s=1, c=COLOR_PALETTE['G'])
			all_s_ys += s_ys
			
			s_x += 1

		plt.ylim(top=max(all_s_ys) + y_padding)

	if multi_bar_labels:
		if len(raw) > 0:
			handles, labels = ax.get_legend_handles_labels()
			bootstrap_patch = Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_PALETTE['G'], markersize=4, label='Results from bootstrapping') 
			handles.append(bootstrap_patch)
			plt.legend(loc=legend_loc, frameon=False, fontsize=7, handles=handles)
		else:
			plt.legend(loc=legend_loc, frameon=False, fontsize=7)

	# Manual tuning
	if subplot_letter == 'd':
		plt.yticks(np.arange(0, 0.055, 0.01), [0, 0.01, 0.02, 0.03, 0.04, 0.05])


def get_gnomad_variants(db, transcript_id, exome=True):
	if exome:
		variants = db.exac.exome_variants.find({'transcripts': transcript_id})
	else:
		variants = db.exac.genome_variants.find({'transcripts': transcript_id})

	variant_positions = set([])
	for variant in variants:
		if variant['filter'] not in VALID_FILTERS:
			continue

		veps = variant['vep_annotations']
		for vep in veps:
			if vep['Feature'] == transcript_id:
				csq = worst_csq_from_csq(vep['Consequence'])
				if csq in VALID_CSQS:
					protein_pos = vep['Protein_position']
					if '-' in protein_pos:
						protein_pos = protein_pos.split('-')[0]
					if is_int(protein_pos):
						variant_positions.add(int(protein_pos))

	return variant_positions


def get_clin_var_variants(db, transcript_id, phenotype_id=''):
	lof_csqs = set(['stop_gained', 'frameshift_variant'])
	clin_vars = db.gevir.clin_vars.find({"canonical_transcripts": transcript_id })
	miss_clin_var_poses = set([])
	lof_clin_var_poses = set([])

	for clin_var in clin_vars:
		clin_sigs = clin_var['clin_sig'].split(', ')
		if 'transcript_consequences' not in clin_var:
			continue

		if phenotype_id and ('OMIM' not in clin_var['phenotype_ids'] or phenotype_id not in clin_var['phenotype_ids']['OMIM']):
			continue
		
		if 'Likely pathogenic' in clin_sigs or 'Pathogenic/Likely pathogenic' in clin_sigs or 'Pathogenic' in clin_sigs:
			for transcript_consequence in clin_var['transcript_consequences']:
				if transcript_consequence['transcript_id'] == transcript_id:
					csq = clin_var['csqs'][transcript_id]

					if csq != 'missense_variant' and csq not in lof_csqs:
						continue

					pos = transcript_consequence['protein_start']
					
					if csq == 'missense_variant':
						miss_clin_var_poses.add(pos)

					if csq in lof_csqs:
						lof_clin_var_poses.add(pos)


	miss_clin_var_poses = list(miss_clin_var_poses)
	miss_clin_var_poses.sort()

	lof_clin_var_poses = list(lof_clin_var_poses)
	lof_clin_var_poses.sort()

	return miss_clin_var_poses, lof_clin_var_poses


def draw_example_gene_variants_subplot(db, subplot_num, subplot_letter='', single_figure=False, row_num=3):
	transcript_data = OrderedDict()
	transcript_data['ENST00000240185'] = ('TARDBP', '612069', 'Amyotrophic lateral sclerosis, with or without Frontotemporal lobar degeneration (AD)')
	transcript_data['ENST00000216124'] = ('ARSA', '250100', 'Metachromatic leukodystrophy (AR)')
	transcript_data['ENST00000571688'] = ('LITAF', '601098', 'Charcot-Marie-Tooth disease, type 1C (AD)')
	transcript_data['ENST00000398339'] = ('TCF4', '610954', 'Pitt-Hopkins syndrome (AD)')

	if single_figure:
		ax = plt.subplot(111)
	else:
		ax = plt.subplot2grid((row_num, 2), (0, 0), colspan=2)
		if subplot_letter:
			ax.text(-0.08, 1.05, subplot_letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='right')
		else:
			ax.text(-0.08, 1.05, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='right')

	y = 0.5
	gene_names = []
	create_legend = True
	for transcript_id, info  in transcript_data.iteritems():
		gene_name, phenotype_id, phenotype_name = info
		gene_names.append(gene_name)
		gene_scores = db.gevir.common_gene_scores.find_one({'_id': transcript_id})

		omim = db.gevir.omim.find_one({'transcript_id': transcript_id})

		phenotypes_num = 0

		print transcript_id, gene_name, phenotypes_num, 'LoF oe:', gene_scores['gnomad_oe_lof_upper'], 'Miss oe:', gene_scores['gnomad_oe_mis_upper'], 'Miss Z:', gene_scores['gnomad_miss_z']
		if omim:
			phenotypes_num = len(omim['phenotypes'])
			for p in omim['phenotypes']:
				print p
		
		ens_transcript = db.gevir.ens_aa_fasta.find_one({'_id': transcript_id})
		cds_aa_len = len(ens_transcript['cds'])

		exome_missense = get_gnomad_variants(db, transcript_id, exome=True)
		genome_missense = get_gnomad_variants(db, transcript_id, exome=False)
		control_missense = exome_missense | genome_missense

		control_missense = list(control_missense)
		control_missense.sort()

		miss_clin_var_poses, lof_clin_var_poses = get_clin_var_variants(db, transcript_id, phenotype_id=phenotype_id)

		control_ys = [y-0.1] * len(control_missense)
		miss_pathogenic_ys = [y] * len(miss_clin_var_poses)
		lof_pathogenic_ys = [y + 0.1] * len(lof_clin_var_poses)

		p = mpatches.Rectangle((0, y-0.25), cds_aa_len, 0.5, fill=True, color=C_LIGHT_GRAY, alpha=0.5, zorder=0)
		ax.add_patch(p)
		p = mpatches.Rectangle((0, y-0.25), cds_aa_len, 0.5, fill=False, color=C_BLACK, alpha=0.5, zorder=1)
		ax.add_patch(p)

		plt.text(x=0, y=y+0.35, s=phenotype_name, size=7)

		if create_legend:
			plt.scatter(control_missense, control_ys, color=COLOR_PALETTE['B'], s=2, label='protein altering (gnomAD)')
			plt.scatter(miss_clin_var_poses, miss_pathogenic_ys, color='black', s=2, label='pathogenic missense (ClinVar)')
			plt.scatter(lof_clin_var_poses, lof_pathogenic_ys, color=COLOR_PALETTE['V'], s=2, label='pathogenic stop gained, frameshift (ClinVar)')
			create_legend = False
		else:
			plt.scatter(control_missense, control_ys, color=COLOR_PALETTE['B'], s=2)
			plt.scatter(miss_clin_var_poses, miss_pathogenic_ys, color='black', s=2)
			plt.scatter(lof_clin_var_poses, lof_pathogenic_ys, color=COLOR_PALETTE['V'], s=2)			

		y += 1

	ys = np.arange(-0.5, len(gene_names) + 1, 1)
	plt.xticks(fontsize=7)
	plt.yticks(ys, [''] + gene_names, fontsize=7)
	plt.ylim(0,y)

	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	ax.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		left=False)        # labels along the bottom edge are off

	legend_elements = [Line2D([0], [0], marker='o', color='w', label='protein altering (gnomAD)', markerfacecolor=COLOR_PALETTE['B'], markersize=7),
					   Line2D([0], [0], marker='o', color='w', label='pathogenic missense (ClinVar)', markerfacecolor='black', markersize=7),
					   Line2D([0], [0], marker='o', color='w', label='pathogenic stop gained,\nframeshift (ClinVar)', markerfacecolor=COLOR_PALETTE['V'], markersize=7),]

	plt.legend(handles=legend_elements, loc='center right', frameon=False, fontsize=7)
	plt.xlabel('Protein length (amino acids)', fontsize=7)

	if single_figure:
		fig = plt.figure(1)
		fig.set_size_inches(6.4, 2.3)
		plt.tight_layout(rect=[0.01, -0.05, 1, 1.1])
		fig_format = 'pdf'
		plt.savefig(FIGURES_FOLDER + 'ClinVar_vs_regions_example_genes.' + fig_format, format=fig_format, dpi=300, transparent=True)


def draw_clin_var_figure(db, length_ranges, groups_region_count, groups_region_miss, miss_percents, miss_per_positions, miss_per_position_cis, miss_per_position_raw,
	                     groups_region_lof, lof_percents, lof_per_positions, lof_per_position_cis, lof_per_position_raw, include_example_genes_subplot=False):

	miss_per_positions_str = []
	for miss_per_position in miss_per_positions:
		miss_per_positions_str.append("%.3f" % round(miss_per_position,3))

	lof_per_positions_str = []
	for lof_per_position in lof_per_positions:
		lof_per_positions_str.append("%.3f" % round(lof_per_position,3))

	multi_bar_labels = ['Missense', 'Stop gained, frameshift']
	subplot_num = 1
	if include_example_genes_subplot:
		row_num = 3
		draw_example_gene_variants_subplot(db, subplot_num, subplot_letter='a', row_num=row_num)
		subplot_num += 2	
		subplot_letters = ['b', 'c', 'd', 'e']
	else:
		row_num = 2
		subplot_letters = ['a', 'b', 'c', 'd']

	draw_region_count_subplot(subplot_num, length_ranges, groups_region_count, groups_region_count, ylabel='Variant Intolerant Regions\n(VIRs, %)', subplot_letter=subplot_letters[0], row_num=row_num)
	subplot_num += 1
	draw_region_count_subplot(subplot_num, length_ranges, [miss_percents, lof_percents], [groups_region_miss, groups_region_lof], ylabel='Pathogenic variants (%)', multi_bar_labels=multi_bar_labels, subplot_letter=subplot_letters[1], row_num=row_num)
	subplot_num += 1
	draw_region_count_subplot(subplot_num, length_ranges, [miss_per_positions, lof_per_positions], [miss_per_positions_str, lof_per_positions_str],
							  cis=[miss_per_position_cis, lof_per_position_cis], raw=[miss_per_position_raw, lof_per_position_raw], 
							  ylabel='Pathogenic variants\nper amino acid', multi_bar_labels=multi_bar_labels, legend_loc='upper left', subplot_letter=subplot_letters[2], row_num=row_num)
	subplot_num += 1

	draw_regions_gerp_whisker_figure(db, subplot_num, subplot_letter=subplot_letters[3], row_num=row_num)

	plt.figure(figsize = (2,2))
	gs1 = gridspec.GridSpec(2, 2)
	gs1.update(wspace=0.05, hspace=0.05)
	fig = plt.figure(1)
	if include_example_genes_subplot:
		fig.set_size_inches(7.5, 6.5)
	else:
		fig.set_size_inches(7, 3.7)

	plt.tight_layout(rect=[0, 0, 1, 1])
	plt.savefig(FIGURES_FOLDER + 'ClinVar_vs_regions.pdf', format='pdf', dpi=300)


def get_transcript_id_to_chrom_dict(db):
	transcript_id_to_chrom = {}
	transcripts = db.exac.transcripts.find({})
	for transcript in transcripts:
		transcript_id_to_chrom[transcript['transcript_id']] = transcript['chrom']
	return transcript_id_to_chrom


def remove_overlapping_boxplot_outliers(boxplot):
	for l in boxplot['fliers']:
		unique_ys = set()
		f_xs = []
		f_ys = []
		xs, ys = l.get_data()
		for i in range(0, len(ys)):
			y = ys[i]
			y = round(y, 2)
			if y not in unique_ys:
				f_xs.append(xs[i])
				f_ys.append(ys[i])
				unique_ys.add(y)
		l.set_data(f_xs, f_ys)


def draw_regions_gerp_whisker_figure(db, subplot_num, subplot_letter='', row_num=3):
	ax = plt.subplot(row_num,2,subplot_num)

	transcript_id_to_chrom = get_transcript_id_to_chrom_dict(db)
	gnomad_transcript_ids = set([])

	if INCLUDE_GNOMAD_OUTLIERS:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True })
	else:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	db.gevir.temp.remove({'_id': 'region_length_gerp_raw_bins'})
	unique_genes = set()
	unique_regions = 0

	bin_names = ['1-5', '6-10', '11-15', '16-20', '21+']
	bin_stats = OrderedDict()
	for bin_name in bin_names:
		bin_stats[bin_name] = []

	regions = db.gevir[REGIONS_COLLECTION].find({"lenght": { "$gte": 1 }, "not_in_cds": False })

	print 'Collecting regions GERP++ data'
	total_lines = regions.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for region in regions:


		chrom = transcript_id_to_chrom[region['transcript_id']]
		coverage_threshold = MIN_AUTOSOMAL_COVERAGE
		if chrom == 'X' or chrom == 'Y':
			coverage_threshold = MIN_XY_COVERAGE

		# count only high covered regions
		if region['exome_coverage'] >= coverage_threshold:
			unique_regions += 1
			unique_genes.add(region['transcript_id'])
			length = region['lenght']
			gerp = region['gerp_mean']

			if length >= 1 and length <= 5:
				bin_stats['1-5'].append(gerp)
			elif length >= 6 and length <= 10:
				bin_stats['6-10'].append(gerp)
			elif length >= 11 and length <= 15:
				bin_stats['11-15'].append(gerp)
			elif length >= 16 and length <= 20:
				bin_stats['16-20'].append(gerp)
			elif length >= 21:
				bin_stats['21+'].append(gerp)
			else:
				print 'unexpected length', length

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)			
	bar.finish()

	print 'unique genes', len(unique_genes)
	print 'unique regions', unique_regions

	# Draw boxplots
	flierprops = dict(marker='o', markerfacecolor=COLOR_PALETTE['B'], markersize=2, markeredgecolor=COLOR_PALETTE['B'])
	
	xs = []
	x = 0
	boxplots = []
	for bin_name, gerps in bin_stats.iteritems():
		print bin_name, '| regions:', len(gerps)
		boxplot = plt.boxplot(gerps, positions=[x], widths=0.8, notch=False, showfliers=True, patch_artist=True, flierprops=flierprops)
		# Enhance the speed of PDF load
		remove_overlapping_boxplot_outliers(boxplot)

		boxplots.append(boxplot)
		x += 1
		xs.append(x)

	# Colour boxplots
	for boxplot in boxplots:
		for patch in boxplot['boxes']:
			patch.set_facecolor(C_LIGHT_GRAY) #COLOR_PALETTE['B']

		for patch in boxplot['medians']:
			patch.set_color(C_BLACK) # 'yellow'
	
	# Plot expected (median) and set axis ticks and labels
	if subplot_letter:
		ax.text(-0.13, 1.11, subplot_letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='right')
	else:
		ax.text(-0.13, 1.11, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='right')

	xs = list(xs)
	normal_xs = [-0.8] + xs + [max(xs) + 0.8]
	plt.xticks(xs, bin_names, fontsize=7)
	plt.yticks([-12, -9, -6, -3, 0, 3, 6], fontsize=7)
	plt.xticks(range(-1, len(xs), 1), [''] + bin_names)
	plt.ylabel('GERP++ (mean)', fontsize=7)
	plt.xlabel('VIR length bins', fontsize=7)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)


def draw_regions_gerp_figure(db, subplot_num, subplot_letter='', row_num=3):
	region_length_gerps_strs = db.gevir.temp.find_one({'_id': 'region_length_gerps'})

	if not region_length_gerps_strs:
		gnomad_transcript_ids = set([])

		if INCLUDE_GNOMAD_OUTLIERS:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True })
		else:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

		for gnomad_gene in gnomad_genes:
			gnomad_transcript_ids.add(gnomad_gene['_id'])

		region_length_gerps = {}
		regions = db.gevir[REGIONS_COLLECTION].find({"exome_coverage": { "$gte": 50.0 }, "lenght": { "$gte": 1 }})

		total_lines = regions.count()
		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for region in regions:
			if region['transcript_id'] not in gnomad_transcript_ids:
				line_number += 1
				continue
			length = region['lenght']
			if length > 20:
				length = 21
			if length not in region_length_gerps:
				region_length_gerps[length] = []
			region_length_gerps[length].append(region['gerp_mean'])

			line_number += 1
			bar.update((line_number + 0.0) / total_lines)
		bar.finish()

		region_length_gerps_means = OrderedDict()
		region_length_gerps_sems = OrderedDict()
		for length in region_length_gerps:
			region_length_gerps_means[str(length)] = np.mean(region_length_gerps[length])
			region_length_gerps_sems[str(length)] = stats.sem(region_length_gerps[length])
		db.gevir.temp.insert({'_id': 'region_length_gerps', 'means': region_length_gerps_means, 'sems': region_length_gerps_sems })

	region_length_gerps_strs = db.gevir.temp.find_one({'_id': 'region_length_gerps'})
	region_length_gerps_means = region_length_gerps_strs['means']
	region_length_gerps_sems = region_length_gerps_strs['sems']

	xs = range(1, 22)
	ys = []
	ax = plt.subplot(row_num,2,subplot_num)

	for x in xs:
		ys.append(region_length_gerps_means[str(x)])

	if subplot_letter:
		ax.text(-0.08, 1.11, subplot_letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='right')
	else:
		ax.text(-0.08, 1.11, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='right')

	plt.errorbar(xs, ys, yerr=region_length_gerps_sems.values(), capsize=4, color=COLOR_PALETTE['B'], fmt='o', markersize=3) # 

	plt.xticks(range(0, 20, 3) + [21], fontsize=7)
	plt.yticks(fontsize=7)
	plt.ylim(0,4)
	ax.set_xticklabels([str(n) for n in range(0, 20, 3)] + ['21+'])

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	plt.ylabel('GERP++ (mean)', fontsize=7)
	plt.xlabel('VIR length', fontsize=7)

#############################################################################
### Figure 3, Supplementary Figures 2 and 4 : GeVIR vs gnomAD constraints ###
#############################################################################

def calculate_cumulative_percentage(ranked_genes, metric_set, precise=True):
	metric_num = float(len(metric_set))
	percentages = []

	# Calculate cumulative percentage for each gene (more accurate, but requires longer time to run)
	if precise:
		genes_so_far = set([])
		for gene in ranked_genes:
			genes_so_far.add(gene)
			percentages.append(len(genes_so_far & metric_set) / metric_num * 100)
	else:
		# Calculate cumulative percentage in bins for each %
		bins = [set(x) for x in np.array_split(ranked_genes, 100)]
		overlap_so_far = 0
		for x in range(0, 100):
			overlap_so_far += len(bins[x] & metric_set) / metric_num * 100
			percentages.append(overlap_so_far)
	return percentages


def calculate_cumulative_f1_percentage(ranked_genes, ad_set, ar_set, precise=True):
	all_ad = len(ad_set)
	f1_percentiles = []

	if precise:
		genes_so_far = set([])
		for gene in ranked_genes:
			genes_so_far.add(gene)
			ad = len(ad_set & genes_so_far)
			ar = len(ar_set & genes_so_far)
			metrics = get_precision_recall_f1(ad, ar, all_ad)
			f1_percentiles.append(metrics[2])
	else:
		ad_so_far = set([])
		ar_so_far = set([])
		bins = [set(x) for x in np.array_split(ranked_genes, 100)]
		for x in range(0, 100):
			ad_so_far |= bins[x] & ad_set
			ar_so_far |= bins[x] & ar_set

			ad = len(ad_so_far)
			ar = len(ar_so_far)
			metrics = get_precision_recall_f1(ad, ar, all_ad)
			f1_percentiles.append(metrics[2])

	return f1_percentiles


def calculate_similiarity_percentage(ranked_genes_1, ranked_genes_2):
	overlap_so_far = 0
	ranked_genes_1_so_far = set([])
	ranked_genes_2_so_far = set([])
	percentages = []

	bins_1 = [set(x) for x in np.array_split(ranked_genes_1, 100)]
	bins_2 = [set(x) for x in np.array_split(ranked_genes_2, 100)]
	for x in range(0, 100):
		ranked_genes_1_so_far |= bins_1[x]
		ranked_genes_2_so_far |= bins_2[x]

		overlap = len(ranked_genes_1_so_far & ranked_genes_2_so_far) / float(len(ranked_genes_1_so_far)) * 100
		percentages.append(overlap)

	return percentages


def draw_cumulative_percentage_subplot(subplot_num, scores_ys, title='', linestyle='-', report_auc=False, show_expected=False, report_peak=False, legend_loc='lower right', legend_order_reverse=True, precise=True):
	ax = plt.subplot(3,3,subplot_num)
	ax.set_title(title, loc='center', fontsize=7)
	ax.text(-0.25, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')

	if precise:
		points_num = len(scores_ys[scores_ys.keys()[0]])
		xs_raw = range(0, points_num)
		points_num = float(points_num)
		xs = []
		for x in xs_raw:
			xs.append(((x + 1) / points_num) * 100)
	else:
		xs = range(1,101)

	if report_peak:
		score_peaks = OrderedDict()
		for score_name, ys in scores_ys.iteritems():
			max_score_y = 0
			max_score_x = 0
			if precise:
				for i in range(0, len(xs)):
					score_y = ys[i]
					x = xs[i]
					if score_y > max_score_y:
						max_score_y = score_y
						max_score_x = x					
			else:
				for x in xs:
					score_y = ys[x-1]
					if score_y > max_score_y:
						max_score_y = score_y
						max_score_x = x
			score_peaks[score_name] = (max_score_x, max_score_y)

	scores_auc = OrderedDict()
	scores_xy = OrderedDict()
	scores_labels = OrderedDict()
	scores_colors = []
	for score_name, ys in scores_ys.iteritems():
		score_auc = sklearn_auc(xs, ys) / 100
		scores_auc[score_name] = sklearn_auc(xs, ys)
		scores_xy[score_name] = (xs, ys)
		if report_auc:
			label = score_name + '\nAUC: %.2f%%' % score_auc
		elif report_peak:
			max_score_x, max_score_y = score_peaks[score_name]
			scores_auc[score_name] = max_score_y
			label = score_name + '\nPeak: %.2f%%' % max_score_y + ' (F1), ' + ('%.2f%%' % max_score_x) + ' (Rank)'
		else:
			label = score_name
		scores_labels[score_name] = label

	if report_auc or report_peak:
		scores_auc = sort_dict_by_values(scores_auc, reverse=legend_order_reverse)

	for score_name in scores_auc:
		xs, ys = scores_xy[score_name]
		plt.plot(xs, ys, color=SCORE_COLORS[score_name], linestyle=linestyle, label=scores_labels[score_name], linewidth=1)
		if report_peak:
			max_score_x, max_score_y = score_peaks[score_name]
			plt.scatter([max_score_x], [max_score_y], s=9, color=SCORE_COLORS[score_name])
			score_peaks
		scores_colors.append(SCORE_COLORS[score_name])

	#if report_auc or report_peak:
	l = plt.legend(loc=legend_loc, frameon=False, fontsize=5, handlelength=1.0)
	for line in l.legendHandles:
		line.set_linewidth(2.0)

	if show_expected:
		plt.plot(xs, xs, '--', label='Expected', color=C_GRAY, linewidth=1)

	plt.xticks(range(0, 101, 10), fontsize=7)
	plt.yticks(fontsize=7)

	ax.set_xticklabels([str(n) for n in range(0, 110, 10)])
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.ylabel('Cumulative percentage (%)', fontsize=7)
	plt.xlabel('Rank (%)', fontsize=7) # Percentiles


class FoldEnrichmentAR():
	def __init__(self):
		self.decile = 0
		self.ar = 0
		self.genes = 0
		self.fe = 0.0
		self.p_value = 0.0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['decile'] = self.decile
		dictionary['ar'] = self.ar
		dictionary['genes'] = self.genes
		dictionary['fe'] = self.fe
		dictionary['p_value'] = self.p_value
		return dictionary


def draw_gaps_vs_miss_z_ad_ar_subplot(subplot_num, scores_data, metric_set, ylabel='', title='', linestyle='-', ax='', for_letter=False):
	print 'AR ANALYSIS'
	font_size = 7

	total_ar = len(metric_set)
	total_genes = set([])
	for score_set in scores_data[scores_data.keys()[0]]:
		total_genes |= score_set
	total_genes = len(total_genes)
	print 'AR', total_ar, 'TOTAL', total_genes

	n_bars = len(scores_data)
	xs = np.arange(len(scores_data.values()[0]))#range(0, n_bars)
	bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(n_bars)
	metric_num = float(len(metric_set))

	if metric_num == 0:
		metric_num = 1

	if not ax:
		ax = plt.subplot(1,3,subplot_num)

	ax.set_title(title, loc='center', fontsize=font_size)

	if for_letter:
		ax.text(-0.11, 1.05, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')
	else:
		ax.text(-0.2, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	top_deciles = set([1,2])
	mid_deciles = set([4,5,6])
	last_deciles = set([8,9,10])

	top_row = ['1, 2 (first 20%)']
	mid_row = ['4, 5, 6 (mid 30%)']
	last_row = ['8, 9, 10 (last 30%)']

	score_fears = OrderedDict()
	ys_scores = []
	for x in range(0, n_bars):
		score_name = scores_data.keys()[x]
		score_fears[score_name] = OrderedDict()
		ys_score = []
		decile = 1

		top_ar = set([])
		top_all = set([])

		mid_ar = set([])
		mid_all = set([])

		last_ar = set([])
		last_all = set([])
		for score_set in scores_data[score_name]:
			n = len(score_set & metric_set)
			p = n * 100 / metric_num
			ys_score.append(p)

			fear = FoldEnrichmentAR()
			fear.decile = decile
			fear.ar = n
			fear.genes = len(score_set)
			fear.fe, fear.p_value = fisher_exact([[fear.ar, fear.genes],[total_ar, total_genes]])
			score_fears[score_name][decile] = fear

			if decile in top_deciles:
				top_ar |= score_set & metric_set
				top_all |= score_set

			if decile in mid_deciles:
				mid_ar |= score_set & metric_set
				mid_all |= score_set

			if decile in last_deciles:
				last_ar |= score_set & metric_set
				last_all |= score_set

			decile += 1

		ys_scores.append(ys_score)

		top_ar_fe, top_ar_p_value = fisher_exact([[len(top_ar), len(top_all)],[total_ar, total_genes]])
		mid_ar_fe, mid_ar_p_value = fisher_exact([[len(mid_ar), len(mid_all)],[total_ar, total_genes]])
		last_ar_fe, last_ar_p_value = fisher_exact([[len(last_ar), len(last_all)],[total_ar, total_genes]])

		top_row += [len(top_ar), len(top_all), top_ar_fe, top_ar_p_value]
		mid_row += [len(mid_ar), len(mid_all), mid_ar_fe, mid_ar_p_value]
		last_row += [len(last_ar), len(last_all), last_ar_fe, last_ar_p_value]
		if for_letter:
			plt.plot(xs, ys_score, linestyle=linestyle, color=SCORE_COLORS[score_name], label=score_name) # .replace('\n', ' ')
			plt.scatter(xs, ys_score, s=5, color=SCORE_COLORS[score_name])
		else:
			plt.plot(xs, ys_score, linestyle=linestyle, color=SCORE_COLORS[score_name], label=score_name, linewidth=1)
			plt.scatter(xs, ys_score, s=3, color=SCORE_COLORS[score_name])			

	xs = list(xs)
	normal_xs = [-0.8] + xs + [max(xs) + 0.8]
	if for_letter:
		plt.plot(normal_xs, [10]*len(normal_xs), '--', label='Expected', color=C_GRAY)
	else:
		plt.plot(normal_xs, [10]*len(normal_xs), '--', label='Expected', color=C_GRAY, linewidth=1)
	plt.xticks(range(0, 10), fontsize=font_size)
	plt.yticks(fontsize=font_size)

	ax.set_xticklabels([str(n) for n in range(10, 110, 10)])

	if for_letter:
		plt.ylabel('Percentage (%) of AR genes', fontsize=font_size)
	else:
		plt.ylabel('Genes from a group (%)', fontsize=font_size)

	if ylabel:
		plt.ylabel(ylabel)
	plt.xlabel('Rank (%)', fontsize=font_size)

	if for_letter:
		l = plt.legend(loc='lower center', frameon=False, fontsize=6)
	else:
		l = plt.legend(loc='lower center', frameon=False, fontsize=5, handlelength=1)
		for line in l.legendHandles:
			line.set_linewidth(2.0)

	# AR Fold Enrichment Report
	headers = ['Decile']
	rows = OrderedDict()
	for score_name in scores_data.keys():
		for col_name in ['AR', 'genes', 'FE', 'p-value']:
			headers.append(score_name + ' ' + col_name)

	rows = OrderedDict()
	for decile in range(1, 11):
		rows[decile] = []
		for score_name, score_fear in score_fears.iteritems():
			score_fear = score_fear[decile]
			rows[decile] += [score_fear.ar, score_fear.genes, score_fear.fe, score_fear.p_value]

	table = [headers]
	for decile, row in rows.iteritems():
		row = [decile] + row
		table.append(row)

	# Add combined deciles stats
	table.append(top_row)
	table.append(mid_row)
	table.append(last_row)

	output_csv = OUTPUT_FOLDER + 'f3_ar_stats.csv'
	write_table_to_csv(table, output_csv)


def draw_gaps_vs_miss_z_score_length_subplot(subplot_num, scores_data, score_sets, title='', ax='', scores_ranked_transcripts={}, for_letter=False):
	font_size = 7
	if not ax:
		ax = plt.subplot(1,3,subplot_num)

	ax.set_title(title, loc='center', fontsize=font_size)

	if for_letter:
		ax.text(-0.13, 1.05, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')
	else:
		ax.text(-0.25, 1.1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	print 'Length Analysis'
	normal_length = score_sets.get_transcripts_median_length(score_sets.all_transcripts)
	print 'Median', normal_length

	# Calculate correlation
	scores_rs = OrderedDict()
	if scores_ranked_transcripts:
		for score_name, ranked_transcripts in scores_ranked_transcripts.iteritems():
			score_ranks = []
			score_lengths = []
			for x in range(0, len(ranked_transcripts)):
				score_ranks.append(x)
				score_lengths.append(score_sets.length_dict[ranked_transcripts[x]])

			scores_rs[score_name] = spearmanr(score_ranks, score_lengths)[0]

	scores_rs = sort_dict_by_values(scores_rs, reverse=True)

	xs = np.arange(len(scores_data.values()[0]))

	scores_colors = []
	for score_name, score_rs in scores_rs.iteritems():
		bins_transcripts = scores_data[score_name]
		ys = []
		for bin_transcripts in bins_transcripts:
			ys.append(score_sets.get_transcripts_median_length(bin_transcripts))

		print score_name, ys[0], float(ys[0]) / normal_length

		label = score_name
		if scores_ranked_transcripts:
			label = score_name + ' (r=%.2f)' % score_rs #scores_rs[score_name]

		if for_letter:
			plt.plot(xs, ys, color=SCORE_COLORS[score_name], label=label)
			plt.scatter(xs, ys, s=5, color=SCORE_COLORS[score_name])
		else:
			plt.plot(xs, ys, color=SCORE_COLORS[score_name], label=label, linewidth=1)
			plt.scatter(xs, ys, s=3, color=SCORE_COLORS[score_name])

		scores_colors.append(SCORE_COLORS[score_name])
		print score_name


	plt.ylabel('Protein length (median)', fontsize=font_size)
	plt.xlabel('Rank (%)', fontsize=font_size)


	xs = list(xs)
	normal_xs = [-0.8] + xs + [max(xs) + 0.8]
	print normal_xs
	if for_letter:
		plt.plot(normal_xs, [normal_length]*len(normal_xs), '--', label='Median (all genes)', color=C_GRAY)
	else:
		plt.plot(normal_xs, [normal_length]*len(normal_xs), '--', label='Median (all genes)', color=C_GRAY, linewidth=1)

	if scores_ranked_transcripts:
		if for_letter:
			l = plt.legend(loc='upper right', frameon=False, fontsize=6)
		else:
			l = plt.legend(loc='upper right', frameon=False, fontsize=5, handlelength=1)
			for line in l.legendHandles:
				line.set_linewidth(2.0)

	plt.xticks(range(0, 10), fontsize=font_size)
	plt.yticks(fontsize=font_size)
	ax.set_xticklabels([str(n) for n in range(10, 110, 10)])


def draw_gaps_vs_gnomad_constraint_scores(db, for_letter=False, clean=True):
	omim_sets = OmimSets(db)
	score_sets = ScoreSets(db, filters={"is_gnomad": True})
	essential_sets = EssentialSets(db)

	ad_set = omim_sets.ad & score_sets.all_transcripts
	ar_set = omim_sets.ar & score_sets.all_transcripts
	null_set = essential_sets.nulls & score_sets.all_transcripts

	mouse_het_lethal_set = essential_sets.mouse_het_lethal & score_sets.all_transcripts
	crispr_essential_set = essential_sets.crispr_essential & score_sets.all_transcripts
	crispr_non_essential_set = essential_sets.crispr_non_essential & score_sets.all_transcripts

	score_ranked_lists = OrderedDict()
	score_ranked_lists[MY_NAME] = score_sets.gevir_dict.keys()
	score_ranked_lists[MY_NO_GERP_NAME] = score_sets.gevir_no_gerp_dict.keys()
	score_ranked_lists[MISS_Z_NAME] = score_sets.gnomad_miss_z_dict.keys()
	score_ranked_lists[MISS_OE_NAME] = score_sets.gnomad_oe_mis_upper_dict.keys()
	score_ranked_lists[LOF_OE_NAME] = score_sets.gnomad_oe_lof_upper_dict.keys()
	score_ranked_lists[COMBINED_NAME] = score_sets.combined_rank_dict.keys()
	# pLI and UNEECON metrics were not used in the evoluation in the original manuscript 
	#score_ranked_lists[PLI_NAME] = score_sets.gnomad_pli_dict.keys()
	#score_ranked_lists[UNEECON_G_NAME] = score_sets.uneecon_g_dict.keys()

	score_names = [MY_NAME, MY_NO_GERP_NAME, MISS_Z_NAME, MISS_OE_NAME, LOF_OE_NAME, COMBINED_NAME] #, PLI_NAME UNEECON_G_NAME

	score_ad_percentages = OrderedDict()
	score_ar_percentages = OrderedDict()
	score_lof_hom_percentages = OrderedDict()
	score_f1_percentages = OrderedDict()
	score_similiarity_percentages = OrderedDict()
	score_decils = OrderedDict()

	score_mouse_het_lethal_percentages = OrderedDict()
	score_crispr_essential_percentages = OrderedDict()
	score_crispr_non_essential_percentages = OrderedDict()


	if INCLUDE_GNOMAD_OUTLIERS:
		temp_data_id = 'gaps_vs_gnomad_constraint_scores'
	else:
		temp_data_id = 'gaps_vs_gnomad_constraint_scores_no_outliers'

	data = db.gevir.temp.find_one({'_id': temp_data_id})
	if clean or not data:
		for score_name in score_names:
			score_ad_percentages[score_name] = calculate_cumulative_percentage(score_ranked_lists[score_name], ad_set)
			score_ar_percentages[score_name] = calculate_cumulative_percentage(score_ranked_lists[score_name], ar_set)
			score_lof_hom_percentages[score_name] = calculate_cumulative_percentage(score_ranked_lists[score_name], null_set)
			score_f1_percentages[score_name] = calculate_cumulative_f1_percentage(score_ranked_lists[score_name], ad_set, ar_set)

			score_mouse_het_lethal_percentages[score_name] = calculate_cumulative_percentage(score_ranked_lists[score_name], mouse_het_lethal_set)
			score_crispr_essential_percentages[score_name] = calculate_cumulative_percentage(score_ranked_lists[score_name], crispr_essential_set)
			score_crispr_non_essential_percentages[score_name] = calculate_cumulative_percentage(score_ranked_lists[score_name], crispr_non_essential_set)
			
			if score_name != MY_NAME:
				print score_name
				score_similiarity_percentages[score_name] = calculate_similiarity_percentage(score_ranked_lists[score_name], score_ranked_lists[MY_NAME])

			score_decils[score_name] = [set(x) for x in np.array_split(score_ranked_lists[score_name], 10)]

		data = OrderedDict()
		data['score_ad_percentages'] = score_ad_percentages
		data['score_ar_percentages'] = score_ar_percentages
		data['score_lof_hom_percentages'] = score_lof_hom_percentages
		data['score_f1_percentages'] = score_f1_percentages
		data['score_mouse_het_lethal_percentages'] = score_mouse_het_lethal_percentages
		data['score_crispr_essential_percentages'] = score_crispr_essential_percentages
		data['score_crispr_non_essential_percentages'] = score_crispr_non_essential_percentages
		data['score_similiarity_percentages'] = score_similiarity_percentages
		data['score_ranked_lists'] = score_ranked_lists
	else:
		data = data['data']
		score_ad_percentages = data['score_ad_percentages']
		score_ar_percentages = data['score_ar_percentages']
		score_lof_hom_percentages = data['score_lof_hom_percentages']
		score_f1_percentages = data['score_f1_percentages']
		score_mouse_het_lethal_percentages = data['score_mouse_het_lethal_percentages']
		score_crispr_essential_percentages = data['score_crispr_essential_percentages']
		score_crispr_non_essential_percentages = data['score_crispr_non_essential_percentages']
		score_similiarity_percentages = data['score_similiarity_percentages']
		score_ranked_lists = data['score_ranked_lists']
		for score_name in score_names:
			score_decils[score_name] = [set(x) for x in np.array_split(score_ranked_lists[score_name], 10)]

	fig = plt.figure()

	subplot_num = 1
	if for_letter:
		ax = plt.subplot(1,2,subplot_num)
		draw_gaps_vs_miss_z_ad_ar_subplot(subplot_num, score_decils, ar_set, ax=ax, for_letter=True)
		subplot_num += 1
		ax = plt.subplot(1,2,subplot_num)
		draw_gaps_vs_miss_z_score_length_subplot(subplot_num, score_decils, score_sets, ax=ax, scores_ranked_transcripts=score_ranked_lists, for_letter=True)

	else:	
		draw_cumulative_percentage_subplot(subplot_num, score_ad_percentages, title='Autosomal Dominant (AD) (n={:,})'.format(len(ad_set)), report_auc=True, show_expected=True)
		subplot_num += 1
		draw_cumulative_percentage_subplot(subplot_num, score_mouse_het_lethal_percentages, title='Mouse het lethal knockout (n={:,})'.format(len(mouse_het_lethal_set)), report_auc=True, show_expected=True)
		subplot_num += 1
		draw_cumulative_percentage_subplot(subplot_num, score_crispr_essential_percentages, title='Cell essential (n={:,})'.format(len(crispr_essential_set)), report_auc=True, show_expected=True)
		subplot_num += 1

		draw_cumulative_percentage_subplot(subplot_num, score_lof_hom_percentages, title='Nulls (n={:,})'.format(len(null_set)), report_auc=True, show_expected=True, legend_order_reverse=False, legend_loc='upper left')
		subplot_num += 1
		draw_cumulative_percentage_subplot(subplot_num, score_crispr_non_essential_percentages, title='Cell non-essential (n={:,})'.format(len(crispr_non_essential_set)), report_auc=True, show_expected=True, legend_order_reverse=False, legend_loc='upper left')
		subplot_num += 1
		ax = plt.subplot(3,3,subplot_num)
		draw_gaps_vs_miss_z_ad_ar_subplot(subplot_num, score_decils, ar_set, title='Autosomal Recessive (AR) (n={:,})'.format(len(ar_set)), ax=ax)
		subplot_num += 1

		draw_cumulative_percentage_subplot(subplot_num, score_f1_percentages, title='AD/AR classification AD F1; All (n={:,})'.format(len(score_sets.all_transcripts)), report_peak=True)
		subplot_num += 1
		draw_cumulative_percentage_subplot(subplot_num, score_similiarity_percentages, title='Similarity with ' + MY_NAME + '; All (n={:,})'.format(len(score_sets.all_transcripts)), report_auc=False)
		subplot_num += 1
		ax = plt.subplot(3,3,subplot_num)
		draw_gaps_vs_miss_z_score_length_subplot(subplot_num, score_decils, score_sets, title='All (n={:,})'.format(len(score_sets.all_transcripts)), ax=ax, scores_ranked_transcripts=score_ranked_lists)

	# Legend
	my_patch = mpatches.Patch(color=SCORE_COLORS[MY_NAME])
	my_no_gerp_patch = mpatches.Patch(color=SCORE_COLORS[MY_NO_GERP_NAME])
	miss_z_patch = mpatches.Patch(color=SCORE_COLORS[MISS_Z_NAME])
	miss_oe_patch = mpatches.Patch(color=SCORE_COLORS[MISS_OE_NAME])
	lof_oe_patch = mpatches.Patch(color=SCORE_COLORS[LOF_OE_NAME])
	combined_patch = mpatches.Patch(color=SCORE_COLORS[COMBINED_NAME])
	#pli_patch = mpatches.Patch(color=SCORE_COLORS[PLI_NAME])
	#uneecon_g_patch = mpatches.Patch(color=SCORE_COLORS[UNEECON_G_NAME])

	patches = OrderedDict()
	patches[MY_NAME] = my_patch
	patches[MY_NO_GERP_NAME] = my_no_gerp_patch
	patches[MISS_Z_NAME] = miss_z_patch
	patches[MISS_OE_NAME] = miss_oe_patch
	patches[LOF_OE_NAME] = lof_oe_patch
	patches[COMBINED_NAME] = combined_patch
	#patches[PLI_NAME] = pli_patch
	#patches[UNEECON_G_NAME] = uneecon_g_patch

	patches['Expected'] = Line2D([0], [0], linestyle='--', color=C_GRAY)
	
	if for_letter:
		fig.set_size_inches(7, 3)
		plt.tight_layout(rect=[0, 0, 1, 1])
		plt.savefig(FIGURES_FOLDER + 'small_GeVIR_vs_gnomAD_metrics_AUC.pdf', format='pdf', dpi=300)
	else:	
		fig.set_size_inches(7, 8)
		plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99])
		plt.savefig(FIGURES_FOLDER + 'GeVIR_vs_gnomAD_metrics_AUC.png', format='png', dpi=300)
	plt.close()

#####################################################
### Figure 4: GeVIR vs LOEUF (lof_oe)  (Top ~15%) ###
#####################################################

def draw_venn_gaps_vs_lof_oe_subplot(db, subplot_num, gaps_set, oe_lof_set, combined_set, set_names, title=''):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(0, 1, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')

	ax.set_title(title, loc='center', fontsize=7)
	out = venn3([gaps_set, oe_lof_set, combined_set], set_names)
	for text in out.set_labels:
		text.set_fontsize(7)

	for text in out.subset_labels:
		if text:
			text.set_text("{:,}".format(int(text.get_text())))
			text.set_fontsize(7)
			text.set_color('white')

	black_subset_indexes = [0, 2]
	for black_subset_index in black_subset_indexes:
		if out.subset_labels[black_subset_index]: # check if subset is not none
			out.subset_labels[black_subset_index].set_color('black')

	patches = ['100', '110', '101', '111', '001', '010', '011']

	for patch in patches:
		subset = out.get_patch_by_id(patch)
		if subset:
			subset.set_edgecolor('none')
			subset.set_alpha(0.8)
			'''
			Alternative Colours
			#627676 - orange + blue
			#789a32 - orange + green
			#268866 - blue + green
			#528a59 - all three

			# Color blind friendly
			'''

			if patch == '100':
				subset.set_color(COLOR_PALETTE['O'])
			elif patch == '010':
				subset.set_color(COLOR_PALETTE['B'])
			elif patch == '001':
				subset.set_color(COLOR_PALETTE['G'])
			if patch == '110': # O + B
				subset.set_color('#71875a')
			if patch == '101': # O + G
				subset.set_color('#599c46')
				subset.set_alpha(0.9)
			if patch == '011': # B + G
				subset.set_color('#008f8d')
				subset.set_alpha(0.9)
			if patch == '111': # B + G + O
				subset.set_color('#4c8b78')
				subset.set_alpha(1)


def draw_evaluation_gaps_vs_lof_oe_subplot(db, subplot_num, gaps_set, oe_lof_set, combined_set, ad_set, ar_set, set_names, title='', y_lim_coef=1.1):
	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.08, 1.05, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')
	ax.set_title(title, loc='center', fontsize=7)
	gaps_metrics = get_precision_recall_f1(len(gaps_set & ad_set), len(gaps_set & ar_set), len(ad_set))
	oe_lof_metrics = get_precision_recall_f1(len(oe_lof_set & ad_set), len(oe_lof_set & ar_set), len(ad_set))
	combined_metrics = get_precision_recall_f1(len(combined_set & ad_set), len(combined_set & ar_set), len(ad_set))

	bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(3)
	xs = np.arange(3)

	plt.bar(xs + x_paddings[0], gaps_metrics, width=bar_width, label=set_names[0], color=COLOR_PALETTE['O'])
	plt.bar(xs + x_paddings[1], oe_lof_metrics, width=bar_width, label=set_names[1], color=COLOR_PALETTE['B'])
	plt.bar(xs + x_paddings[2], combined_metrics, width=bar_width, label=set_names[2], color=COLOR_PALETTE['G'])

	max_y = max(gaps_metrics + oe_lof_metrics + combined_metrics)
	y_padding = max_y / 40.0

	# Add text annotation
	for x in range(0, len(gaps_metrics)):
		p_text = "%.2f" % gaps_metrics[x] 
		plt.text(x=x + x_paddings[0], y=gaps_metrics[x] + y_padding, s=p_text, size=7, color='black', ha='center')

	for x in range(0, len(oe_lof_metrics)):
		p_text = "%.2f" % oe_lof_metrics[x] 
		plt.text(x=x + x_paddings[1], y=oe_lof_metrics[x] + y_padding, s=p_text, size=7, color='black', ha='center')

	for x in range(0, len(combined_metrics)):
		p_text = "%.2f" % combined_metrics[x] 
		plt.text(x=x + x_paddings[2], y=combined_metrics[x] + y_padding, s=p_text, size=7, color='black', ha='center')

	plt.ylim(top=max_y * y_lim_coef)
	plt.yticks(fontsize=7)
	ax.set_xticks(xs)
	ax.set_xticklabels(['Precision', 'Recall', 'F1 score'], fontsize=7)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	plt.tick_params(axis='x',          # changes apply to the x-axis
				    which='both',      # both major and minor ticks are affected
				    bottom=False,      # ticks along the bottom edge are off
				    top=False,         # ticks along the top edge are off
				    labelbottom=True) # labels along the bottom edge are off

	plt.legend(loc='upper center', frameon=False, fontsize=7)


def draw_3_scores_venn_plot(db, score_3_sets, all_transcripts, top_num, top_percent, fig_name, performance_plot_y_lim_coef=1.1):
	omim_sets = OmimSets(db)

	ad_set = omim_sets.ad & all_transcripts
	ar_set = omim_sets.ar & all_transcripts

	set_1, set_2, set_3 = score_3_sets.values()
	set_names = score_3_sets.keys()

	fig = plt.figure(figsize = (2,2))

	subplot_num = 1

	draw_venn_gaps_vs_lof_oe_subplot(db, subplot_num, set_1, set_2, set_3, set_names)
	subplot_num += 1

	draw_venn_gaps_vs_lof_oe_subplot(db, subplot_num, set_1 & omim_sets.ad, set_2 & omim_sets.ad, set_3 & omim_sets.ad, set_names)
	subplot_num += 1

	draw_venn_gaps_vs_lof_oe_subplot(db, subplot_num, set_1 & omim_sets.ar, set_2 & omim_sets.ar, set_3 & omim_sets.ar, set_names)
	subplot_num += 1

	draw_evaluation_gaps_vs_lof_oe_subplot(db, subplot_num, set_1, set_2, set_3, ad_set, ar_set, set_names, y_lim_coef=performance_plot_y_lim_coef)

	plt.ylabel('Percentage (%)', fontsize=7)
	
	fig.set_size_inches(7, 4.4)
	plt.tight_layout(rect=[-0.05, 0, 1, 1])
	plt.savefig(FIGURES_FOLDER + fig_name, format='pdf', dpi=300)
	plt.close()


def draw_gaps_vs_lof_oe(db, full_report=False):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	highly_intolerant_genes = db.gevir.common_gene_scores.find({ "is_gnomad": True, "gnomad_oe_lof_upper": { "$lt": 0.35 } })
	top_num = highly_intolerant_genes.count()
	top_percent = int(round(float(top_num) * 100 / len(score_sets.all_transcripts)))

	score_3_sets = OrderedDict()
	score_3_sets[MY_NAME] = set(score_sets.gevir_dict.keys()[:top_num])
	score_3_sets[LOF_OE_NAME] = set(score_sets.gnomad_oe_lof_upper_dict.keys()[:top_num])
	score_3_sets[COMBINED_NAME] = set(score_sets.combined_rank_dict.keys()[:top_num])

	draw_3_scores_venn_plot(db, score_3_sets, score_sets.all_transcripts, top_num, top_percent, 'GeVIR_vs_oe_LoF_top.pdf')
	if full_report:
		report_top_important_genes(db, top_set=score_3_sets[MY_NAME] & score_3_sets[LOF_OE_NAME])


def get_gnomad_pli_gte_09_genes(db):
	transcript_ids = set()
	genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "pLI": { "$gt": 0.9 } })
	for gene in genes:
		transcript_ids.add(gene['transcript'])
	return transcript_ids


def draw_gaps_vs_lof_oe_vs_pli(db):
	pli_set = get_gnomad_pli_gte_09_genes(db)
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	top_num = len(pli_set)
	top_percent = int(round(float(top_num) * 100 / len(score_sets.all_transcripts)))

	score_3_sets = OrderedDict()
	score_3_sets[MY_NAME] = set(score_sets.gevir_dict.keys()[:top_num])
	score_3_sets[LOF_OE_NAME] = set(score_sets.gnomad_oe_lof_upper_dict.keys()[:top_num])
	score_3_sets[PLI_NAME] = pli_set

	lof_common_num = len(score_3_sets[LOF_OE_NAME] & score_3_sets[PLI_NAME])
	print 'LOEUF & pLI', lof_common_num, '/', len(score_3_sets[LOF_OE_NAME]), ';', lof_common_num / float(len(score_3_sets[LOF_OE_NAME]))

	draw_3_scores_venn_plot(db, score_3_sets, score_sets.all_transcripts, top_num, top_percent, 'GeVIR_vs_oe_LoF_vs_pLI_top.png', performance_plot_y_lim_coef=1.3)


############################################
### Extended Data Figure 1: GeVIR vs CCR ###
############################################

class GeVIRvsCCRs():
	def __init__(self):
		self.gevir_ads = []
		self.gevir_ars = []
		self.ccrs_ads = []
		self.ccrs_ars = []
		self.common_ads = []
		self.common_ars = []
		self.gevir_f1s = []
		self.ccrs_f1s = []
		self.gevir_no_gerp_ads = []
		self.gevir_no_gerp_ars = []
		self.gevir_no_gerp_f1s = []

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['gevir_ads'] = self.gevir_ads
		dictionary['gevir_ars'] = self.gevir_ars
		dictionary['ccrs_ads'] = self.ccrs_ads
		dictionary['ccrs_ars'] = self.ccrs_ars
		dictionary['common_ads'] = self.common_ads
		dictionary['common_ars'] = self.common_ars
		dictionary['gevir_f1s'] = self.gevir_f1s
		dictionary['ccrs_f1s'] = self.ccrs_f1s

		dictionary['gevir_no_gerp_ads'] = self.gevir_no_gerp_ads
		dictionary['gevir_no_gerp_ars'] = self.gevir_no_gerp_ars
		dictionary['gevir_no_gerp_f1s'] = self.gevir_no_gerp_f1s

		return dictionary

def get_f1_peak(f1_list):
	index = 0
	max_f1 = 0
	for x in range(0, len(f1_list)):
		f1 = f1_list[x]
		if f1 > max_f1:
			max_f1 = f1
			index = x
	return index, max_f1


def draw_gevir_vs_ccrs_subplot(subplot_num, gevir_ys, gevir_no_gerp_ys, ccrs_ys, title='', x_label='', y_label='', report_peak=False, all_genes_num=1):
	gevir_colour = SCORE_COLORS[MY_NAME]
	gevir_no_gerp_colour = SCORE_COLORS[MY_NO_GERP_NAME]
	ccrs_colour = COLOR_PALETTE['B']

	ax = plt.subplot(2,2,subplot_num)
	ax.text(-0.12, 1.08, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')
	ax.set_title(title, loc='center', fontsize=7)
	xs = range(1, len(gevir_ys) + 1)

	if report_peak:
		gevir_peak_gene, gevir_max_f1 = get_f1_peak(gevir_ys)
		gevir_peak_gene_percentile = (gevir_peak_gene * 100) / float(all_genes_num)
		gevir_label = MY_NAME +' Peak: %.2f%%' % gevir_max_f1 + ' (F1),\nat ' + str(gevir_peak_gene) + ', %.2f%%' % gevir_peak_gene_percentile + ' (Rank)'
		gevir_no_gerp_peak_gene, gevir_no_gerp_max_f1 = get_f1_peak(gevir_no_gerp_ys)
		gevir_no_gerp_peak_gene_percentile = (gevir_no_gerp_peak_gene * 100) / float(all_genes_num)
		gevir_no_gerp_label = MY_NO_GERP_NAME.replace('\n', ' ') + '\nPeak: %.2f%%' % gevir_no_gerp_max_f1 + ' (F1),\nat ' + str(gevir_no_gerp_peak_gene) + ', %.2f%%' % gevir_no_gerp_peak_gene_percentile + ' (Rank)'
		ccr_peak_gene, ccr_max_f1 = get_f1_peak(ccrs_ys)
		ccr_peak_gene_percentile = (ccr_peak_gene * 100) / float(all_genes_num)
		ccr_label = CCRS_NAME + ' Peak: %.2f%%' % ccr_max_f1 + ' (F1),\nat ' + str(ccr_peak_gene) + ', %.2f%%' % ccr_peak_gene_percentile + ' (Rank)'
		plt.plot(xs, gevir_ys, color=gevir_colour, label=gevir_label)
		plt.plot(xs, gevir_no_gerp_ys, color=gevir_no_gerp_colour, label=gevir_no_gerp_label)
		plt.plot(xs, ccrs_ys, color=ccrs_colour, label=ccr_label)
		plt.scatter([gevir_peak_gene], [gevir_max_f1], color=gevir_colour, s=15)
		plt.scatter([gevir_no_gerp_peak_gene], [gevir_no_gerp_max_f1], color=gevir_no_gerp_colour, s=15)
		plt.scatter([ccr_peak_gene], [ccr_max_f1], color=ccrs_colour, s=15)
	else:
		plt.plot(xs, gevir_ys, color=gevir_colour, label=MY_NAME)
		plt.plot(xs, gevir_no_gerp_ys, color=gevir_no_gerp_colour, label=MY_NO_GERP_NAME)
		plt.plot(xs, ccrs_ys, color=ccrs_colour, label=CCRS_NAME)

	plt.xticks(fontsize=7)
	plt.yticks(fontsize=7)
	plt.ylabel(y_label, fontsize=7)
	plt.xlabel(x_label, fontsize=7)

	plt.xticks(range(0, 7001, 1000), ['{:,.0f}'.format(n) for n in range(0, 7001, 1000)])

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	scores_colors = [gevir_colour, gevir_no_gerp_colour, ccrs_colour]
	l = plt.legend(loc='lower right', frameon=False, fontsize=6)


def draw_gevir_vs_ccrs_length_subplot(subplot_num, scores_data, score_sets, title='', ax='', scores_ranked_transcripts={}):
	gevir_colour = SCORE_COLORS[MY_NAME]
	gevir_no_gerp_colour = SCORE_COLORS[MY_NO_GERP_NAME]
	ccrs_colour = SCORE_COLORS[CCRS_NAME]

	if not ax:
		ax = plt.subplot(2,2,subplot_num)

	ax.text(-0.14, 1.08, SUBPLOT_LETTERS[subplot_num], transform=ax.transAxes, fontsize=7, fontweight='bold', va='top', ha='right')
	ax.set_title(title, loc='center', fontsize=7)

	# Calculate correlation
	scores_rs = OrderedDict()
	if scores_ranked_transcripts:
		for score_name, ranked_transcripts in scores_ranked_transcripts.iteritems():
			score_ranks = []
			score_lengths = []
			for x in range(0, len(ranked_transcripts)):
				score_ranks.append(x)
				score_lengths.append(score_sets.length_dict[ranked_transcripts[x]])

			scores_rs[score_name] = spearmanr(score_ranks, score_lengths)[0]


	bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(3)
	
	score_n = 0

	score_boxplots = OrderedDict()
	for score_name, bins_transcripts in scores_data.iteritems():
		x = 0
		bins_transcripts = scores_data[score_name]
		for bin_transcripts in bins_transcripts:
			transcript_lengths = score_sets.get_transcripts_length(bin_transcripts)
			
			# Report median and 25/75 quartiles
			y = np.median(transcript_lengths)
			y_err_lower = np.percentile(transcript_lengths, 25)
			y_err_upper = np.percentile(transcript_lengths, 75)
			print score_name, x, y_err_lower, y_err_upper

			boxplot = plt.boxplot(transcript_lengths, positions=[x + x_paddings[score_n]], widths=bar_width*0.8, notch=True, showfliers=False, patch_artist=True)
			
			if x == 0:
				score_boxplots[score_name] = [boxplot]
			else:
				score_boxplots[score_name].append(boxplot)

			x += 1
		score_n += 1

	# Colour boxplots
	for score_name, boxplots in score_boxplots.iteritems():
		for boxplot in boxplots:
			for patch in boxplot['boxes']:
				patch.set_facecolor(SCORE_COLORS[score_name])

			for patch in boxplot['medians']:
				patch.set_color('yellow')

	# Plot expected (median) and set axis ticks and labels
	xs = np.arange(len(scores_data.values()[0]))
	xs = list(xs)
	normal_xs = [-0.8] + xs + [max(xs) + 0.8]
	plt.xticks(range(0, 7), fontsize=7)
	plt.yticks(fontsize=7)
	plt.xticks(range(-1, 7, 1), [''] + ['{:,.0f}'.format(n) for n in range(1000, 7001, 1000)])
	plt.ylabel('Protein length', fontsize=7)
	plt.xlabel('Rank (genes)', fontsize=7)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	normal_length = score_sets.get_transcripts_median_length(score_sets.all_transcripts)

	# Draw secondary legend (gene score spearman correlations)
	scores_order = [MY_NAME, MY_NO_GERP_NAME, CCRS_NAME]
	my_patch = mpatches.Patch(color=SCORE_COLORS[MY_NAME])
	my_no_gerp_patch = mpatches.Patch(color=SCORE_COLORS[MY_NO_GERP_NAME])
	ccrs_patch = mpatches.Patch(color=SCORE_COLORS[CCRS_NAME])
	patches = OrderedDict()
	patches[MY_NAME] = my_patch
	patches[MY_NO_GERP_NAME] = my_no_gerp_patch
	patches[CCRS_NAME] = ccrs_patch

	spearmanr_patches = OrderedDict()
	scores_colors = []

	for score_name in scores_order:
		score_rs = scores_rs[score_name]
		label = score_name.replace('\n', ' ') + ' (r=%.2f)' % score_rs
		spearmanr_patches[label] = patches[score_name]
		scores_colors.append(SCORE_COLORS[score_name])

	l = ax.legend(spearmanr_patches.values(), spearmanr_patches.keys(), loc='upper right', frameon=False, fontsize=6)
	ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))


def draw_gevir_vs_ccrs(db):
	fig = plt.figure()
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	omim_sets = OmimSets(db)
	ccrs_dict = {}
	all_ccrs_transcripts = []
	ccr_genes = db.gevir.ccr_genes.find({})
	for ccr_gene in ccr_genes:
		ccr_transcript_id = ccr_gene['transcript_id']
		if ccr_transcript_id:
			all_ccrs_transcripts.append(ccr_transcript_id)
			if ccr_gene['ccrs_gte_95'] > 0:
				ccrs_dict[ccr_transcript_id] = ccr_gene['ccrs_gte_95']

	gevir_transcripts = score_sets.gevir_dict.keys()
	gevir_no_gerp_transcripts = score_sets.gevir_no_gerp_dict.keys()
	all_common_transcripts = len(set(gevir_transcripts) & set(all_ccrs_transcripts))
	omim_ads = score_sets.all_transcripts & set(all_ccrs_transcripts) & omim_sets.ad
	omim_ars = score_sets.all_transcripts & set(all_ccrs_transcripts) & omim_sets.ar
	omim_ads_num = len(omim_ads)
	omim_ars_num = len(omim_ars)

	print 'All Common Genes', all_common_transcripts
	print 'OMIM ad', omim_ads_num

	ccrs_dict = sort_dict_by_values(ccrs_dict, reverse=True)

	ccrs_transcripts = ccrs_dict.keys()
	print 'Top Genes limited by CCR>=95', len(set(gevir_transcripts) & set(ccrs_transcripts))
	top_common_genes = len(set(gevir_transcripts) & set(ccrs_transcripts))

	gevir_ad_num = 0
	gevir_ar_num = 0
	gevir_no_gerp_ad_num = 0
	gevir_no_gerp_ar_num = 0
	ccrs_ad_num = 0
	ccrs_ar_num = 0

	ranked_ccrs_transcripts = []
	ranked_gevir_transcripts = []
	ranked_gevir_no_gerp_transcripts = []

	gevir_vs_ccrs = GeVIRvsCCRs()
	for x in range(0, len(ccrs_transcripts)):
		ccr_transcript = ccrs_transcripts[x]
		if ccr_transcript not in gevir_transcripts:
			continue
		gevir_transcript = gevir_transcripts[x]
		gevir_no_gerp_transcript = gevir_no_gerp_transcripts[x]

		if gevir_transcript in omim_ads:
			gevir_ad_num += 1
		if gevir_transcript in omim_ars:
			gevir_ar_num += 1

		if gevir_no_gerp_transcript in omim_ads:
			gevir_no_gerp_ad_num += 1
		if gevir_no_gerp_transcript in omim_ars:
			gevir_no_gerp_ar_num += 1

		if ccr_transcript in omim_ads:
			ccrs_ad_num += 1
		if ccr_transcript in omim_ars:
			ccrs_ar_num += 1

		gevir_vs_ccrs.gevir_ads.append(gevir_ad_num)
		gevir_vs_ccrs.gevir_ars.append(gevir_ar_num)
		gevir_vs_ccrs.gevir_no_gerp_ads.append(gevir_no_gerp_ad_num)
		gevir_vs_ccrs.gevir_no_gerp_ars.append(gevir_no_gerp_ar_num)
		gevir_vs_ccrs.ccrs_ads.append(ccrs_ad_num)
		gevir_vs_ccrs.ccrs_ars.append(ccrs_ar_num)
		gevir_vs_ccrs.gevir_f1s.append(get_precision_recall_f1(gevir_ad_num, gevir_ar_num, omim_ads_num)[2])
		gevir_vs_ccrs.gevir_no_gerp_f1s.append(get_precision_recall_f1(gevir_no_gerp_ad_num, gevir_no_gerp_ar_num, omim_ads_num)[2])
		gevir_vs_ccrs.ccrs_f1s.append(get_precision_recall_f1(ccrs_ad_num, ccrs_ar_num, omim_ads_num)[2])

		ranked_ccrs_transcripts.append(ccr_transcript)
		ranked_gevir_transcripts.append(gevir_transcript)
		ranked_gevir_no_gerp_transcripts.append(gevir_no_gerp_transcript)

	scores_ranked_transcripts = OrderedDict()
	scores_ranked_transcripts[CCRS_NAME] = ranked_ccrs_transcripts
	scores_ranked_transcripts[MY_NAME] = ranked_gevir_transcripts
	scores_ranked_transcripts[MY_NO_GERP_NAME] = ranked_gevir_no_gerp_transcripts

	print 'Top Common Genes', len(set(ranked_ccrs_transcripts) & set(ranked_gevir_transcripts))

	score_decils = OrderedDict()
	score_names = [MY_NAME, MY_NO_GERP_NAME, CCRS_NAME]
	for score_name in score_names:
		score_decils[score_name] = [set(x) for x in np.array_split(scores_ranked_transcripts[score_name], 10)]
		for decile in score_decils[score_name]:
			print len(decile)

	score_bins_7 = OrderedDict()
	for score_name in score_names:
		score_bins_7[score_name] = [set(x) for x in np.array_split(scores_ranked_transcripts[score_name], 7)]
		for score_bin in score_bins_7[score_name]:
			print len(score_bin)

	subplot_num = 1
	draw_gevir_vs_ccrs_subplot(subplot_num, gevir_vs_ccrs.gevir_ads, gevir_vs_ccrs.gevir_no_gerp_ads, gevir_vs_ccrs.ccrs_ads, x_label='Rank (genes)', y_label='Cumulative number of AD genes')
	subplot_num += 1
	draw_gevir_vs_ccrs_subplot(subplot_num, gevir_vs_ccrs.gevir_ars, gevir_vs_ccrs.gevir_no_gerp_ars, gevir_vs_ccrs.ccrs_ars, x_label='Rank (genes)', y_label='Cumulative number of AR genes')
	subplot_num += 1
	draw_gevir_vs_ccrs_subplot(subplot_num, gevir_vs_ccrs.gevir_f1s, gevir_vs_ccrs.gevir_no_gerp_f1s, gevir_vs_ccrs.ccrs_f1s, x_label='Rank (genes)', y_label='Cumulative percentage (%)', report_peak=True, all_genes_num=all_common_transcripts)
	subplot_num += 1
	draw_gevir_vs_ccrs_length_subplot(subplot_num, score_bins_7, score_sets, ax='', scores_ranked_transcripts=scores_ranked_transcripts)

	gevir_patch = mpatches.Patch(color=SCORE_COLORS[MY_NAME])
	gevir_no_gerp_patch = mpatches.Patch(color=SCORE_COLORS[MY_NO_GERP_NAME])
	ccrs_patch = mpatches.Patch(color=COLOR_PALETTE['B'])

	for x in range(0, len(gevir_vs_ccrs.gevir_ads)):
		diff = gevir_vs_ccrs.gevir_ads[x] - gevir_vs_ccrs.ccrs_ads[x]
		if diff >= 15:
			print 'GeVIR - CCR AD count difference reached 15 at', str(x), 'genes'
			break 

	patches = OrderedDict()
	patches[MY_NAME] = gevir_patch
	patches[MY_NO_GERP_NAME] = gevir_no_gerp_patch
	patches['CCRS'] = ccrs_patch
	
	fig.set_size_inches(7, 4.5)
	plt.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.1)
	plt.savefig(FIGURES_FOLDER + 'GeVIR_vs_CCRS.png', format='png', dpi=300)
	plt.close()


#####################################################################################################
### Supplementary Figure 5: Known AD and AR renes enrichement among GeVIR, LOEUF and VIRLoF ranks ###
#####################################################################################################

class GeneScoresFE():
	def __init__(self):
		self.gevir_fe_ad = []
		self.gevir_fe_ar = []
		self.oe_lof_fe_ad = []
		self.oe_lof_fe_ar = []
		self.gevir_and_oe_lof_fe_ad = []
		self.gevir_and_oe_lof_fe_ar = []
		self.gevir_fe_ad_p = []
		self.gevir_fe_ar_p = []
		self.oe_lof_fe_ad_p = []
		self.oe_lof_fe_ar_p = []
		self.gevir_and_oe_lof_fe_ad_p = []
		self.gevir_and_oe_lof_fe_ar_p = []

def read_gene_scores_fe(db):
	gene_scores_fe = GeneScoresFE()
	p_value_threshold = 0.00001

	# Creates web gene scores database if it does not exist
	web_gene_scores_exists = db.gevir.web_gene_scores.find_one({})
	if not web_gene_scores_exists:
		if INCLUDE_GNOMAD_OUTLIERS:
			export_gene_scores(db, enrichment_offset=AD_AR_ENRICHMENT_OFFSET)
		else:
			export_gene_scores(db, enrichment_offset=AD_AR_ENRICHMENT_OFFSET_NO_OUTLIERS)

	gene_scores = db.gevir.web_gene_scores.find({})
	for gene in gene_scores:
		gevir_percentile = gene['gevir_percentile']
		oe_lof_percentile = gene['loeuf_percentile']
		gevir_and_oe_lof_percentile = gene['virlof_percentile']

		# GeVIR
		if gene['gevir_ad_p'] < p_value_threshold:
			gene_scores_fe.gevir_fe_ad_p.append((gevir_percentile, gene['gevir_ad_enrichment']))
		else:
			gene_scores_fe.gevir_fe_ad.append((gevir_percentile, gene['gevir_ad_enrichment']))

		if gene['gevir_ar_p'] < p_value_threshold:
			gene_scores_fe.gevir_fe_ar_p.append((gevir_percentile, gene['gevir_ar_enrichment']))
		else:
			gene_scores_fe.gevir_fe_ar.append((gevir_percentile, gene['gevir_ar_enrichment']))

		# LOEUF (oe LoF)
		if gene['loeuf_ad_p'] < p_value_threshold:
			gene_scores_fe.oe_lof_fe_ad_p.append((oe_lof_percentile, gene['loeuf_ad_enrichment']))
		else:
			gene_scores_fe.oe_lof_fe_ad.append((oe_lof_percentile, gene['loeuf_ad_enrichment']))

		if gene['loeuf_ar_p'] < p_value_threshold:
			gene_scores_fe.oe_lof_fe_ar_p.append((oe_lof_percentile, gene['loeuf_ar_enrichment']))
		else:
			gene_scores_fe.oe_lof_fe_ar.append((oe_lof_percentile, gene['loeuf_ar_enrichment']))

		# GeVIR + LOEUF = VIRLoF 
		if gene['virlof_ad_p'] < p_value_threshold:
			gene_scores_fe.gevir_and_oe_lof_fe_ad_p.append((gevir_and_oe_lof_percentile, gene['virlof_ad_enrichment']))
		else:
			gene_scores_fe.gevir_and_oe_lof_fe_ad.append((gevir_and_oe_lof_percentile, gene['virlof_ad_enrichment']))

		if gene['virlof_ar_p'] < p_value_threshold:
			gene_scores_fe.gevir_and_oe_lof_fe_ar_p.append((gevir_and_oe_lof_percentile, gene['virlof_ar_enrichment']))
		else:
			gene_scores_fe.gevir_and_oe_lof_fe_ar.append((gevir_and_oe_lof_percentile, gene['virlof_ar_enrichment']))

	return gene_scores_fe


def draw_gene_scores_fold_enrcihment_scale_subplot(ax, color_blind_friendly=False):
	width = 2

	f1 = 0.33
	f2 = 0.66
	f3 = 1.5
	f4 = 3
	f5 = 6

	if color_blind_friendly:
		# http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=5
		ax.bar(1, f1, width, color='#2c7bb6')
		ax.bar(1, f2, width, bottom=f1, color='#abd9e9')
		ax.bar(1, f3, width, bottom=f2, color='#ffffbf')
		ax.bar(1, f4, width, bottom=f3, color='#fdae61')
		ax.bar(1, f5, width, bottom=f4, color='#d7191c')
	else:
		ax.bar(1, f1, width, color=C_DARK_GREEN)
		ax.bar(1, f2, width, bottom=f1, color=C_LIGHT_GREEN)
		ax.bar(1, f3, width, bottom=f2, color=C_YELLOW)
		ax.bar(1, f4, width, bottom=f3, color=C_ORANGE)
		ax.bar(1, f5, width, bottom=f4, color=C_RED)

	ax.set_ylim(0,5)
	ax.set_xlim(0,1)
	#ax.plot([0,1], [1,1], '--', color=C_GRAY)

	ax.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=False) # labels along the bottom edge are off
	ax.tick_params(
		axis='y',         
		labelsize=7)

	yticks = [0.33, 0.66, 1, 1.5, 2, 3, 4, 5]
	ax.set_yticks(yticks)
	ax.set_ylabel("Fold enrichment", fontsize=7)
	ax.set_yticklabels(yticks)


def draw_gene_scores_ad_fold_enrcihment(ax, gene_scores_fe):
	gevir_fe_ad_p_x = list(x[0] for x in gene_scores_fe.gevir_fe_ad_p)
	gevir_fe_ad_p_y = list(x[1] for x in gene_scores_fe.gevir_fe_ad_p)

	gevir_fe_ad_x = list(x[0] for x in gene_scores_fe.gevir_fe_ad)
	gevir_fe_ad_y = list(x[1] for x in gene_scores_fe.gevir_fe_ad)


	oe_lof_fe_ad_p_x = list(x[0] for x in gene_scores_fe.oe_lof_fe_ad_p)
	oe_lof_fe_ad_p_y = list(x[1] for x in gene_scores_fe.oe_lof_fe_ad_p)

	oe_lof_fe_ad_x = list(x[0] for x in gene_scores_fe.oe_lof_fe_ad)
	oe_lof_fe_ad_y = list(x[1] for x in gene_scores_fe.oe_lof_fe_ad)


	gevir_and_oe_lof_fe_ad_p_x = list(x[0] for x in gene_scores_fe.gevir_and_oe_lof_fe_ad_p)
	gevir_and_oe_lof_fe_ad_p_y = list(x[1] for x in gene_scores_fe.gevir_and_oe_lof_fe_ad_p)

	gevir_and_oe_lof_fe_ad_x = list(x[0] for x in gene_scores_fe.gevir_and_oe_lof_fe_ad)
	gevir_and_oe_lof_fe_ad_y = list(x[1] for x in gene_scores_fe.gevir_and_oe_lof_fe_ad)

	ax.set_ylim(0,5)
	ax.scatter(gevir_fe_ad_p_x, gevir_fe_ad_p_y, s=1, color=SCORE_COLORS[MY_NAME])
	ax.scatter(gevir_fe_ad_x, gevir_fe_ad_y, s=1, color=SCORE_COLORS[MY_NAME])

	ax.scatter(oe_lof_fe_ad_p_x, oe_lof_fe_ad_p_y, s=1, color=SCORE_COLORS[LOF_OE_NAME])
	ax.scatter(oe_lof_fe_ad_x, oe_lof_fe_ad_y, s=1, color=SCORE_COLORS[LOF_OE_NAME])

	ax.scatter(gevir_and_oe_lof_fe_ad_p_x, gevir_and_oe_lof_fe_ad_p_y, s=1, color=SCORE_COLORS[COMBINED_NAME])
	ax.scatter(gevir_and_oe_lof_fe_ad_x, gevir_and_oe_lof_fe_ad_y, s=1, color=SCORE_COLORS[COMBINED_NAME])

	ax.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		right=False,       # ticks along the right edge are off
		left=False,        # ticks along the left edge are off
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=False,
		direction='in') # labels along the bottom edge are off
	ax.tick_params(
		axis='x',         
		labelsize=7)
	plt.setp(ax.get_yticklabels(), visible=False)

	upper_significance_y = max(gevir_fe_ad_y + oe_lof_fe_ad_y + gevir_and_oe_lof_fe_ad_y)
	lower_significance_y = min(gevir_fe_ad_y + oe_lof_fe_ad_y + gevir_and_oe_lof_fe_ad_y)

	ax.plot([0,100], [1,1], '--', color=C_GRAY)
	ax.plot([0,100], [upper_significance_y,upper_significance_y], '-.', color=COLOR_PALETTE['R'])
	ax.plot([0,100], [lower_significance_y,lower_significance_y], '-.', color=COLOR_PALETTE['R'])

	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	ax.set_title('a', loc='left', fontweight='bold', fontsize=7)
	ax.set_title("Autosomal Dominant (AD)", fontsize=7)
	ax.set_xlabel("Rank (%)", fontsize=7)

	# Legend
	gevir_patch = mpatches.Patch(color=SCORE_COLORS[MY_NAME])
	oe_lof_patch = mpatches.Patch(color=SCORE_COLORS[LOF_OE_NAME])
	gevir_and_oe_lof_patch = mpatches.Patch(color=SCORE_COLORS[COMBINED_NAME])

	patches = OrderedDict()
	patches[MY_NAME] = gevir_patch
	patches[LOF_OE_NAME] = oe_lof_patch
	patches[COMBINED_NAME] = gevir_and_oe_lof_patch	
	patches['Expected'] = Line2D([0], [0], linestyle='--', color=C_GRAY)
	patches['Significance threshold\n(p-value < 1e-5)'] = Line2D([0], [0], linestyle='-.', color=COLOR_PALETTE['R'])
	
	plt.legend(patches.values(), patches.keys(), loc='upper right', frameon=False, ncol=1, fontsize=7)


def draw_gene_scores_ar_fold_enrcihment(ax, gene_scores_fe):
	gevir_fe_ar_p_x = list(x[0] for x in gene_scores_fe.gevir_fe_ar_p)
	gevir_fe_ar_p_y = list(x[1] for x in gene_scores_fe.gevir_fe_ar_p)

	gevir_fe_ar_x = list(x[0] for x in gene_scores_fe.gevir_fe_ar)
	gevir_fe_ar_y = list(x[1] for x in gene_scores_fe.gevir_fe_ar)


	oe_lof_fe_ar_p_x = list(x[0] for x in gene_scores_fe.oe_lof_fe_ar_p)
	oe_lof_fe_ar_p_y = list(x[1] for x in gene_scores_fe.oe_lof_fe_ar_p)

	oe_lof_fe_ar_x = list(x[0] for x in gene_scores_fe.oe_lof_fe_ar)
	oe_lof_fe_ar_y = list(x[1] for x in gene_scores_fe.oe_lof_fe_ar)


	gevir_and_oe_lof_fe_ar_p_x = list(x[0] for x in gene_scores_fe.gevir_and_oe_lof_fe_ar_p)
	gevir_and_oe_lof_fe_ar_p_y = list(x[1] for x in gene_scores_fe.gevir_and_oe_lof_fe_ar_p)

	gevir_and_oe_lof_fe_ar_x = list(x[0] for x in gene_scores_fe.gevir_and_oe_lof_fe_ar)
	gevir_and_oe_lof_fe_ar_y = list(x[1] for x in gene_scores_fe.gevir_and_oe_lof_fe_ar)

	ax.set_ylim(0,5)
	ax.scatter(gevir_fe_ar_p_x, gevir_fe_ar_p_y, s=1, color=SCORE_COLORS[MY_NAME])
	ax.scatter(gevir_fe_ar_x, gevir_fe_ar_y, s=1, color=SCORE_COLORS[MY_NAME])

	ax.scatter(oe_lof_fe_ar_p_x, oe_lof_fe_ar_p_y, s=1, color=SCORE_COLORS[LOF_OE_NAME])
	ax.scatter(oe_lof_fe_ar_x, oe_lof_fe_ar_y, s=1, color=SCORE_COLORS[LOF_OE_NAME])

	ax.scatter(gevir_and_oe_lof_fe_ar_p_x, gevir_and_oe_lof_fe_ar_p_y, s=1, color=SCORE_COLORS[COMBINED_NAME])
	ax.scatter(gevir_and_oe_lof_fe_ar_x, gevir_and_oe_lof_fe_ar_y, s=1, color=SCORE_COLORS[COMBINED_NAME])


	ax.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		right=False,       # ticks along the right edge are off
		left=True,         # ticks along the left edge are off
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=False,
		direction='in') # labels along the bottom edge are off

	ax.tick_params(
		axis='both',         
		labelsize=7)
	plt.setp(ax.get_yticklabels(), visible=False)

	upper_significance_y = max(gevir_fe_ar_y + oe_lof_fe_ar_y + gevir_and_oe_lof_fe_ar_y)
	lower_significance_y = min(gevir_fe_ar_y + oe_lof_fe_ar_y + gevir_and_oe_lof_fe_ar_y)

	ax.plot([0,100], [1,1], '--', color=C_GRAY)
	ax.plot([0,100], [upper_significance_y,upper_significance_y], '-.', color=COLOR_PALETTE['R'])
	ax.plot([0,100], [lower_significance_y,lower_significance_y], '-.', color=COLOR_PALETTE['R'])

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	ax.set_title('b', loc='left', fontweight='bold', fontsize=7)
	ax.set_title("Autosomal Recessive (AR)", fontsize=7)
	ax.set_xlabel("Rank (%)", fontsize=7)


def draw_web_gene_scores(db, color_blind_friendly=False):
	fig = plt.figure()
	gene_scores_fe = read_gene_scores_fe(db)

	gs = gridspec.GridSpec(1, 3, width_ratios=[1, 20, 20])
	gs.update(wspace=0, hspace=0.1)
	ax0 = plt.subplot(gs[0])
	draw_gene_scores_fold_enrcihment_scale_subplot(ax0, color_blind_friendly=color_blind_friendly)

	ax1 = plt.subplot(gs[1], sharey=ax0)
	draw_gene_scores_ad_fold_enrcihment(ax1, gene_scores_fe)

	ax2 = plt.subplot(gs[2], sharey=ax0)
	draw_gene_scores_ar_fold_enrcihment(ax2, gene_scores_fe)

	'''
	# Alternative Legend
	gevir_patch = mpatches.Patch(color=SCORE_COLORS[MY_NAME]) # color=COLOR_PALETTE['R']
	oe_lof_patch = mpatches.Patch(color=SCORE_COLORS[LOF_OE_NAME]) # color=COLOR_PALETTE['G']
	gevir_and_oe_lof_patch = mpatches.Patch(color=SCORE_COLORS[COMBINED_NAME]) # color=COLOR_PALETTE['P']

	patches = OrderedDict()
	patches[MY_NAME] = gevir_patch
	patches[LOF_OE_NAME] = oe_lof_patch
	patches[COMBINED_NAME] = gevir_and_oe_lof_patch	
	patches['Expected'] = Line2D([0], [0], linestyle='--', color=C_GRAY)
	patches['Significance threshold (p-value < 1e-5)'] = Line2D([0], [0], linestyle='-.', color=COLOR_PALETTE['R'])
	
	fig.legend(patches.values(), patches.keys(), 'lower center', frameon=False, ncol=5, fontsize=7)
	'''
	fig.set_size_inches(7, 4)
	plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.125) # bottom=0.16
	plt.savefig(FIGURES_FOLDER + 'WEB_percentiles_ad_ar_fold_enrichemnt.png', format='png', dpi=300)
	#plt.show()
	plt.close()

#######################################################################################################
### Supplementary Figure 3: GeVIR and gnomAD gene constraint length analysis with notched boxplots  ###
#######################################################################################################

def draw_gene_scores_length(db):
	fig = plt.figure()

	# Gene score traqnscripts (sorted)
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	scores_ranked_transcripts = OrderedDict()
	scores_ranked_transcripts[MY_NAME] = score_sets.gevir_dict.keys()
	scores_ranked_transcripts[MY_NO_GERP_NAME] = score_sets.gevir_no_gerp_dict.keys()
	scores_ranked_transcripts[MISS_Z_NAME] = score_sets.gnomad_miss_z_dict.keys()
	scores_ranked_transcripts[MISS_OE_NAME] = score_sets.gnomad_oe_mis_upper_dict.keys()
	scores_ranked_transcripts[LOF_OE_NAME] = score_sets.gnomad_oe_lof_upper_dict.keys()
	scores_ranked_transcripts[COMBINED_NAME] = score_sets.combined_rank_dict.keys()

	# Order of gene scores on the plot
	score_names = [MY_NAME, MY_NO_GERP_NAME, MISS_Z_NAME, MISS_OE_NAME, LOF_OE_NAME, COMBINED_NAME]

	# Split gene score 
	scores_data = OrderedDict()
	for score_name in score_names:
		scores_data[score_name] = [set(x) for x in np.array_split(scores_ranked_transcripts[score_name], 10)]

	print 'Length Analysis'
	normal_length = score_sets.get_transcripts_median_length(score_sets.all_transcripts)

	print 'Median (all genes)', normal_length
	print 'score_name | bin | lower quartile (25) | upper quartile (75)'

	# Calculate correlation
	scores_rs = OrderedDict()
	if scores_ranked_transcripts:
		for score_name, ranked_transcripts in scores_ranked_transcripts.iteritems():
			score_ranks = []
			score_lengths = []
			for x in range(0, len(ranked_transcripts)):
				score_ranks.append(x)
				score_lengths.append(score_sets.length_dict[ranked_transcripts[x]])

			scores_rs[score_name] = spearmanr(score_ranks, score_lengths)[0]

	scores_rs = sort_dict_by_values(scores_rs, reverse=True)


	bar_width, x_paddings = bar_chart_get_bar_width_and_x_paddings(6)
	
	score_n = 0

	score_boxplots = OrderedDict()
	for score_name, bins_transcripts in scores_data.iteritems():
		x = 0
		bins_transcripts = scores_data[score_name]
		for bin_transcripts in bins_transcripts:
			transcript_lengths = score_sets.get_transcripts_length(bin_transcripts)
			
			# Report median and 25/75 quartiles
			y = np.median(transcript_lengths)
			y_err_lower = np.percentile(transcript_lengths, 25)
			y_err_upper = np.percentile(transcript_lengths, 75)
			print score_name, x, y_err_lower, y_err_upper

			boxplot = plt.boxplot(transcript_lengths, positions=[x + x_paddings[score_n]], widths=bar_width*0.8, notch=True, showfliers=False, patch_artist=True)
			
			if x == 0:
				score_boxplots[score_name] = [boxplot]
			else:
				score_boxplots[score_name].append(boxplot)

			x += 1
		score_n += 1

	# Colour boxplots
	for score_name, boxplots in score_boxplots.iteritems():
		for boxplot in boxplots:
			for patch in boxplot['boxes']:
				patch.set_facecolor(SCORE_COLORS[score_name])

			for patch in boxplot['medians']:
				patch.set_color('yellow')

	# Plot expected (median) and set axis ticks and labels
	xs = np.arange(len(scores_data.values()[0]))
	xs = list(xs)
	normal_xs = [-0.8] + xs + [max(xs) + 0.8]
	plt.plot(normal_xs, [normal_length]*len(normal_xs), '--', color=C_GRAY, zorder=10)
	plt.xticks(range(0, 10), fontsize=7)
	plt.yticks(fontsize=7)
	plt.xticks(range(-1, 10, 1), [''] + [str(n) for n in range(10, 110, 10)])
	plt.ylabel('Protein length', fontsize=7)
	plt.xlabel('Rank (%)', fontsize=7)

	# Draw main legend (gene score names)
	my_patch = mpatches.Patch(color=SCORE_COLORS[MY_NAME])
	my_no_gerp_patch = mpatches.Patch(color=SCORE_COLORS[MY_NO_GERP_NAME])
	miss_z_patch = mpatches.Patch(color=SCORE_COLORS[MISS_Z_NAME])
	miss_oe_patch = mpatches.Patch(color=SCORE_COLORS[MISS_OE_NAME])
	lof_oe_patch = mpatches.Patch(color=SCORE_COLORS[LOF_OE_NAME])
	combined_patch = mpatches.Patch(color=SCORE_COLORS[COMBINED_NAME])

	'''
	# Alternative legend; Combined with coloured Spearman correlation ranks
	patches = OrderedDict()
	patches[MY_NAME] = my_patch
	patches[MY_NO_GERP_NAME] = my_no_gerp_patch
	patches[MISS_Z_NAME] = miss_z_patch
	patches[MISS_OE_NAME] = miss_oe_patch
	patches[LOF_OE_NAME] = lof_oe_patch
	patches[COMBINED_NAME] = combined_patch
	patches['Median (all genes)'] = Line2D([0], [0], linestyle='--', color=C_GRAY)
	
	fig.legend(patches.values(), patches.keys(), 'lower center', frameon=False, ncol=7, fontsize=7)
	'''
	# Draw secondary legend (gene score spearman correlations)
	#spearmanr_patches = OrderedDict()
	#scores_colors = []

	patches = OrderedDict()
	for score_name, score_rs in scores_rs.iteritems():
		patch = mpatches.Patch(color=SCORE_COLORS[score_name])
		label = score_name + ' (r=%.2f)' % score_rs
		patches[label] = patch
		'''
		label = 'Spearman r=%.2f' % score_rs
		spearmanr_patches[label] = patches[score_name]
		scores_colors.append(SCORE_COLORS[score_name])
		'''
	patches['Median (all genes)'] = Line2D([0], [0], linestyle='--', color=C_GRAY)

	l = fig.legend(patches.values(), patches.keys(), loc=(0.76,0.7), frameon=False, fontsize=7)

	'''
	# Alternative legend style, colour subplot legends
	legend_auc_texts = l.get_texts()
	legend_auc_lines = l.get_patches()

	for x in range(0, len(legend_auc_texts)):
		legend_auc_texts[x].set_color(scores_colors[x])
		legend_auc_lines[x].set_alpha(0)
	'''

	ax = fig.axes[0]
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

	# Save resulted figure
	fig.set_size_inches(7, 4)
	plt.tight_layout(rect=[0, 0, 1, 1])
	plt.savefig(FIGURES_FOLDER + 'gene_lengths.png', format='png', dpi=300)
	plt.close()

###################################################################################################
### Supplementary Figure 3: Correlation between gene length and known pathogenic variant types  ###
###################################################################################################

def analyse_clin_var_and_gene_length(db):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	all_set = score_sets.all_transcripts

	omim_sets = OmimSets(db)
	ad_set = set(omim_sets.ad & all_set)
	ar_set = set(omim_sets.ar & all_set)
	transcript_ids = ad_set | ar_set

	lof_only_lengths = []
	miss_only_lengths = []
	lof_and_miss_lengths = []

	all_lengths = []
	for transcript_id in all_set:
		ens_protein = db.gevir.ens_aa_fasta.find_one({'_id': transcript_id})
		length = len(ens_protein['cds']) - 1 # do not count stop codon
		all_lengths.append(length)

	for transcript_id in transcript_ids:
		ens_protein = db.gevir.ens_aa_fasta.find_one({'_id': transcript_id})
		length = len(ens_protein['cds']) - 1 # do not count stop codon

		clinvar_gene = db.gevir.clin_var_genes.find_one({'canonical_transcript': transcript_id})
		if clinvar_gene:
			lof_num = clinvar_gene['pathogenic']['lof']
			miss_num = clinvar_gene['pathogenic']['miss_and_indels']

			if lof_num > 0 and miss_num == 0:
				lof_only_lengths.append(length)
			elif lof_num == 0 and miss_num > 0:
				miss_only_lengths.append(length)
			elif lof_num > 0 and miss_num > 0:
				lof_and_miss_lengths.append(length)

	print 'AD or AR genes:', len(transcript_ids)
	print 'AD or AR genes with Pathogenic variants:', len(lof_only_lengths + miss_only_lengths + lof_and_miss_lengths)
	bins = np.linspace(0, 9000, 90)

	fig = plt.figure()
	plt.hist(lof_only_lengths, bins=bins, normed=True, alpha=0.5, label='Only LoF (N=' + str(len(lof_only_lengths)) + ')')
	plt.hist(lof_and_miss_lengths, bins=bins, normed=True, alpha=0.5, label='LoF and Missense/INDELs (N=' + str(len(lof_and_miss_lengths)) + ')')
	plt.hist(miss_only_lengths, bins=bins, normed=True, alpha=0.5, label='Only Missense/INDELs Only (N=' + str(len(miss_only_lengths)) + ')')
	plt.yticks(fontsize=14)
	plt.xticks(np.arange(0, 9000, 500), fontsize=14)
	plt.legend(loc='upper right', frameon=False, fontsize=14)
	plt.xlabel('Gene Length (amino acids)', fontsize=14)
	plt.ylabel('AD or AR Genes (normalised)', fontsize=14)
	
	fig.set_size_inches(14, 10)
	plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99])
	plt.savefig(FIGURES_FOLDER + 'clin_var_and_gene_length.png', format='png', dpi=300)
	plt.close()

	lof_only_vs_miss_only_stats, lof_only_vs_miss_only_p_value = stats.ttest_ind(lof_only_lengths, miss_only_lengths, equal_var=False)
	lof_and_miss_vs_miss_only_stats, lof_and_miss_vs_miss_only_p_value = stats.ttest_ind(lof_and_miss_lengths, miss_only_lengths, equal_var=False)
	miss_only_vs_miss_only_stats, miss_only_vs_miss_only_p_value = stats.ttest_ind(miss_only_lengths, miss_only_lengths, equal_var=False)
	all_vs_miss_only_stats, all_vs_miss_only_p_value = stats.ttest_ind(all_lengths, miss_only_lengths, equal_var=False)

	headers = ['Category', 'Gene Number', 'Mean Length', 'Median Length', "Length comparison with Miss/INDELs (Welch's t-test)"]
	table_data = [headers, 
				['Only LoF', len(lof_only_lengths), '{:.1f}'.format(np.mean(lof_only_lengths)), np.median(lof_only_lengths), "{:.2E}".format(lof_only_vs_miss_only_p_value)],
				['LoF and Missense/INDELs', len(lof_and_miss_lengths), '{:.1f}'.format(np.mean(lof_and_miss_lengths)), np.median(lof_and_miss_lengths), "{:.2E}".format(lof_and_miss_vs_miss_only_p_value)],
				['Only Missense/INDELs', len(miss_only_lengths), '{:.1f}'.format(np.mean(miss_only_lengths)), np.median(miss_only_lengths), "{:.2E}".format(miss_only_vs_miss_only_p_value)],
				['All', len(all_lengths), '{:.1f}'.format(np.mean(all_lengths)), np.median(all_lengths), "{:.2E}".format(all_vs_miss_only_p_value)],
	]

	table = AsciiTable(table_data)
	print table.table

	output_csv = OUTPUT_FOLDER + 'clinvar_and_gene_length.csv'
	write_table_to_csv(table_data, output_csv)


#############
### DAVID ###
#############

def export_gene_list_for_david(db, gene_list='GeVIR'):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	print gene_list, 'All genes:', len(score_sets.all_transcripts)
	highly_intolerant_genes = db.gevir.common_gene_scores.find({ "is_gnomad": True, "gnomad_oe_lof_upper": { "$lt": 0.35 } })

	top_num = highly_intolerant_genes.count()
	gevir_set = set(score_sets.gevir_dict.keys()[:top_num])
	oe_lof_set = set(score_sets.gnomad_oe_lof_upper_dict.keys()[:top_num])

	if gene_list == 'GeVIR':
		result_set = gevir_set - oe_lof_set
	elif gene_list == 'LOEUF':
		result_set = oe_lof_set - gevir_set
	else:
		result_set = oe_lof_set & gevir_set

	headers = ['gene_id', 'transcript_id', 'gene_name']
	table = [headers]
	for transcript_id in result_set:
		gene = db.exac.genes.find_one({'canonical_transcript': transcript_id})
		row = [gene['gene_id'], transcript_id, gene['gene_name']]
		table.append(row)

	if gene_list == 'GeVIR':
		report_name = 'david_gene_list_gevir.csv'
	elif gene_list == 'LOEUF':
		report_name = 'david_gene_list_loeuf.csv'
	else:
		report_name = 'david_gene_list_both.csv'

	output_csv = OUTPUT_FOLDER + report_name
	write_table_to_csv(table, output_csv)

	# Report Statistical enrichemnts:
	essential_sets = EssentialSets(db)

	mouse_het_lethal_set = set(score_sets.all_transcripts) & essential_sets.mouse_het_lethal
	crispr_essential_set = set(score_sets.all_transcripts) & essential_sets.crispr_essential

	transcripts_all = len(score_sets.all_transcripts)
	mouse_het_lethal_all = len(mouse_het_lethal_set)
	crispr_essential_all = len(crispr_essential_set)

	top_gevir = len(result_set)
	top_gevir_mouse_het_lethal = len(result_set & mouse_het_lethal_set)
	top_gevir_crispr_essential = len(result_set & crispr_essential_set)

	mouse_het_lethal_table = [[top_gevir_mouse_het_lethal, top_gevir],[mouse_het_lethal_all, transcripts_all]]
	cell_essential_table = [[top_gevir_crispr_essential, top_gevir],[crispr_essential_all, transcripts_all]]
	mouse_het_lethal_fold_enrichemnt, mouse_het_lethal_p = fisher_exact(mouse_het_lethal_table)
	crispr_essential_fold_enrichemnt, crispr_essential_p = fisher_exact(cell_essential_table)

	print 'Mouse Het Lethal', mouse_het_lethal_fold_enrichemnt, '(fold-enrichment)', float_to_sci_str(mouse_het_lethal_p), '(p-value)'
	print report_2_x_2_test_table_as_obs_exp_str(mouse_het_lethal_table)
	print 'Cell Essential', crispr_essential_fold_enrichemnt, '(fold-enrichment)', float_to_sci_str(crispr_essential_p), '(p-value)'
	print report_2_x_2_test_table_as_obs_exp_str(cell_essential_table)


def get_gene_id_to_length_dict(db):
	genes = db.exac.genes.find({})
	gene_id_to_length = {}
	for gene in genes:
		if 'canonical_transcript' not in gene:
			continue
		transcript_id = gene['canonical_transcript']
		gene_protein_sequence = db.gevir.ens_aa_fasta.find_one({'_id': transcript_id})
		if not gene_protein_sequence:
			continue

		length = len(gene_protein_sequence['cds']) - 1
		gene_id_to_length[gene['gene_id']] = length
	return gene_id_to_length


def merge_david_reports_and_add_gene_length(db):
	DAVID_FOLDER = '/home/niab/cs/phd/phd_projects/gevir/tables/david_gevir_vs_loeuf/david_original/'
	OUTPUT_DAVID_FOLDER = '/home/niab/cs/phd/phd_projects/gevir/tables/david_gevir_vs_loeuf/'

	GEVIR_DAVID_REPORTS = [
		'gevir_bp_all.tsv',
		'gevir_cc_all.tsv',
		'gevir_kegg.tsv',
		'gevir_mf_all.tsv',
	]

	LOEUF_DAVID_REPORTS = [
		'loeuf_bp_all.tsv',
		'loeuf_cc_all.tsv',
		'loeuf_kegg.tsv',
		'loeuf_mf_all.tsv',
	]

	GEVIR_AND_LOEUF_DAVID_REPORTS = [
		'both_bp_all.tsv',
		'both_cc_all.tsv',
		'both_kegg.tsv',
		'both_mf_all.tsv',
	]

	gene_id_to_length = get_gene_id_to_length_dict(db)

	table = []
	for report_name in GEVIR_DAVID_REPORTS + LOEUF_DAVID_REPORTS + GEVIR_AND_LOEUF_DAVID_REPORTS:
		print report_name
		file_name = DAVID_FOLDER + report_name
		report = CsvReader(file_name, delimiter='\t')
		if len(table) == 0:
			headers = ['gene_list'] + report.headers + ['mean_length', 'median_length', 'standard_deviation']
			table.append(headers)

		for document in report.data:
			if report_name in GEVIR_DAVID_REPORTS:
				gene_list = 'Only GeVIR'
			elif report_name in LOEUF_DAVID_REPORTS:
				gene_list = 'Only LOEUF'
			else:
				gene_list = 'GeVIR & LOEUF'

			gene_ids = document['Genes']
			gene_ids = gene_ids.split(', ')
			gene_lengths = []
			for gene_id in gene_ids:
				gene_lengths.append(gene_id_to_length[gene_id])

			mean_length = np.mean(gene_lengths)
			median_length = np.median(gene_lengths)
			std_length = np.std(gene_lengths)
			row = [gene_list] + document.values() + [mean_length, median_length, std_length]
			table.append(row)

	output_csv = OUTPUT_DAVID_FOLDER + 'david_gevir_vs_loeuf.csv'
	write_table_to_csv(table, output_csv)


def get_group_ad_ar_enrichment(group_name, group_set, ad_set, ar_set, all_set):
	ad_num = len(ad_set)
	ar_num = len(ar_set)
	all_num = len(all_set)
	group_num = len(group_set)
	group_ad_num = len(group_set & ad_set)
	group_ar_num = len(group_set & ar_set)
	group_ad_proportion = proportion_to_percents_str(float(group_ad_num) / group_num)
	group_ar_proportion = proportion_to_percents_str(float(group_ar_num) / group_num)
	group_ad_fold_enrichemnt, group_ad_p_value = fisher_exact([[group_ad_num, group_num], [ad_num, all_num]])
	group_ar_fold_enrichemnt, group_ar_p_value = fisher_exact([[group_ar_num, group_num], [ar_num, all_num]])

	return [group_name,
			group_num,
			'{} ({}%)'.format(group_ad_num, group_ad_proportion),
			'{:.2f}'.format(group_ad_fold_enrichemnt),
			float_to_sci_str(group_ad_p_value), 
			
			'{} ({}%)'.format(group_ar_num, group_ar_proportion),
			'{:.2f}'.format(group_ar_fold_enrichemnt),
			float_to_sci_str(group_ar_p_value), 
			]

#############################################################
### Supplementary Figure 1 Shows VIRs GERP++ distribution ###
#############################################################

def draw_region_gerp_weights_hist(db, clean=False):
	
	if clean:
		transcript_ids = []
		transcript_ids_xy = []

		if INCLUDE_GNOMAD_OUTLIERS:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True }) # "no_issues": True, 
		else:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

		for gnomad_gene in gnomad_genes:
			if gnomad_gene['chrom'] == 'X' or gnomad_gene['chrom'] == 'Y':
				transcript_ids_xy.append(gnomad_gene['_id'])
			else:
				transcript_ids.append(gnomad_gene['_id'])

		zero_length_regions = 0
		gerps = []
		genes_gerps = []
		db.gevir.temp_gene_region_gerps.drop()

		total_lines = len(transcript_ids)
		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for transcript_id in transcript_ids:
			gene_gerps = []
			regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_AUTOSOMAL_COVERAGE }, 'not_in_cds': False, "lenght": { "$gte": 1 } })
			for region in regions:
				gerps.append(region['gerp_mean'])
				gene_gerps.append(region['gerp_mean'])
			genes_gerps.append({'_id': transcript_id, 'region_gerps': gene_gerps })

			regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_AUTOSOMAL_COVERAGE }, 'not_in_cds': False, "lenght": 0 })
			zero_length_regions += regions.count()
			line_number += 1
			bar.update((line_number + 0.0) / total_lines)
		bar.finish()

		total_lines = len(transcript_ids_xy)
		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()	
		for transcript_id in transcript_ids_xy:
			gene_gerps = []
			regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_XY_COVERAGE }, 'not_in_cds': False, "lenght": { "$gte": 1 } })
			for region in regions:
				gerps.append(region['gerp_mean'])
				gene_gerps.append(region['gerp_mean'])
			genes_gerps.append({'_id': transcript_id, 'region_gerps': gene_gerps })

			regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_XY_COVERAGE }, 'not_in_cds': False, "lenght": 0 })
			zero_length_regions += regions.count()
			line_number += 1
			bar.update((line_number + 0.0) / total_lines)
		bar.finish()

		db.gevir.temp_gene_region_gerps.insert_many(genes_gerps)
		print 'Zero length regions:', zero_length_regions

	gerps = []
	genes_gerps = db.gevir.temp_gene_region_gerps.find({})
	for gene_gerps in genes_gerps:
		gerps += gene_gerps['region_gerps']

	print len(gerps)
	print 'Regions mean GERP++ CI 99%'
	print np.percentile(gerps, 1)
	print np.percentile(gerps, 99)
	print min(gerps)

	bins = np.linspace(-13, 7, 100)
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
	plt.hist(gerps, bins=bins, color=COLOR_PALETTE['B'])
	plt.yticks(fontsize=7)
	plt.xticks(np.arange(-13, 7, 1), fontsize=7)
	plt.xlabel('Mean GERP++', fontsize=7)
	plt.ylabel('Varinat Intolerant Regions (VIRs)', fontsize=7)
	fig.set_size_inches(7, 4)
	plt.tight_layout(rect=[0, 0, 1, 1])
	plt.savefig(FIGURES_FOLDER + 'region_gerps.png', format='png', dpi=300)
	plt.close()	

###################
### EXTRA STATS ###
###################

def compare_essential_and_omim_genes(db):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	essential_sets = EssentialSets(db)
	omim_sets = OmimSets(db)
	all_set = score_sets.all_transcripts
	ad_set = set(omim_sets.ad & all_set)
	ar_set = set(omim_sets.ar & all_set)

	mouse_het_lethal_set = essential_sets.mouse_het_lethal & all_set
	cell_essential_set = essential_sets.crispr_essential & all_set
	cell_non_essential = essential_sets.crispr_non_essential & all_set
	null_set = essential_sets.nulls & all_set

	headers = ['Group Name', 'Genes', 'AD Genes', 'AD Fold Enrichment', 'AD p-value', 'AR Genes', 'AR Fold Enrichment', 'AR p-value']

	mouse_het_lethal = get_group_ad_ar_enrichment('Mouse het lethal', mouse_het_lethal_set, ad_set, ar_set, all_set)
	cell_essential = get_group_ad_ar_enrichment('Cell essential', cell_essential_set, ad_set, ar_set, all_set)
	cell_non_essential = get_group_ad_ar_enrichment('Cell non-essential', cell_non_essential, ad_set, ar_set, all_set)
	nulls = get_group_ad_ar_enrichment('Null', null_set, ad_set, ar_set, all_set)

	table = [headers, mouse_het_lethal, cell_essential, cell_non_essential, nulls]
	output_csv = OUTPUT_FOLDER + 'essential_and_omim_genes.csv'
	write_table_to_csv(table, output_csv)


def report_number_of_well_covered_regions(db):
	print 'Number of regions with high mean coverage.'
	transcript_ids = []
	transcript_ids_xy = []

	if INCLUDE_GNOMAD_OUTLIERS:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True }) # "no_issues": True, 
	else:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	for gnomad_gene in gnomad_genes:
		if gnomad_gene['chrom'] == 'X' or gnomad_gene['chrom'] == 'Y':
			transcript_ids_xy.append(gnomad_gene['_id'])
		else:
			transcript_ids.append(gnomad_gene['_id'])

	ex_min = 8
	ex_max = 12
	ex_regions_length_small = 0
	ex_regions_length_21_plus = 0

	autosomal_regions_all = 0
	autosomal_regions_pass = 0
	autosomal_regions_pass_0 = 0

	for transcript_id in transcript_ids:
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, 'not_in_cds': False})
		autosomal_regions_all += regions.count()
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_AUTOSOMAL_COVERAGE }, 'not_in_cds': False})
		autosomal_regions_pass += regions.count()
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_AUTOSOMAL_COVERAGE }, 'not_in_cds': False, "lenght": 0 })
		autosomal_regions_pass_0 += regions.count()

		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_AUTOSOMAL_COVERAGE }, 'not_in_cds': False, "$and": [ { "lenght": { "$gte": ex_min } }, { "lenght": { "$lte": ex_max } } ] })
		ex_regions_length_small += regions.count()
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_AUTOSOMAL_COVERAGE }, 'not_in_cds': False, "lenght": { "$gt": 20 } })
		ex_regions_length_21_plus += regions.count()

	print 'All autosomal regions:', autosomal_regions_all
	print 'Autosomal regions with coverage >=' + str(MIN_AUTOSOMAL_COVERAGE), autosomal_regions_pass
	print 'Autosomal proportion:', autosomal_regions_pass / float(autosomal_regions_all)

	xy_regions_all = 0
	xy_regions_pass = 0
	xy_regions_pass_0 = 0
	for transcript_id in transcript_ids_xy:
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, 'not_in_cds': False})
		xy_regions_all += regions.count()
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_XY_COVERAGE }, 'not_in_cds': False})
		xy_regions_pass += regions.count()
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_XY_COVERAGE }, 'not_in_cds': False, "lenght": 0 })
		xy_regions_pass_0 += regions.count()

		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_XY_COVERAGE }, 'not_in_cds': False, "$and": [ { "lenght": { "$gte": ex_min } }, { "lenght": { "$lte": ex_max } } ] })
		ex_regions_length_small += regions.count()
		regions = db.gevir[REGIONS_COLLECTION].find({'transcript_id': transcript_id, "exome_coverage": { "$gte": MIN_XY_COVERAGE }, 'not_in_cds': False, "lenght": { "$gt": 20 } })
		ex_regions_length_21_plus += regions.count()

	print 'All XY regions:', xy_regions_all	
	print 'XY regions with coverage >=' + str(MIN_XY_COVERAGE), xy_regions_pass
	print 'XY proportion:', xy_regions_pass / float(xy_regions_all)

	print 'Length small regions (' + str(ex_min) + '-' + str(ex_max) + '):', ex_regions_length_small
	print 'Length >20 regions:', ex_regions_length_21_plus

	print '0 Length regions:', autosomal_regions_pass_0 + xy_regions_pass_0


def get_omim_known_transcript_ids(db):
	known_transcript_ids = set([])
	omim_genes = CsvReader(OMIM_TSV, delimiter='\t')
	for document in omim_genes.data:
		gene_name_38 = document['Approved Symbol'].upper()
		gene_id_38 = document['Ensembl Gene ID']
		chrom = document['Chromosome'][3:]
		phenotypes = document['Phenotypes']
		exac_gene_by_id = db.exac.genes.find_one({'gene_id': gene_id_38, 'chrom': chrom})
		exac_gene_by_name = db.exac.genes.find_one({'gene_name_upper': gene_name_38, 'chrom': chrom})
		transcript_id = ''
		if exac_gene_by_id and 'canonical_transcript' in exac_gene_by_id:
			transcript_id = exac_gene_by_id['canonical_transcript']
		elif exac_gene_by_name and 'canonical_transcript' in exac_gene_by_name:
			transcript_id = exac_gene_by_name['canonical_transcript']

		if transcript_id and phenotypes:
			known_transcript_ids.add(transcript_id)
	return known_transcript_ids


def report_top_important_genes(db, top_set=None):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	essential_sets = EssentialSets(db)
	all_set = score_sets.all_transcripts

	print len(score_sets.all_transcripts)
	mouse_het_lethal_set = essential_sets.mouse_het_lethal & all_set

	known_transcripts = get_omim_known_transcript_ids(db)
	omim_sets = OmimSets(db)
	ad_set = set(omim_sets.ad & all_set)
	ar_set = set(omim_sets.ar & all_set)
	
	if not top_set:
		top_set = set([])
		genes = db.gevir.web_gene_scores.find({ "virlof_percentile": { "$lte": 10.0 } })
		for gene in genes:
			top_set.add(gene['canonical_transcript'])

		print 'Top 10% genes', len(top_set)
	else:
		print 'Top genes: {}'.format(len(top_set))

	print 'AD', len(top_set & ad_set), len(ad_set), '{:.2f}'.format(float(len(top_set & ad_set)) * 100 / len(ad_set))
	print 'AR', len(top_set & ar_set), len(ar_set), '{:.2f}'.format(float(len(top_set & ar_set)) * 100 / len(ar_set))


	top_unknown_genes_set = top_set - known_transcripts
	oth_unknown_genes_set = all_set - top_set - known_transcripts

	top_unknown_genes = len(top_unknown_genes_set)
	oth_unknown_genes = len(oth_unknown_genes_set)

	ad_fold_enrichemnt, ad_p_value = fisher_exact([[len(top_set & ad_set), len(top_set & ar_set)], [len(ad_set), len(ar_set)]])
	print 'AD Fold enrichemnt', '{:.2f}'.format(ad_fold_enrichemnt)
	print 'AD p-value', "{:.2E}".format(ad_p_value)
	print 'Top unknown genes', top_unknown_genes, '{:.2f}'.format(float(top_unknown_genes) * 100 / len(top_set))
	print 'Other unknown genes', oth_unknown_genes

	print 'Mouse Het Lethal in top unknown', len(mouse_het_lethal_set & top_unknown_genes_set), top_unknown_genes
	print 'Mouse Het Lethal in oth unknown', len(mouse_het_lethal_set & oth_unknown_genes_set), oth_unknown_genes
	mouse_het_lethal_fold_enrichemnt, mouse_het_lethal_p_value = fisher_exact([[len(mouse_het_lethal_set & top_unknown_genes_set), top_unknown_genes], [len(mouse_het_lethal_set & oth_unknown_genes_set), oth_unknown_genes]])

	print 'Mouse Het Lethal Fold-enrichment', mouse_het_lethal_fold_enrichemnt
	print 'Mouse Het Lethal p-value', mouse_het_lethal_p_value


def report_gene_constraint_metrics_percentiles(db, transcript_id):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	all_genes = float(len(score_sets.all_transcripts))
	gevir_list = score_sets.gevir_dict.keys()
	miss_z_list = score_sets.gnomad_miss_z_dict.keys()
	oe_mis_upper_list = score_sets.gnomad_oe_mis_upper_dict.keys()
	oe_lof_upper_lsit = score_sets.gnomad_oe_lof_upper_dict.keys()

	gevir_rank = '{:.2f}'.format((gevir_list.index(transcript_id) * 100 / all_genes))
	miss_z_rank = '{:.2f}'.format((miss_z_list.index(transcript_id) * 100 / all_genes))
	oe_mis_upper_rank = '{:.2f}'.format((oe_mis_upper_list.index(transcript_id) * 100 / all_genes))
	oe_lof_upper_rank = '{:.2f}'.format((oe_lof_upper_lsit.index(transcript_id) * 100 / all_genes))
	print 'GeVIR rank', gevir_rank
	print 'Miss z rank', miss_z_rank
	print 'MOEUF rank', oe_mis_upper_rank
	print 'LOEUF rank', oe_lof_upper_rank


def get_transcript_gerp_scores(db, transcript_id):
	gerps = OrderedDict()
	exons = db.exac.exons.find({ "transcript_id": transcript_id, "feature_type": "CDS" }).sort([("xstart", pymongo.ASCENDING)])
	exon_gerps = []
	for exon in exons:
		start = exon['start']
		stop = exon['stop']	
		poses = db.gerp[exon['chrom']].find({ '$and': [ { "_id": { '$gte': start } }, { "_id": { '$lte': stop } } ] })
		for pos in poses:
			exon_gerps.append(pos['RS'])

	print 'MEAN GERP++', np.mean(exon_gerps)


def check_gevir_gerp_impact(db):
	score_sets = ScoreSets(db, filters={"is_gnomad": True })
	essential_sets = EssentialSets(db)
	all_set = score_sets.all_transcripts
	omim_sets = OmimSets(db)
	ad_set = set(omim_sets.ad & all_set)
	ar_set = set(omim_sets.ar & all_set)
	cell_non_essential_set = essential_sets.crispr_non_essential

	genes = db.gevir.common_gene_scores.find({})
	top_gevir = set([])
	shifted_genes = set([])
	for gene in genes:
		transcript_id = gene['_id']
		gevir = gene['gevir_percentile']
		gevir_no_gerp = gene['gevir_no_gerp_percentile']
		if gevir < 30:
			top_gevir.add(transcript_id)
		if gevir_no_gerp < 30 and gevir > 70:
			shifted_genes.add(transcript_id)
	print len(shifted_genes), len(shifted_genes & ad_set), len(shifted_genes & ar_set)
	print len(top_gevir), len(ad_set & top_gevir)
	print len(cell_non_essential_set & shifted_genes)


# Recreates common_gene_scores and supports original order of items
def recreate_common_gene_scores(db):
	db.gevir.common_gene_scores.drop()
	gene_scores = db.gevir.common_gene_scores_19361.find({})

	for gene_score in gene_scores:
		db.gevir.common_gene_scores.insert(gene_score)


######################################################################
### Extra Figure for GeVIR webite: colourful fold enrichment scale ###
######################################################################

def draw_gene_scores_fold_enrcihment_scale_for_web():
	fig = plt.figure()
	width = 2

	f1 = 5
	f2 = 8
	f3 = 11
	f4 = 14
	f5 = 19
	ax = plt.subplot(111)

	ax.barh(1, f1, width, color=C_DARK_GREEN)
	ax.barh(1, f2, width, left=f1, color=C_LIGHT_GREEN)
	ax.barh(1, f3, width, left=f2, color=C_YELLOW)
	ax.barh(1, f4, width, left=f3, color=C_ORANGE)
	ax.barh(1, f5, width, left=f4, color=C_RED)

	ax.set_ylim(0,0.5)
	ax.set_xlim(0,19)

	ax.tick_params(
		axis='x',
		labelsize=14)
	ax.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=False) # labels along the bottom edge are off


	xticks = [5, 8, 9.5, 11, 14]
	xtick_labels = ['0.33', '0.66', '1', '1.5', '3']
	ax.set_xticks(xticks)
	ax.set_xticklabels(xtick_labels)

	fig.set_size_inches(8, 1.3) # 5.5
	plt.yticks([])
	plt.title('Fold Enrichment Color Codes', fontsize=16)
	plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99])
	plt.savefig(FIGURES_FOLDER + 'genes_fe_barchart.png', format='png', dpi=150)
	plt.close()


def main():
	# Fix somehow to be able to run all figure methods one by one without corrupting
	db = MongoDB()

	# Uncomment the functions to produce the figures/tables
	# Some functions requires some store temporary data in the database to have "clean" attribute, which can be set to false
	# after the first run, to save the 

	### Figures ###

	# Figure 1a
	draw_example_gene_variants_subplot(db, 0, single_figure=True)

	# Figure 2, Supplementary Table 1
	analyse_clin_var(db, clean_temp_data=False, include_example_genes_subplot=False)

	# Figure 3, Supplementary Table 3
	draw_gaps_vs_gnomad_constraint_scores(db, clean=True, for_letter=True)

	# Supplementary Figures 2 (contains results for Table 1) and 4 
	draw_gaps_vs_gnomad_constraint_scores(db, clean=True)

	# Supplementary Figure 1 Shows VIRs GERP++ distribution
	draw_region_gerp_weights_hist(db, clean=True)

	# Supplementary Figure 3 (Gene lengths boxplots)
	draw_gene_scores_length(db)

	# Supplementary Figures pLI [Figure not included, but stats are discussed in the supplementary note]
	#draw_gaps_vs_lof_oe_vs_pli(db)

	# Figure 4
	draw_gaps_vs_lof_oe(db)

	# Extended Data Figure 1
	draw_gevir_vs_ccrs(db)

	# Supplementary Figure 5
	draw_web_gene_scores(db, color_blind_friendly=True) 

	# Figure for website (fold enrichment scale)
	#draw_gene_scores_fold_enrcihment_scale_for_web() 

	### Supplementary Tables ###

	# Reports enrichment of AD and AR genes in cell essential, non-essential and Null datasets
	# Creates Supplementary Table 4
	compare_essential_and_omim_genes(db)

	# Exports ~15% top most intolerant genes for DAVID functional enrichment analysis
	export_gene_list_for_david(db, gene_list='GeVIR')
	export_gene_list_for_david(db, gene_list='LOEUF')
	export_gene_list_for_david(db, gene_list='BOTH')

	# Requires manually compaund reports exported from DAVID online tool:
	# https://david.ncifcrf.gov/
	# Supplementary Tables 5 (manually picked) and 6
	# merge_david_reports_and_add_gene_length(db)

	### Report stats in the terminal ###

	# Used for testing purposes and to report percentiles of LITAF gene (ENST00000571688)
	#report_gene_constraint_metrics_percentiles(db , 'ENST00000571688')

	# Reports VIRLoF 10% most intolerant gene stats for discussion
	#report_top_important_genes(db)

	# Reports well covered region stats for methods
	#report_number_of_well_covered_regions(db)

	# Report genes which ranks shifted significantly without GERP++
	#check_gevir_gerp_impact(db)

	### Not used in the original manuscript ###
	#count_gnomad_miss_variants_in_omim_genes(db)
	#calculate_gnomad_pathogenic_missense_variants_stats(db)
	#analyse_clin_var_and_gene_length(db)
	#draw_web_gene_scores(db, color_blind_friendly=False) 


if __name__ == "__main__":
	sys.exit(main())

'''
# DAVID; stats from export_gene_list_for_david
------------------------------------------------------------------
GeVIR All genes: 19361
Mouse Het Lethal 1.97021503104 (fold-enrichment) 3.13E-5 (p-value)
Observed 52 out of 1317 | Expected 388 out of 19361
Cell Essential 2.21732054775 (fold-enrichment) 5.07E-11 (p-value)
Observed 100 out of 1317 | Expected 663 out of 19361

LOEUF All genes: 19361
Mouse Het Lethal 2.23543628522 (fold-enrichment) 2.82E-7 (p-value)
Observed 59 out of 1317 | Expected 388 out of 19361
Cell Essential 1.37473873961 (fold-enrichment) 2.52E-2 (p-value)
Observed 62 out of 1317 | Expected 663 out of 19361


BOTH All genes: 19361
Mouse Het Lethal 3.61114690722 (fold-enrichment) 7.42E-27 (p-value)
Observed 121 out of 1672 | Expected 388 out of 19361
Cell Essential 2.82939119704 (fold-enrichment) 3.07E-25 (p-value)
Observed 162 out of 1672 | Expected 663 out of 19361
'''