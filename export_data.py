import sys
import csv
import pymongo
import math
import progressbar
import numpy as np
from collections import OrderedDict
from scipy.stats import fisher_exact, rankdata
from common import MongoDB, calculate_percentiles, sort_dict_by_values, write_table_to_csv
from gevir import get_transcript_amino_acids_xpositions
from gnomad_utils import xpos_to_pos
from gene_sets import OmimSets, EssentialSets

OUTPUT_FOLDER = './tables/'

MIN_AUTOSOMAL_COVERAGE = 50
MIN_XY_COVERAGE = 45

# 5% offset = 18,352 / 20 ~= 918
# 5% offset = 19,361 / 20 ~= 968
AD_AR_ENRICHMENT_OFFSET = 968
AD_AR_ENRICHMENT_OFFSET_NO_OUTLIERS = 918

def export_gene_identifiers(db):
	transcript_ids = []
	gevir_genes = db.gevir.gevir_scores.find({})

	for gevir_gene in gevir_genes:
		transcript_ids.append(gevir_gene['_id'])

	headers = ['gene_identifier_upper', 'canonical_transcript']
	table = [headers]
	for transcript_id in transcript_ids:
		exac_gene = db.exac.genes.find_one({'canonical_transcript': transcript_id})
		table.append([transcript_id, transcript_id])
		table.append([exac_gene['gene_id'], transcript_id])
		table.append([exac_gene['gene_name'].upper(), transcript_id])

		if 'other_names' in exac_gene:
			for gene_name in exac_gene['other_names']:
				table.append([gene_name.upper(), transcript_id])

	output_csv = OUTPUT_FOLDER + 'gene_identifiers.csv'
	write_table_to_csv(table, output_csv)


class GeneScore():
	def __init__(self):
		self.gene_name = ''
		self.gene_id = ''
		self.canonical_transcript = ''
		self.approved_symbol = ''
		self.length = 0

		self.gevir_percentile = 0.0
		self.oe_lof_percentile = 0.0
		self.gevir_and_oe_lof_percentile = 0.0

		# Start of statistical comparison raw data
		# Examined group stats
		self.gevir_ex_group_range = ''
		self.oe_lof_ex_group_range = ''
		self.gevir_and_oe_lof_ex_group_range = ''

		self.gevir_ex_group_ad_num = 0
		self.oe_lof_ex_group_ad_num = 0
		self.gevir_and_oe_lof_ex_group_ad_num = 0

		self.gevir_ex_group_ar_num = 0
		self.oe_lof_ex_group_ar_num = 0
		self.gevir_and_oe_lof_ex_group_ar_num = 0

		self.gevir_ex_group_all_num = 0
		self.oe_lof_ex_group_all_num = 0
		self.gevir_and_oe_lof_ex_group_all_num = 0

		# Reference (i.e. all other genes) group stats
		self.ref_group_ad_num = 0
		self.ref_group_ar_num = 0
		self.ref_group_all_num = 0

		# End of statistical comparison raw data

		self.gevir_ad_enrichment = 0.0
		self.oe_lof_ad_enrichment = 0.0
		self.gevir_and_oe_lof_ad_enrichment = 0.0

		self.gevir_ar_enrichment = 0.0
		self.oe_lof_ar_enrichment = 0.0
		self.gevir_and_oe_lof_ar_enrichment = 0.0

		self.gevir_null_enrichment = 0.0
		self.oe_lof_null_enrichment = 0.0
		self.gevir_and_oe_lof_null_enrichment = 0.0

		self.gevir_ad_p = 1.0
		self.oe_lof_ad_p = 1.0
		self.gevir_and_oe_lof_ad_p = 1.0

		self.gevir_ar_p = 1.0
		self.oe_lof_ar_p = 1.0
		self.gevir_and_oe_lof_ar_p = 1.0

		self.gevir_null_p = 1.0
		self.oe_lof_null_p = 1.0
		self.gevir_and_oe_lof_null_p = 1.0

		self.include_gene_groups = False
		self.ad = 'N'
		self.ar = 'N'
		self.null = 'N'
		self.cell_essential = 'N'
		self.cell_non_essential = 'N'
		self.mouse_het_lethal = 'N'
		self.ccrs_gte_95 = ''

		self.gnomad_outlier = 'N'
		self.gnomad_gene_issues = ''

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['gene_name'] = self.gene_name
		dictionary['gene_id'] = self.gene_id
		dictionary['canonical_transcript'] = self.canonical_transcript
		dictionary['approved_symbol'] = self.approved_symbol
		dictionary['protein_length'] = self.length

		dictionary['gevir_percentile'] = self.gevir_percentile
		dictionary['loeuf_percentile'] = self.oe_lof_percentile
		dictionary['virlof_percentile'] = self.gevir_and_oe_lof_percentile

		# Raw stats number

		dictionary['gevir_examined_range (ex_r)'] = self.gevir_ex_group_range
		dictionary['loeuf_examined_range (ex_r)'] = self.oe_lof_ex_group_range
		dictionary['virlof_examined_range (ex_r)'] = self.gevir_and_oe_lof_ex_group_range

		dictionary['gevir_ex_r_ad_num'] = self.gevir_ex_group_ad_num
		dictionary['loeuf_ex_r_ad_num'] = self.oe_lof_ex_group_ad_num
		dictionary['virlof_ex_r_ad_num'] = self.gevir_and_oe_lof_ex_group_ad_num

		dictionary['gevir_ex_r_ar_num'] = self.gevir_ex_group_ar_num
		dictionary['loeuf_ex_r_ar_num'] = self.oe_lof_ex_group_ar_num
		dictionary['virlof_ex_r_ar_num'] = self.gevir_and_oe_lof_ex_group_ar_num

		dictionary['gevir_ex_r_total'] = self.gevir_ex_group_all_num
		dictionary['loeuf_ex_r_total'] = self.oe_lof_ex_group_all_num
		dictionary['virlof_ex_r_total'] = self.gevir_and_oe_lof_ex_group_all_num

		dictionary['all_genes_ad_num'] = self.ref_group_ad_num
		dictionary['all_genes_ar_num'] = self.ref_group_ar_num
		dictionary['all_genes_total'] = self.ref_group_all_num


		dictionary['gevir_ad_enrichment'] = self.gevir_ad_enrichment
		dictionary['loeuf_ad_enrichment'] = self.oe_lof_ad_enrichment
		dictionary['virlof_ad_enrichment'] = self.gevir_and_oe_lof_ad_enrichment

		dictionary['gevir_ar_enrichment'] = self.gevir_ar_enrichment
		dictionary['loeuf_ar_enrichment'] = self.oe_lof_ar_enrichment
		dictionary['virlof_ar_enrichment'] = self.gevir_and_oe_lof_ar_enrichment

		dictionary['gevir_ad_p'] = self.gevir_ad_p
		dictionary['loeuf_ad_p'] = self.oe_lof_ad_p
		dictionary['virlof_ad_p'] = self.gevir_and_oe_lof_ad_p

		dictionary['gevir_ar_p'] = self.gevir_ar_p
		dictionary['loeuf_ar_p'] = self.oe_lof_ar_p
		dictionary['virlof_ar_p'] = self.gevir_and_oe_lof_ar_p

		if self.include_gene_groups == True:
			dictionary['ad_group'] = self.ad
			dictionary['ar_group'] = self.ar
			dictionary['null_group'] = self.null
			dictionary['cell_essential_group'] = self.cell_essential
			dictionary['cell_non_essential_group'] = self.cell_non_essential
			dictionary['mouse_het_lethal_group'] = self.mouse_het_lethal
			dictionary['ccrs_gte_95'] = self.ccrs_gte_95

		dictionary['gnomad_outlier'] = self.gnomad_outlier
		dictionary['gnomad_gene_issues'] = self.gnomad_gene_issues

		return dictionary


def get_transcript_ccrs(db):
	transcripts_ccrs = {}
	ccr_genes = db.gevir.ccr_genes.find({})
	for ccr_gene in ccr_genes:
		transcripts_ccrs[ccr_gene['transcript_id']] = ccr_gene['ccrs_gte_95']
	return transcripts_ccrs


def get_transcript_lengths(db):
	transcript_lengths = {}
	ens_proteins = db.gevir.ens_aa_fasta.find({})
	for ens_protein in ens_proteins:
		length = len(ens_protein['cds']) - 1 # do not count stop codon
		transcript_lengths[ens_protein['_id']] = length
	return transcript_lengths


def get_transcript_id_to_gene_id_dict(db):
	transcript_id_to_gene_id = {}
	genes = db.exac.genes.find({})
	for gene in genes:
		if 'canonical_transcript' not in gene:
			continue
		transcript_id_to_gene_id[gene['canonical_transcript']] = gene['gene_id']
	return transcript_id_to_gene_id


def get_gene_id_to_hugo_name(db):
	gene_id_to_hugo_name = {}
	hugo_genes = db.gevir.hugo.find({ "ensembl_gene_id": {"$exists": True }, "symbol": { "$exists": True }, "status": "Approved" })
	for hugo_gene in hugo_genes:
		gene_id = hugo_gene['ensembl_gene_id']
		gene_name = hugo_gene['symbol']
		gene_id_to_hugo_name[gene_id] = gene_name
	return gene_id_to_hugo_name


def export_gene_scores(db, enrichment_offset=1000, include_gene_groups=False):
	gene_scores = {}
	gevir_genes = db.gevir.gevir_scores.find({})
	omim_sets = OmimSets(db)
	essential_sets = EssentialSets(db)

	transcript_ccrs = get_transcript_ccrs(db)
	transcript_lengths = get_transcript_lengths(db)
	transcript_id_to_gene_id = get_transcript_id_to_gene_id_dict(db)
	gene_id_to_hugo_name = get_gene_id_to_hugo_name(db)

	for gevir_gene in gevir_genes:
		gene_score = GeneScore()
		transcript_id = gevir_gene['_id']
		gene_score.gene_name = gevir_gene['gene_name']
		gene_score.canonical_transcript = gevir_gene['_id']
		gene_id = transcript_id_to_gene_id[transcript_id]
		gene_score.gene_id = gene_id
		
		if gene_id in gene_id_to_hugo_name:
			gene_score.approved_symbol = gene_id_to_hugo_name[gene_id]

		gene_score.length = transcript_lengths[transcript_id]

		if include_gene_groups:
			gene_score.include_gene_groups = True
			if transcript_id in omim_sets.ad:
				gene_score.ad = 'Y'
			if transcript_id in omim_sets.ar:
				gene_score.ar = 'Y'
			if transcript_id in essential_sets.nulls:
				gene_score.null = 'Y'
			if transcript_id in essential_sets.crispr_essential:
				gene_score.cell_essential = 'Y'
			if transcript_id in essential_sets.crispr_non_essential:
				gene_score.cell_non_essential = 'Y'
			if transcript_id in essential_sets.mouse_het_lethal:
				gene_score.mouse_het_lethal = 'Y'

			if transcript_id in transcript_ccrs:
				gene_score.ccrs_gte_95 = transcript_ccrs[transcript_id]

		gnomad_gene = db.gevir.gnomad_scores.find_one({'_id': transcript_id})
		if not gnomad_gene['no_issues']:
			gene_score.gnomad_outlier = 'Y'
			gene_score.gnomad_gene_issues = ', '.join(gnomad_gene['gene_issues'])

		gene_scores[gevir_gene['_id']] = gene_score

	# Read gene scores
	gevir_scores = {}
	oe_lof_scores = {}
	gevir_and_oe_lof_scores = {}
	transcript_ids = gene_scores.keys()
	for transcript_id in transcript_ids:
		common_gene_scores = db.gevir.common_gene_scores.find_one({'_id': transcript_id})
		gevir_scores[transcript_id] = common_gene_scores['gevir_score']
		oe_lof_scores[transcript_id] = common_gene_scores['gnomad_oe_lof_upper']
		gevir_and_oe_lof_scores[transcript_id] = common_gene_scores['combined_rank']

	# Sort genes (most intolerant -> least intolerant)
	gevir_scores = sort_dict_by_values(gevir_scores, reverse=True)
	oe_lof_scores = sort_dict_by_values(oe_lof_scores)
	gevir_and_oe_lof_scores = sort_dict_by_values(gevir_and_oe_lof_scores)

	# Calculate percentiles (low percentile = most intolerant)
	gevir_percentiles = calculate_percentiles(gevir_scores.keys())
	oe_lof_percentiles = calculate_percentiles(oe_lof_scores.keys())
	gevir_and_oe_lof_percentiles = calculate_percentiles(gevir_and_oe_lof_scores.keys())

	for transcript_id, gene_score in gene_scores.iteritems():
		gene_score.gevir_percentile = gevir_percentiles[transcript_id]
		gene_score.oe_lof_percentile = oe_lof_percentiles[transcript_id]
		gene_score.gevir_and_oe_lof_percentile = gevir_and_oe_lof_percentiles[transcript_id]

	# Get sets of OMIM AD & AR genes
	omim_ad_set = set([])
	omim_ar_set = set([])

	omim_genes = db.gevir.omim.find({})
	for omim_gene in omim_genes:
		inheritance = omim_gene['inheritance']
		transcript_id = omim_gene['transcript_id']

		if inheritance == 'AD':
			omim_ad_set.add(transcript_id)
		if inheritance == 'AR':
			omim_ar_set.add(transcript_id)

	# Consider only genes which have scores
	omim_ad_set = omim_ad_set & set(gene_scores.keys())
	omim_ar_set = omim_ar_set & set(gene_scores.keys())

	ad_num = len(omim_ad_set)
	ar_num = len(omim_ar_set)
	all_num = len(gene_scores.keys())

	gevir_transcripts = gevir_percentiles.keys()
	oe_lof_transcripts = oe_lof_percentiles.keys()
	gevir_and_oe_lof_transcripts = gevir_and_oe_lof_percentiles.keys()

	gene_num = len(gene_scores)

	temp_data = {}
	total_lines = gene_num
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for x in range(0,gene_num):
		start = x - enrichment_offset
		stop = x + enrichment_offset
		if start <= 0:
			start = 0
			stop += 1
		elif stop >= gene_num:
			stop = gene_num
		else:
			start -= 1

		# Add 1 to start to include first gene which has 0 index in the list
		examined_range = '{:.2f}-{:.2f}%'.format(float(start + 1) * 100 / gene_num, float(stop) * 100 / gene_num)


		gevir_window = set(gevir_transcripts[start:stop])
		gevir_window_num = len(gevir_window)

		gevir_transcript = gevir_transcripts[x]
		gevir_ad_num = len(omim_ad_set & gevir_window)
		gevir_ar_num = len(omim_ar_set & gevir_window)

		gene_scores[gevir_transcript].gevir_ex_group_range = examined_range
		gene_scores[gevir_transcript].gevir_ad_enrichment, gene_scores[gevir_transcript].gevir_ad_p = fisher_exact([[gevir_ad_num, gevir_window_num],[ad_num, all_num]])
		gene_scores[gevir_transcript].gevir_ar_enrichment, gene_scores[gevir_transcript].gevir_ar_p = fisher_exact([[gevir_ar_num, gevir_window_num],[ar_num, all_num]])

		gene_scores[gevir_transcript].gevir_ex_group_ad_num = gevir_ad_num
		gene_scores[gevir_transcript].gevir_ex_group_ar_num = gevir_ar_num
		gene_scores[gevir_transcript].gevir_ex_group_all_num = gevir_window_num

		# Note reference group stats the same for all ranking lists
		gene_scores[gevir_transcript].ref_group_ad_num = ad_num
		gene_scores[gevir_transcript].ref_group_ar_num = ar_num
		gene_scores[gevir_transcript].ref_group_all_num = all_num


		oe_lof_window = set(oe_lof_transcripts[start:stop])
		oe_lof_window_num = len(oe_lof_window)
		oe_lof_transcript = oe_lof_transcripts[x]
		oe_lof_ad_num = len(omim_ad_set & oe_lof_window)
		oe_lof_ar_num = len(omim_ar_set & oe_lof_window)

		gene_scores[oe_lof_transcript].oe_lof_ex_group_range = examined_range
		gene_scores[oe_lof_transcript].oe_lof_ad_enrichment, gene_scores[oe_lof_transcript].oe_lof_ad_p = fisher_exact([[oe_lof_ad_num, oe_lof_window_num],[ad_num, all_num]])
		gene_scores[oe_lof_transcript].oe_lof_ar_enrichment, gene_scores[oe_lof_transcript].oe_lof_ar_p = fisher_exact([[oe_lof_ar_num, oe_lof_window_num],[ar_num, all_num]])

		gene_scores[oe_lof_transcript].oe_lof_ex_group_ad_num = oe_lof_ad_num
		gene_scores[oe_lof_transcript].oe_lof_ex_group_ar_num = oe_lof_ar_num
		gene_scores[oe_lof_transcript].oe_lof_ex_group_all_num = oe_lof_window_num


		gevir_and_oe_lof_window = set(gevir_and_oe_lof_transcripts[start:stop])
		gevir_and_oe_lof_window_num = len(gevir_and_oe_lof_window)
		gevir_and_oe_lof_transcript = gevir_and_oe_lof_transcripts[x]
		gevir_and_oe_lof_ad_num = len(omim_ad_set & gevir_and_oe_lof_window)
		gevir_and_oe_lof_ar_num = len(omim_ar_set & gevir_and_oe_lof_window)

		gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ex_group_range = examined_range
		gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ad_enrichment, gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ad_p = fisher_exact([[gevir_and_oe_lof_ad_num, gevir_and_oe_lof_window_num],[ad_num, all_num]])
		gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ar_enrichment, gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ar_p = fisher_exact([[gevir_and_oe_lof_ar_num, gevir_and_oe_lof_window_num],[ar_num, all_num]])

		gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ex_group_ad_num = gevir_and_oe_lof_ad_num
		gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ex_group_ar_num = gevir_and_oe_lof_ar_num
		gene_scores[gevir_and_oe_lof_transcript].gevir_and_oe_lof_ex_group_all_num = gevir_and_oe_lof_window_num

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	# Write results to csv and database
	ex_gene_score = GeneScore()
	ex_gene_score.include_gene_groups = include_gene_groups
	ex_gene_score = ex_gene_score.get_dictionary()
	headers = ex_gene_score.keys()
	table = [headers]
	db.gevir.web_gene_scores.drop()
	for gene_score in gene_scores.values():
		gene_score = gene_score.get_dictionary()
		table.append(gene_score.values())
		db.gevir.web_gene_scores.insert(gene_score)

	output_csv = OUTPUT_FOLDER + 'gene_scores.csv'
	write_table_to_csv(table, output_csv)


class RegionData():
	def __init__(self):
		self.region_id = ''
		self.transcript_id = ''
		self.chrom = ''
		self.strand = ''
		self.start = 0
		self.stop = 0
		self.xstart = 0
		self.xstop = 0
		self.codon_xstart = 0
		self.codon_xstop = 0
		self.region_flag = ''
		self.length = 0
		self.exome_coverage = 0.0
		self.gerp_mean = 0.0
		self.percentile = 0.0
		self.start_variant = ''
		self.start_variant_csq = ''
		self.stop_variant = ''
		self.stop_variant_csq = ''

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['region_id'] = self.region_id
		dictionary['transcript_id'] = self.transcript_id
		dictionary['chrom'] = self.chrom
		dictionary['strand'] = self.strand
		dictionary['start'] = self.start
		dictionary['stop'] = self.stop
		dictionary['xstart'] = self.xstart
		dictionary['xstop'] = self.xstop
		dictionary['codon_xstart'] = self.codon_xstart
		dictionary['codon_xstop'] = self.codon_xstop
		dictionary['region_flag'] = self.region_flag
		dictionary['length'] = self.length
		dictionary['exome_coverage'] = self.exome_coverage
		dictionary['gerp_mean'] = self.gerp_mean
		dictionary['percentile'] = self.percentile
		dictionary['start_variant'] = self.start_variant
		dictionary['start_variant_csq'] = self.start_variant_csq
		dictionary['stop_variant'] = self.stop_variant
		dictionary['stop_variant_csq'] = self.stop_variant_csq
		return dictionary


def export_regions(db):
	db.gevir.web_regions.drop()
	transcript_ids = []
	gevir_genes = db.gevir.gevir_scores.find({})
	for gevir_gene in gevir_genes:
		transcript_ids.append(gevir_gene['_id'])

	# Get region weights calculated for each length (X) based on how common are regions of X length or longer
	autosomal_region_weights = db.gevir.region_weights.find_one({'_id': 'autosomal'})
	autosomal_region_weights = autosomal_region_weights['region_size_weights']

	xy_region_weights = db.gevir.region_weights.find_one({'_id': 'xy'})
	xy_region_weights = xy_region_weights['region_size_weights']

	transcript_chroms = {}
	transcript_strands = {}

	region_ids = []
	region_scores = []

	# Loop over PASS regions to get data for percentiles calculation
	#transcript_ids = ['ENST00000263100', 'ENST00000302118']
	total_lines = len(transcript_ids)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id in transcript_ids:
		exac_transcript = db.exac.transcripts.find_one({'transcript_id': transcript_id})
		chrom = exac_transcript['chrom']
		strand = exac_transcript['strand']
		transcript_chroms[transcript_id] = chrom
		transcript_strands[transcript_id] = strand

		# Note: We use different coverage thresholds for autosomes and xy chromosomes
		is_xy = False
		coverage_threshold = MIN_AUTOSOMAL_COVERAGE

		if chrom == 'X' or chrom == 'Y':
			 is_xy = True
			 coverage_threshold = MIN_XY_COVERAGE
		
		regions = db.gevir.variant_regions.find({'transcript_id': transcript_id, 'exome_coverage': { '$gte': coverage_threshold } })
		for region in regions:
			region_id = str(region['xstart']) + '-' + str(region['xstop'])
			region_ids.append(region_id)

			if is_xy:
				region_weight = xy_region_weights[str(region['lenght'])]
			else:
				region_weight = autosomal_region_weights[str(region['lenght'])]

			region_weight = region_weight['weight']
			gerp = region['gerp_mean']
			# Note: this is added to avoid severe penalties for fractional gerp scores
			if gerp >= 0 and gerp < 1:
				gerp = 1
			if gerp > -1 and gerp < 0:
				gerp = -1

			region_score = region_weight * gerp
			region_scores.append(region_score)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	# rank, reverse ranking and convert it to percentiles: 
	region_ranks = rankdata(region_scores, method='min')
	region_percentiles = region_ranks / float(len(region_scores))
	min_percentile = min(region_percentiles)
	reversed_region_percentiles = []
	for region_percentile in region_percentiles:
		reversed_percentile = ((1 - region_percentile) + min_percentile) * 100
		reversed_region_percentiles.append(reversed_percentile)

	region_id_percentiles = OrderedDict()
	for x in range(0, len(region_scores)):
		region_id = region_ids[x]
		region_percentile = reversed_region_percentiles[x]
		region_id_percentiles[region_id] = region_percentile

	# Loop over regions again and collect all data
	total_lines = len(transcript_ids)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id in transcript_ids:
		transcript_amino_acids_xpositions = get_transcript_amino_acids_xpositions(db, transcript_id)

		transcript_regions_data = []
		regions = db.gevir.variant_regions.find({'transcript_id': transcript_id})
		for region in regions:
			region_id = str(region['xstart']) + '-' + str(region['xstop'])
			chrom = transcript_chroms[transcript_id]
			strand = transcript_strands[transcript_id]

			region_data = RegionData()

			region_data.region_id = region_id
			region_data.transcript_id = transcript_id
			region_data.chrom = chrom
			region_data.strand = strand

			region_data.start = xpos_to_pos(region['xstart'])
			region_data.stop = xpos_to_pos(region['xstop'])
			region_data.xstart = region['xstart']
			region_data.xstop = region['xstop']

			# Here we get start and stop protein codons xpositions, so that 
			# the actual unafected amino acids will be between codon_xstart and codon_xstop (but not equal to).
			# Note that we reverse start and stop variants if strand is negative.
			if strand == '+':
				start_pos = region['start_variant']['protein_pos']
				stop_pos = region['stop_variant']['protein_pos']
				region_data.codon_xstart = transcript_amino_acids_xpositions[start_pos][1]
				region_data.codon_xstop = transcript_amino_acids_xpositions[stop_pos][0]
			else:
				start_pos = region['stop_variant']['protein_pos']
				stop_pos = region['start_variant']['protein_pos']
				region_data.codon_xstart = transcript_amino_acids_xpositions[start_pos][0]
				region_data.codon_xstop = transcript_amino_acids_xpositions[stop_pos][1]

			# PASS - good coverage (50 autosome, 45 XY), functional variant csq (e.g. LoF and Missense)
			# LOW_COVERAGE - coverage less tham 50 for autosomes or 45 for X & Y chromosomes
			# INCONSISTENT - at least one variant has a feunctional csq, but position is not in coding sequence (exons) according to ExAC
			region_data.region_flag = 'PASS'

			coverage_threshold = MIN_AUTOSOMAL_COVERAGE
			if chrom == 'X' or chrom == 'Y':
				coverage_threshold = MIN_XY_COVERAGE

			if region['exome_coverage'] < coverage_threshold:
				region_data.region_flag = 'LOW_COVERAGE'

			if region['not_in_cds']:
				region_data.region_flag = 'INCONSISTENT'

			region_data.length = region['lenght']
			region_data.exome_coverage = region['exome_coverage']
			region_data.gerp_mean = region['gerp_mean']

			# LOW_COVERAGE and INCONSISTENT regions which were not used in percentile calculation
			# are assigned the maximum percentile (least interesting) - 100
			if region_id in region_id_percentiles:
				region_data.percentile = region_id_percentiles[region_id]
			else:
				region_data.percentile = 100.0

			region_data.start_variant = region['start_variant']['variant_id']
			region_data.start_variant_csq = region['start_variant']['csq']
			region_data.stop_variant = region['stop_variant']['variant_id']
			region_data.stop_variant_csq = region['stop_variant']['csq']

			if region_data.start_variant_csq == 'fake_start':
				region_data.start_variant_csq = 'START CODON'

			if region_data.stop_variant_csq == 'fake_stop':
				region_data.stop_variant_csq = 'STOP CODON'

			region_data = region_data.get_dictionary()
			transcript_regions_data.append(region_data)

		db.gevir.web_regions.insert_many(transcript_regions_data)
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	# Write results to csv and database
	print 'Exporting data to csv...'
	output_csv = OUTPUT_FOLDER + 'regions.csv'
	output_file = open(output_csv,'w+')
	writer = csv.writer(output_file)

	ex_region_data = RegionData()
	ex_region_data = ex_region_data.get_dictionary()
	headers = ex_region_data.keys()
	writer.writerow(headers)

	web_regions = db.gevir.web_regions.find({})
	total_lines = web_regions.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for web_region in web_regions:
		del web_region['_id']
		row = web_region.values()
		writer.writerow(row)
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()
	output_file.close()


class CodingRegion():
	def __init__(self):
		self.transcript_id = ''
		self.xstart = ''
		self.xstop = ''

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['transcript_id'] = self.transcript_id
		dictionary['xstart'] = self.xstart
		dictionary['xstop'] = self.xstop
		return dictionary


def export_coding_regions(db):
	coding_regions = []
	transcript_ids = []
	gevir_genes = db.gevir.gevir_scores.find({})
	for gevir_gene in gevir_genes:
		transcript_id = gevir_gene['_id']
		exons = db.exac.exons.find({'transcript_id': transcript_id, 'feature_type': 'CDS'})
		for exon in exons:
			coding_region = CodingRegion()
			coding_region.transcript_id = exon['transcript_id']
			if exon['xstart'] <= exon['xstop']:
				coding_region.xstart = exon['xstart']
				coding_region.xstop = exon['xstop']
			else:
				coding_region.xstart = exon['xstop']
				coding_region.xstop = exon['xstart']
			coding_regions.append(coding_region)

	ex_coding_region = CodingRegion()
	ex_coding_region = ex_coding_region.get_dictionary()
	headers = ex_coding_region.keys()
	table = [headers]
	for coding_region in coding_regions:
		coding_region = coding_region.get_dictionary()
		table.append(coding_region.values())

	output_csv = OUTPUT_FOLDER + 'coding_regions.csv'
	write_table_to_csv(table, output_csv)


def main():
	db = MongoDB()
	# Creates main supplementary table with GeVIR, LOEUF and VIRLoF scores
	# is required to run "draw_web_gene_scores" method in figures.py
	# enrichment_offset is used to calculate AD/AR fold-enrichment
	# 5% offset = 18,352 / 20 ~= 918
	# 5% offset = 19,361 / 20 ~= 968
	# Supplementary Table 2
	export_gene_scores(db, enrichment_offset=AD_AR_ENRICHMENT_OFFSET, include_gene_groups=True)

	# Exports datasets required to build the website
	#export_gene_identifiers(db)
	#export_regions(db)
	#export_coding_regions(db)

if __name__ == "__main__":
	sys.exit(main())