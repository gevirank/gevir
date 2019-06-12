import os
import sys
import progressbar
import pymongo
import csv
import numpy as np
from collections import OrderedDict
from common import MongoDB, is_int, sort_dict_by_values, write_table_to_csv
from gnomad_utils import worst_csq_from_csq, xpos_to_pos

#################
### CONSTANTS ###
#################

VALID_CSQS = set(['stop_gained',
				  'frameshift_variant',
				  'stop_lost',
				  'start_lost',
				  'inframe_insertion',
				  'inframe_deletion',
				  'missense_variant',])

OUTPUT_FOLDER = './tables/'

REGIONS_COLLECTION = 'variant_regions'
GEVIR_SCORES_COLLECTION = 'gevir_scores'

VALID_FILTERS = set(['PASS', 'SEGDUP', 'LCR'])

MIN_AUTOSOMAL_COVERAGE = 50
MIN_XY_COVERAGE = 45

INCLUDE_GNOMAD_OUTLIERS = True

'''
IMPORTANT!
1) This script has to be run AFTER data_import script which creates custom collections used here.
2) This script REQUIRES following collections from gnomAD database:
   exome_coverage, exome_variants, genome_variants, exons, genes, transcripts.
'''

##########################
### REGIONS COLLECTION ###
##########################


class Variant():
	def __init__(self):
		self.variant_id = ''
		self.csq = ''
		self.xpos = 0
		self.protein_pos = 0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['variant_id'] = self.variant_id
		dictionary['csq'] = self.csq
		dictionary['xpos'] = self.xpos
		dictionary['protein_pos'] = self.protein_pos
		return dictionary


class Region():
	def __init__(self):
		self.transcript_id = ''
		self.lenght = 0
		self.exome_coverage = 0.0
		self.gerp_mean = 0.0
		self.gerp_std = 0.0
		self.start_variant = {}
		self.stop_variant = {}
		self.is_fake = False
		self.not_in_cds = False

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['transcript_id'] = self.transcript_id
		dictionary['lenght'] = self.lenght
		dictionary['exome_coverage'] = self.exome_coverage
		dictionary['gerp_mean'] = self.gerp_mean
		dictionary['gerp_std'] = self.gerp_std
		if self.start_variant.xpos < self.stop_variant.xpos:
			dictionary['xstart'] = self.start_variant.xpos
			dictionary['xstop'] = self.stop_variant.xpos
		else:
			dictionary['xstart'] = self.stop_variant.xpos
			dictionary['xstop'] = self.start_variant.xpos

		dictionary['start_variant'] = self.start_variant.get_dictionary()
		dictionary['stop_variant'] = self.stop_variant.get_dictionary()
		dictionary['is_fake'] = self.is_fake
		dictionary['not_in_cds'] = self.not_in_cds
		return dictionary


def read_variants(variants, transcript_id, transcript_vars):
	for variant in variants:
		var_filter = variant['filter']
		# Skip non-pass quality variants
		if var_filter not in VALID_FILTERS:
			continue

		variant_id = variant['variant_id']
		veps = variant['vep_annotations']
		for vep in veps:
			if vep['Feature'] == transcript_id:
				csq = worst_csq_from_csq(vep['Consequence'])
				# Analyse only non-synonymous protein coding variants (except splicing)
				if csq in VALID_CSQS:
					var = Variant()
					var.variant_id = variant_id
					var.csq = csq
					var.xpos = variant['xpos']
					protein_pos = vep['Protein_position']

					# For indels and frameshifts use start position
					if '-' in protein_pos:
						protein_pos = protein_pos.split('-')[0]
					# Check that protein position can be converted to integer
					if is_int(protein_pos):
						var.protein_pos = int(protein_pos)
						transcript_vars[int(protein_pos)] = var
	return transcript_vars


def get_transcript_amino_acids_xpositions(db, transcript_id):
	exons = db.exac.exons.find({'transcript_id': transcript_id, 'feature_type': 'CDS'})
	
	strand = '+'
	exon_xstarts = []
	exon_xstops = []
	for exon in exons:
		if exon['strand'] == '-':
			strand = '-'
		exon_xstarts.append(exon['xstart'])
		exon_xstops.append(exon['xstop'])

	exon_xposes = []
	for exon_n in range(0, len(exon_xstarts)):
		exon_xstart = exon_xstarts[exon_n]
		exon_xstop = exon_xstops[exon_n]

		exon_xposes += range(exon_xstart-1, exon_xstop)

	# Last exon does not include stop codon which has to be added separately
	# Strand changes order of xpositions (+:Ascending;-:Descending)
	if strand == '+':
		exon_xposes.sort()
		stop_codon = range(exon_xposes[-1:][0] + 1, exon_xposes[-1:][0] + 4)
		exon_xposes += stop_codon
	else:
		exon_xposes.sort(reverse=True)
		stop_codon = range(exon_xposes[-1:][0] - 3, exon_xposes[-1:][0])
		stop_codon.sort(reverse=True)
		exon_xposes += stop_codon

	codon_xposes = [exon_xposes[i:i + 3] for i in xrange(0, len(exon_xposes), 3)]

	# calcualte each xstart and xstop for each transcript protein position (codon)
	codons = {}
	for x in range(0, len(codon_xposes)):
		codons[x + 1] = (codon_xposes[x][0], codon_xposes[x][2])
		#codons[x + 1] = (codon_xposes[x][0], codon_xposes[x][1], codon_xposes[x][2])

	return codons


def get_coverage(db, region, coverage_collection, transcript_amino_acids_xpositions):
	chrom = region.start_variant.variant_id.split('-')[0]
	'''
	Coverage is calculated between start and stop variant xpositions and includes variant xpositions (e.g. greater/lower or equal).
	We include variants xpositions to know coverage even of regions with 0 length (2 adjacent variants). 
	GERP score is calculated between start and stop xpositions of AMINO ACIDS between variants and excludes affected amino acids.
	We excluded all xpositions of amino acids in codons affected by variants for GERP score calculation because
	they might be less conservative and might bias results.
	Therefore GERP is calculated only for regions with length in amino acids > 0
	'''

	# calculate first and last xpos for each amino acid (codon) in a transcript
	#amino_acids_xpositions = get_transcript_amino_acids_xpositions(db, region.transcript_id)

	# xposes of start and stop codon between variants
	start_inside_codon_xposes = transcript_amino_acids_xpositions[region.start_variant.protein_pos + 1]
	stop_inside_codon_xposes = transcript_amino_acids_xpositions[region.stop_variant.protein_pos - 1]

	# start xpos is lower than stop xpos if strand is +, but values have to be swapped if strand is -
	if region.start_variant.xpos < region.stop_variant.xpos:
		start_xpos = region.start_variant.xpos
		stop_xpos = region.stop_variant.xpos

		gerp_start_xpos = start_inside_codon_xposes[0]
		gerp_stop_xpos = stop_inside_codon_xposes[1]
	else:
		start_xpos = region.stop_variant.xpos
		stop_xpos = region.start_variant.xpos

		gerp_start_xpos = stop_inside_codon_xposes[1]
		gerp_stop_xpos = start_inside_codon_xposes[0]

	# region length in nucleotide bases
	x_length = stop_xpos - start_xpos
	
	mean_coverages = []
	gerp_scores = []
	not_in_cds = False

	# Check if both variants are located in the same exon.
	# region.lenght is in amino acids so has to be multiplied by 3 to be compared with length in xpositions (nucleotides)
	# 5 is a maximum possible difference between 2 nucleotides in adjacent codons (e.g. 1 and 6)
	if x_length <= (5 + region.lenght * 3):
		exome_coverages = db.exac[coverage_collection].find({ "$and": [ { "xpos": { "$gte": start_xpos } }, { "xpos": { "$lte": stop_xpos } } ] })
		for exome_coverage in exome_coverages:
			mean_coverages.append(exome_coverage['mean'])

		if region.lenght > 0:
			gerp_start = xpos_to_pos(gerp_start_xpos)
			gerp_stop = xpos_to_pos(gerp_stop_xpos)

			gerp_data = db.gerp[chrom].find({ '$and': [ { "_id": { '$gte': gerp_start } }, { "_id": { '$lte': gerp_stop } } ] })
			for gerp_pos in gerp_data:
				gerp_scores.append(gerp_pos['RS'])
	else:
		# Find exons of start and stop variants
		start_exon = db.exac.exons.find_one({"feature_type": "CDS", "transcript_id": region.transcript_id, "xstart": { '$lte': start_xpos }, "xstop": { '$gte': start_xpos }})
		stop_exon = db.exac.exons.find_one({"feature_type": "CDS", "transcript_id": region.transcript_id, "xstart": { '$lte': stop_xpos }, "xstop": { '$gte': stop_xpos }})
		
		if start_exon and stop_exon:
			# Add mean coverages from start variant xposition to end of an exon with start variant
			exome_coverages = db.exac[coverage_collection].find({ "$and": [ { "xpos": { "$gte": start_xpos } }, { "xpos": { "$lte": start_exon['xstop'] } } ] })
			for exome_coverage in exome_coverages:
				mean_coverages.append(exome_coverage['mean'])

			# Add GERP scores if variants are not in adjacent codons
			if region.lenght > 0:
				gerp_start = xpos_to_pos(gerp_start_xpos)
				gerp_stop = xpos_to_pos(start_exon['xstop'])
				gerp_data = db.gerp[chrom].find({ '$and': [ { "_id": { '$gte': gerp_start } }, { "_id": { '$lte': gerp_stop } } ] })
				for gerp_pos in gerp_data:
					gerp_scores.append(gerp_pos['RS'])
			
			# Add mean coverages from start of an exon with stop variant to stop variant xposition
			exome_coverages = db.exac[coverage_collection].find({ "$and": [ { "xpos": { "$gte": stop_exon['xstart'] } }, { "xpos": { "$lte": stop_xpos } } ] })
			for exome_coverage in exome_coverages:
				mean_coverages.append(exome_coverage['mean'])

			# Add GERP scores if variants are not in adjacent codons
			if region.lenght > 0:
				gerp_start = xpos_to_pos(stop_exon['xstart'])
				gerp_stop = xpos_to_pos(gerp_stop_xpos)
				gerp_data = db.gerp[chrom].find({ '$and': [ { "_id": { '$gte': gerp_start } }, { "_id": { '$lte': gerp_stop } } ] })
				for gerp_pos in gerp_data:
					gerp_scores.append(gerp_pos['RS'])

			# In case there are exons without variants between exons with start and stop variants 
			# (i.e. whole exon without non-synonymous variants, see SERHL2 gene for example)
			# we add all their mean coverages and GERP scores
			middle_exons = db.exac.exons.find({"feature_type": "CDS", "transcript_id": region.transcript_id, "xstart": { '$gte': start_exon['xstop'] }, "xstop": { '$lte': stop_exon['xstart'] }})
			for middle_exon in middle_exons:
				exome_coverages = db.exac[coverage_collection].find({ "$and": [ { "xpos": { "$gte": middle_exon['xstart'] } }, { "xpos": { "$lte": middle_exon['xstop'] } } ] })
				for exome_coverage in exome_coverages:
					mean_coverages.append(exome_coverage['mean'])

				# Add GERP scores if variants are not in adjacent codons
				if region.lenght > 0:
					gerp_start = xpos_to_pos(middle_exon['xstart'])
					gerp_stop = xpos_to_pos(middle_exon['xstop'])
					gerp_data = db.gerp[chrom].find({ '$and': [ { "_id": { '$gte': gerp_start } }, { "_id": { '$lte': gerp_stop } } ] })
					for gerp_pos in gerp_data:
						gerp_scores.append(gerp_pos['RS'])
		else:
			# Start, stop or both variants were not in CDS exons according to gnomAD database
			# Such regions will have 0 coverage and gerp scores, consequently not used in future analysis
			#print region.transcript_id, region.start_variant.variant_id, region.start_variant.protein_pos, region.stop_variant.variant_id, region.stop_variant.protein_pos
			not_in_cds = True

	mean_coverage = 0.0
	mean_gerp = 0.0
	std_gerp = 0.0

	if len(mean_coverages) > 0:
		# Calculate mean coverage of ALL positions mean coverages (e.g. mean of all individual coverages)
		mean_coverage = np.mean(mean_coverages)

	if len(gerp_scores) > 0:
		mean_gerp = np.mean(gerp_scores)
		std_gerp = np.std(gerp_scores)

	return mean_coverage, mean_gerp, std_gerp, not_in_cds


def create_fake_start_stop_variants(db, transcript_id):
	exons = db.exac.exons.find({"transcript_id": transcript_id, "feature_type": "CDS"})
	gene = db.exac.genes.find_one({"canonical_transcript": transcript_id})

	strand = gene['strand']

	'''
	# uncomment if script is used for transcripts with no coding exons in gnomAD
	if exons.count() == 0:
		return {}
	'''

	length = 0
	start_xposes = []
	start_poses = []

	stop_xposes = []
	stop_poses = []

	chrom = ''
	for exon in exons:
		length += exon['stop'] - exon['start'] + 1 # +1 to include each exon start nucleotide
		chrom = exon['chrom']
		start_xposes.append(exon['xstart'])
		start_poses.append(exon['start'])

		stop_xposes.append(exon['xstop'])
		stop_poses.append(exon['stop'])

	# Swap start and stop for transcripts with negative (-) strand
	if strand == '+':
		start_pos = min(start_poses)
		start_xpos = min(start_xposes)

		stop_pos = max(stop_poses)
		stop_xpos = max(stop_xposes)
	else:
		start_pos = max(stop_poses)
		start_xpos = max(stop_xposes)

		stop_pos = min(start_poses)
		stop_xpos = min(start_xposes)

	# Calculate length of a transcript in amino acids, which is protein position of stop codon
	length = length / 3 + 1 # +1 to include stop codon

	start_variant = Variant()
	start_variant.variant_id = chrom + '-' + str(start_pos) + '-START' 
	start_variant.csq = 'fake_start'
	start_variant.xpos = start_xpos
	start_variant.protein_pos = 1

	stop_variant = Variant()
	stop_variant.variant_id = chrom + '-' + str(stop_pos) + '-STOP' 
	stop_variant.csq = 'fake_stop'
	stop_variant.xpos = stop_xpos
	stop_variant.protein_pos = length

	fake_vars = {1: start_variant, length: stop_variant}
	return fake_vars


def calculate_transcript_regions(db, transcript_id):
	transcript_vars = {}

	# Read exome and genome variants into transcript_vars dictionary.
	# Note that variants at the same protein positions replace each other,
	# but it does not matter for region analysis.
	variants = db.exac.exome_variants.find({'transcripts': transcript_id})
	read_variants(variants, transcript_id, transcript_vars)
	variants = db.exac.genome_variants.find({'transcripts': transcript_id})
	read_variants(variants, transcript_id, transcript_vars)

	# Add fake start and stop variants to take into account 
	# potential regions without variations at the beginning and end of a transcript.
	# For example, if first variant in a transcript is at 15th amino acid, then
	# we'll analyse it as a 13 amino acids long region between fake start and valid variant.
	fake_vars = create_fake_start_stop_variants(db, transcript_id)
	for key, value in fake_vars.iteritems():
		transcript_vars[key] = value

	positions = transcript_vars.keys()
	positions.sort()

	# calculate first and last xpos for each amino acid (codon) in a transcript
	transcript_amino_acids_xpositions = get_transcript_amino_acids_xpositions(db, transcript_id)

	for x in range(0, len(positions) - 1):
		region = Region()
		region.transcript_id = transcript_id
		region.lenght = positions[x + 1] - positions[x] - 1 # Do not count amino acids affected by variants
		region.start_variant = transcript_vars[positions[x]]
		region.stop_variant = transcript_vars[positions[x + 1]]

		exome_coverage, mean_gerp, std_gerp, not_in_cds = get_coverage(db, region, 'exome_coverage', transcript_amino_acids_xpositions)

		region.exome_coverage = exome_coverage
		region.gerp_mean = mean_gerp
		region.gerp_std = std_gerp
		region.not_in_cds = not_in_cds
		if region.start_variant.csq == 'fake_start' or region.stop_variant.csq == 'fake_stop':
			region.is_fake = True

		db.gevir[REGIONS_COLLECTION].insert(region.get_dictionary())


def calculate_all_transcripts_regions(db):
	db.gevir[REGIONS_COLLECTION].drop()

	# We get only canonical transcripts with no reported sequencing issues in gnomAD, 
	# which starts start codon, ends with stop codon and which length in nucleotides is divisible by 3

	transcript_ids = []

	if INCLUDE_GNOMAD_OUTLIERS:
		gnomad_scores = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True })
	else:
		gnomad_scores = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	for gnomad_score in gnomad_scores:
		transcript_ids.append(gnomad_score['_id'])

	print 'Transcripts Number:', len(transcript_ids)
	total_lines = len(transcript_ids)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id in transcript_ids:
		calculate_transcript_regions(db, transcript_id)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	db.gevir[REGIONS_COLLECTION].create_index([('transcript_id', pymongo.ASCENDING)], name='transcript_id_1')
	db.gevir[REGIONS_COLLECTION].create_index([('lenght', pymongo.ASCENDING)], name='lenght_1')
	db.gevir[REGIONS_COLLECTION].create_index([('exome_coverage', pymongo.ASCENDING)], name='exome_coverage_1')


####################
### GEVIR SCORES ###
####################


class GevirGene():
	def __init__(self):
		self.transcript_id = ''
		self.gene_id = ''
		self.gene_name = ''
		self.chrom = ''
		self.gevir_score = 0
		self.gevir_rank = 0
		self.gevir_percentile = 0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['_id'] = self.transcript_id
		dictionary['gene_id'] = self.gene_id
		dictionary['gene_name'] = self.gene_name
		dictionary['chrom'] = self.chrom
		dictionary['gevir_score'] = self.gevir_score
		dictionary['gevir_rank'] = self.gevir_rank
		dictionary['gevir_percentile'] = self.gevir_percentile
		return dictionary


class RegionSizeWeight():
	def __init__(self):
		self.size = 0
		self.count = 0
		self.count_longer_regions = 0
		self.frequency = 0
		self.weight = 0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['length'] = self.size
		dictionary['count'] = self.count
		dictionary['count_longer_regions'] = self.count_longer_regions
		dictionary['frequency'] = self.frequency
		dictionary['weight'] = self.weight
		return dictionary


def calculate_region_weights(db, transcript_ids, min_mean_coverage, weights_report_name='', xy=False):
	# Count number of regions with sufficient coverage in length bins 
	region_sizes = {}
	all_regions = 0

	total_lines = len(transcript_ids)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id in transcript_ids:
		regions = db.gevir[REGIONS_COLLECTION].find({ "transcript_id": transcript_id, "exome_coverage": { "$gte": min_mean_coverage } })
		for region in regions:
			length = region['lenght']
			all_regions += 1
			if length in region_sizes:
				region_sizes[length] += 1
			else:
				region_sizes[length] = 1

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	# Sort all OBSERVED region sizes (e.g. 0, 1, 2, 3... 21, 22... 65, 68, 74 etc)
	sizes = region_sizes.keys()
	sizes.sort()

	# Calculate region size (length) X weights as 
	# 1 divided by proportion of regions with size < X out of all regions.
	# So region size 0 will have minimal weight 1 and all other region sizes will have larger weights.
	if weights_report_name:
		ex_region_weight = RegionSizeWeight()
		ex_region_weight = ex_region_weight.get_dictionary()
		#headers = ['size', 'count', 'count_longer_regions', 'frequency', 'weight']
		headers = ex_region_weight.keys()
		table = [headers]

	region_size_weights = OrderedDict()
	sum_regions_below = 0
	region_weights = OrderedDict()
	for size in sizes:
		weight = float(all_regions) / (all_regions - sum_regions_below)
		region_weights[size] = weight
		#region_weights[size] = 1 / ((all_regions - sum_regions_below) / float(all_regions))
		region_size_weight = RegionSizeWeight()
		region_size_weight.size = size
		region_size_weight.count = region_sizes[size]
		region_size_weight.count_longer_regions = all_regions - sum_regions_below
		region_size_weight.frequency = (all_regions - sum_regions_below) / float(all_regions)
		#region_size_weight.weight = 1 / ((all_regions - sum_regions_below) / float(all_regions))
		region_size_weight.weight = weight
		region_size_weight = region_size_weight.get_dictionary()
		region_size_weights[str(size)] = region_size_weight

		if weights_report_name:
			'''
			row = [size, region_sizes[size], all_regions - sum_regions_below, 
			       (all_regions - sum_regions_below) / float(all_regions),
			       1 / ((all_regions - sum_regions_below) / float(all_regions))]
			'''
			row = region_size_weight.values()
			table.append(row)

		sum_regions_below += region_sizes[size]

	if xy:
		db.gevir.region_weights.remove({'_id': 'xy'})
		db.gevir.region_weights.insert({'_id': 'xy', 'region_size_weights': region_size_weights})
	else:
		db.gevir.region_weights.remove({'_id': 'autosomal'})
		db.gevir.region_weights.insert({'_id': 'autosomal', 'region_size_weights': region_size_weights})

	if weights_report_name:
		output_csv = OUTPUT_FOLDER + weights_report_name
		write_table_to_csv(table, output_csv)
	return region_weights


def calculate_transcripts_scores(db, transcript_ids, min_mean_coverage, transcript_scores, weights_report_name='', xy=False):
	# calculate region weights (longer regions are rarely observed and have larger weights)
	region_weights = calculate_region_weights(db, transcript_ids, min_mean_coverage, weights_report_name=weights_report_name, xy=xy)
	
	total_lines = len(transcript_ids)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id in transcript_ids:
		# Get all regions in a transcript with sufficient coverage
		regions = db.gevir[REGIONS_COLLECTION].find({ "transcript_id": transcript_id, "exome_coverage": { "$gte": min_mean_coverage } })
		# Get all regions in a transcript
		regions_count = db.gevir[REGIONS_COLLECTION].find({ "transcript_id": transcript_id })
		regions_count = float(regions_count.count())

		'''
		Sum all region weights multiplied by GERP (can be positive or negative) in a gene
		Values between 1 and -1 are rounded to closest integer (i.e. 1 or -1)

		Note: GeVIR score can be negative if there are many non-coservative gene regions
		      Such regions might contain non-pass quality variants (e.g. MUC2)
		      and therefore have opposite impact on GeVIR score
		'''

		gevir_score = 0
		for region in regions:
			length = region['lenght']
			gerp = region['gerp_mean']
			
			if gerp < 1 and gerp > 0:
				gerp = 1
			elif gerp > -1 and gerp < 0:
				gerp = -1
			
			gevir_score += region_weights[length] * gerp

		# Normalise by number of regions in a transcript INCLUDING low covered regions
		# Longer transcript normally have more regions in general and 
		# have greater chance to contain long region, thus have to be penelised more


		gevir_score = gevir_score / regions_count

		transcript_scores[transcript_id] = gevir_score
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	return transcript_scores


def create_gevir_scores(db):
	transcript_ids = []
	transcript_ids_xy = []

	if INCLUDE_GNOMAD_OUTLIERS:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True })
	else:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	for gnomad_gene in gnomad_genes:
		if gnomad_gene['chrom'] == 'X' or gnomad_gene['chrom'] == 'Y':
			transcript_ids_xy.append(gnomad_gene['_id'])
		else:
			transcript_ids.append(gnomad_gene['_id'])

	transcript_scores = OrderedDict()
	print 'Autosomal Transcripts:', len(transcript_ids)
	calculate_transcripts_scores(db, transcript_ids, MIN_AUTOSOMAL_COVERAGE, transcript_scores, weights_report_name='autosomal_region_weights_table.csv', xy=False)
	print 'XY Transcripts:', len(transcript_ids_xy)
	calculate_transcripts_scores(db, transcript_ids_xy, MIN_XY_COVERAGE, transcript_scores, weights_report_name='xy_region_weights_table.csv', xy=True)
	# Sort transcripts by Gevir score (descending), used to calculate percentiles (low percentile = high intolerance)
	transcript_scores = sort_dict_by_values(transcript_scores, reverse=True)
	
	rank = 1
	max_rank = float(len(transcript_scores))

	db.gevir[GEVIR_SCORES_COLLECTION].drop()

	total_lines = len(transcript_scores)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for transcript_id, gevir_score in transcript_scores.iteritems():
		exac_gene = db.exac.genes.find_one({'canonical_transcript': transcript_id})
		gevir_gene = GevirGene()
		gevir_gene.transcript_id = transcript_id
		gevir_gene.gene_id = exac_gene['gene_id']
		gevir_gene.gene_name = exac_gene['gene_name']
		gevir_gene.chrom = exac_gene['chrom']
		gevir_gene.gevir_score = gevir_score
		gevir_gene.gevir_rank = rank
		gevir_gene.gevir_percentile = (rank / max_rank) * 100
		rank += 1

		db.gevir[GEVIR_SCORES_COLLECTION].insert(gevir_gene.get_dictionary())

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	db.gevir[GEVIR_SCORES_COLLECTION].create_index([('gene_name', pymongo.ASCENDING)], name='gene_name_1')
	db.gevir[GEVIR_SCORES_COLLECTION].create_index([('gene_id', pymongo.ASCENDING)], name='gene_id_1')


def main():
	db = MongoDB()
	# Create "variant_regions" collection from gnomAD variants data
	#calculate_all_transcripts_regions(db)

	# Create GeVIR scores, requires "variant_regions" collection
	# If regions collection was created with INCLUDE_GNOMAD_OUTLIERS=True,
	# only this method has to be rerun with INCLUDE_GNOMAD_OUTLIERS=False 
	# to create gene scores for a dataset without outliers
	# Creates Supplementary Tables S6, S7
	#create_gevir_scores(db)


if __name__ == "__main__":
	sys.exit(main())