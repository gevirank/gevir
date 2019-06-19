import os
import sys
import numpy as np
import math
from collections import OrderedDict
from common import MongoDB, is_int, sort_dict_by_values
from common import calculate_percentiles, is_float

INCLUDE_GNOMAD_OUTLIERS = True

########################
### DATA PREPARATION ###
########################

class GeneScores():
	"""Document in "common_gene_scores" collection, used to store all gene scores"""
	def __init__(self):
		self.transcript_id = ''
		self.gene_name = ''
		self.length = 0
		self.gevir_score = ''
		self.gevir_percentile = ''

		self.gnomad_miss_z = ''
		self.gnomad_oe_mis_upper = ''
		self.gnomad_oe_lof_upper = ''

		self.combined_rank = ''

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['_id'] = self.transcript_id
		dictionary['gene_name'] = self.gene_name
		dictionary['length'] = self.length
		dictionary['gevir_score'] = self.gevir_score
		dictionary['gevir_percentile'] = self.gevir_percentile

		dictionary['gnomad_miss_z'] = self.gnomad_miss_z
		dictionary['gnomad_oe_mis_upper'] = self.gnomad_oe_mis_upper
		dictionary['gnomad_oe_lof_upper'] = self.gnomad_oe_lof_upper

		dictionary['combined_rank'] = self.combined_rank

		if self.gnomad_miss_z != '':
			dictionary['is_gnomad'] = True
		else:
			dictionary['is_gnomad'] = False

		return dictionary


def replace_nan_or_string_values_with_max_in_dictionary(dictionary, reverse=False):
	non_numeric_keys = []
	numeric_values = []
	for key, value in dictionary.iteritems():
		if is_float(value):
			if math.isnan(value):
				non_numeric_keys.append(key)
			else:
				numeric_values.append(value)
		else:
			non_numeric_keys.append(key)

	non_numeric_values_replacer = 0
	if reverse:
		non_numeric_values_replacer = min(numeric_values) - 1
	else:
		non_numeric_values_replacer = max(numeric_values) + 1		

	for non_numeric_key in non_numeric_keys:
		dictionary[non_numeric_key] = non_numeric_values_replacer

	return dictionary


def create_common_gene_scores(db):
	print 'Reading gnomAD Constraints data...'
	gnomad_miss_z = {}
	gnomad_oe_mis_upper = {}
	gnomad_oe_lof_upper = {}

	if INCLUDE_GNOMAD_OUTLIERS:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True }) # "no_issues": True
	else:
		gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	for gnomad_gene in gnomad_genes:
		transcript_id = gnomad_gene['_id']
		
		gnomad_miss_z[transcript_id] = gnomad_gene['mis_z']
		gnomad_oe_mis_upper[transcript_id] = gnomad_gene['oe_mis_upper']
		gnomad_oe_lof_upper[transcript_id] = gnomad_gene['oe_lof_upper']

	# Missense scores contain 1 NaN value; oe_lof (LOEUF) contains 482 genes with NA values.
	# We replace this values with the least siginificant values (lowest for Missense z-score, highers for LOEUF)
	# Just in case for future run this analysis on oe missense (MOEUF) too
	gnomad_miss_z = replace_nan_or_string_values_with_max_in_dictionary(gnomad_miss_z, reverse=True)
	gnomad_oe_mis_upper = replace_nan_or_string_values_with_max_in_dictionary(gnomad_oe_mis_upper, reverse=False)
	gnomad_oe_lof_upper = replace_nan_or_string_values_with_max_in_dictionary(gnomad_oe_lof_upper, reverse=False)

	print 'Reading gnomAD GeVIR data...'
	gevir = {}
	gevir_genes = db.gevir.gevir_scores.find({})
	for gnomad_gene in gevir_genes:
		transcript_id = gnomad_gene['_id']
		score = gnomad_gene['gevir_score']
		gevir[transcript_id] = score

	sorted_gevir = sort_dict_by_values(gevir, reverse=True)
	gevir_percentiles = calculate_percentiles(sorted_gevir.keys(), reverse=False)

	print 'Merging and importing data...'

	any_transcripts = set(gevir.keys() + gnomad_miss_z.keys())

	genes_scores = {}
	for transcript_id in any_transcripts:
		gene = db.exac.genes.find_one({'canonical_transcript': transcript_id})
		gene_name = gene['gene_name']
		gene_protein_sequence = db.gevir.ens_aa_fasta.find_one({'_id': transcript_id})
		length = len(gene_protein_sequence['cds']) - 1 # Note: -1 to exclude stop codon

		gene_scores = GeneScores()
		gene_scores.transcript_id = transcript_id
		gene_scores.gene_name = gene_name
		gene_scores.length = length

		genes_scores[transcript_id] = gene_scores

	for transcript_id, score in gnomad_miss_z.iteritems():
		genes_scores[transcript_id].gnomad_miss_z = score

	for transcript_id, score in gevir.iteritems():
		genes_scores[transcript_id].gevir_score = score
		genes_scores[transcript_id].gevir_percentile = gevir_percentiles[transcript_id]

	for transcript_id, score in gnomad_oe_mis_upper.iteritems():
		genes_scores[transcript_id].gnomad_oe_mis_upper = score

	for transcript_id, score in gnomad_oe_lof_upper.iteritems():
		genes_scores[transcript_id].gnomad_oe_lof_upper = score

	# Combined rank is a resorted sum of gnomAD lof observe/expected upper confidence interval and GeVIR ranks.

	common_transcripts = set(gnomad_oe_lof_upper.keys()) & set(gevir.keys())
	
	# Note that intolerant genes should have LOW observed/expected score, but HIGH GeVIR score.
	sorted_gnomad_oe_lof_upper = sort_dict_by_values(gnomad_oe_lof_upper, reverse=False)
	sorted_gevir = sort_dict_by_values(gevir, reverse=True)

	ranked_gevir = OrderedDict()
	x = 1
	for transcript_id in sorted_gevir.keys():
		if transcript_id in common_transcripts:
			ranked_gevir[transcript_id] = x
			x += 1

	ranked_oe_lof_upper = OrderedDict()
	x = 1
	for transcript_id in sorted_gnomad_oe_lof_upper.keys():
		if transcript_id in common_transcripts:
			ranked_oe_lof_upper[transcript_id] = x
			x += 1

	combined_ranks = OrderedDict()
	for transcript_id, gevir_rank in ranked_gevir.iteritems():
		combined_ranks[transcript_id] = gevir_rank + ranked_oe_lof_upper[transcript_id]

	combined_ranks = sort_dict_by_values(combined_ranks, reverse=False)

	x = 1
	for transcript_id in combined_ranks:
		genes_scores[transcript_id].combined_rank = x
		x += 1

	db.gevir.common_gene_scores.drop()
	for gene_scores in genes_scores.values():
		db.gevir.common_gene_scores.insert(gene_scores.get_dictionary())


class ScoreSets():
	def __init__(self, db, filters={}):
		self.all_transcripts = set([])
		self.length_dict = {}
		self.gevir_dict = {}

		self.gnomad_miss_z_dict = {}
		self.gnomad_oe_mis_upper_dict = {}
		self.gnomad_oe_lof_upper_dict = {}

		self.combined_rank_dict = {}

		genes_scores = db.gevir.common_gene_scores.find(filters)

		for gene_scores in genes_scores:
			transcript_id = gene_scores['_id']
			self.all_transcripts.add(transcript_id)
			self.length_dict[transcript_id] = gene_scores['length']
			self.gevir_dict[transcript_id] = gene_scores['gevir_score']
			self.gnomad_miss_z_dict[transcript_id] = gene_scores['gnomad_miss_z']
			self.gnomad_oe_mis_upper_dict[transcript_id] = gene_scores['gnomad_oe_mis_upper']
			self.gnomad_oe_lof_upper_dict[transcript_id] = gene_scores['gnomad_oe_lof_upper']
			self.combined_rank_dict[transcript_id] = gene_scores['combined_rank']

		self.gevir_dict = sort_dict_by_values(self.gevir_dict, reverse=True)
		self.gnomad_miss_z_dict = sort_dict_by_values(self.gnomad_miss_z_dict, reverse=True)
		self.gnomad_oe_mis_upper_dict = sort_dict_by_values(self.gnomad_oe_mis_upper_dict, reverse=False)
		self.gnomad_oe_lof_upper_dict = sort_dict_by_values(self.gnomad_oe_lof_upper_dict, reverse=False)
		self.combined_rank_dict = sort_dict_by_values(self.combined_rank_dict, reverse=False)		

	def get_transcripts_median_length(self, transcript_ids):
		transcripts_length = []
		for transcript_id in transcript_ids:
			transcripts_length.append(self.length_dict[transcript_id])
		return np.median(transcripts_length)

	def get_transcripts_mean_length(self, transcript_ids):
		transcripts_length = []
		for transcript_id in transcript_ids:
			transcripts_length.append(self.length_dict[transcript_id])
		return np.mean(transcripts_length)

	def get_transcripts_length(self, transcript_ids):
		transcripts_length = []
		for transcript_id in transcript_ids:
			transcripts_length.append(self.length_dict[transcript_id])
		return transcripts_length


class OmimSets():
	def __init__(self, db):
		self.ad = set([])
		self.ar = set([])

		# OMIM
		self.ad_and_ar = set([])
		self.xld = set([])
		self.xlr = set([])
		self.xld_and_xlr = set([])

		omim_genes = db.gevir.omim.find({})

		for omim_gene in omim_genes:
			inheritance = omim_gene['inheritance']
			transcript_id = omim_gene['transcript_id']

			if inheritance == 'AD':
				self.ad.add(transcript_id)
			if inheritance == 'AR':
				self.ar.add(transcript_id)
			if 'AD' in inheritance and 'AR' in inheritance:
				self.ad_and_ar.add(transcript_id)

			if inheritance == 'XLD':
				self.xld.add(transcript_id)
			if inheritance == 'XLR':
				self.xlr.add(transcript_id)
			if 'XLD' in inheritance and 'XLR' in inheritance:
				self.xld_and_xlr.add(transcript_id)

	
class EssentialSets():
	def __init__(self, db, nulls_not_in_omim=True):
		self.haploid_essential = set([])
		self.mouse_essential = set([])
		self.nulls = set([])

		# Nulls
		lof_tolerant_genes = db.gevir.mac_arthur_gene_lists.find_one({'_id': 'lof_tolerant'})
		self.nulls = set(lof_tolerant_genes['transcript_ids'])

		mouse_het_lethal = db.gevir.mouse_het_lethal.find_one({'_id': 'transcript_ids'})
		self.mouse_het_lethal= set(mouse_het_lethal['data'])
		crispr_essential = db.gevir.mac_arthur_gene_lists.find_one({'_id': 'crispr_essential'})
		self.crispr_essential = set(crispr_essential['transcript_ids'])
		crispr_non_essential = db.gevir.mac_arthur_gene_lists.find_one({'_id': 'crispr_non_essential'})
		self.crispr_non_essential = set(crispr_non_essential['transcript_ids'])


def main():
	# Fix somehow to be able to run all figure methods one by one without corrupting
	db = MongoDB()
	create_common_gene_scores(db)

if __name__ == "__main__":
	sys.exit(main())