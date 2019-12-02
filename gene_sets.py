import os
import sys
import numpy as np
import math
from collections import OrderedDict
from common import MongoDB, is_int, sort_dict_by_values
from common import calculate_percentiles, is_float

INCLUDE_GNOMAD_OUTLIERS = True
EXCLUDE_XY_CHORMS = False

# NOTE: UNEECON scores were not evaluated in the original GeVIR manuscript,
# since UNEECON preprint was published at the late stage of our manuscript review.
# The code which analyse UNEECON scores is commented (ignored).

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
		self.gevir_no_gerp_score = ''
		self.gevir_no_gerp_percentile = ''

		self.gnomad_pli = ''
		self.gnomad_miss_z = ''
		self.gnomad_oe_mis_upper = ''
		self.gnomad_oe_lof_upper = ''

		self.combined_rank = ''
		#self.uneecon_g = '' 

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['_id'] = self.transcript_id
		dictionary['gene_name'] = self.gene_name
		dictionary['length'] = self.length
		dictionary['gevir_score'] = self.gevir_score
		dictionary['gevir_percentile'] = self.gevir_percentile
		dictionary['gevir_no_gerp_score'] = self.gevir_no_gerp_score
		dictionary['gevir_no_gerp_percentile'] = self.gevir_no_gerp_percentile

		dictionary['gnomad_pli'] = self.gnomad_pli
		dictionary['gnomad_miss_z'] = self.gnomad_miss_z
		dictionary['gnomad_oe_mis_upper'] = self.gnomad_oe_mis_upper
		dictionary['gnomad_oe_lof_upper'] = self.gnomad_oe_lof_upper

		dictionary['combined_rank'] = self.combined_rank

		#dictionary['uneecon_g'] = self.uneecon_g

		if self.gnomad_miss_z != '':
			dictionary['is_gnomad'] = True
		else:
			dictionary['is_gnomad'] = False

		'''
		if self.uneecon_g != '':
			dictionary['is_uneecon'] = True
		else:
			dictionary['is_uneecon'] = False
		'''

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
	gnomad_pli = {}
	gnomad_miss_z = {}
	gnomad_oe_mis_upper = {}
	gnomad_oe_lof_upper = {}


	if INCLUDE_GNOMAD_OUTLIERS:
		if EXCLUDE_XY_CHORMS:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "$and": [ { "chrom": { "$ne": "X" } }, { "chrom": { "$ne": "Y" } } ] })
		else:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True })
	else:
		if EXCLUDE_XY_CHORMS:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True, "$and": [ { "chrom": { "$ne": "X" } }, { "chrom": { "$ne": "Y" } } ] })
		else:
			gnomad_genes = db.gevir.gnomad_scores.find({ "canonical": True, "valid_transcript": True, "no_issues": True })

	for gnomad_gene in gnomad_genes:
		transcript_id = gnomad_gene['_id']
		
		gnomad_pli[transcript_id] = gnomad_gene['pLI']
		gnomad_miss_z[transcript_id] = gnomad_gene['mis_z']
		gnomad_oe_mis_upper[transcript_id] = gnomad_gene['oe_mis_upper']
		gnomad_oe_lof_upper[transcript_id] = gnomad_gene['oe_lof_upper']

	# Missense scores contain 1 NaN value; oe_lof (LOEUF) and pLI contains 482 genes with NA values.
	# We replace this values with the least siginificant values (lowest for pLI and Missense z-score, highers for LOEUF and MOEUF)
	# Just in case for future run this analysis on oe missense (MOEUF) too
	gnomad_pli = replace_nan_or_string_values_with_max_in_dictionary(gnomad_pli, reverse=True)
	gnomad_miss_z = replace_nan_or_string_values_with_max_in_dictionary(gnomad_miss_z, reverse=True)
	gnomad_oe_mis_upper = replace_nan_or_string_values_with_max_in_dictionary(gnomad_oe_mis_upper, reverse=False)
	gnomad_oe_lof_upper = replace_nan_or_string_values_with_max_in_dictionary(gnomad_oe_lof_upper, reverse=False)

	print 'Reading gnomAD GeVIR data...'
	gevir = {}
	gevir_genes = db.gevir.gevir_scores.find({})
	for gnomad_gene in gevir_genes:
		if EXCLUDE_XY_CHORMS and (gnomad_gene['chrom'] == 'X' or gnomad_gene['chrom'] == 'Y'):
			continue
		transcript_id = gnomad_gene['_id']
		score = gnomad_gene['gevir_score']
		gevir[transcript_id] = score

	sorted_gevir = sort_dict_by_values(gevir, reverse=True)
	gevir_percentiles = calculate_percentiles(sorted_gevir.keys(), reverse=False)

	print 'Reading gnomAD GeVIR NO GERP++ data...'
	gevir_no_gerp = {}
	gevir_no_gerp_genes = db.gevir.gevir_no_gerp_scores.find({})
	for gnomad_gene in gevir_no_gerp_genes:
		if EXCLUDE_XY_CHORMS and (gnomad_gene['chrom'] == 'X' or gnomad_gene['chrom'] == 'Y'):
			continue
		transcript_id = gnomad_gene['_id']
		score = gnomad_gene['gevir_score']
		gevir_no_gerp[transcript_id] = score

	sorted_gevir_no_gerp = sort_dict_by_values(gevir_no_gerp, reverse=True)
	gevir_no_gerp_percentiles = calculate_percentiles(sorted_gevir_no_gerp.keys(), reverse=False)

	'''
	print 'Reading Uneecon data...'
	uneecon_g = {}
	uneecon_genes = db.gevir.uneecon_genes.find({})
	for uneecon_gene in uneecon_genes:
		transcript_id = uneecon_gene['_id']
		score = uneecon_gene['uneecon_g']
		uneecon_g[transcript_id] = score		

	#uneecon_g = sort_dict_by_values(uneecon_g, reverse=True)
	'''
	print 'Merging and importing data...'

	any_transcripts = set(gevir.keys() + gnomad_miss_z.keys()) # + uneecon_g.keys()) # TODO remove uneecon_g !!!

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
	'''
	for transcript_id, score in gnomad_pli.iteritems():
		genes_scores[transcript_id].gnomad_pli = score
	'''
	for transcript_id, score in gnomad_miss_z.iteritems():
		genes_scores[transcript_id].gnomad_miss_z = score

	for transcript_id, score in gevir.iteritems():
		genes_scores[transcript_id].gevir_score = score
		genes_scores[transcript_id].gevir_percentile = gevir_percentiles[transcript_id]

	for transcript_id, score in gevir_no_gerp.iteritems():
		genes_scores[transcript_id].gevir_no_gerp_score = score
		genes_scores[transcript_id].gevir_no_gerp_percentile = gevir_no_gerp_percentiles[transcript_id]

	for transcript_id, score in gnomad_oe_mis_upper.iteritems():
		genes_scores[transcript_id].gnomad_oe_mis_upper = score

	for transcript_id, score in gnomad_oe_lof_upper.iteritems():
		genes_scores[transcript_id].gnomad_oe_lof_upper = score

	#for transcript_id, score in uneecon_g.iteritems():
	#	genes_scores[transcript_id].uneecon_g = score

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
		self.gevir_percentile_dict = {}
		self.gevir_no_gerp_dict = {}
		self.gevir_no_gerp_percentile_dict = {}
		self.gnomad_pli_dict = {}
		self.gnomad_miss_z_dict = {}
		self.gnomad_oe_mis_upper_dict = {}
		self.gnomad_oe_lof_upper_dict = {}

		self.combined_rank_dict = {}
		#self.uneecon_g_dict = {}

		genes_scores = db.gevir.common_gene_scores.find(filters)

		for gene_scores in genes_scores:
			transcript_id = gene_scores['_id']
			self.all_transcripts.add(transcript_id)
			self.length_dict[transcript_id] = gene_scores['length']
			self.gevir_dict[transcript_id] = gene_scores['gevir_score']
			self.gevir_percentile_dict[transcript_id] = gene_scores['gevir_percentile']
			self.gevir_no_gerp_dict[transcript_id] = gene_scores['gevir_no_gerp_score']
			self.gevir_no_gerp_percentile_dict[transcript_id] = gene_scores['gevir_no_gerp_percentile']
			self.gnomad_pli_dict[transcript_id] = gene_scores['gnomad_pli']
			self.gnomad_miss_z_dict[transcript_id] = gene_scores['gnomad_miss_z']
			self.gnomad_oe_mis_upper_dict[transcript_id] = gene_scores['gnomad_oe_mis_upper']
			self.gnomad_oe_lof_upper_dict[transcript_id] = gene_scores['gnomad_oe_lof_upper']
			self.combined_rank_dict[transcript_id] = gene_scores['combined_rank']
			#self.uneecon_g_dict[transcript_id] = gene_scores['uneecon_g']

		self.gevir_dict = sort_dict_by_values(self.gevir_dict, reverse=True)
		self.gevir_percentile_dict = sort_dict_by_values(self.gevir_percentile_dict, reverse=True)
		self.gevir_no_gerp_dict = sort_dict_by_values(self.gevir_no_gerp_dict, reverse=True)
		self.gevir_no_gerp_percentile_dict = sort_dict_by_values(self.gevir_no_gerp_percentile_dict, reverse=True)
		self.gnomad_pli_dict = sort_dict_by_values(self.gnomad_pli_dict, reverse=True)
		self.gnomad_miss_z_dict = sort_dict_by_values(self.gnomad_miss_z_dict, reverse=True)
		self.gnomad_oe_mis_upper_dict = sort_dict_by_values(self.gnomad_oe_mis_upper_dict, reverse=False)
		self.gnomad_oe_lof_upper_dict = sort_dict_by_values(self.gnomad_oe_lof_upper_dict, reverse=False)
		self.combined_rank_dict = sort_dict_by_values(self.combined_rank_dict, reverse=False)
		#self.uneecon_g_dict = sort_dict_by_values(self.uneecon_g_dict, reverse=True)

		self.gnomad_pli_percentile_dict = calculate_percentiles(self.gnomad_pli_dict.keys())
		self.gnomad_miss_z_percentile_dict = calculate_percentiles(self.gnomad_miss_z_dict.keys())
		self.gnomad_oe_mis_upper_percentile_dict = calculate_percentiles(self.gnomad_oe_mis_upper_dict.keys())
		self.gnomad_oe_lof_upper_percentile_dict = calculate_percentiles(self.gnomad_oe_lof_upper_dict.keys())
		self.combined_rank_percentile_dict = calculate_percentiles(self.combined_rank_dict.keys())
		#self.uneecon_g_percentile_dict = calculate_percentiles(self.uneecon_g_dict.keys())

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
	db = MongoDB()
	create_common_gene_scores(db)

if __name__ == "__main__":
	sys.exit(main())