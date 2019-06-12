import os
import sys
import progressbar
import pymongo
import csv
import pronto
import json
from collections import OrderedDict
from intermine.webservice import Service
from common import MongoDB, get_gene_canonical_transcript_by_name, write_table_to_csv, file_len
from csv_reader import CsvReader
from gnomad_utils import worst_csq_from_list, get_minimal_representation

SOURCE_DATA_FOLDER = './source_data/'

# The GERP++ scores
# http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_scores.tar.gz
GERP_FOLDER = SOURCE_DATA_FOLDER + 'hg19.GERP_scores/' 

# The Ensembl coding and peptide sequences
# http://grch37.ensembl.org/biomart/martview/e19a3419814f357a70908603f15f7bae?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.ensembl_gene_id|hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.coding&FILTERS=&VISIBLEPANEL=resultspanel
ENS_CDS_FASTA = SOURCE_DATA_FOLDER + 'ens_cds_fasta.txt'
# http://grch37.ensembl.org/biomart/martview/e19a3419814f357a70908603f15f7bae?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.ensembl_gene_id|hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.peptide&FILTERS=&VISIBLEPANEL=resultspanel
ENS_AA_FASTA = SOURCE_DATA_FOLDER + 'ens_protein_fasta.txt'

# The gnomAD constrained metrics
# https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
GNOMAD_SCORES_TSV = SOURCE_DATA_FOLDER + 'gnomad_constraint_original.tsv' # gnomad_constraint.tsv gnomad.v2.1.1.lof_metrics.by_transcript.txt

# The OMIM genemap2.txt file can be found, !!!AFTER REGISTRATION!!!, at
# https://omim.org/downloads/
# IMPORTANT: Comment lines at the start and end of the file have to be removed, 
# except column names ('# Chromosome' has to be renamed to 'Chromosome').
OMIM_TSV = SOURCE_DATA_FOLDER + 'genemap2.txt'

# The ClinVar datasets
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
CLIN_VARS_TSV = SOURCE_DATA_FOLDER + '/clin_var/' + 'variant_summary.txt'
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt
CLIN_VAR_CITES_TSV = SOURCE_DATA_FOLDER + '/clin_var/' + 'var_citations.txt'

# The Constraned Coding Regions (CCRs) datasets
# https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz
CCRS_AUTOSOMAL_TSV = SOURCE_DATA_FOLDER + 'ccrs.autosomes.v2.20180420.bed'
# https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz
CCRS_X_TSV = SOURCE_DATA_FOLDER + 'ccrs.xchrom.v2.20180420.bed'

# The cell essential, non-essential and Loss-of-Function tolerant (i.e. Null) datasets obtained from MacArthur repository
# https://github.com/macarthur-lab/gene_lists/blob/master/lists/homozygous_lof_tolerant_twohit.tsv
MAC_ARTHUR_LOF_TOLERANT = SOURCE_DATA_FOLDER + 'ma_lof_hom_tolerant.tsv'
# https://github.com/macarthur-lab/gene_lists/blob/master/lists/CEGv2_subset_universe.tsv
MAC_ARTHUR_CRISPR_ESSENTIAL = SOURCE_DATA_FOLDER + 'ma_crispr_essential.tsv'
# https://github.com/macarthur-lab/gene_lists/blob/master/lists/NEGv1_subset_universe.tsv
MAC_ARTHUR_CRISPR_NON_ESSENTIAL = SOURCE_DATA_FOLDER + 'ma_crispr_non_essential.tsv'

# The human-mouse ortholog mapping file can be found at 
# http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
MOUSE_TO_HUMAN_TSV = SOURCE_DATA_FOLDER + 'mouse_to_human.tsv'

'''
The mouse heterozygous lethal genes can be obtained from http://mousemine.org/ by quiring the database with following constraints:
path="OntologyAnnotation.ontologyTerm" type="MPTerm";
path="OntologyAnnotation.subject" type="SequenceFeature";
path="OntologyAnnotation.evidence.baseAnnotations.subject" type="Genotype";
path="OntologyAnnotation.evidence.baseAnnotations.subject.zygosity" op="=" value="ht" code="B";
path="OntologyAnnotation.ontologyTerm.name" op="CONTAINS" value="lethal"
Alternatively, check mousemine mehod to generate the file.
'''
MOUSE_HET_LETHAL_TSV = SOURCE_DATA_FOLDER + 'mouse_het_lethal_original.tsv'


############
### GERP ###
############

def import_chrom(db, file_name):
	input_file = open(GERP_FOLDER + file_name, 'rt')
	reader = csv.reader(input_file, delimiter='\t')

	pos = 1
	chrom = file_name.split('.', 1)[0][3:]

	db.gerp[chrom].drop()

	bulk = db.gerp[chrom].initialize_unordered_bulk_op()

	total_lines = len(open(GERP_FOLDER + file_name).readlines())
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for row in reader:
		bulk.insert({'_id': pos, 'RS': float(row[1])})
		if (pos % 10000 == 0):
			bulk.execute()
			bulk = db.gerp[chrom].initialize_unordered_bulk_op()		
		pos += 1
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	if (pos % 10000 != 0):
		bulk.execute()


def import_gerp(db):
	for file_name in os.listdir(GERP_FOLDER):
		print file_name
		import_chrom(db, file_name)


#########################
### Ensembl CDS FASTA ###
#########################

def import_ens_cds_fasta(db, fasta_file, collection_name):
	db.gevir[collection_name].drop()
	with open(fasta_file) as f:
		rows = f.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	rows = [x.strip() for x in rows]

	transcripts = []
	cds = ''
	transcript_id = ''
	for row in rows:
		if row and row[0] == '>':
			gene_id, new_transcript_id = row.split('|')
			if not transcript_id:
				transcript_id = new_transcript_id
			if cds:
				transcripts.append({'_id' : transcript_id,'cds' : cds})
				transcript_id = new_transcript_id
				cds = ''
		else:
			cds += row

	db.gevir[collection_name].insert_many(transcripts)


#####################
### gnomAD Scores ###
#####################

# Note: gnomAD scores were updated since original GeVIR analysis (22/10/18).
#       New file contains 'constraint_flag' column instead of 'gene_issues'
#       with slightly different format.
#       Overall the trends remains the same, last checked with scores at 11/06/19.
def import_gnomad_scores(db, new_gnomad_file=True):
	gnomad_scores = CsvReader(GNOMAD_SCORES_TSV, delimiter='\t')

	int_keys = ['obs_lof', 'obs_mis', 'obs_syn']
	float_keys = ['exp_lof', 'oe_lof', 'oe_lof_lower', 'oe_lof_upper', 'exp_mis', 'oe_mis',
				  'oe_mis_lower', 'oe_mis_upper', 'exp_syn', 'oe_syn', 'oe_syn_lower', 
				  'oe_syn_upper', 'lof_z', 'mis_z', 'syn_z', 'pLI', 'pRec', 'pNull']

	for document in gnomad_scores.data:
		document['_id'] = document['transcript']

		if new_gnomad_file:
			gene_issues = document['constraint_flag']
		else:
			gene_issues = document['gene_issues']
		
		gene_issues_list = []

		for int_key in int_keys:
			if document[int_key] != 'NA':
				document[int_key] = int(document[int_key])

		for float_key in float_keys:
			if document[float_key] != 'NA':
				document[float_key] = float(document[float_key])

		if new_gnomad_file:
			if gene_issues != '':
				gene_issues = gene_issues.split('|')
				for gene_issue in gene_issues:
					gene_issues_list.append(gene_issue)
		else:
			if gene_issues != '[]':
				gene_issues = gene_issues[1:-1].split(',')
				for gene_issue in gene_issues:
					gene_issues_list.append(gene_issue[1:-1])

		if document['canonical'] == 'true':
			document['canonical'] = True
		else:
			document['canonical'] = False

		document['gene_issues'] = gene_issues_list

		if gene_issues_list:
			document['no_issues'] = False
		else:
			document['no_issues'] = True

		cds = db.gevir.ens_cds_fasta.find_one({'_id': document['transcript']})
		cds = cds['cds']

		valid_transcript = True
		transcript_issues = []

		if len(cds) % 3 != 0:
			valid_transcript = False
			transcript_issues.append('not_divisible_by_3')

		stop_codons = set(['TAA', 'TAG', 'TGA'])
		if cds[-3:] not in stop_codons:
			valid_transcript = False
			transcript_issues.append('stop_missing')

		if cds[:3] != 'ATG':
			valid_transcript = False
			transcript_issues.append('start_missing')

		document['valid_transcript'] = valid_transcript
		document['transcript_issues'] = transcript_issues

		if not cds:
			print 'not found', document['transcript']

		exac_transcript = db.exac.transcripts.find_one({'transcript_id': document['transcript']})
		document['chrom'] = exac_transcript['chrom']

	gnomad_scores.import_to_db(db.gevir, 'gnomad_scores')
	db.gevir.gnomad_scores.create_index([('gene', pymongo.ASCENDING)], name='gene_1')
	db.gevir.gnomad_scores.create_index([('canonical', pymongo.ASCENDING)], name='canonical_1')


############
### OMIM ###
############

class OmimGene():
	def __init__(self):
		self.gene_name = ''
		self.gene_id = ''
		self.transcript_id = ''
		self.mim_number = ''
		self.build_persistent = False
		self.gene_name_38 = ''
		self.gene_id_38 = ''
		self.chrom = ''
		self.start = 0
		self.stop = 0
		self.inheritance = ''
		self.susceptibility = False
		self.phenotypes_ids = []
		self.phenotypes_inheritance = {}
		self.phenotypes = []

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['gene_name'] = self.gene_name
		dictionary['gene_id'] = self.gene_id
		dictionary['transcript_id'] = self.transcript_id
		dictionary['mim_number'] = self.mim_number
		dictionary['build_persistent'] = self.build_persistent
		dictionary['gene_name_38'] = self.gene_name_38
		dictionary['gene_id_38'] = self.gene_id_38
		dictionary['chrom'] = self.chrom
		dictionary['start'] = self.start
		dictionary['stop'] = self.stop
		dictionary['inheritance'] = self.inheritance
		dictionary['susceptibility'] = self.susceptibility
		dictionary['phenotypes_ids'] = self.phenotypes_ids
		dictionary['phenotypes_inheritance'] = self.phenotypes_inheritance
		dictionary['phenotypes'] = self.phenotypes
		return dictionary


def import_omim(db):
	db.gevir.omim.drop()

	omim_genes = CsvReader(OMIM_TSV, delimiter='\t')

	total_lines = len(omim_genes.data)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for document in omim_genes.data:
		phenotypes = document['Phenotypes']
		if phenotypes == '':
			line_number += 1
			continue

		inheritance = set([])
		phenotypes_ids = []
		phenotypes_inheritance = {}
		phenotypes = phenotypes.split(';')
		susceptibility = False
		for phenotype in phenotypes:
			
			# OMIM FAQ: A question mark, "?", before the disease name indicates an unconfirmed or possibly spurious mapping.
			if '?' in phenotype:
				continue

			# OMIM FAQ: Brackets, "[ ]", indicate "nondiseases," mainly genetic variations that lead to 
			# apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
			if '[' in phenotype and ']' in phenotype:
				continue

			
			# OMIM FAQ: Braces, "{ }", indicate mutations that contribute to susceptibility to 
			# multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).
			if '{' in phenotype and '}' in phenotype and '(3)' in phenotype:
				susceptibility = True
				#continue

			phenotype = phenotype.strip()

			# 3 is a confirmed phenotype
			if '(3)' in phenotype:
				phenotype_id = ''
				phen_strs = phenotype.split(', ')
				for phen_sub_str in phen_strs:
					if '(3)' in phen_sub_str:
						phenotype_id_str = phen_sub_str.split(' ')[0]
						try:
							int(phenotype_id_str)
							phenotype_id = phenotype_id_str
							phenotypes_ids.append(phenotype_id)
						except ValueError:
							pass


				known_inheritance = False
				# Make sure that all characters are in the same case
				lower_phenotype = phenotype.lower()
				# Finally check inheritance pattern.
				# There might be multiple diseases associated with the same gene

				if '{' in phenotype and '}' in phenotype:
					inheritance.add('')
				else:
					if 'autosomal recessive' in lower_phenotype:
						inheritance.add('AR')
						phenotypes_inheritance[phenotype_id] = 'AR'
						known_inheritance = True
					if 'x-linked recessive' in lower_phenotype:
						inheritance.add('XLR')
						phenotypes_inheritance[phenotype_id] = 'XLR'
						known_inheritance = True
					if 'autosomal dominant' in lower_phenotype:
						inheritance.add('AD')
						phenotypes_inheritance[phenotype_id] = 'AD'
						known_inheritance = True
					if 'x-linked dominant' in lower_phenotype:
						inheritance.add('XLD')
						phenotypes_inheritance[phenotype_id] = 'XLD'
						known_inheritance = True
					if not known_inheritance:
						inheritance.add('unknown')
						phenotypes_inheritance[phenotype_id] = 'unknown'
					

		if len(inheritance) > 0 or susceptibility:
			omim_gene = OmimGene()
			omim_gene.mim_number = document['Mim Number']
			omim_gene.gene_name_38 = document['Approved Symbol'].upper()
			omim_gene.gene_id_38 = document['Ensembl Gene ID']
			omim_gene.chrom = document['Chromosome'][3:]
			omim_gene.phenotypes = phenotypes
			if len(inheritance) > 1 and '' in inheritance:
				inheritance.remove('')
			inheritance = list(inheritance)
			inheritance.sort()
			omim_gene.inheritance = ', '.join(inheritance)
			omim_gene.susceptibility = susceptibility
			omim_gene.phenotypes_ids = phenotypes_ids
			omim_gene.phenotypes_inheritance = phenotypes_inheritance
			omim_gene.build_persistent = False

			insert_to_db = False
			exac_gene_by_id = db.exac.genes.find_one({'gene_id': omim_gene.gene_id_38, 'chrom': omim_gene.chrom})
			exac_gene_by_name = db.exac.genes.find_one({'gene_name_upper': omim_gene.gene_name_38, 'chrom': omim_gene.chrom})
			if exac_gene_by_id and exac_gene_by_name and exac_gene_by_id['gene_id'] == exac_gene_by_name['gene_id']:
				if 'canonical_transcript' in exac_gene_by_id:
					omim_gene.transcript_id = exac_gene_by_id['canonical_transcript']
					omim_gene.gene_name = omim_gene.gene_name_38
					omim_gene.gene_id = omim_gene.gene_id_38
					omim_gene.start = exac_gene_by_id['start']
					omim_gene.stop = exac_gene_by_id['stop']
					omim_gene.build_persistent = True
					insert_to_db = True
			elif exac_gene_by_id:
				if 'canonical_transcript' in exac_gene_by_id:
					omim_gene.transcript_id = exac_gene_by_id['canonical_transcript']
					omim_gene.gene_name = exac_gene_by_id['gene_name_upper']
					omim_gene.gene_id = omim_gene.gene_id_38
					omim_gene.start = exac_gene_by_id['start']
					omim_gene.stop = exac_gene_by_id['stop']
					insert_to_db = True
			elif exac_gene_by_name:
				if 'canonical_transcript' in exac_gene_by_name:
					omim_gene.transcript_id = exac_gene_by_name['canonical_transcript']
					omim_gene.gene_name = omim_gene.gene_name_38
					omim_gene.gene_id = exac_gene_by_name['gene_id']
					omim_gene.start = exac_gene_by_name['start']
					omim_gene.stop = exac_gene_by_name['stop']
					insert_to_db = True
			if insert_to_db:
				db.gevir.omim.insert(omim_gene.get_dictionary())

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	db.gevir.omim.create_index([('gene_name', pymongo.ASCENDING)], name='gene_name_1')
	db.gevir.omim.create_index([('gene_id', pymongo.ASCENDING)], name='gene_id_1')
	db.gevir.omim.create_index([('transcript_id', pymongo.ASCENDING)], name='transcript_id_1')


###############
### ClinVar ###
###############

class ClinVar():
	def __init__(self):
		self.allele_id = ''
		self.rsid = ''
		self.var_class = ''
		self.chrom = ''
		self.start = 0
		self.stop = 0
		self.ref = ''
		self.alt = ''
		self.gene_name = ''
		self.clin_sig = ''
		self.clin_sig_num = -1
		self.submitters_num = 0
		self.submitters_categories = 0
		self.review_status = ''
		self.phenotype_ids = OrderedDict()
		self.phenotype_list = []
		self.cites = []


	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['allele_id'] = self.allele_id
		dictionary['rsid'] = self.rsid
		dictionary['var_class'] = self.var_class
		dictionary['chrom'] = self.chrom
		dictionary['start'] = self.start
		dictionary['stop'] = self.stop
		dictionary['ref'] = self.ref
		dictionary['alt'] = self.alt
		dictionary['gene_name'] = self.gene_name
		dictionary['clin_sig'] = self.clin_sig
		dictionary['clin_sig_num'] = self.clin_sig_num
		dictionary['submitters_num'] = self.submitters_num
		dictionary['submitters_categories'] = self.submitters_categories
		dictionary['review_status'] = self.review_status
		dictionary['phenotype_ids'] = self.phenotype_ids
		dictionary['phenotype_list'] = self.phenotype_list
		dictionary['cites'] = self.cites
		dictionary['cites_count'] = len(self.cites)
		return dictionary


def read_clin_var_cites():
	clin_var_cites = {}
	clin_vars = CsvReader(CLIN_VAR_CITES_TSV, delimiter='\t')

	total_lines = len(clin_vars.data)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for document in clin_vars.data:
		allele_id = document['#AlleleID']
		citation_id = document['citation_id']

		if allele_id in clin_var_cites:
			clin_var_cites[allele_id].add(citation_id)
		else:
			clin_var_cites[allele_id] = set([citation_id])

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()
	return clin_var_cites


def import_clin_var_raw(db):
	clin_var_cites = read_clin_var_cites()

	db.gevir.clin_vars.drop()

	clin_vars = CsvReader(CLIN_VARS_TSV, delimiter='\t')

	total_lines = len(clin_vars.data)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for document in clin_vars.data:
		if document['Assembly'] != 'GRCh37':
			continue
			
		clin_var = ClinVar()
		clin_var.allele_id = document['#AlleleID']
		clin_var.rsid = document['RS# (dbSNP)']
		if clin_var.rsid != '':
			clin_var.rsid = 'rs' + str(clin_var.rsid)

		clin_var.chrom = document['Chromosome']
		clin_var.start = int(document['Start'])
		clin_var.stop = int(document['Stop'])
		if clin_var.start == clin_var.stop:
			clin_var.var_class = 'SNV'
		else:
			clin_var.var_class = 'INDEL'

		clin_var.ref = document['ReferenceAllele']
		clin_var.alt = document['AlternateAllele']
		clin_var.gene_name = document['GeneSymbol']
		clin_var.clin_sig = document['ClinicalSignificance']
		clin_var.clin_sig_num = int(document['ClinSigSimple'])


		clin_var.submitters_num = int(document['NumberSubmitters'])
		clin_var.submitters_categories = int(document['SubmitterCategories'])

		clin_var.review_status = document['ReviewStatus']

		phenotype_ids = document['PhenotypeIDS'].split(';')
		for phenotype_id in phenotype_ids:
			if ',' in phenotype_id:
				phenotypes = phenotype_id.split(',')
			else:
				phenotypes = [phenotype_id]

			for phenotype in phenotypes:
				if ':' in phenotype:
					key, value = phenotype.split(':', 1)
					if key not in clin_var.phenotype_ids:
						clin_var.phenotype_ids[key] = [value]
					else:
						clin_var.phenotype_ids[key].append(value)

		clin_var.phenotype_list = document['PhenotypeList'].split(';')

		if clin_var.allele_id in clin_var_cites:
			clin_var.cites = list(clin_var_cites[clin_var.allele_id])

		db.gevir.clin_vars.insert(clin_var.get_dictionary())
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	db.gevir.clin_vars.create_index([('chrom', pymongo.ASCENDING)], name='chrom_1')
	db.gevir.clin_vars.create_index([('start', pymongo.ASCENDING)], name='start_1')
	db.gevir.clin_vars.create_index([('stop', pymongo.ASCENDING)], name='stop_1')
	db.gevir.clin_vars.create_index([('ref', pymongo.ASCENDING)], name='ref_1')
	db.gevir.clin_vars.create_index([('alt', pymongo.ASCENDING)], name='alt_1')
	db.gevir.clin_vars.create_index([('allele_id', pymongo.ASCENDING)], name='allele_id_1')
	db.gevir.clin_vars.create_index([('gene_name', pymongo.ASCENDING)], name='gene_name_1')


class EnsemblVar():
	def __init__(self):
		self.chrom = '.'
		self.start = '.'
		self.end = '.'
		self.allele = '.'
		self.strand = '.'
		self.id = '.'

	def get_row(self):
		row = [self.chrom, self.start, self.end, self.allele, self.strand, self.id]
		return row


def export_clin_var_ensembl(db):
	ens_var_data = []
	clin_vars = db.gevir.clin_vars.find({})

	total_lines = clin_vars.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for clin_var in clin_vars:
		ens_var = EnsemblVar()
		ens_var.chrom = clin_var['chrom']

		if clin_var['ref'] != '-':
			ens_var.start = clin_var['start']
			ens_var.end = clin_var['stop']			
		else:
			ens_var.start = clin_var['stop'] + 1
			ens_var.end = clin_var['stop']

		ens_var.allele = clin_var['ref'] + '/' + clin_var['alt']
		ens_var.strand = '+'
		ens_var.id = '_'.join([clin_var['chrom'], str(clin_var['start']), clin_var['ref'], clin_var['alt']])

		row = ens_var.get_row()
		ens_var_data.append(row)

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	output_csv = SOURCE_DATA_FOLDER + '/clin_var/clin_var_ensembl.txt'
	write_table_to_csv(ens_var_data, output_csv, delimiter=' ')


def import_clin_var_temp_json(db):
	db.gevir.clin_var_vep.drop()

	with open(SOURCE_DATA_FOLDER + '/clin_var/clin_var_vep.json', 'r') as f:
		table = [json.loads(line) for line in f]

	for row in table:
		db.gevir.clin_var_vep.insert(row)

	db.gevir.clin_var_vep.create_index([('id', pymongo.ASCENDING)], name='id_1')


def update_clin_vars_with_canonical_transcripts(db):
	gene_name_transcript = {}
	exac_genes = db.exac.genes.find({})
	for exac_gene in exac_genes:
		if 'canonical_transcript' in exac_gene:
			gene_name_transcript[exac_gene['gene_name']] = exac_gene['canonical_transcript']

	clin_vars = db.gevir.clin_vars.find({})

	total_lines = clin_vars.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for clin_var in clin_vars:
		transcripts = set([])
		gene_name = clin_var['gene_name']
		if ';' in gene_name:
			gene_names = gene_name.split(';')
			for gene_name in gene_names:
				if gene_name in gene_name_transcript:
					transcripts.add(gene_name_transcript[gene_name])
				else:
					exac_gene = db.exac.genes.find_one({'other_names': gene_name})
					if exac_gene and 'canonical_transcript' in exac_gene:
						transcripts.add(exac_gene['canonical_transcript'])

		if gene_name in gene_name_transcript:
			transcripts.add(gene_name_transcript[gene_name])
		else:
			exac_gene = db.exac.genes.find_one({'other_names': gene_name})
			if exac_gene and 'canonical_transcript' in exac_gene:
				transcripts.add(exac_gene['canonical_transcript'])

		db.gevir.clin_vars.update_one({ 'allele_id': clin_var['allele_id'] },  { '$set': {'canonical_transcripts': list(transcripts) } })
		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()


# IMPORTANT: genes without canonical transcripts in gnomAD were ignored
def update_clin_vars_with_csqs(db):
	clin_vars = db.gevir.clin_vars.find({})

	total_lines = clin_vars.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for clin_var in clin_vars:
		if 'canonical_transcripts' not in clin_var:
			continue

		clin_var_id = '_'.join([clin_var['chrom'], str(clin_var['start']), clin_var['ref'], clin_var['alt']])
		clin_var_vep = db.gevir.clin_var_vep.find_one({'id': clin_var_id})
		if clin_var_vep and 'transcript_consequences' in clin_var_vep:
			db.gevir.clin_vars.update_one({ 'allele_id': clin_var['allele_id'] },  { '$set': {'transcript_consequences': clin_var_vep['transcript_consequences']} })

			csqs = {}
			for vep in clin_var_vep['transcript_consequences']:
				if vep['transcript_id'] in clin_var['canonical_transcripts']:
					csqs[vep['transcript_id']] = worst_csq_from_list(vep['consequence_terms'])

		db.gevir.clin_vars.update_one({ 'allele_id': clin_var['allele_id'] },  { '$set': {'csqs': csqs} })

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()


class ClinVarTypes():
	def __init__(self):
		self.splice_acceptor_donor_variant = 0
		self.stop_gained = 0
		self.frameshift_variant = 0
		self.stop_lost = 0
		self.start_lost = 0
		self.inframe_insertion = 0
		self.inframe_deletion = 0
		self.missense_variant = 0
		self.splice_region_variant = 0
		self.synonymous_variant = 0

	def count_csq(self, csq):
		if csq == 'splice_acceptor_variant' or csq == 'splice_donor_variant':
			self.splice_acceptor_donor_variant += 1
		elif csq == 'stop_gained':
			self.stop_gained += 1
		elif csq == 'frameshift_variant':
			self.frameshift_variant += 1		
		elif csq == 'stop_lost':
			self.stop_lost += 1	
		elif csq == 'start_lost':
			self.start_lost += 1	
		elif csq == 'inframe_insertion':
			self.inframe_insertion += 1
		elif csq == 'inframe_deletion':
			self.inframe_deletion += 1
		elif csq == 'missense_variant':
			self.missense_variant += 1
		elif csq == 'splice_region_variant':
			self.splice_region_variant += 1
		elif csq == 'synonymous_variant':
			self.synonymous_variant += 1

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['splice_acceptor_donor_variant'] = self.splice_acceptor_donor_variant
		dictionary['stop_gained'] = self.stop_gained
		dictionary['frameshift_variant'] = self.frameshift_variant
		dictionary['stop_lost'] = self.stop_lost
		dictionary['start_lost'] = self.start_lost
		dictionary['inframe_insertion'] = self.inframe_insertion
		dictionary['inframe_deletion'] = self.inframe_deletion
		dictionary['missense_variant'] = self.missense_variant
		dictionary['splice_region_variant'] = self.splice_region_variant
		dictionary['synonymous_variant'] = self.synonymous_variant
		dictionary['lof'] = self.splice_acceptor_donor_variant + self.stop_gained + self.frameshift_variant
		dictionary['miss_and_indels'] = self.stop_lost + self.start_lost + self.inframe_insertion + self.inframe_deletion + self.missense_variant
		dictionary['syn'] = self.splice_region_variant + self.synonymous_variant
		return dictionary


class ClinVarGene():
	def __init__(self):
		self.gene_name = ''
		self.canonical_transcript = ''
		self.omim_inheritance = ''
		self.pathogenic = ClinVarTypes()
		self.pathogenic_conflict = ClinVarTypes()
		self.uncertain_significance = ClinVarTypes()
		self.cites_pathogenic = set([])
		self.cites_pathogenic_conflict = set([])
		self.cites_uncertain_significance = set([])

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['gene_name'] = self.gene_name
		dictionary['canonical_transcript'] = self.canonical_transcript
		dictionary['omim_inheritance'] = self.omim_inheritance
		dictionary['pathogenic'] = self.pathogenic.get_dictionary()
		dictionary['pathogenic_conflict'] = self.pathogenic_conflict.get_dictionary()
		dictionary['uncertain_significance'] = self.uncertain_significance.get_dictionary()
		dictionary['cites_pathogenic'] = len(self.cites_pathogenic)
		dictionary['cites_pathogenic_conflict'] = len(self.cites_pathogenic_conflict)
		dictionary['cites_uncertain_significance'] = len(self.cites_uncertain_significance)
		return dictionary


def create_clin_var_genes(db):
	db.gevir.clin_var_genes.drop()
	omim_genes = db.gevir.omim.find({})

	pass_csqs = set(['splice_acceptor_variant',
					'splice_donor_variant',
					'stop_gained',
					'frameshift_variant',
					'stop_lost',
					'start_lost',
					'inframe_insertion',
					'inframe_deletion',
					'missense_variant',
					'splice_region_variant',
					'synonymous_variant'])

	total_lines = omim_genes.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for omim_gene in omim_genes:
		inheritance = omim_gene['inheritance']
		if inheritance == '': # not susceptability only
			line_number += 1
			continue

		gene_name = omim_gene['gene_name']
		transcript_id = omim_gene['transcript_id']

		clin_var_gene = ClinVarGene()
		clin_var_gene.gene_name = gene_name
		clin_var_gene.canonical_transcript = transcript_id
		clin_var_gene.omim_inheritance = inheritance

		clin_vars = db.gevir.clin_vars.find({ "gene_name": gene_name })

		count = 0
		for clin_var in clin_vars:
			if 'csqs' not in clin_var or transcript_id not in clin_var['csqs']:
				continue

			csq = clin_var['csqs'][transcript_id]

			if csq not in pass_csqs:
				continue

			clin_sigs = clin_var['clin_sig'].split(', ')

			if 'Likely pathogenic' in clin_sigs or 'Pathogenic/Likely pathogenic' in clin_sigs or 'Pathogenic' in clin_sigs:
				clin_var_gene.pathogenic.count_csq(csq)
				count += 1
				for cite in clin_var['cites']:
					clin_var_gene.cites_pathogenic.add(cite)
			elif 'Conflicting interpretations of pathogenicity' in clin_sigs:
				clin_var_gene.pathogenic_conflict.count_csq(csq)
				count += 1
				for cite in clin_var['cites']:
					clin_var_gene.cites_pathogenic_conflict.add(cite)
			elif 'Uncertain significance' in clin_sigs:
				clin_var_gene.uncertain_significance.count_csq(csq)
				count += 1
				for cite in clin_var['cites']:
					clin_var_gene.cites_uncertain_significance.add(cite)

		if count > 0:
			db.gevir.clin_var_genes.insert(clin_var_gene.get_dictionary())

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()


def import_clin_var(db):
	import_clin_var_raw(db)
	export_clin_var_ensembl(db)

	'''
	"clin_var_ensembl.txt" obtained from previous step has to ve reannotated with ensembl VEP
	to create json file for the next step.
	Example command to run Ensembl VEP to reannotate ClinVar files
	./vep --offline --everything --fork 6 --cache --force_overwrite --json  
	-i "./source_data/clin_var/clin_var_ensembl.txt" 
	-o "./source_data/clin_var/clin_var_vep.json" 
	'''

	import_clin_var_temp_json(db)
	update_clin_vars_with_canonical_transcripts(db)
	update_clin_vars_with_csqs(db)
	create_clin_var_genes(db)


##########################################
### Conservative Coding RegionS (CCRS) ###
##########################################

class CCR():
	def __init__(self):
		self.chrom = ''
		self.start = 0
		self.end = 0
		self.ccr_pct = 0.0
		self.gene = ''
		self.ranges = ''
		self.varflag = ''
		self.syn_density = 0.0
		self.cpg = 0.0
		self.cov_score = 0.0
		self.resid = 0.0
		self.resid_pctile = 0.0
		self.unique_key = 0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['chrom'] = self.chrom
		dictionary['start'] = self.start
		dictionary['end'] = self.end
		dictionary['ccr_pct'] = self.ccr_pct
		dictionary['gene'] = self.gene
		dictionary['ranges'] = self.ranges
		dictionary['varflag'] = self.varflag
		dictionary['syn_density'] = self.syn_density
		dictionary['cpg'] = self.cpg
		dictionary['cov_score'] = self.cov_score
		dictionary['resid'] = self.resid
		dictionary['resid_pctile'] = self.resid_pctile
		dictionary['unique_key'] = self.unique_key
		return dictionary


def import_ccrs_file(db, ccrs_bed):
	input_file = open(ccrs_bed, 'rt')
	reader = csv.reader(input_file, delimiter='\t')
	headers = next(reader)

	chrom_index = headers.index('#chrom')
	start_index = headers.index('start')
	end_index = headers.index('end')
	ccr_pct_index = headers.index('ccr_pct')
	gene_index = headers.index('gene')
	ranges_index = headers.index('ranges')
	varflag_index = headers.index('varflag')
	syn_density_index = headers.index('syn_density')
	cpg_index = headers.index('cpg')
	cov_score_index = headers.index('cov_score')
	resid_index = headers.index('resid')
	resid_pctile_index = headers.index('resid_pctile')
	unique_key_index = headers.index('unique_key')

	bulk = db.gevir.ccrs.initialize_unordered_bulk_op()
	total_lines = file_len(ccrs_bed)
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for row in reader:
		ccr = CCR()
		ccr.chrom = row[chrom_index]
		ccr.start = int(row[start_index])
		ccr.end = int(row[end_index])
		ccr.ccr_pct = float(row[ccr_pct_index])
		ccr.gene = row[gene_index]
		ccr.ranges = row[ranges_index]
		ccr.varflag = row[varflag_index]
		ccr.syn_density = float(row[syn_density_index])
		ccr.cpg = float(row[cpg_index])
		ccr.cov_score = float(row[cov_score_index])
		ccr.resid = float(row[resid_index])
		ccr.resid_pctile = float(row[resid_pctile_index])
		ccr.unique_key = int(row[unique_key_index])
		ccr = ccr.get_dictionary()

		bulk.insert(ccr)
		if (line_number % 10000 == 0):
			bulk.execute()
			bulk = db.gevir.ccrs.initialize_unordered_bulk_op()

		line_number += 1	
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	if (line_number % 10000 != 0):
		bulk.execute()


def import_ccrs(db):
	db.gevir.ccrs.drop()
	import_ccrs_file(db, CCRS_AUTOSOMAL_TSV)
	import_ccrs_file(db, CCRS_X_TSV)

	db.gevir.ccrs.create_index([('ccr_pct', pymongo.ASCENDING)], name='ccr_pct_1')
	db.gevir.ccrs.create_index([('gene', pymongo.ASCENDING)], name='gene_1')
	db.gevir.ccrs.create_index([('varflag', pymongo.ASCENDING)], name='varflag_1')


class GeneCCRs():
	def __init__(self):
		self.transcript_id = ''
		self.gene_name = ''
		self.ccrs_gte_95 = 0
		self.ccrs_gte_96 = 0
		self.ccrs_gte_97 = 0
		self.ccrs_gte_98 = 0
		self.ccrs_gte_99 = 0

	def get_dictionary(self):
		dictionary = OrderedDict()
		dictionary['transcript_id'] = self.transcript_id
		dictionary['gene_name'] = self.gene_name
		dictionary['ccrs_gte_95'] = self.ccrs_gte_95
		dictionary['ccrs_gte_96'] = self.ccrs_gte_96
		dictionary['ccrs_gte_97'] = self.ccrs_gte_97
		dictionary['ccrs_gte_98'] = self.ccrs_gte_98
		dictionary['ccrs_gte_99'] = self.ccrs_gte_99
		return dictionary


def count_gene_ccrs(db):
	genes_ccrs = {}
	ccrs = db.gevir.ccrs.find({})

	total_lines = ccrs.count()
	line_number = 0
	bar = progressbar.ProgressBar(maxval=1.0).start()
	for ccr in ccrs:
		gene_name = ccr['gene']
		if gene_name not in genes_ccrs:
			genes_ccrs[gene_name] = GeneCCRs()
			genes_ccrs[gene_name].gene_name = gene_name

		ccr_pct = ccr['ccr_pct']
		if ccr_pct >= 95:
			genes_ccrs[gene_name].ccrs_gte_95 += 1
		if ccr_pct >= 96:
			genes_ccrs[gene_name].ccrs_gte_96 += 1
		if ccr_pct >= 97:
			genes_ccrs[gene_name].ccrs_gte_97 += 1
		if ccr_pct >= 98:
			genes_ccrs[gene_name].ccrs_gte_98 += 1
		if ccr_pct >= 99:
			genes_ccrs[gene_name].ccrs_gte_99 += 1

		line_number += 1
		bar.update((line_number + 0.0) / total_lines)
	bar.finish()

	for gene_name in genes_ccrs:
		transcript_id = get_gene_canonical_transcript_by_name(db, gene_name)
		genes_ccrs[gene_name].transcript_id = transcript_id

	db.gevir.ccr_genes.drop()
	for gene_ccrs in genes_ccrs.values():
		gene_ccrs = gene_ccrs.get_dictionary()
		db.gevir.ccr_genes.insert(gene_ccrs)


#############################
### Mac Arthur gene lists ###
#############################

def read_mac_arthur_gene_list(db, gene_list_name, gene_list_tsv):
	with open(gene_list_tsv) as f:
		genes = f.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	gene_names = [gene.strip() for gene in genes]

	transcript_ids = []
	for gene_name in gene_names:
		transcript_id = get_gene_canonical_transcript_by_name(db, gene_name)
		if transcript_id:
			transcript_ids.append(transcript_id)

	db.gevir.mac_arthur_gene_lists.remove({'_id': gene_list_name})
	db.gevir.mac_arthur_gene_lists.insert({'_id': gene_list_name, 'gene_names': gene_names, 'transcript_ids': transcript_ids})


def import_mac_arthur_gene_lists(db):
	read_mac_arthur_gene_list(db, 'lof_tolerant', MAC_ARTHUR_LOF_TOLERANT)
	read_mac_arthur_gene_list(db, 'crispr_essential', MAC_ARTHUR_CRISPR_ESSENTIAL)
	read_mac_arthur_gene_list(db, 'crispr_non_essential', MAC_ARTHUR_CRISPR_NON_ESSENTIAL)

########################
### Mouse Het Lethal ###
########################

def query_mousemine_to_create_mouse_het_lethal_knockout_genes():
	service = Service("http://www.mousemine.org/mousemine/service")
	query = service.new_query("OntologyAnnotation")
	query.add_constraint("ontologyTerm", "MPTerm")
	query.add_constraint("subject", "SequenceFeature")
	query.add_constraint("evidence.baseAnnotations.subject", "Genotype")
	query.add_view(
	    "subject.primaryIdentifier", "subject.symbol",
	    "evidence.baseAnnotations.subject.symbol",
	    "evidence.baseAnnotations.subject.background.name",
	    "evidence.baseAnnotations.subject.zygosity", "ontologyTerm.identifier",
	    "ontologyTerm.name"
	)
	query.add_sort_order("OntologyAnnotation.subject.symbol", "ASC")
	query.add_constraint("evidence.baseAnnotations.subject.zygosity", "=", "ht", code = "B")
	query.add_constraint("ontologyTerm.name", "CONTAINS", "lethal", code = "A")

	headers = ['Subject Primary Identifier',
	           'Ontology Annotation Subject . Symbol',
	           'Base Annotations Subject . Symbol',
	           'Subject Background',
	           'Subject Zygosity',
	           'Ontology Annotation Ontology Term . Identifier',
	           'Ontology Annotation Term Name',]
	table = [headers]

	for query_row in query.rows():
		row = [query_row["subject.primaryIdentifier"], 
		       query_row["subject.symbol"],
		       query_row["evidence.baseAnnotations.subject.symbol"],
		       query_row["evidence.baseAnnotations.subject.background.name"],
		       query_row["evidence.baseAnnotations.subject.zygosity"],
		       query_row["ontologyTerm.identifier"],
		       query_row["ontologyTerm.name"],
		      ]
		table.append(row)

	output_csv = SOURCE_DATA_FOLDER + 'mouse_het_lethal.tsv'
	write_table_to_csv(table, output_csv, delimiter='\t')


def import_mouse_het_lethal_knockout_genes(db):
	mouse_to_human_gene_names = {}
	mouse_to_human = CsvReader(MOUSE_TO_HUMAN_TSV, delimiter='\t')
	for document in mouse_to_human.data:
		mouse_to_human_gene_names[document['mouse_gene_name']] = document['human_gene_name']

	gene_names_het_lethal = set([])
	gene_transcripts_het_lethal = set([])

	mouse_abnormal_survival = CsvReader(MOUSE_HET_LETHAL_TSV, delimiter='\t')
	for document in mouse_abnormal_survival.data:
		mouse_gene_name = document['Ontology Annotation Subject . Symbol']
		mouse_zygosity = document['Subject Zygosity']
		mouse_phenotype = document['Ontology Annotation Term Name']
		# Ignore genes for which there are no human homologs
		if mouse_gene_name not in mouse_to_human_gene_names:
			continue

		human_gene_name = mouse_to_human_gene_names[mouse_gene_name]
		transcript_id = get_gene_canonical_transcript_by_name(db, human_gene_name)

		# Ignore genes which names were not found in ExAC
		if not transcript_id:
			continue

		exac_gene = db.exac.genes.find_one({'canonical_transcript': transcript_id})
		gene_name = exac_gene['gene_name']
		gene_names_het_lethal.add(gene_name)
		gene_transcripts_het_lethal.add(transcript_id)

	db.gevir.mouse_het_lethal.drop()
	db.gevir.mouse_het_lethal.insert({'_id': 'gene_names', 'data': list(gene_names_het_lethal)})
	db.gevir.mouse_het_lethal.insert({'_id':'transcript_ids', 'data': list(gene_transcripts_het_lethal)})


def main():
	db = MongoDB()
	# Creates additional database (gerp) with gerp score for each chromosomal position
	# IMPORTANT: this operation might require quite a lot of time to run (progress for each chromosome will be displayed)
	# Final database size should be ~36.6 GB (storage size occupied on disk)
	#import_gerp(db)

	# ENS CDS FASTA
	#import_ens_cds_fasta(db, ENS_CDS_FASTA, 'ens_cds_fasta')
	#import_ens_cds_fasta(db, ENS_AA_FASTA, 'ens_aa_fasta')

	# gnomAD scores (22/10/18)
	#import_gnomad_scores(db, new_gnomad_file=False)

	# OMIM data donwloaded (11/13/18)
	#import_omim(db)

	# ClinVar data donwloaded (21/08/18)
	#import_clin_var(db)

	# Conservative Coding RegionS
	# https://s3.us-east-2.amazonaws.com/ccrs/ccr.html 
	#import_ccrs(db)
	#count_gene_ccrs(db)

	# Mac Arthur Datasets
	# https://github.com/macarthur-lab/gene_lists
	#import_mac_arthur_gene_lists(db)

	# Mouse Het Lethal Knockout Genes (5/02/19)
	# http://www.mousemine.org/mousemine/templates.do
	# Mammalian phenotypes (MP terms) --> Mouse genes and models
	# Search for *lethal*
	# Alternatively, use following method to query mousemine database:
	# query_mousemine_to_create_mouse_het_lethal_knockout_genes()

	# import_mouse_het_lethal_knockout_genes(db)

if __name__ == "__main__":
	sys.exit(main())