import csv
import pymongo
import numpy as np
import math
from collections import OrderedDict
from decimal import Decimal


#################
### CONSTANTS ###
#################
DB_HOST = 'localhost'
DB_PORT = 27017
DB_NAME_EXAC = 'exac'
DB_NAME_GERP = 'gerp'
DB_NAME_GEVIR = 'gevir'
DB_NAME_REFERENCE_GENOME = 'reference_genome'
		

class MongoDB():
	"""Database Client."""
	def __init__(self):
		client = pymongo.MongoClient(host=DB_HOST, port=DB_PORT, document_class=OrderedDict)
		self.exac = client[DB_NAME_EXAC]
		self.gerp = client[DB_NAME_GERP]
		self.gevir = client[DB_NAME_GEVIR]
		self.ref = client[DB_NAME_REFERENCE_GENOME]


def file_len(fname):
	"""Calculate length of a file."""
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1


def is_float(x):
	"""Check if value (e.g. string) can be converted to float."""
	try:
		a = float(x)
	except ValueError:
		return False
	else:
		return True


def is_int(x):
	"""Check if value (e.g. string) can be converted to integer."""
	try:
		a = float(x)
		b = int(a)
	except ValueError:
		return False
	else:
		return a == b


def calculate_percentiles(ranked_list, reverse=False):
	"""Return list of percetiles based on number of elements in input list."""
	percentiles = OrderedDict()
	max_num = len(ranked_list)

	percentile = 0.0
	for x in range(0, max_num):
		if reverse:
			percentile = (1 - float(x + 1) / max_num) * 100
		else:
			percentile = float(x + 1) / max_num * 100
		percentiles[ranked_list[x]] = percentile

	return percentiles


def write_table_to_csv(table, output_csv, delimiter=','):
	"""Write table (list of lists) to csv."""
	output_file = open(output_csv,'w+')
	writer = csv.writer(output_file, delimiter=delimiter)

	for row in table:
		writer.writerow(row)

	output_file.close()


def sort_dict_by_values(dictionary, reverse=False):
	"""Return dictionary sorted by values."""
	sorted_tuples = sorted(dictionary.items(), key=lambda x: x[1], reverse=reverse)
	result = OrderedDict()
	for x in range(0, len(sorted_tuples)):
		result[sorted_tuples[x][0]] = sorted_tuples[x][1]
	return result


def get_gene_canonical_transcript_by_name(db, gene_name):
	""" Return gene canonical transcript if gene name is in exac database.
	
	Taken from exac lookups.py (modified to return canonical transcript).
	"""
	
	gene = db.exac.genes.find_one({'gene_name': gene_name}, projection={'_id': False})
	if not gene:
		gene = db.exac.genes.find_one({'other_names': gene_name}, projection={'_id': False})
	if gene and 'canonical_transcript' in gene:
		return gene['canonical_transcript']
	else:
		return ''

def float_to_sci_str(num):
	return "{:.2E}".format(Decimal(num))


def proportion_to_percents_str(proportion):
	return "{0:.1f}".format(proportion*100)

# Code taken from this tutorial:
# https://machinelearningmastery.com/calculate-bootstrap-confidence-intervals-machine-learning-results-python/
def calculate_confidence_interval(stats, alpha=0.95):
	p = ((1.0-alpha)/2.0) * 100
	lower = max(0.0, np.percentile(stats, p))
	p = (alpha+((1.0-alpha)/2.0)) * 100
	upper = min(1.0, np.percentile(stats, p))
	return lower, upper


def report_2_x_2_test_table_as_obs_exp_str(table):
	return 'Observed ' + str(table[0][0]) + ' out of ' + str(table[0][1])  + ' | ' + 'Expected ' + str(table[1][0])  + ' out of ' + str(table[1][1])
