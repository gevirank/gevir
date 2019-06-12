import os
import sys
import csv
import progressbar
from collections import OrderedDict
import time

from common import is_float, is_int

csv.field_size_limit(sys.maxsize)

class CsvReader():
	"""Reads csv data from the given file.
	
	CSV is parsed by provided delimiter, default is ','.
	By default, also tries to convert column values to int or float.
	Convertion occures only when ALL values in the column can be converted.
	For example, values in column with chromosome number are not converted
	because X and Y are strings. This behaviour can be disabled by setting
	auto_detect_types = False.
	"""
	def __init__(self, path_to_csv, delimiter=',', auto_detect_types=True, key_column_name=''):
		self.path_to_csv = path_to_csv
		self.delimiter = delimiter
		self.auto_detect_types = auto_detect_types
		self.row_num = 0
		self.column_num = 0
		self.headers = []
		self.data = []
		self.data_dict = {}
		self._read_csv(key_column_name)


	def import_to_db(self, db, collection, remove_existing=True):
		print 'Importing data to "' + collection + '" collection'
		if remove_existing:
			db[collection].drop()

		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for document in self.data:
			db[collection].insert(document)
			line_number += 1
			bar.update((line_number + 0.0) / self.row_num)
		bar.finish()

		
	def _read_csv(self, key_column_name):
		path, file = os.path.split(self.path_to_csv)
		print 'Reading file "' + file + '"...'
		input_file = open(self.path_to_csv, 'rt')
		reader = csv.reader(input_file, delimiter=self.delimiter)

		self.headers = next(reader)
		self.column_num = len(self.headers)
		self.row_num = len(open(self.path_to_csv).readlines())

		if self.auto_detect_types:
			int_columns = [True] * len(self.headers)
			float_columns = [True] * len(self.headers)

		line_number = 0
		bar = progressbar.ProgressBar(maxval=1.0).start()
		for row in reader:

			for x in range(0, self.column_num):
				if self.auto_detect_types:
					if int_columns[x]:
						int_columns[x] = is_int(row[x])
					if float_columns[x]:
						float_columns[x] = is_float(row[x])

			row_dict = self.__row_to_dict(row)
			self.data.append(row_dict)

			if key_column_name:
				self.data_dict[row_dict[key_column_name]] = row_dict

			line_number += 1
			bar.update((line_number + 0.0) / self.row_num)
		bar.finish()

		if self.auto_detect_types:
			self.__format_data(int_columns, float_columns)


	def __format_data(self, int_columns, float_columns):
		int_headers = set()
		float_headers = set()

		for x in range(0, self.column_num):
			if int_columns[x]:
				int_headers.add(self.headers[x])
			elif float_columns[x]:
				float_headers.add(self.headers[x])

		for row_dict in self.data:
			for key, value in row_dict.iteritems():
				if key in int_headers:
					row_dict[key] = int(float(value))
				elif key in float_headers:
					row_dict[key] = float(value)


	def __row_to_dict(self, row):
		dictionary = OrderedDict()
		for x in range(0, len(self.headers)):
			dictionary[self.headers[x]] = row[x]
		return dictionary