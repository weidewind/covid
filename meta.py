import pandas as pd
import re
import datetime
import optparse
import os.path


options_rpn_meta ='/export/home/popova/workspace/covid/data/raw/17_06_2021_meta.csv'
options_meta ='/export/home/popova/workspace/covid/data/raw/metadata.tsv'
options_genbank_meta ='/export/home/popova/workspace/covid/data/raw/public-latest.metadata.tsv'



def get_location_dict(colnum):
	print("Parsing metas..")
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'Место забора']]
	meta.columns = ['seq_id', 'location']
	if colnum == 1:
		meta['country'] = 'Russia'
	elif colnum == 2:
		meta['country'] = meta['location'].str.split('/', expand=True)[1].str.strip()
		meta['country'] = meta['country'].fillna(meta['location'].str.split('/', expand=True)[0].str.strip())
		meta['country'] = meta['country'].fillna("unknown")
	print(meta.head())
	gismeta = pd.read_csv(options_meta, sep="\t")[['Accession ID', 'Location']]
	gismeta.columns = ['seq_id', 'location']
	gismeta['country'] = gismeta['location'].str.split('/', expand=True)[colnum].str.strip()
	gismeta['country'] = gismeta['country'].fillna("unknown")
	print(gismeta.head())
	genmeta = pd.read_csv(options_genbank_meta, sep="\t")[['strain', 'country']]
	genmeta.columns = ['seq_id', 'location']
	if colnum == 1:
		genmeta['country'] = genmeta['location'].fillna("unknown")
	elif colnum == 2:
		genmeta['country'] = "unknown"
	print(genmeta.head())
	meta = pd.concat([meta, gismeta, genmeta])
	meta_dict = dict(zip(meta['seq_id'], meta['country']))
	return(meta_dict)


def get_country_dict():
	return(get_location_dict(colnum=1))


def get_region_dict():
	return(get_location_dict(colnum=2))


def clean_dates(date):
	date = date.split(' ')[0] # get rid of hours and minutes
	try:
		date = re.sub('-XX', '', date)
	except:
		print("date_" + str(date) + "_date")
		raise
	if (len(date) < 8): # day or month unknown
		date = "unknown"
	else:
		try:
			date = datetime.datetime.strptime(date, '%Y-%m-%d')
			date = datetime.datetime.strftime(date, '%Y-%m-%d')
		except ValueError:
			try:
				date = datetime.datetime.strptime(date, '%d.%m.%Y')
				date = datetime.datetime.strftime(date, '%Y-%m-%d')
			except ValueError:
				try:
					date = datetime.datetime.strptime(date, '%d.%m.%y')
					date = datetime.datetime.strftime(date, '%Y-%m-%d')
				except ValueError:
						print("Incorrect data format, should be YYYY-MM-DD [HH:MM:SS] or at least dd.mm.yy or dd.mm.yyyy : " + date)
						date = "unknown"
						raise
	return(date)


def get_date_dict(dates_file):
	if os.path.isfile(dates_file):
		dates = pd.read_csv(dates_file, sep="\t", names=['seq_id', 'date'])
		dates_dict = dict(zip(dates['seq_id'], dates['date']))
		return(dates_dict)
	else:
		print("Parsing metas..")
		meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'Дата забора']]
		meta.columns = ['seq_id', 'date']
		gismeta = pd.read_csv(options_meta, sep="\t")[['Accession ID', 'Collection date']]
		gismeta.columns = ['seq_id', 'date']
		print("gisaid head: ")
		print(gismeta.head())
		genmeta = pd.read_csv(options_genbank_meta, sep="\t")[['strain', 'date']]
		genmeta.columns = ['seq_id', 'date']
		print("genbank head: ")
		print(genmeta.head())
		meta = pd.concat([meta, gismeta, genmeta])
		try:
			print("Cleaning dates..")
			meta['date'] = meta['date'].apply(clean_dates)
		except:
			print("Something went wrong when I was trying to clean these dates in meta.py: " + [str(d) for d in meta['date']])
			raise
		meta_dict = dict(zip(meta['seq_id'], meta['date']))
		with open(dates_file, "w") as out:
			for strain, date in meta_dict.items():
				out.write(strain + "\t" + date + "\n")
		return(meta_dict)


def get_gisaid_duplicates():
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'gisaid_id']]
	meta.columns = ['seq_id', 'gisaid_id']
	return(dict(zip(meta['seq_id'], meta['gisaid_id'])))