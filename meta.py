import pandas as pd
import re
import datetime
#import optparse

#parser = optparse.OptionParser()
#parser.add_option('-r', '--rpn_meta', help='', type='str')
#parser.add_option('-m', '--meta', help='', type='str')

#options, args = parser.parse_args()

options_rpn_meta="/export/home/popova/workspace/covid/data/raw/17_06_2021_meta.csv"
options_meta="/export/home/popova/workspace/covid/data/raw/metadata.tsv"
options_genbank_meta="/export/home/popova/workspace/covid/data/raw/public-latest.metadata.tsv"


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
	try:
		date = re.sub('-XX', '', date)
	except:
		print("date_" + str(date) + "_date")
		raise
	if (len(date) < 10):
		date = "unknown"
	else:
		try:
			datetime.datetime.strptime(date, '%Y-%m-%d')
		except ValueError:
			try:
				datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S')
			except ValueError:
					print("Incorrect data format, should be YYYY-MM-DD [HH:MM:SS]: " + date)
					date = "unknown"
	return(date)


def get_date_dict():
	print("Parsing metas..")
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'Дата забора']]
	meta.columns = ['seq_id', 'date']
	gismeta = pd.read_csv(options_meta, sep="\t")[['Accession ID', 'Collection date']]
	gismeta.columns = ['seq_id', 'date']
	print(gismeta.head())
	genmeta = pd.read_csv(options_genbank_meta, sep="\t")[['strain', 'date']]
	genmeta.columns = ['seq_id', 'date']
	print(genmeta.head())
	meta = pd.concat([meta, gismeta, genmeta])
	try:
		meta['date'] = meta['date'].apply(clean_dates)
	except:
		print([str(d) for d in meta['date']])
		raise
	meta_dict = dict(zip(meta['seq_id'], meta['date']))
	return(meta_dict)


def get_gisaid_duplicates():
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'gisaid_id']]
	meta.columns = ['seq_id', 'gisaid_id']
	return(dict(zip(meta['seq_id'], meta['gisaid_id'])))