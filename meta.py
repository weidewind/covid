import pandas as pd
options_rpn_meta = "/export/home/popova/workspace/covid/data/raw/21_05_2021_meta.csv"
options_meta = "/export/home/popova/workspace/covid/data/raw/metadata_2021_05_31.tsv"


def get_country_dict():
	print("Parsing metas..")
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер']]
	meta.columns = ['seq_id']
	meta['country'] = 'Russia'
	gismeta = pd.read_csv(options_meta, sep="\t")[['Accession ID', 'Location']]
	gismeta.columns = ['seq_id', 'location']
	gismeta['country'] = gismeta['location'].str.split('/',expand=True)[1].str.strip()
	print(gismeta.head())
	meta = pd.concat([meta, gismeta])
	meta_dict = dict(zip(meta['seq_id'], meta['country']))
	return(meta_dict)


def get_date_dict():
	print("Parsing metas..")
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'Дата забора']]
	meta.columns = ['seq_id', 'date']
	gismeta = pd.read_csv(options_meta, sep="\t")[['Accession ID', 'Collection date']]
	gismeta.columns = ['seq_id', 'date']
	print(gismeta.head())
	meta = pd.concat([meta, gismeta])
	meta_dict = dict(zip(meta['seq_id'], meta['date']))
	return(meta_dict)


def get_gisaid_duplicates():
	meta = pd.read_csv(options_rpn_meta, sep="\t")[['Внутренний номер', 'gisaid_id']]
	meta.columns = ['seq_id', 'gisaid_id']
	return(dict(zip(meta['seq_id'], meta['gisaid_id'])))