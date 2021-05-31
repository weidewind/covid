import argparse
import pandas as pd
import csv
#from transliterate import translit
import googletrans


parser = argparse.ArgumentParser(description='Collect the results of global analysis into xlsx')
parser.add_argument("--pango_file", dest='pango_file', type=str, help="full path to the alignment file")
parser.add_argument("--meta_file", dest='meta_file', type=str, help="full path to the metadata file")
parser.add_argument("--gisaid_meta_file", dest='gisaid_meta_file', type=str, help="full path to the metadata file")
parser.add_argument("--output", dest='output', type=str, help="full path to the output file")
parser.add_argument("--mut_file", dest='mut_file', type=str, help="full path to the file with mutations of interest (in genomic coordinates, one mutation  per line)")
args = parser.parse_args()

pango = pd.read_csv(args.pango_file, sep = ",").iloc[:, 0:3]
pango.rename(columns = {'taxon': 'seq_id'}, inplace = True)
print(pango[1:3])

muts = pd.read_csv(args.mut_file, sep = ",")
print(muts[1:3])

meta = pd.read_excel(args.meta_file)[['Внутренний номер', 'Место забора']]
meta.columns = ['seq_id', 'region']
translator = googletrans.Translator()
#meta['region'] = meta['region'].map(lambda a: translator.translate(a, src='ru', dest='en'))
print(meta['region'].head().map(lambda a: translator.translate(a, src='ru', dest='en')))
#meta['region'] = meta['region'].map(lambda a: translit(a, "ru", reversed=True))

gismeta = pd.read_csv(args.gisaid_meta_file, sep = "\t")[['strain', 'location']]
gismeta.columns = ['seq_id', 'region']
print(gismeta.head())

meta = pd.concat([meta, gismeta])
print (meta.head())



merged = pd.merge(pango, muts, on='seq_id')
merged = pd.merge(meta, merged, how='right', on='seq_id')


print(merged.head())

def summ_exact(myseries):
	mylist = myseries.tolist()
	summ = 0
	for val in mylist:
		if val == 1 or val == "1":
			summ += 1
	if (summ == 0):
		summ = ""
	return(summ)

agg_over_muts = {mut: summ_exact for mut in merged.head().iloc[:,4:]}
out = merged.groupby(['lineage', 'region'], dropna=False).agg(agg_over_muts)
total = merged.groupby(['lineage', 'region'], dropna=False).count().iloc[:,:1]
total.columns = ['total']
print(total)
print(out)

out = pd.merge(out, total, on=['lineage','region'])


#print(out)
merged.to_csv(args.output, index = False)
out.to_csv(args.output + ".summary")