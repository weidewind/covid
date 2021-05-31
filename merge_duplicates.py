import pandas as pd
import csv

duplicates = {}

output = "/export/home/popova/workspace/covid/data/munched/gennady/merged_far20.duplicates"

clusters = "/export/home/popova/workspace/covid/data/munched/gennady/nonrus_clusters_far20.txt"
justdupl = "/export/home/popova/workspace/covid/data/munched/gennady/duplicates_21.03.txt"
rusdupl = "/export/home/popova/workspace/covid/data/munched/gennady/rus.duplicates"
nonrusdupl = "/export/home/popova/workspace/covid/data/munched/gennady/nonrus.duplicates"
rus_to_nonrus_dupl = "/export/home/popova/workspace/covid/data/munched/gennady/rus_to_nonrus.duplicates"

meta = "/export/home/popova/workspace/covid/data/russian/meta_all.xlsx"

with open(justdupl, "r") as rd:
	for line in rd:
			splitter = line.strip().split(';')
			duplicates[splitter[0]] = splitter[1:]

rusduplicates = {}
with open(rusdupl, "r") as rd:
	for line in rd:
			splitter = line.strip().split(';')
			rusduplicates[splitter[0]] = splitter[1:]

with open(rus_to_nonrus_dupl, "r") as rnd:
	for line in rnd:
		splitter = line.strip().split('\t')
		russ = splitter[0].split(';')
		for r in russ:
			if r in rusduplicates:
				duplist = splitter[1].split(';')
				oldlist = []
				if r in duplicates:
					oldlist = duplicates[r]
				for d in duplist:
					if d in duplicates:
						prev = duplicates.pop(d)
						oldlist.extend(prev)
				oldlist.extend(duplist)
				duplicates[r] = oldlist
				break


for f in [rusdupl, nonrusdupl, clusters]:
	with open(f, "r") as rd:
		for line in rd:
			splitter = line.strip().split(';')
			duplist = splitter[1:]
			oldlist = []
			if splitter[0] in duplicates:
				oldlist = duplicates[splitter[0]]
			for d in duplist:
				if d in duplicates:
					prev = duplicates.pop(d)
					oldlist.extend(prev)
			oldlist.extend(duplist)
			duplicates[splitter[0]] = oldlist

meta = pd.read_excel(meta)[['Внутренний номер', 'gisaid_id']]
meta.columns = ['seq_id', 'gisaid_id']
for id in duplicates:
	if id in meta['gisaid_id']:
		print("Warning! duplicate of russian strain, " + id + ", is one of the centroids!")

gisaid_dupl = meta['gisaid_id'].tolist()
with open (output, 'w') as out:
	for key, value in duplicates.items():
		clean_value = [i for i in value if i not in gisaid_dupl]
		out.write(key+"\t"+";".join(clean_value)+"\n")
