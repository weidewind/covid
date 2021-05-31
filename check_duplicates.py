import pandas as pd
from meta import get_country_dict

duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_rus_to_rus_dups.tsv"
rtn_duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_rus_to_nonrus_dups.tsv"
ntn_duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_nonrus_to_nonrus_dups.tsv"


meta_dict = get_country_dict()


print("Parsing rus duplicates..")
with open(duplicates_file, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		for s in splitter:
			dups = s.split(';')
			for d in dups:
				if meta_dict.get(d, "NA") != "Russia":
					print(d + " is not russian! " + meta_dict.get(d, "NA"))


print("Parsing rus to nonrus duplicates..")
with open(rtn_duplicates_file, "r") as dfile:
	for line in dfile:
		rus, nonrus = line.strip().split('\t')
		rdups = rus.split(';')
		for rd in rdups:
			if meta_dict.get(rd, "NA") != "Russia":
				print(rd + " is not russian! " + meta_dict.get(rd, "NA"))
		ndups = nonrus.split(';')
		for nd in ndups:
			if meta_dict.get(nd, "NA") == "Russia":
				print(rd + " is  russian! " + meta_dict.get(nd, "NA"))


print("Parsing nonrus duplicates..")
with open(ntn_duplicates_file, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		for s in splitter:
			dups = s.split(';')
			for d in dups:
				if meta_dict.get(d, "NA") == "Russia":
					print(d + " is russian! " + meta_dict.get(d, "NA"))