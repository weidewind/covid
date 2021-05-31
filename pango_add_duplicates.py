import pandas as pd

pango_file = "/export/home/popova/workspace/covid/data/munched/gennady/may25/rus_unique.pangolined"
duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may25/21_05_2021_rus_to_rus_dups.tsv"
output = "/export/home/popova/workspace/covid/data/munched/gennady/may25/rus.pangolined.withduplicates"

pangodict = {}
with open(output, "w") as out:
	with open(pango_file, "r") as pango:
		for l in pango:
			out.write(l)
			splitter = l.strip().split(",")
			pangodict[splitter[0]] = ",".join(splitter[1:])

	with open(duplicates_file, "r") as duplicates:
		for l in duplicates:
			splitter = l.strip().split("\t")
			str_in_tree = splitter[0]
			dups = splitter[1].split(";")
			for d in dups:
				out.write(d + "," + pangodict[str_in_tree] + "\n")



