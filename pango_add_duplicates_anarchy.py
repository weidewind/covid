import pandas as pd
import optparse

# Used in pango_stat.sh
#pango_file = "/export/home/popova/workspace/covid/data/munched/gennady/may25/rus_unique.pangolined"
#duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may25/21_05_2021_rus_to_rus_dups.tsv"
#output = "/export/home/popova/workspace/covid/data/munched/gennady/may25/rus.pangolined.withduplicates"

parser = optparse.OptionParser()
parser.add_option('-p', '--pango_file', help='pango_file', type='str')
parser.add_option('-d', '--duplicates_file', help='rus_to_rus duplicates ', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

pangodict = {}
with open(options.output, "w") as out:
	with open(options.pango_file, "r") as pango:
		for l in pango:
			out.write(l)
			splitter = l.strip().split(",")
			pangodict[splitter[0]] = ",".join(splitter[1:])

	with open(options.duplicates_file, "r") as duplicates:
		for l in duplicates:
			splitter = l.strip().split("\t")
			str_in_tree = splitter[0]
			dups = splitter[1].split(";")
			dups.append(str_in_tree)
			for d in dups:
				if d in pangodict:
					dups.remove(d)
					for s in dups:
						out.write(s + "," + pangodict[d] + "\n")
					break



