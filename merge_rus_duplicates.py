from meta import get_country_dict, get_gisaid_duplicates
import pandas as pd
import optparse

# add duplicates of russian strains from --duplicates to --rus_duplicates
# !!! unluckily, also add any nonrussian strains with russian duplicates 
# (=> the resulting file can contain nonrussian strains at the start of the line)

parser = optparse.OptionParser()
parser.add_option('-r', '--rtr_duplicates', help='file with rus duplicates of rus strains', type='str')
parser.add_option('-n', '--rtn_duplicates', help='file with nonrus duplicates of rus strains', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

meta_dict = get_country_dict()

dup_list = {}

print("Parsing --rtr_duplicates..")
with open(options.rtr_duplicates, "r") as rus_dupls:
	for line in rus_dupls:
		rus, dupls = line.strip().split('\t')
		dupls_arr = dupls.split(';')
		if rus in dup_list:
			dup_list[rus].extend(dupls_arr)
		else:
			dup_list[rus] = dupls_arr


print("Parsing --rtn_duplicates..")
with open(options.rtn_duplicates, "r") as nonrus_dupls:
	for line in nonrus_dupls:
		rus, dupls = line.strip().split('\t')
		dupls_arr = dupls.split(';')
		if rus in dup_list:
			dup_list[rus].extend(dupls_arr)
		else:
			dup_list[rus] = dupls_arr



all_rus = list(dup_list.keys())
all_rus.extend(dup_list.values())

gisaid_dupl = get_gisaid_duplicates()

for ruskey in gisaid_dupl.keys():
	if ruskey in all_rus:
		if gisaid_dupl[ruskey] in all_rus:
			print("Oh noes! both rus and gisaid id for the same strain are found in rus files!")


print("Writing new duplicates list..")
with open(options.output, "w") as out:
	for k,v in dup_list.items():
		#clean_v = [i for i in v if i not in gisaid_dupl]
		out.write(k + "\t" + ";".join(v) + "\n")
