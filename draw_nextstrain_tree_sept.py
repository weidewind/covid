import optparse
import pandas as pd
import numpy as np
import re
from collections import Counter
from meta import get_region_dict
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace
from tree_utils import add_data_from_data_and_duplicates_files
from datetime import datetime
import webcolors
from meta import get_date_dict
import random


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--dates', help='leaf_dates.csv', type='str')
parser.add_option('-s', '--states', help='leaf_states.csv', type='str')
parser.add_option('-i', '--important_lineages', help='lineages_of_interest', type='str')
parser.add_option('-r', '--rusmeta', help='', type='str')
parser.add_option('-l', '--latest_date', help='', type='str')
parser.add_option('-n', '--gismeta', help='', type='str')
parser.add_option('-p', '--pangolined', help='', type='str')
parser.add_option('-m', '--mode', help='c for circular or r for rectangular', type='str')
parser.add_option('-e', '--print_delta', action = "store_true", dest = "print_delta", default = False, help='should i also print delta tree separately?')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()


def find_color_generalized(nname, default_color):
	lin = pango_dict[nname]
	splitter = lin.split(".")
	for i in reversed(range(1,len(splitter)+1)):
		upperlin = ".".join(splitter[0:i])
		if upperlin in pango_color_dict:
			return pango_color_dict[upperlin]
	return(default_color)
	#c = pango_color_dict.get(lin, default_color)
	#lin_up = ".".join(lin.split(".")[:-1])
	#c = pango_color_dict.get(lin, pango_color_dict.get(lin_up, default_color))
	#print(lin_up)
	#return(c)


def find_faded_color(nname, mycolor, default_color, blending_color):
	c = mycolor
	if c == default_color:
		return (c)
	strdate = dates_dict.get(nname, "unknown")
	if strdate[:1] not in [str(i) for i in list(range(0, 9))]:
		return (blend(c, blending_color))
	else:
		if len(strdate) > 10:
			d = strdate.split(" ")[0]
			date =  datetime.strptime(d, "%d.%m.%y") 
		else:
			date = datetime.strptime(strdate, "%Y-%m-%d")

		#now = datetime.now()
		#diff = now - date
		
		diff = latest_date - date
		mydiff = round((diff.days-20)/30)
		i = 0
		while i < mydiff:
			c = blend(c, blending_color)
			i += 1
		#print(",".join([str(mydiff), c]))
		return(c)


def blend(c1, c2):
	if c1[0] != "#":
		c1 = str(webcolors.name_to_hex(c1))
	if c2[0] != "#":
		c2 = str(webcolors.name_to_hex(c2))
	#print(c1+" "+c2)
	r = int((int(("0x"+c1[1:3]),16)+int(("0x"+c2[1:3]),16))/0x2)
	g = int((int(("0x"+c1[3:5]),16)+int(("0x"+c2[3:5]),16))/0x2)
	b = int((int(("0x"+c1[5:]),16)+int(("0x"+c2[5:]),16))/0x2)
	return("#"+str(hex(r))[2:].zfill(2)+str(hex(g))[2:].zfill(2)+str(hex(b)[2:].zfill(2)))


def color_tree(t, size):
	for node in t.traverse():
		nstyle = NodeStyle()
		nstyle["size"] = 0
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		nstyle["vt_line_color"] = "#aaaaaa"
		nstyle["hz_line_color"] = "#aaaaaa"
		if node.is_leaf():
			col = blend(find_color_generalized(node.name, "#aaaaaa"),"#aaaaaa")
			nstyle["vt_line_color"] = col
			nstyle["hz_line_color"] = col
			#if node.name == "EPI_ISL_2533828" or pango_dict[node.name] == "B.1.617.2":
			#	node.add_face(TextFace(node.name + " " + pango_dict[node.name], fgcolor="#000000",fsize = 50), column=0)
		if country_dict.get(node.name, "unknown") == "Russia":
			if node.is_leaf():
				mycol = find_color_generalized(node.name, "#777777")
				col = find_faded_color(node.name, mycol, "#777777", "#000000")
				nstyle["fgcolor"] = col
				nstyle["size"] = 4
				#if pango_dict[node.name] in ["B.1.1.523", "AT.1"]:
				#	node.add_face(TextFace(node.name + " " + pango_dict[node.name], fgcolor="#000000",fsize = size), column=0)
		node.set_style(nstyle)


def label_lineages(t, size):
	leaves = t.get_leaves()
	for lin in pango_color_dict:
		lin_nodes = [l for l in leaves if pango_dict[l.name] == lin and country_dict.get(l.name, "unknown") == "Russia"]
		if lin_nodes:
			lnodes = random.sample(lin_nodes, 1)
			for n in lnodes:
				n.add_face(TextFace(lin, fgcolor="#000000", fsize = size), column=0, position="aligned")



def rus_maxdate():
	alldates = [datetime.strptime(strdate, "%Y-%m-%d") for strain, strdate in dates_dict.items() if strdate != "unknown" and country_dict.get(strain, "unknown") == "Russia"]
	return(max(alldates))

# def maxdate():
# 	#dates_num = pd.to_numeric(pd.Series(dates_dict.values()))
# 	#max = pd.to_datetime(np.max(dates_num)).strftime("%Y-%m-%d")
# 	maxdate = datetime.strptime(max(dates_dict.values()), "%d.%m.%y")
# 	return(maxdate)


print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing pangolineages..")
pango = pd.read_csv(options.pangolined, sep=",")
pango_dict = dict(zip(pango['taxon'], pango['lineage']))
print(dict(list(pango_dict.items())[0:5]))


print("Parsing countries..")
countries = pd.read_csv(options.states, sep="\t", names=['seq_id', 'state', 'region'])
country_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(country_dict.items())[0:5]))


# print("Parsing countries..")
# rusmeta = pd.read_csv(options.rusmeta, sep="\t")[['Внутренний номер', 'Место забора. Страна']]
# rusmeta.columns = ['seq_id', 'location']
# rusmeta['location'] = 'Russia'
# gismeta = pd.read_csv(options.gismeta, sep="\t")[['Accession ID', 'Location']]
# gismeta.columns = ['seq_id', 'location']
# gismeta['location'] = gismeta['location'].fillna("unknown")
# print(gismeta.head())
# meta = pd.concat([rusmeta, gismeta])
# country_dict = dict(zip(meta['seq_id'], meta['location']))
# print(dict(list(country_dict.items())[0:5]))


print("Parsing dates..")
# rusmeta = pd.read_csv(options.rusmeta, sep="\t")[['Внутренний номер', 'Дата забора']]
# rusmeta.columns = ['seq_id', 'date']
# dates_dict = dict(zip(rusmeta['seq_id'], rusmeta['date']))
dates_dict = get_date_dict(options.dates)
print(dict(list(dates_dict.items())[0:5]))
if not options.latest_date:
	latest_date = rus_maxdate()
else:
	latest_date = datetime.strptime(options.latest_date, '%Y-%m-%d')
#latest_date = datetime.strptime("2021-09-06", "%Y-%m-%d")
print("latest date: " + latest_date.strftime("%Y-%m-%d"))

print("Parsing colors..")
#B.1.1.451->B.1.1.523
# pango_color_dict = {"B.1.1.7":"blue", "B.1.351":"orange", "B.1.617.2":"crimson", 
# 					"B.1.1.523":"violet", "B.1.1.317":"lawngreen", "B.1.525":"#abdda4",
# 					"A.23.1":"#9e0142", "AT.1":"green", "P.1":"#fee08b", "B.1.1.28.1":"#fee08b","AY":"dodgerblue",
# 					"AY.4":"gold","AY.12":"aqua"}

colors = pd.read_csv(options.important_lineages, sep=",", dtype=str)[["lineage", "color"]]
pango_color_dict = dict(zip(colors["lineage"], colors["color"]))
print(pango_color_dict)





# print("Parsing regions..")
# regions_dict = get_region_dict()
# print(dict(list(regions_dict.items())[0:5]))

# print("Parsing dates..")
# dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
# dates_dict = dict(zip(dates['seq_id'], dates['date']))
# print(dict(list(dates_dict.items())[0:5]))


print("Styling tree..")

ts = TreeStyle()
ts.mode = options.mode
#ts.root_opening_factor = 1
ts.show_leaf_name = False
if options.mode == "r":
	ts.scale = 100
	size = 50
else:
	size = 200 #600
for lin in pango_color_dict:
	ts.legend.add_face(CircleFace(size, pango_color_dict[lin]), column=0)
	ts.legend.add_face(TextFace(lin, fsize = size), column=1)
ts.legend.add_face(TextFace("last date: " + latest_date.strftime("%Y-%m-%d"), fsize = size), column=1) 
ts.allow_face_overlap = True
#ts.layout_fn = multiple_country_layout


color_tree(tree, size)
label_lineages(tree, size)


tree.ladderize()
tree.render(options.output + ".svg", tree_style=ts, dpi=72, w=10000, units="mm")
tree.render(options.output + ".pdf", tree_style=ts, dpi=72, w=10000, units="mm")

if options.print_delta:
	print("Looking for delta..")
	leaves = tree.get_leaves()
	delta_nodes = [l for l in leaves if len(l.name)>1 and (pango_dict[l.name] == "B.1.617.2" or pango_dict[l.name][:2] == "AY")]
	print("Found "+str(len(delta_nodes))+" delta strains, looking for the lca..")
	#lca = tree.get_common_ancestor(delta_nodes, random.sample(delta_nodes,int(len(delta_nodes)/10))) 
	lca = delta_nodes[0].up.up
	print("lca of part of deltas is " + lca.name)
	print("Writing delta tree..")
	deltafile = options.output + "_delta.nwk"
	lca.write( format=1, outfile=deltafile, format_root_node=True)
	lca_tree = Tree(deltafile, format=1)

	color_tree(lca_tree, size)

	lca_tree.ladderize()

	ts = TreeStyle()
	ts.mode = "r"
	ts.show_leaf_name = False
	ts.scale = 100
	size = 50

	for lin in pango_color_dict:
		ts.legend.add_face(CircleFace(size, pango_color_dict[lin]), column=0)
		ts.legend.add_face(TextFace(lin, fsize = size), column=1)
	ts.legend.add_face(TextFace("last date: " + latest_date.strftime("%Y-%m-%d"), fsize = size), column=1)
	ts.mode = "r"
	lca_tree.render(options.output + "_delta.pdf", tree_style=ts, dpi=72, w=10000, units="mm")
	lca_tree.render(options.output + "_delta.svg", tree_style=ts, dpi=72, w=10000, units="mm")









