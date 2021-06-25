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



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-r', '--rusmeta', help='', type='str')
parser.add_option('-n', '--nextmeta', help='', type='str')
parser.add_option('-p', '--pangolined', help='', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()


def find_color(nname, default_color):
	lin = pango_dict[nname]
	c = pango_color_dict.get(lin, default_color)
	#lin_up = ".".join(lin.split(".")[:-1])
	#c = pango_color_dict.get(lin, pango_color_dict.get(lin_up, default_color))
	#print(lin_up)
	return(c)


def find_faded_color(nname, default_color, blending_color):
	c = find_color(nname, default_color)
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
		
		diff = maxdate() - date
		mydiff = round((diff.days-20)/30)
		i = 0
		while i < mydiff:
			c = blend(c, blending_color)
			i += 1
		print(",".join([str(mydiff), c]))
		return(c)


def blend(c1, c2):
	if c1[0] != "#":
		c1 = str(webcolors.name_to_hex(c1))
	if c2[0] != "#":
		c2 = str(webcolors.name_to_hex(c2))
	print(c1+" "+c2)
	r = int((int(("0x"+c1[1:3]),16)+int(("0x"+c2[1:3]),16))/0x2)
	g = int((int(("0x"+c1[3:5]),16)+int(("0x"+c2[3:5]),16))/0x2)
	b = int((int(("0x"+c1[5:]),16)+int(("0x"+c2[5:]),16))/0x2)
	return("#"+str(hex(r))[2:].zfill(2)+str(hex(g))[2:].zfill(2)+str(hex(b)[2:].zfill(2)))


def color_tree(t):
	for node in t.traverse():
		nstyle = NodeStyle()
		nstyle["size"] = 0
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		nstyle["vt_line_color"] = "#aaaaaa"
		nstyle["hz_line_color"] = "#aaaaaa"
		if node.is_leaf():
			col = blend(find_color(node.name, "#aaaaaa"),"#aaaaaa")
			nstyle["vt_line_color"] = col
			nstyle["hz_line_color"] = col
			#if node.name == "EPI_ISL_2533828" or pango_dict[node.name] == "B.1.617.2":
			#	node.add_face(TextFace(node.name + " " + pango_dict[node.name], fgcolor="#000000",fsize = 50), column=0)
		if country_dict.get(node.name, "unknown") == "Russia":
			if node.is_leaf():
				col = find_faded_color(node.name, "#777777", "#000000")
				nstyle["fgcolor"] = col
				nstyle["size"] = 4
				if pango_dict[node.name] in ["B.1.1.523", "AT.1"]:
					node.add_face(TextFace(node.name + " " + pango_dict[node.name], fgcolor="#000000",fsize = 50), column=0)
		node.set_style(nstyle)


def maxdate():
	return(datetime.strptime("2021-06-03", "%Y-%m-%d"))

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
rusmeta = pd.read_csv(options.rusmeta, sep="\t")[['Внутренний номер', 'Место забора']]
rusmeta.columns = ['seq_id', 'location']
rusmeta['location'] = 'Russia'
nextmeta = pd.read_csv(options.nextmeta, sep="\t")[['gisaid_epi_isl', 'Location']]
nextmeta.columns = ['seq_id', 'location']
nextmeta['location'] = nextmeta['location'].fillna("unknown")
print(nextmeta.head())
meta = pd.concat([rusmeta, nextmeta])
country_dict = dict(zip(meta['seq_id'], meta['location']))
print(dict(list(country_dict.items())[0:5]))


print("Parsing dates..")
rusmeta = pd.read_csv(options.rusmeta, sep="\t")[['Внутренний номер', 'Дата забора']]
rusmeta.columns = ['seq_id', 'date']
dates_dict = dict(zip(rusmeta['seq_id'], rusmeta['date']))
print(dict(list(dates_dict.items())[0:5]))
 
#B.1.1.451->B.1.1.523
pango_color_dict = {"B.1.1.7":"blue", "B.1.351":"orange", "B.1.617.2":"red", "AT.1":"green", "B.1.1.523":"violet", "B.1.1.317":"lawngreen"}




# print("Parsing regions..")
# regions_dict = get_region_dict()
# print(dict(list(regions_dict.items())[0:5]))

# print("Parsing dates..")
# dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
# dates_dict = dict(zip(dates['seq_id'], dates['date']))
# print(dict(list(dates_dict.items())[0:5]))



print("Drawing trees..")

print("Styling tree..")
color_tree(tree)


ts = TreeStyle()
#ts.mode = "c"
#ts.root_opening_factor = 1
ts.show_leaf_name = False
for lin in pango_color_dict:
	ts.legend.add_face(CircleFace(30, pango_color_dict[lin]), column=0)
	ts.legend.add_face(TextFace(lin, fsize = 50), column=1)
ts.legend.add_face(TextFace("last date: " + maxdate().strftime("%Y-%m-%d"), fsize = 50), column=0) 
#ts.layout_fn = multiple_country_layout
ts.scale = 100
tree.ladderize()
tree.render(options.output + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
tree.render(options.output + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")






