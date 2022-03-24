import optparse
import pandas as pd
import re
from collections import Counter
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, TextFace, TreeNode
from datetime import datetime, timedelta
from meta import get_date_dict
import math



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-u', '--states', help='states file', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-d', '--districts', help='file with ids and districts, output from meta_to_states.py', type='str')
parser.add_option('-a', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('--origins', help='file with ids and travel history, output from meta_to_origins.py', type='str')
parser.add_option('-l', '--lineage', help='', type='str')
parser.add_option('-o', '--output', help='output', type='str')
parser.add_option('-s', '--steps', help='go n nodes up before printing', type='int')
parser.add_option('-m', '--maxdist', help='max number of branches between the entry and the closest nodes', type='int')
parser.add_option('-g', '--pangolined', help='rus.pangolined.withduplicates', type='str')
parser.add_option('--merged_data', help='30jan_pangolin_3.1.20_middle_data.csv', type='str')
#parser.add_option('-e', '--entries_file', help='transmission_lineages.withduplicates.out', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing pangolineages..")
pango = pd.read_csv(options.pangolined, sep=",")
pango_dict = dict(zip(pango['taxon'], pango['lineage']))
print(dict(list(pango_dict.items())[0:5]))

print("Parsing countries..")
countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
country_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(country_dict.items())[0:5]))

print("Parsing districts..")
districts = pd.read_csv(options.districts, sep="\t", names=['seq_id', 'region'])
district_dict = dict(zip(districts['seq_id'], districts['region']))
print(dict(list(district_dict.items())[0:5]))

print("Parsing states..")
states = pd.read_csv(options.states, sep="\s", names=['seq_id', 'foreign_prob', 'rus_prob'])
state_dict = dict(zip(states['seq_id'], states['rus_prob']))
# for k, v in state_dict.items():
#     if len(k)>len("_middle"):
#         print("k[-7:]")
#         print(k[-7:])
#     if len(k)>len("_middle") and k[-7:] == "_middle":
#         print("k and k[:len(k)-7]")
#         print(k)
#         print(k[:len(k)-7])
#         state_dict[k] = state_dict[k[:len(k)-7]]
print(dict(list(state_dict.items())[0:5]))

print("Parsing dates..")
dates_dict = get_date_dict(options.dates)
print(dict(list(dates_dict.items())[0:5]))

print("Parsing origins..")
origins = pd.read_csv(options.origins, sep="\t", names=['seq_id', 'origin'])
origin_dict = dict(zip(origins['seq_id'], origins['origin']))
print(dict(list(origin_dict.items())[0:5]))


# distance in number of branches
def collectClosestSearch(entry_node, maxdist):
    output_nodes = []
    tnode = entry_node
    output_nodes = collectClosest(entry_node, 0, maxdist, output_nodes)
    while(not tnode.is_root() and get_distance_to_parent(child=entry_node, parent=tnode) <= maxdist):
        siss = tnode.get_sisters()
        for s in siss:
            #print("sis " + s.name)
            curdist = get_distance_to_parent(child=entry_node, parent=tnode) + tnode.dist + s.dist
            output_nodes = collectClosest(s, curdist, maxdist, output_nodes)
        tnode = tnode.up
    print("Finished searching for the closest.")
    return(output_nodes)


# distance in number of branches
def collectClosest(node, curdist, maxdist, output_nodes):  # curdepth - current distance to the leaf in question; mindist - best distance so far
    #print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(maxdist)]))
    if node.is_leaf() and curdist <= maxdist:
        output_nodes.append(node.name)
    elif curdist <= maxdist:
        for ch in node.children:
            chdist = curdist + ch.dist
            if chdist <= maxdist: ## add or output_nodes does not contain dates
                output_nodes = collectClosest(ch, chdist, maxdist, output_nodes)
        #print("After traversing all children of node " + node.name + ":")
        #print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_nodes", ",".join([n.name for n in output_nodes])]))
    return(output_nodes)


# distance in number of branches
def get_distance_to_parent(child, parent):
    temp = child
    dist = 0
    while not temp.name == parent.name:
        dist += 1
        temp = temp.up
    return(dist)





def filter_lin(row):
    splitter = re.split(r'[;:]', row['lineage_stat'])
    if options.lineage in splitter:
        return True
    else:
        return False


print("Parsing entry to lineage data..")
merged =  pd.read_csv(options.merged_data, sep="\t")
print(merged[0:5])
entries = merged[merged.apply(filter_lin, axis=1)]["entry"].tolist()
print("Entries:")
print(entries)
mrca = entries[0]
for e in entries:
    print("e:")
    print(e)
    mrca = tree.get_common_ancestor(mrca, e)
    print("mrca " + mrca.name )
print("mrca " + mrca.name )
#closest_nodes = collectClosestSearch(mrca, options.maxdist)
closest_nodes = []
for e in entries:
	closest_nodes.extend(collectClosestSearch(tree.search_nodes(name=e)[0], options.maxdist))
print("closest_nodes: ")
print(closest_nodes)


def collapsed_leaf(node):
    if node.is_leaf():
    	return True
    if hasattr(node, "closest") and node.closest:
        return False
    if hasattr(node, "entry_inside") and node.entry_inside:
    	return False
    if options.lineage not in node.lineages.keys(): 
        return True
    if len(node.lineages.keys()) == 1 and len(node.countries.keys()) < 3:
       return True
    if "Russia" not in node.countries:
    	return True
    else:
        return False


def detach_node(node):
    if collapsed_leaf(node) and not (hasattr(node, "closest") and node.closest) and not (hasattr(node, "entry_outside") and node.entry_outside): 
        return True
#    if node in closest_nodes:
#        return False
#    size = sum(node.countries.values())
 #   if collapsed_leaf(node) and size < 5 and not node.mindate < thres_date and not (len(node.countries.keys()) == 1 and "Russia" in node.countries): # узел сколлапсирован, он маленький, в нем нет ранних данных ии что??
 #       return True
 #   else:
 #       return False
    else:
        return False





def dict_to_string(dict):
    return ";".join(k + ":" + str(v) for k,v in dict.items())


def dict_from_string(string):
    d = {}
    for lin in string.split(";"):
        splitter = lin.split(":")
        d[splitter[0]] = splitter[1]

i = 0
tnode = mrca
while i < options.steps and not tnode.is_root():
    tnode = tnode.up
    i = i + 1

print("upper node: ")
print(tnode.name)


for n in tnode.traverse("postorder"):
    if n.is_leaf():
        co = country_dict[n.name]
        lin = pango_dict.get(n.name, "unknown")
        da = datetime.strptime(dates_dict[n.name], "%Y-%m-%d") if not dates_dict[n.name] == "unknown" else None
        entry_inside = True if n.name in entries else False	
        closest = True if n.name in closest_nodes else False
        n.add_features(countries={co: 1}, mindate=da, maxdate=da, lineages={lin: 1}, entry_inside=entry_inside, closest=closest)
        if co == "Russia":
            di = district_dict[n.name]
            n.add_features(districts={di: 1})
    else:
        countries = {}
        lineages = {}
        districts = {}
        mindate = None #datetime.now() - timedelta(days=100 * 365)
        maxdate = None #datetime.now() + timedelta(days=100 * 365)
        entry_inside = True if n.name in entries else False
        closest = False
        for ch in n.get_children():
            for k, v in ch.countries.items():
                countries[k] = countries.get(k, 0) + v
            for k, v in ch.lineages.items():
                lineages[k] = lineages.get(k, 0) + v
            if mindate is None:
                mindate = ch.mindate
            if maxdate is None:
                maxdate = ch.maxdate
            if ch.mindate is not None and mindate > ch.mindate:
                mindate = ch.mindate
            if ch.maxdate is not None and maxdate < ch.maxdate:
                maxdate = ch.maxdate
            if hasattr(ch, "districts"):
                for k, v in ch.districts.items():
                    districts[k] = districts.get(k, 0) + v
            if hasattr(ch, "entry_inside") and ch.entry_inside:
            	entry_inside = True
            if hasattr(ch, "closest") and ch.closest:
                closest = True
        n.add_features(countries=countries, mindate=mindate, 
                        maxdate=maxdate, lineages=lineages, 
                        districts=districts, entry_inside=entry_inside, closest=closest)

entry_nodes = [tnode.search_nodes(name=e)[0] for e in entries]
print("entry nodes:")
print(entry_nodes)
mind = min([e.mindate for e in entry_nodes if e.mindate is not None])
maxd = max([e.maxdate for e in entry_nodes if e.maxdate is not None])
thres_date = mind + timedelta(days=int((maxd - mind).days/10))
#thres_date = mrca.mindate + timedelta(days=int((mrca.maxdate - mrca.mindate).days/10))
print("days:")
print(str(int((maxd - mind).days/10)))
#tree.write(outfile=options.output + "nwk", features=["mintime", "lineages"])
for e in entry_nodes:
    for c in e.traverse("preorder"):
        c.add_features(entry_outside=True)

#tnode.write(outfile=options.output + ".nwk", is_leaf_fn=collapsed_leaf, format=9, features=["mindate"], format_root_node=True)
#newtree = Tree(options.output + ".nwk", format=9)
#newtree.get_tree_root().name = tnode.name

newtree = tnode.detach()
newtree.write(outfile=options.output + "/" + options.lineage + "_full.nwk", format=8, features=["mindate","closest"], format_root_node=True)

for n in newtree.traverse("postorder"):
    #print("collapsing and detaching: node " + n.name)
    # if n.name.startswith("264770_middle"):
    #     print("Debugging 264770_middle:")
    #     print(n.name)
    #     print(n)
    #     print("collapsed?")
    #     print(collapsed_leaf(n))
    #     print("detached?")
    #     print(detach_node(n))
    if detach_node(n):
    	n.detach()
    elif collapsed_leaf(n):
        for ch in n.get_children():
            ch.detach()
    else:
        dupls = [ch for ch in n.get_children() if ch.is_leaf() and  ch.dist == 0]
        if len(dupls) > 3:

            countries = {}
            lineages = {}
            districts = {}
            mindate = None #datetime.now() - timedelta(days=100 * 365)
            maxdate = None #datetime.now() + timedelta(days=100 * 365)
            entry_inside = True if n.name in entries else False
            closest = False
            for ch in dupls:
                for k, v in ch.countries.items():
                    countries[k] = countries.get(k, 0) + v
                for k, v in ch.lineages.items():
                    lineages[k] = lineages.get(k, 0) + v
                if mindate is None:
                    mindate = ch.mindate
                if maxdate is None:
                    maxdate = ch.maxdate
                if ch.mindate is not None and mindate > ch.mindate:
                    mindate = ch.mindate
                if ch.maxdate is not None and maxdate < ch.maxdate:
                    maxdate = ch.maxdate
                if hasattr(ch, "districts"):
                    for k, v in ch.districts.items():
                        districts[k] = districts.get(k, 0) + v
                if hasattr(ch, "entry_inside") and ch.entry_inside:
                    entry_inside = True
                if hasattr(ch, "closest") and ch.closest:
                    closest = True
            dupls_node = n.add_child(name=n.name + "_identical")
            dupls_node.add_features(countries=countries, mindate=mindate, 
                            maxdate=maxdate, lineages=lineages, 
                            districts=districts, entry_inside=entry_inside, closest=closest, duplicates=True)
            for ch in dupls:
                ch.detach()


newtree.write(outfile=options.output + "/" + options.lineage + ".nwk", format=8, features=["mindate"], format_root_node=True)

def shade_of_grey(prob):
    if prob < 0.1:
        prob = 0.1
    t = str(hex(int(256*(1-prob)))[2:].zfill(2))
    return("#"+t+t+t)


def shade_of_green(lineage, lin_dict):
    prop = lin_dict.get(lineage, 0)/sum(lin_dict.values())
    t = str(hex(int(48 + prop * (255 - 48)))[2:].zfill(2))
    return("#30"+t+"30")


def style_tree(t):
    for node in t.traverse():
        nstyle = NodeStyle()
        if len(node.name)>len("_middle") and node.name[-7:] == "_middle":
            name = node.name[:len(node.name)-7]
            print("my name is")
            print(name)
        else:
            name = node.name
        state_prob = state_dict.get(name, node.countries.get("Russia",0)/sum(node.countries.values()))
        nstyle["fgcolor"] = shade_of_grey(state_prob)
        nstyle["vt_line_width"] = 2
        nstyle["hz_line_width"] = 2
        nstyle["vt_line_color"] = shade_of_grey(state_prob)
        nstyle["hz_line_color"] = shade_of_grey(state_prob)
        text = ""
        size = sum(node.lineages.values())
        if node.name in entries:
            nstyle["fgcolor"] = shade_of_green(options.lineage, node.lineages)
            nstyle["size"] = 2+2*math.log(size,2)
            text = node.name
            node.add_face(TextFace(node.name, fgcolor="#000000", fsize=20), column=0) #fsize = size
        if node.name == mrca.name:
            node.add_face(TextFace("MRCA " + node.name, fgcolor="#000000"), column=0) 
        if node.is_leaf():
            nstyle["size"] = 2+2*math.log(size,2)
            nstyle["fgcolor"] = shade_of_green(options.lineage, node.lineages)
            if size > 1:
            	text = text + " " + str(size)
            if node.name in origin_dict:
                text = text + " " + origin_dict[node.name]
            if node.mindate is not None:
                mindate = datetime.strftime(node.mindate, "%Y-%m-%d")
            else:
                mindate = ""
            if node.mindate is not None and node.mindate < thres_date: # datetime.strptime("2021-03-10", "%Y-%m-%d"):
            	text = text + " " + mindate
            if not (len(node.lineages.keys()) == 1 and options.lineage in node.lineages):
                text = text + " " + ";".join(k + ":" + str(v) for k,v in node.lineages.items())
            if not (len(node.countries.keys()) == 1 and "Russia" in node.countries):
                text = text + " " + ";".join(k + ":" + str(v) for k,v in node.countries.items())
            if hasattr(node, "duplicates") and node.duplicates:
                nstyle["shape"] = "square"

  #          if state_dict[node.name] == 1:
   #                 text = text + " " + district_dict[node.name]
   #         else:
   #                 text = text + " " + country_dict[node.name]

            fgc = "black"
            if node.mindate is not None and node.mindate < thres_date: #datetime.strptime("2021-03-10", "%Y-%m-%d"):
                fgc = "orange"
            if node.name in origin_dict:
                fgc = "red"
                
            node.add_face(TextFace(text,  fgcolor=fgc), column=0) #fsize = size
        node.set_style(nstyle)


ts = TreeStyle()
#ts.mode = options.mode
#ts.root_opening_factor = 1
ts.show_leaf_name = False
#if options.mode == "r":
#    ts.scale = 10

ts.allow_face_overlap = True
style_tree(newtree)


newtree.ladderize()
newtree.render(options.output + "/" + options.lineage + ".svg", tree_style=ts, dpi=72, w=10000, units="mm")
newtree.render(options.output + "/" + options.lineage + ".pdf", tree_style=ts, dpi=72, w=10000, units="mm")