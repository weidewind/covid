import optparse
from ete3 import Tree
import pandas as pd

parser = optparse.OptionParser()
parser.add_option('-v', '--variants_file', help='taxon_mutations', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()
print(options)

variants = pd.read_csv(options.variants_file, sep = ",")
unique_entries = list(set(variants["entry"]))

df = pd.DataFrame(columns=["entry","mut","nonmut"])
for e in unique_entries:
	df = df.append({
		"entry": e,
		"mut": variants[(variants["entry"] == e) & (variants["is_double_mut"] == 1)]["taxon"].count(),
		"nonmut":  variants[(variants["entry"] == e) & (variants["is_double_mut"] == 0)]["taxon"].count()
		}, ignore_index=True)
df.sort_values(by='mut', ascending=False, inplace=True)
df.to_csv(options.output, index = False)
