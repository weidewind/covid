import optparse
import matplotlib
import pandas as pd
import re



parser = optparse.OptionParser()
parser.add_option('-m', '--melted_foreign', type="str", help='as in /export/home/popova/workspace/covid/output/variant_analysis/3oct_opt_min_newpango/melted_foreign')
parser.add_option('-e', '--entries', type="str", help='pair of entries to count the intersection of closest foreign strains')
options, args = parser.parse_args()
print(options)

melted = pd.read_csv(options.melted_foreign, sep = ";", dtype = {"entry": str})
e1,e2 = options.entries.split(",")
print(e1 + " " + e2)
e1_f = melted[(melted["entry"] == e1)]["strain"]
print(e1_f)
e2_f = melted[(melted["entry"] == e2)]["strain"]
print(str(len(list(set(e1_f) & set(e2_f)))))

