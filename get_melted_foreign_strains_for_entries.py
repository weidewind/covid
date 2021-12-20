import pandas as pd
import re
import datetime
import optparse
import os.path
from ete3 import Tree

parser = optparse.OptionParser()
parser.add_option('-d','--dates_stats', help='dates-stats', type='str')
parser.add_option('-v', '--variants_file', help='', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

dates_stats = pd.read_csv(options.states, sep="\t", names=['seq_id', 'state', 'region'])