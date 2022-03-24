import pandas as pd
import optparse
import datetime


parser = optparse.OptionParser()
parser.add_option('-s', '--states', help='leaf states file', type='str')
parser.add_option('-d', '--dates', help='leaf dates file', type='str')
parser.add_option('-o', '--output', help='output file', type='str')

options, args = parser.parse_args()

dates = pd.read_csv(options.dates, sep="\t")
dates.columns = ['id', 'date']
known = dates.loc[dates['date'] != "unknown"]
unknown_number = dates.loc[dates['date'] == "unknown"].shape[0]
known['date'] = pd.to_datetime(known['date'])

print(known['date'][:10])
ax = known['date'].hist(bins = 100, figsize=(10, 8))
ax.set_title("date known:" + str(known.shape[0]) + ", date unknown:"+ str(unknown_number) + ", min:" + str(known['date'].min()) + " max:" + str(known['date'].max()))
fig = ax.get_figure()
fig.savefig(options.output)

