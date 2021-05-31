from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
#from Bio.Alphabet import Gapped
#from Bio.Alphabet import IUPAC
#from Bio.Alphabet import ambigious_dna
import argparse
import re
import math

parser = argparse.ArgumentParser(description='Collect the results of global analysis into xlsx')
parser.add_argument("--fasta", dest='fasta', type=str, help="full path to the alignment file")
parser.add_argument("--output", dest='output', type=str, help="full path to the output file")
parser.add_argument("--mutations", dest='mutations', type=str, help="full path to the file with mutations of interest (in genomic coordinates, one mutation  per line)")
parser.add_argument("--ref_ind", dest='ref_ind', type=int, default=0, help="reference sequence index, default 0")
#parser.add_argument("--translate", action='store_true', dest='translate',default=True, help="dna? (prot otherwise)")
args = parser.parse_args()

# todo: class Mutation? parse, tostr..
# todo: a==A


# alignment ind -> site number in seq
def seq_numbering(seq):
	numbering = []
	num = 0
	for letter in list(refseq):
		if (letter != "-"):
			num += 1
		numbering.append(num)
	return(numbering)


def call_mutations(seq, refseq, ref_num):
	mut_list = []
	ins_open = -1
	del_open = -1
	for ind, letter in enumerate(list(seq)): 
		ref_letter = refseq[ind]
		ref_ind = ref_num[ind]
		if (letter != ref_letter):
			if (letter == "-" and del_open == -1):
				del_open = ind
			if (ref_letter == "-" and ins_open == -1):
				ins_open = ind
			if (letter not in ["-", "n", "N"] and ref_letter != "-"):  # todo: process all ambigious letters correctly
				mut_list.append(ref_letter + str(ref_num[ind]) + letter)
		if (del_open != -1 and letter != "-"):
			mut_list.append(del_tostring(ref_num[del_open], ref_num[ind], str(refseq[del_open:ind])))
			del_open = -1
		if (ins_open != -1 and ref_letter != "-"):
			mut_list.append(ins_tostring(ref_num[ins_open], str(seq[ins_open:ind])))
			ins_open = -1
	return(mut_list)


def del_tostring(del_open, del_close, delseq):
	string = "del" + str(del_open)
	if (del_close - 1 > del_open):
		string += "_" + str(del_close - 1)
	trimmed_seq = delseq.replace("-", "")
	string += trimmed_seq
	return(string)


def ins_tostring(ins_open, insseq):
	string = "ins" + str(ins_open)
	trimmed_seq = insseq.replace("-", "")
	string += trimmed_seq
	return(string)

#  todo: add search for overlapping or close indels?
def find_close_variants(mut, sorted_mut_list):
	site = re.match(r"[a-zA-Z]([0-9]+)[a-zA-Z]", mut)
	if site:
		s = site.group(1)
		for m in sorted_mut_list:
			match = re.match(rf"[a-zA-Z]{s}[a-zA-Z]", m)
			if match:
				return match.group(0)

# same subst or indel at the same position and, in case of del, of the same length
def has_exact_match(mut, sorted_mut_list):
	site = re.match(r"[a-zA-Z][0-9]+[a-zA-Z]", mut)
	if site:
		if mut in sorted_mut_list:
			return(True)
	else:
		indel = re.match(r"([a-zA-Z]+)([0-9_]+)([a-zA-Z]*)", mut)
		if not indel:
			print ("invalid mutation format "+ mut)
		itype = indel.group(1)
		iloc = indel.group(2)
		iseq = indel.group(3)
		for m in sorted_mut_list:
			match = re.match(rf"{itype}{iloc}([a-zA-Z]*)", m)
			if match:
				return(True)
	return(False)


align = AlignIO.read(args.fasta, "fasta")

refseq = align[args.ref_ind].seq
ref_numbering = seq_numbering(refseq)
#print(refseq[ref_numbering.index(21563):ref_numbering.index(25385)])
#if (args.translate):
#	align = align[:, ref_numbering.index(21563):ref_numbering.index(25385)]
#	align_prot = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna, "-"))   # Gapped(IUPAC.unambiguous_dna, "-")
#	for rec in align:
		#try:
#		align_prot.add_sequence(rec.id, rec.seq.translate(gap = "-"))
		#except:
			#pass
			#print(rec.id +"\n" + rec.seq + "\n")
#	align = align_prot
#	refseq = align[args.ref_ind].seq
#	ref_numbering = seq_numbering(refseq)
# вырезать, траслировать кусок, создать нормлаьный рпнпрот
#print("ref numbering " + " ".join([str(n) for n in ref_numbering]))

#  todo: convert from nucl to prot automatically
with open(args.mutations) as mf:
	risky_muts = [line.strip() for line in mf]

		
with open(args.output, "w") as out:
	with open(args.output + ".table.csv", "w") as table:
		table.write("seq_id," + ",".join(risky_muts) + "\n")
		for a in align:
			table_str = a.id + ","
			mut_list = call_mutations(a.seq, refseq, ref_numbering)
			out.write(a.id + "\t" + ",".join(mut_list) + "\n")
			for m in risky_muts:
				if has_exact_match(m, mut_list):
					table_str += str(1)
				else:
					close_var = find_close_variants(m, mut_list)  # todo in case of indels, there may be multiple close variaants
					if close_var:
						table_str += close_var
				table_str += ","
			table.write(table_str[:-1] + "\n")


