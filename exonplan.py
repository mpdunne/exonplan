#usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#   Copyright Michael Dunne 2018
#
#   www.github.com/mpdunne/exonplan
#
########################################

########################################################
################ Safely import packages ################
########################################################

errors = []
libnames = ["csv", "re", "sys", "math", "copy", "string", "argparse"]

for libname in libnames:
    try:
        lib = __import__(libname)
    except ImportError as e:
        errors.append(e)
    else:
        globals()[libname] = lib

try:
        from Bio import SeqIO
        from Bio.Seq import Seq
except ImportError as e:
        errors.append(e)

if errors:
        errmsg = "Missing modules :(\nThe following module errors need to be resolved before running OMGene:"
        for error in errors: errmsg += "\n-- " + str(error)
        sys.exit(errmsg)

########################################################
###################### The Meat ########################
########################################################


def cleanGtf(gtflines):
	return sorted(deDup([a for a in gtflines if a[2].lower()=="cds"]), key = lambda x: (math.pow(-1,x[6]=="-"))*int(x[3]))

def deDup(dlist):
	nlist = []
	for i in dlist:
		if i not in nlist: nlist.append(i)
	return nlist

def readCsv(path_file, ignoreBlank = True):
	with open(path_file, "r") as f:
		data = [line for line in csv.reader(f, delimiter="\t") if ''.join(line).strip()] if ignoreBlank else list(csv.reader(f, delimiter="\t"))
	return data

def flattenAlignment(alignedseqs, gtfs, l1, l2, l3, starts, ends):
	aa=[l1,l2]
	res = []; resex = []
	for s in alignedseqs:
		if not s.id in gtfs:
			sys.exit("Error: sequence " + s.id + " is not present in the gtf file. Make sure transcript_id tag matches the sequence ID.")
		gtf = cleanGtf(gtfs[s.id])
		if not gtf: continue
		seqExpanded = string.join([a*3 for a in str(s.seq)], "")
		pos = 0; flat = ""; totalLen = 0
		numlines = len(gtf)
		strand = gtf[0][6]
		for i, line in enumerate(gtf):
			cdslen = int(line[4]) - int(line[3]) + 1
			totalLen += cdslen
			if totalLen > len(seqExpanded):
				sys.exit("Error: sequence " + s.id + " does not match gtf coordinates. The sequence is too short.\nDoes the aligned sequence include stop codons?")
			# The stop codon won't be included in the alignment.
			if i == numlines - 1:
				# partial codons are possible...
				cdslen = cdslen - (((totalLen-1) % 3)+1)
			localpos = 0
			localres = []
			while localpos < cdslen:
				if seqExpanded[pos] != "-":
					localres.append(pos)
					flat+=aa[i%2]
					localpos += 1
				else:
					flat += "-"
				pos += 1
		t = copy.deepcopy(s)
		t.seq = Seq(flat)
		resex.append(t)
		outstr = ""
		for qpos in range(len(s)):
			pos=3*qpos
			triplet =flat[pos:pos+3]
			outstr +=  l1 if triplet == l1*3 else (l2 if triplet == l2*3 else ("-" if triplet == "---" else l3))
		outstr = re.sub(r""+l3+"(["+l3+"-]*$)", r"", outstr)
		outstr = outstr + (len(s) - len(outstr))*"-"
		if ends:
			outstr = re.sub(r"(-*)E", r"E\1", outstr + "E")
		if starts:
			outstr = re.sub(r"^(-*)["+l1+l2+"]", r"\1S", outstr)
		t = copy.deepcopy(s)
		t.seq = Seq(outstr)
		res.append(t)
	return res, resex

def checkSequences(seqs, gtfs, factor):
	errors = []
	for tid in gtfs:
		s = seqs[tid]
		g = cleanGtf(gtfs[tid])
		totalLen = 0
		for i, line in enumerate(g):
			cdslen = int(line[4]) - int(line[3]) + 1
			totalLen += cdslen
		faclen = len(str(s.seq).replace("-",""))*factor
		if not totalLen in [faclen, faclen + factor]:
			errors += [s.id]
	if errors:
		err = ""
		for e in errors: err+="\n"+e
		sys.exit("Error: The following sequence lengths do not match their GTF sequences:" + err)


########################################################
##################### Entry code #######################
########################################################

# Accepts an alignment and a gtf file, and spits out the flattened alignment
parser = argparse.ArgumentParser(description="Run ExonPlan")
parser.add_argument("-a", "--alignment", metavar="alignment", dest="AL", default="")
parser.add_argument("-g", "--gtf", metavar="gtffile", dest="GT", default="")
parser.add_argument("-l1", dest="L1", default="A")
parser.add_argument("-l2", dest="L2", default="G")
parser.add_argument("-l3", dest="L3", default="Q")
parser.add_argument("--outfmt", dest="OUTF", default="0")
parser.add_argument("-s", "--starts", dest="ST", default="yes")
parser.add_argument("-e", "--ends", dest="EN", default="yes")

# Parse the arguments...
args = parser.parse_args()

alnin  = args.AL
gtfin  = args.GT
l1     = args.L1
l2     = args.L2
l3     = args.L3
starts = args.ST
ends   = args.EN
outfmt = args.OUTF

# Sanity check for some of the options...
if not alnin or not gtfin:
	sys.exit("Please input a FASTA alignment file (-a) and a GTF file (-g)")

if any(not re.match(r"[A-Za-z]", i) for i in [l1,l2,l3]):
	sys.exit("Error: Symbol choice must be an alphabetic character (a-z or A-Z)")

if not re.match(r"[012]", outfmt):
	sys.exit("Error: Output format (--outf) must be either 0,1 or 2.")
outfmt = int(outfmt)

if not (starts in ["yes", "no"] and ends in ["yes", "no"]):
	sys.exit("Error: -s and -e options must be \"yes\" or \"no\".")

# Read the GTF and FASTA files and double-check that they at least look reasonably sound.

seqs   = list(SeqIO.parse(alnin, "fasta"))
sids   = [s.id for s in seqs]
seqs_d = dict((s.id, s) for s in seqs)

gtfs   = readCsv(gtfin)
gtfs_d={}
for line in gtfs:
	tid = re.sub(r".*transcript_id[= ]\"?([^\"]*)\"?.*", r"\1", line[8])
	if tid in sids:
		gtfs_d[tid] = gtfs_d.get(tid, []) + [line]

# There should be exactly one FASTA tid for each GTF tid.

if not len(set(sids)) == len(sids):
	sys.exit("Error: FASTA file contains duplicate sequence IDs")

if not sorted(gtfs_d) == sorted(sids):
	sys.exit("Error: Sequence IDs in FASTA alignment file do not match those in\n	 the GTF file. Check and try again :)")

# Each sequence should match the length of its gtf entry.
checkSequences(seqs_d, gtfs_d, 3)

# Go!
res, resex = flattenAlignment(seqs, gtfs_d, l1, l2, l3, starts == "yes", ends == "yes")

# Output according to output format choice
if outfmt == 0:
	for i in res: print ">"+i.id+".flat\n"+str(i.seq)
elif outfmt == 1:
	for i in res: print ">"+i.id+"\n" + str(seqs_d[i.id].seq)
	for i in res: print ">"+i.id+".flat\n"+str(i.seq)
elif outfmt == 2:
	for i in res:
		print ">"+i.id+"\n" + str(seqs_d[i.id].seq)
		print ">"+i.id+".flat\n"+str(i.seq)
