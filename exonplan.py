#usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import csv
import re
import math
import sys
import copy
import string
from Bio import SeqIO
from Bio.Seq import Seq
import math
#for l in w: print l, w[limport random

def cleanGtf(gtflines):
        return sorted(deDup([a for a in gtflines if a[2].lower()=="cds"]), key = lambda x: (math.pow(-1,x[6]=="-"))*int(x[3]))

def deDup(dlist):
        nlist = []
        for i in dlist:
                if i not in nlist:
                        nlist.append(i)
        return nlist

def readCsv(path_file, ignoreBlank = True):
        with open(path_file, "r") as f:
                data = [line for line in csv.reader(f, delimiter="\t") if ''.join(line).strip()] if ignoreBlank else list(csv.reader(f, delimiter="\t"))
        return data

def flattenAlignment(alignedseqs, gtfs):
        aa=["A","G"]
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
                        outstr +=  "A" if triplet == "AAA" else ("G" if triplet == "GGG" else ("-" if triplet == "---" else "Q"))
                outstr = re.sub(r"Q([Q-]*$)", r"", outstr)
                outstr = outstr + (len(s) - len(outstr))*"-"
                outstr = re.sub(r"(-*)E", r"E\1", outstr + "E")
                outstr = re.sub(r"^(-*)[AG]", r"\1S", outstr)
		t = copy.deepcopy(s)
		t.seq = Seq(outstr)
		res.append(t)
        return res, resex

# Accepts an alignment and a gtf file, and spits out the flattened alignment
args  = sys.argv[1:]
if not len(args) > 1:
	print "Please input an alignment file (fasta) followed by a gtf file"
alnin = args[0]
gtfin = args[1]

gtfs = readCsv(gtfin)
seqs = list(SeqIO.parse(alnin, "fasta"))
sids = [s.id for s in seqs]

gtfs_d={}
for line in gtfs:
	tid = re.sub(r".*transcript_id[= ]\"?([^\"]*)\"?.*", r"\1", line[8])
	if tid in sids:
		gtfs_d[tid] = gtfs_d.get(tid, []) + [line]

res, resex = flattenAlignment(seqs, gtfs_d)
for i in res:
	print ">"+i.id+".flat\n"+str(i.seq)
