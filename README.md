# Exon plan - add exon data to your alignments

What does exonplan do?
==========
Exonplan is a short python script that generates flattened versions of inputted alignment data, indicating exon junctions.

https://github.com/mpdunne/exonplan

Usage
=====

Exonplan takes as input a fasta alignment file and a set of gtfs containing the coordinates for the genes listed in the alingment file. IDs used in the alignment file must match `transcript_id` tags in the gtf gile. To use exonplan, simple call the following in the terminal:

```
python exonplan.py path_to_aln.fa path_to_gtf.gtf
```

