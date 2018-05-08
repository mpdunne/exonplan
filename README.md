![alt text](http://empede.co.uk/imgrepos/exonplan_head.png? "ExonPlan - exon maps for amino acid alignments")

What does exonplan do?
==========
Exonplan is a simple python script that elucidates the exon structures of genes in amino acid alignments.

https://github.com/mpdunne/exonplan

Usage
=====

ExonPlan takes a GTF gene coordinate file and an amino acid fasta alignment file, and returns a flattened version, with each exon alternately represented by a single character. This allows the exon structure of genes in an amino acid alignment to be seen relative to each other.

To run ExonPlan, type:

```
python exonplan.py -a your_alignment.fa -g your_gtf.gtf
```

Each entry in the alignment must be represented by a single entry in the GTF file, and vice versa (by ```transcript_id``` tags in the case of the GTFs). The GTF coordinates used are those marked as ```CDS``` (or ```cds```), and you must make sure the inputted coordinate lengths match those of the sequences. Sequence and GTF mismatches will be caught by ExonPlan and flagged up. Bear in mind also that stop codons are often removed by alignment software: you may wish to replace "```*```" symbols with "```X```" before alignment. Terminal stop codons do not affect ExonPlan.

Example files, as well as an example result, are included in the ```example``` directory.

Exonplan is a simple program that has not been thoroughly tested - if something breaks, raise a support request!

Options
=======

ExonPlan uses three symbols to indicate exon location. By default, these are alternately ```A``` and ```G```, and ```Q``` for amino acids that bridge two exons. These can be changed using the ```-l1```, ```-l2```, and ```-l3``` options respectively.

ExonPlan assumes that the first codon in a sequence is a start, and the final codon is a stop. This can be turned off using ```--starts no``` and ```--ends no``` respectively.

There are three output options, specified using ```--outfmt i```, where ```i``` is either: ```0```: output plan only; ```1```: output alignment followed by plan; ```2```: output alignment interleaved with plan.

Example
=======
![alt text](http://empede.co.uk/imgrepos/exonplan_eg.png? "ExonPlan - example alignment")
