Short read toolbox.

A collection of scripts for working with genomic data.

**fasta2nuccomp.pl**  - Reads in a fasta file, such as a genome, calculates length, A/T content, G/C content, N content, IUPAC polymorphic content (non A,C,G,T or N) as well as other nucleotides (nucleotides not included in the previous categories).  This information is written to a text file (comma delimited, with header) which is then read into R for plotting.  If you do not have R installed the plotting step will fail.  But you can read the text file into your favarite plotting package and use it as you wish.

*EOF*
