# HmmScan
Two scripts and example data are added, which can be used to parse hmmscan --domtblout output to retrieve non overlapping profile hits for proteins based on their e value.

Files:

Seqq.txt file contains protein sequences in fasta format, which were used as example data
Hit.txt file contains an example output format of hmmscan --domtblout for protein sequences in Seqq.txt file
HitGitHub.py script is used to parse hmmscan --domtblout output and create HitS.txt file (protein name "\t" protein length "\t" short profile annotation "\t" profile accession "\t" hit e value "\t" hit region from - to (tells the region where on the query protein the profile hit was alligned)), which can be used by OverlapGitHub.py script to check overlaps of profile hits in proteins and select nonoverlapping profile hits (or overlapping at most as indicated in overlap_limit variable). It creates a file NoOverlappingHits.txt with output format:
protein "\t" protein length "\t" short profile annotation "\t" profile accession "\t" e value "\t" region of the protein where profile hit |@| next profile..... (if protein had several nonoverlapping profile hits).
