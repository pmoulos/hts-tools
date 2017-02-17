# hts-tools
A legacy Perl library for basic analysis of NGS ChIP-Seq data

## Summary
This library comprises some legacy ChIP-Seq data analysis scripts assembled and
organized into a module. Most are outdated and their functionality covered by 
much faster and modern tools such as [BEDTools](https://github.com/arq5x/bedtools2) 
etc. In addition, none of these (apart from track conversion tools, aimed to 
convert among various UCSC Genome Browser visualization formats) supports the
standard for years now BAM format. 

However, some of these legacy tools may still prove useful for some functionalities,
such as the BED interesection utility (```Intersect.pm```, wrapped by ```intersectbed.pl```
and ```multisectbed.pl``` scripts). This specific one, although orders of magnitude 
slower than ```bedtools intersect```, offers some convenience not present in
bedtools, such as carrying all columns at all operations (very useful if someone
has a file with ChIP-Seq peaks, accompanied by several additional peak metrics),
and several other goodies (present in the wrapper script's help). ```multisectbed.pl```
also automatically created all possible outputs and Venn diagrams. The module
includes even a pure Perl implemented parallel read counter, ideal for romantics!

Tools that are worthy looking to are ```Assign.pm``` (wrapped by ```hyperassignpeaks.pl```),
```Fetch.pm```, ```MotifScan.pm``` and ```Normalize.pm```. You will find most
useful also ```Normalize.pm``` (wrapped by ```normalize_bedgraph.pl``` and
```normalize_bed.pl```) for BED and BedGraph files normalization (the latter being
especially useful for RNA-Seq tracks). Most have wrapper scripts which facilitate
their usage and at the same time they are Object Oriented which allows their
incorporation in your personal scripts. Each module and wrapper script has
extensive help and documentation.
