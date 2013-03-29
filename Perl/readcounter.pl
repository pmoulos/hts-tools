#!/usr/bin/perl

# All documentation moves to HTS::Tools::Count. This now a wrapper for HTS::Tools::Count. Use perldoc
# HTS::Tools::Count to get instructions.
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)

use strict;
use Getopt::Long;

use lib '/media/HD4/Fleming/hts-tools/HTS-Tools/lib';
use HTS::Tools;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "readcounter.pl";
our @input; # Input file(s)
our $regionfile;	# Genomic regions file
our $sort = 0;	 # Sort files
our $percentage = 0.95; # Percentage of tag inside genomic region
our $lscore; # Linear scoring to determine if tag partially inside genomic region included
our $escore; # Exponential scoring to determine if tag partially inside genomic region included
our $c = 3; # Default constant to be used in the exponential scoring scheme
our $countsmall;	# Count tags in genomic regions that are smaller than the tag length?
our $splitper; # Split genomic regions in sub-regions of xbps and count individually
our $nbins; # Split genomic regions in sub-regions of xbps and count individually
our $ncore = 1; # Number of cores for parallel calculations
our $stats; # Return some statistics about tag distributions per area
our $keeporder; # Keep the order of input region file
our $output; # The output file, can be auto
our $source; # Source in case of downloading data
our $splicing; # Splicing in case of downloading data from UCSC or RefSeq
our $silent = 0; # Display verbose messages
our $help = 0; # Help?

&check_inputs;

my $tool = HTS::Tools->new({
		"tool" => "count",
		"silent" => $silent,
		"params" => {
			"input" => \@input,
			"region" => $regionfile,
			"source" => $source,
			"splicing" => $splicing,
			"sort" => $sort,
			"percent" => $percentage,
			"lscore" => $lscore,
			"escore" => $escore,
			"constant" => $c,
			"small" => $countsmall,
			"split" => $splitper,
			"nbins" => $nbins,
			"stats" => $stats,
			"ncore" => $ncore,
			"keeporder" => $keeporder,
			"output" => $output
		}
});
$tool->run;

# Just parse parameters, checking is now performed by the module
sub check_inputs
{
    GetOptions(
		"input|i=s{,}" => \@input,
		"region|r=s" => \$regionfile,
		"source|u=s" => \$source,
		"splicing|x=s" => \$splicing,
		"sort|a" => \$sort,
		"percent|p=f" => \$percentage,
		"lscore|l" => \$lscore,
		"escore|e" => \$escore,
		"constant|c=f" => \$c,
		"small|m" => \$countsmall,
		"split|t=i" => \$splitper,
		"nbins|b=i" => \$nbins,
		"stats|z" => \$stats,
		"ncore|n=i" => \$ncore,
		"keeporder|k" => \$keeporder,
		"output|o=s" => \$output,
		"silent|s" => \$silent,
		"help|h" => \$help
	);
    if ($help)
    {
    	&program_usage;
    	exit;
    }
}

# Print the program usage
sub program_usage
{
	my $usagetext = << "END";

A perl program to join reads in bed file(s) with their genomic regions
in a genome and count presences of reads (uniquely) in each gene. It
simplifies the operation on genomic intervals of Galaxy and solw by
performing all operations locally and unified. Additionally, it
addresses the problem of partial overlaping between a read and a
genomic regions by 3 possible ways: an overlaping coefficient, a
linear probabilitistic score and an exponential probabilistic score.
Please see program description in the header of the script for further
details and file format descriptions.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file1 [file2, file3, ..., filen] --region file [OPTIONS]

--- Required ---
  --input|i  file(s)  Input file(s). The input bed files can be 3 or 6
  			column bedfiles.
  --region|r file Genomic regions file. The first 3 columns should
  			be of the same structure as in bed files. The 4th column
  			should contain a UNIQUE identifier for each region (e.g.
  			Ensembl IDs). The rest columns can contain additional data
  			about the regions. Instead of a local file, it can be one
  			of the following, which will automatically download and use
  			genomic annotations from the latest version of --source:
  			"human-gene" for homo sapiens gene co-ordinates
			"human-exon" for homo sapiens exon co-ordinates
			"human-5utr" for homo sapiens 5'UTR co-ordinates
			"human-3utr" for homo sapiens 3'UTR co-ordinates
			"human-cds" for homo sapiens CDS co-ordinates
			"mouse-gene" for mus musculus gene co-ordinates
			"mouse-exon" for mus musculus exon co-ordinates
			"mouse-5utr" for mus musculus 5'UTR co-ordinates
			"mouse-3utr" for mus musculus 3'UTR co-ordinates
			"mouse-cds" for mus musculus CDS co-ordinates
			"rat-gene" for rattus norvegicus gene co-ordinates
			"rat-exon" for rattus norvegicus exon co-ordinates
			"rat-5utr" for rattus norvegicus 5'UTR co-ordinates
			"rat-3utr" for rattus norvegicus 3'UTR co-ordinates
			"rat-cds" for rattus norvegicus CDS co-ordinates
			"fly-gene" for drosophila melanogaster gene co-ordinates
			"fly-exon" for drosophila melanogaster exon co-ordinates
			"fly-5utr" for drosophila melanogaster 5'UTR co-ordinates
			"fly-3utr" for drosophila melanogaster 3'UTR co-ordinates
			"fly-cds" for drosophila melanogaster CDS co-ordinates
			"zebrafish-gene" for danio rerio gene co-ordinates
			"zebrafish-exon" for danio rerio exon co-ordinates
			"zebrafish-5utr" for danio rerio 5'UTR co-ordinates
			"zebrafish-3utr" for danio rerio 3'UTR co-ordinates
			"zebrafish-cds" for danio rerio CDS co-ordinates
--- Optional ---
  --source|u		Use this option to set the online data source in
			the case of selecting one of the prefefined region templates
			with --region. Can be one of "ucsc", "refseq" or "ensembl".
			Defaults to "ensembl".
  --sort|a		Use this option if you wish to sort the input files
			first. This is not obligatory as the new implementation
			using the interval trees structure does not require sorting,
			however, it will run faster if the input files are already
			sorted. If you wish to sort the files outside the script,
			in Linux systems, the command should look like this:
				sort -k1,1 -k2g,2 inputfile > outputfile 
			In Windows systems, it would be better to use a spreadsheet
			program like Excel to perform sorting by columns. Both the
			input bed files and the region file should be sorted. If
			the structure of the region file is as described (1st column:
			chromosome, 2nd column: start 3rd column: end, 4th column:
			unique ID etc.) the sorting command is exactly the same as
			in common bed files (3 or 6 columns). Keep in mind that
			in future versions, the --sort parameter will be deprecated.
  --percent|p		Use this option to provide the overlaping
  			coefficient according to which tags that partially fall
  			outside provided genomic regions will be assigned to
  			genomic regions. Defaults to 0.95 and is the default
  			algorithm for assigning tags that fall partially inside
  			provided genomic regions.
  --lscore|l		Use this option to use the linear probabilistic
  			scoring scheme to address the problem of partial tag
  			overlap with the provided genomic regions. Please see
  			header of the script file for description. Defaults to 0
  			(not use).
  --escore|e		Use this option to use the exponential probabilistic
  			scoring scheme to address the problem of partial tag
  			overlap with the provided genomic regions. Please see
  			header of the script file for description. Defaults to 0
  			(not use).
  --constant|c		Use this option to provide the constant for the
  			exponential scoring scheme (see description of --escore
  			option). Defaults to 3.
  --small|m		Use this option if you wish to take into consideration
  			genomic regions in your genomic regions file which are
  			smaller than the tag length (rare but may occur, especially
  			with customized genomic regions). Defaults to 0 (not use).
  --split|t		Use this option if you wish to further split your genomic
  			regions in smaller areas. In this case, tags will be counted
  			per area and the distribution of counts per area will be
  			returned as a column in the output file. After --split,
  			you can provide the length of sub-areas. The default is
  			1000 and splits the regions per 1kb, unless the regions are
  			smaller. In this case, the area will consist of only one
  			sub-area of length equal to that of the area. Note also
  			that when using this option, tags that are found inside
  			sub-areas are not assigned to those sub-areas based on
  			scoring schemes (options --percent, --lscore and --escore)
  			but tags are assigned based on the location of their center.
  --nbins|b		Use this option if you wish to further split your
			genomic regions in a predefined number (nbins) of smaller
			areas. The same things as in --split apply.
  --stats|z		Use this option to also return basic statistics of
			counts in the windows used returned by using --split.
  --ncore|n		If the machine has multicore processor(s) and the
			module Parallel::ForkManager is installed, you can use
			parallel processing. Default is 1 and can go up to 12.
  --keeporder|k		Use this parameter if you want to force the lines
			of the output counts table to be in the same order (e.g.
			sorted per chromosome or gene name) as the region file. This
			is accomplished through the use of the module Tie::IxHash::Easy
			which must be present in your machine. If the module is not
			present, --keeporder is deactivated. Keep in mind that
			maintaining the order requires slighlty more memory during
			runtime.
  --output|o		A file to write the output to. If "auto", then it
			generates an automatic filename in the folder where the
			input files are. If not provided, output is written to STDOUT..
  --silent|s		Use this option if you want to turn informative
  			messages off.
  --help|h		Display this help text.

The program returns bed file(s) that contain the provided genomic regions
with any other information together with read counts per region.

END
	print $usagetext;
	exit;
}
