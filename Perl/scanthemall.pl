#!/usr/bin/perl

# All documentation moves to HTS::Tools::Motifscan. This now a wrapper for HTS::Tools::Assign. Use perldoc
# HTS::Tools::Motifscan to get instructions.
#
# Author : Panagiotis Moulos (pmoulos@fleming.gr)

# TODO List
# Incorporate support for MotifLocator
# Put an option to get background from same chromosome
# Automate background generation (somewhen...)

use strict;
use Getopt::Long;

use lib '/media/HD4/Fleming/hts-tools/HTS-Tools/lib';
use HTS::Tools;
 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "scanthemall.pl";
our @seqfile; # The peak-region-whatever file(s) (FASTA format)
our $motfile; # A file with PWMs 
our $backfile; # A possible background regions file to choose from
our @cntfile; # Original peak files, containing the centers in case of peak ids in fasta
our @cntcol; # Columns with peak IDs and peak centers in cntfile
our @range;	# Range of cutoffs to use with p53Scan
our $scanner; # Which scanner to use
our $fpr; # Allowed false positive rate to be used with p53Scan
our $times; # How many times larger background?
our $length; # Length for background sequences
our $besthit; # How many best matches to return (default 1)
our @output; # gff? bed? stats? report?
our $unistats; # In stats files, return unique numbers?
our $justscan; # Just scan given only one threshold, no background
our $log; # Keep log?
our $silent; # Display verbose messages
our $help; # Help?

# Check inputs
&check_inputs;

my $tool = HTS::Tools->new({
		"tool" => "motifscan",
		"log" => $log,
		"silent" => $silent,
		"params" => {
			"input" => \@seqfile,
			"motif" => $motfile,
			"background" => $backfile,
			"scanner" => $scanner,
			"center" => \@cntfile,
			"colext" => \@cntcol,
			"range" => \@range,
			"fpr" => $fpr,
			"times" => $times,
			"length" => $length,
			"besthit" => $besthit,
			"output" => \@output,
			"uniquestats" => $unistats,
			"justscan" => $justscan
		}
});
$tool->run;

# Just parse parameters, checking is now performed by the module
sub check_inputs
{
	GetOptions("input|i=s{,}" => \@seqfile,
    		   "motif|m=s" => \$motfile,
    		   "background|b=s" => \$backfile,
    		   "scanner|a=s" => \$scanner,
    		   "center|c=s{,}" => \@cntfile,
    		   "colext|x=i{,}" => \@cntcol,
    		   "range|r=s{,}" => \@range,
    		   "fpr|p=f" => \$fpr,
    		   "times|n=i" => \$times,
    		   "length|w=i" => \$length,
			   "besthit|e=i" => \$besthit,
    		   "output|o=s{,}" => \@output,
    		   "uniquestats|u" => \$unistats,
    		   "justscan|j=f" => \$justscan,
    		   "log|l" => \$log,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&program_usage;
    	exit;
    }
}

sub programUsage
{
	my $usagetext = << "END";

$scriptname
This script is a wrapper and significance threshold optimizer for several 
sequence motif scanners. At the present point, the scanners supported are 
pwmscan from the GimmeMotifs suite (van Heeringen and Veenstra, 2010) and 
MotifScanner (Thijs et al., 2001). Thus, the installation of these 3rd party 
tools in the system that uses HTS::Tools::Motifscan is necessary (or at 
least of one of them). Then, the module runs a workflow which which massively 
scan the input fasta files (e.g. a set of sequences corresponding to ChIP-Seq 
peaks) for the motifs given in the motif file in the form of positional weight 
matrices (e.g. motifs that have been characterized as significant during a de 
novo motif search in the same set of peaks). This workflow determines an 
optimal significance threshold based on a set of random background sequences 
selected from the background fasta file provided and on the false postive rate 
(fpr) provided. The script provides three basic outputs: a simple file with 
the number of hits of each motif in each input file, gff files with the hits 
and bed files with the hits, including also a significance score of the scan 
hit in the 5th column of the bed file.

Author : Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname --input peak_fasta_file(s) --region regfile --background backfile [OPTIONS]

--- Required ---
  --input|i  file(s)	A set of input FASTA file(s) that will 
			be used for matching of the motifs contained in --motif
			file. Each of these files will be scanned using the 
			scanner defined by --scanner and the hits will be 
			reported either as a simple tab-delimited stats file, 
			a gff file, or a bed file under conditions (see below).
  --motif|m  file	A file containing several motifs in the Position 
			Weight Matrix (PWM) format. The file will be decomposed
			in several motifs and each motif will be matched against 
			each input FASTA file.
--- Optional ---
  --background|b		A FASTA file containing background sequences 
			(must be at least the same length as the input sequences)
			that will be used for random sampling to  define a score 
			cutoff for motif matching significance based on --fpr. It
			is required if the option --justscan is not activated.
  --scanner|a		A scanning algorithm to use. Currently, two 
			algorithms are supported: the pwmscan algorithm (van 
			Heeringen and Veenstra, 2010) and MotifScanner (Thijs 
			et al., 2001).
  --range|r		A range of cutoffs that will be used to determine 
			the final matching score cutoff corresponding to I<fpr>.
			The range can be given in the form a:b or a:x:b where x 
			is an increment step. If the format a:b is chosen, the 
			default increment is 1. However, in the latest versions 
			of pwmscan, the matching score is normalized to one, so 
			this notation will be eventually deprecated.
  --fpr|p		The desired False Positive Rate (FPR), that is the 
			percentage of motif matches found in the background
			sequences. Defaults to 0.05.
  --times|n		How many times should the background set of sequences 
			that will be used for the determination of the cutoff
			score from I<range>, be larger than the input set? It 
			defaults to 10, which means that if an input FASTA file
			contains 100 sequences, the background will contain 1000 
			random sequences.
  --length|w		The length of the background sequences. Generally, 
			it should be larger than the length of input sequences.
			It defaults to 400.
  --output|o		The output types that the user wished to get. It 
			can be one or more of "stats" for a simple delimited
			file containing the hits for each motif and input file, 
			"gff" for the gff output of pwmscan or a related file 
			from MotifScanner, containing the actual hits and 
			positions in the input sequences, or "bed" for an
			output BED file containing the motif matches locations 
			and a score to be used for coloring. The "bed" output is 
			available only if a set of peak files is given so as to 
			determine the relative location of the match inside the 
			peak and construct proper bed lines.
			Example: --output stats gff
  --besthit|e		The number of best hits to be retrieved when --scanner
			is "pwmscan". Defaults to 1.
  --uniquestats|e		If the number of besthits is greater than 1,
			the specifying I<uniquestats> to 1 will cause the output
			"stats" file to contain unique hits. Like this you can 
			avoid paradoxes like having more hits than input FASTA 
			sequences. However, in some situations you might actually 
			want to retrieve multiple hits.
  --justscan|j		Set this to 1, to just scan the input sequences 
			for the motifs without defining an FPR. Useful if you
			have a very limited number of input sequences (e.g. just 
			one).
  --center|c		A set of genomic regions (e.g. peaks) with the
			SAME IDs as the input FASTA sequences (or a superset of
			these) which will be used to assign peak regions to FASTA
			files in order to determine the proper coordinates for the
			generation of BED output. It is optional, however, if BED
			output is requested, these files are not given, and the 
			FASTA sequence length does not correspond to the length
			that can be extracted by the FASTA ID, the coordinates
			will be inaccurate.
  --colext|x		A vector of length 3, containing the column numbers
			of peak ID and peak summit, and the length of the (possible)
			extension upstream and downstream of the peak summit. For
			example --colext 4 5 75.
  --log|l		Output a log file. It can be a file name or empty for
				auto-generation.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.
	
The main output of the program is up to nine files with information on gene-peak
association.

END
	print $usagetext;
	exit;
}
