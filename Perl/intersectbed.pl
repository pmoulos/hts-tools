#!/usr/bin/perl

# All documentation moves to HTS::Tools::Intersect. This now a wrapper for HTS::Tools::Intersect. Use perldoc
# HTS::Tools::Intersect to get instructions.
#
# Author      : Panagiotis Moulos (pmoulos@fleming.gr)

use strict;
use Getopt::Long;

use HTS::Tools;
 			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "intersectbed.pl";
our $fileA; # Input bedfile A 
our $fileB; # Input bedfile B
our $sort; # Sort input files? (Necessary for binary search)
our @percent; # Overlap percentage (not the default)
our $any; # Any overlap (default)
our @extend; # Extend peak modes and use this instead of whole width
our $mode;	# The column where peak mode is in case of extend or distance calculation
our $autoxtend; # Extend peaks based on half of median peak length, peak length column
our $both; # When overlaping with percent, overlap based on both tags or the first?
our $exact; # Percentage intersections exactly as the UCSC Table Browser?
our $maxud; # Number of regions to report in non-overlapping regions
our $agap; # Acceptable gap to say that two peaks from two conditions are close
our $keeporder; # Keep order of input files (Tie::IxHash::Easy)
our $reportonce; # Report overlapping regions only one time
our @out; # Output files
our $multi; # Called from outside to multiply intersect?
our $dryrun; # Dry-run, display stats only
our $waitbar; # Do not use waitbar
our $silent; # Display verbose messages
our $help;	# Help?

# Check inputs
&check_inputs;

my $tool = HTS::Tools->new({
		"tool" => "intersect",
		"silent" => $silent,
		"params" => {
			"inputA" => $fileA,
			"inputB" => $fileB,
			"sort" => $sort,
			"percent" => \@percent,
			"any" => $any,
			"extend" => \@extend,
			"mode" => $mode,
			"autoextend" => $autoxtend,
			"both" => $both,
			"exact" => $exact,
			"maxud" => $maxud,
			"gap" => $agap,
			"output" => \@out,
			"multi" => $multi,
			"keeporder" => $keeporder,
			"reportonce" => $reportonce,
			"dryrun" => $dryrun,
			"waitbar" => $waitbar
	}
});
$tool->run;

# Just parse parameters, checking is now performed by the module
sub check_inputs
{
    my $stop;
    GetOptions("inputA|a=s" => \$fileA,
    		   "inputB|b=s" => \$fileB,
    		   "sort|r" => \$sort,
    		   "percent|p=f{,}" => \@percent,
    		   "any|y" => \$any,
    		   "extend|e=i{,}" => \@extend,
    		   "mode|m=i" => \$mode,
    		   "autoextend|x" => \$autoxtend,
    		   "both|t" => \$both,
    		   "exact|c" => \$exact,
    		   "maxud|n=i" => \$maxud,
    		   "gap|g=i" => \$agap,
    		   "output|o=s{,}" => \@out,
    		   "multi|m" => \$multi,
    		   "keeporder|d" => \$keeporder,
    		   "reportonce|u" => \$reportonce,
    		   "dryrun|z" => \$dryrun,
    		   "waitbar|w" => \$waitbar,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    if ($help)
    {
    	&program_usage;
    	exit;
    }
}

sub program_usage 
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";
	
$scriptname
A pure perl module to calculate intersection between bed files (can contain as many
additional columns as desired, as long as the first 3 are chromosome, start, end). It
requires as input only the bed files and returns the overlap regions between the first
and the second file. However, if fileA and fileB are intersected, the final file containing
name OVERLAP...fileA output will  contain regions from fileA that overlap with regions of
fileB and the final file containing name OVERLAP...B will contain those regions from fileB.
The _ONLY files contain regions only in fileA or only in fileB, when performing  the intersection
of fileA and fileB. The actions performed are similar to those of the UCSC Table Browser, 
and BEDTools intersecting functions, at similar speed and more convenience! This module can be
safely used when the user wants to maintain several precalculated statistics, e.g. the average
number of reads per specified window size under a peak region, something which is not currently
possible using BEDTools. Intersection can be done at specific overlap percentages or any overlap
(like Table Browser, default, 1bp). In addition, one can specify an array of percentages so that
various overlaps and overlap percentages can be calculated in a batch mode to answer questions
like 'is there an overlap saturation between my two peak sets?' Extension from a point (e.g. the
peak mode) towards both directions is also possible as long as a column with this point (the peak
'summit') is given in both files and the appropriate option is used. The user has the option to
retrieve only certain of the four available output file types. Stats are also displayed and the
distances between overlapping and limited non-overlapping regions can be returned. See the help
below for further details.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --inputA fileA --inputB fileB [OPTIONS]

--- Required ---
  --inputA|a  file  First input file
  --inputB|b  file  Second input file
--- Optional ---
  --sort|r		Use this option to sort the input files first. This is
			not necessary but the module runs slighlty faster if the files
			are sorted beforehand.
  --percent|p		The overlap percentage between the two files. It  
			should be a value between 0 and 100. Alternatively, the
			overlap percentage can be a series of values. In this case
			the program will run as a batch procedure, calculating the
			number of overlapping regions for the determined output types
			for each given percentage, creating and writing as output the
			distribution of the number of overlapping regions. In this
			case, no other output files are produced. Instead of a series
			of numbers, it can be an element in the format a:b where a is
			the starting percentage and b the ending. A list between a and
			b sith interval 1 will be autogeneratied. This options does 
			not work with --any option. Examples are --percent 50, 
			--percent 10 20 30 40 50, --percent 1:100.
  --any|y		Use this switch for returning tags with any overlap
			between fileA and fileB. When used, --percent is ignored.
  --extend|e		Use this option to supply the program the upstream and
			downstream extension length. It should be two numbers, the first
			the number of bps to extend upstream of peak center and the second
			to extend downstream. Use this option in combination with --mode
			option. For example --extend 100 200.
  --autoextend|x	Use this switch to calculate half of the median peak 
			length which will then be used for extending from peak mode. Use
			this switch in combination with --mode option to provide the column
			with the peak mode.
  --mode|m		Use this option to supply the program the column in both
			fileA and fileB that contains the peak mode or any point from
			which extension to left and right will start. E.g. --mode 4.
			This parameter has to be provided when --extend and/or
			--autoextend are used. Also when overpairs is chosen as an
			output format.
  --both|t		Normally, when overlaping fileA with fileB at a percentage
			level, if the tag of the fileA overlaps with tag(s) from fileB,
			this percentage will be calculated based on tag from fileA.
			That is, if the required percentage is p and the overlap is d,
			then if d < p*length(tagA) then there is no overlap. However, if
			tagB is smaller than tagA it could d > p*length(tagB), so there
			is p percent overlap between tagA and tagB when using the length
			of tagB. This option allows for such comparisons by checking
			both cases and returning a positive overlap if one of the above
			two is true. This swicth is useless when used with --extend
			option and ignored.
  --exact|c		Use this switch to calculate overlaps based on exact 
  			genomic locations rather than when a region is "totally"
  			included inside the other, when performing overlaps based
  			on percentages. For example if overlaping fileA with fileB,
  			if a region in fileA starts at 300 and ends at 800 and the
  			overlaping region in fileB starts at 500 and ends at 700,
  			if using the --exact swicth, a 50% percent overlap will not
  			return this region as positive because (700-500)<0.5*(800-300)
  			even if tagB is totally included in tagA. This switch forces
  			the program to calculate percent overlaps exactly as the 
  			UCSC Table Browser. Can also be used in combination with the
  			--both switch.
  --maxud|n		Use this option to define the maximum number of
			non-overlapping regions that will be reported using the
			"nonpairs" option of I<output> parameter, within an
			upstream/downstream region defined by the --gap input
			parameter.
  --gap|g		Use this option to retrieve distances between non-
  			overlaping (according to specified criteria) regions 
  			that is not larger than the number specified with the
  			--gap option. Can be useful with onlyA and nonpairs
  			output options. The default gap is set to 0 so the option
  			is not used. If you wish to use that option you should
  			provide a positive integer, for example --gap 10000.
  --output|o		Use this option to determine which intersection output
			filetypes you wish to retrieve.	Possible choices are:
   				"overlapA" for retrieving regions from fileA 
				that overlap with fileB.
   				"overlapB" for retrieving those regions that 
				overlap was found with fileA regions when using 
				fileA as first input (ATTENTION! NOT REGIONS FROM
				fileB THAT OVERLAP WITH fileA).
				"onlyA" for retrieving regions from fileA that
				DO NOT overlap with regions in fileB.
				"onlyB" for retrieving those regions that overlap
				was NOT found with fileA when using fileA as first
				input (ATTENTION! SAME BEHAVIOUR AS "overlapB"
				choice).
				"overpairs" for file with concatenated regions
				from fileA and fileB (similar to the BEDPE format
				of BEDTools) with additional columns that contain
				statistics about the distances of regions found
				to be overlapping between fileA and fileB (the
				statistics depend on the --extend, --autoextend,
				--percent and --mode options. For example, if --mode
				and --extend are given, the distances are from the
				centers of the regions while in any other case,
				all the distances from the start, end and center
				of regions are reported.
				"nonpairs" for a file with concatenated regions
				from fileA and fileB (similar to the BEDPE format
				of BEDTools), containing the distances from
				non-overlapping regions. If the non-overlapping
				regions fall in a range of --gap upstream and
				downstream (where --gap must be specified). The
				maximum number of non-overlapping regions is
				defined by the --maxud parameter.
  --keeporder|d		Use this parameter if you want to force
			the lines of the output files to be in the same order
			(e.g. sorted per chromosome or gene name) as the input
			files. This is accomplished through the module
			Tie::IxHash::Easy which must be present in your
			machine. If the module is not present, the --keeporder
			option is deactivated. Keep in mind that maintaining
			the order requires slighlty more memory during runtime.
  --reportonce|u	Use this option to report only once regions
			from one file that may overlap multiple times from
			regions in the other file (this happens for example
			with BEDTools). Such situations may arise when for
			example there is a broad peak in fileA which is
			split in two peaks in fileB.
  --dryrun|z		Use this option if you wish to do a "dry-run",
			that is just display statistics about chosen overlaps and
			not write any output files.
  --waitbar|w		Use this option if you wish to display a simple  
			progress bar while selecting and printing the final files.   
			For small files it is probably useless as the program finishes
			very quickly.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.

Usage examples:

perl intersectbed.pl --inputA A.bed --inputB B.bed --sort --any --output overlapA 
onlyA --waitbar

perl intersectbed.pl --inputA A.bed --inputB B.bed --percent 50 --mode 4 --extend 100 --both
--output overlapA onlyA onlyB --waitbar
	
The main output of the program is up to four files in BED format containing also
any additional data columns.

END
	print $usagetext;
	exit;
}
