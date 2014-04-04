#!/usr/bin/perl

# All documentation moves to HTS::Tools::Normalize and HTS::Tools::Normalize::Bed.
# This now a wrapper for HTS::Tools::Normalize and HTS::Tools::Normalize::Bed.
# Use perldoc HTS::Tools::Normalize or HTS::Tools::Normalize::Bed to get instructions.
#
# Author : Panagiotis Moulos (pmoulos@fleming.gr)

use strict;
use Getopt::Long;

use HTS::Tools;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "normalize_bed.pl";
our @infile;
our $savrem;
our $rand;
our $sort;
our $log;
our $silent;
our $help;

# Check inputs
&check_inputs;

my $tool = HTS::Tools->new({
        "tool" => "normalize",
        "log" => $log,
        "silent" => $silent,
        "params" => {
            "input" => \@infile,
            "type" => "bed",
            "sumto" => $rand,
            "savrem" => $savrem,
            "sort" => $sort
        }
});
$tool->run;


# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@infile,
           "savrem|t" => \$savrem,
           "rand|r=s" => \$rand,
           "sort|a" => \$sort,
           "silent|s" => \$silent,
           "help|h" => \$help
    );
    # Check if the required arguments are set
    if ($help)
    {
        &program_usage;
        exit;
    }
}

# Print the program usage
sub program_usage
{
    # The look sucks here but it is actually good in the command line
    my $usagetext = << "END";

$scriptname
A perl program to normalize across several samples from Solexa runs
and compensate for the different number of sequence reads for each run.
Input files should be bed format files generated after the genome wide
mapping step. Can run on Linux or Windows.
Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file1 [file2, file3, ..., filen] [OPTIONS]

--- Required ---
  --input|i  file(s)  Input file(s)
--- Optional ---
  --rand|r      The method of random removal of sequences. Use
            "generate" to generate random numbers between 1..#tags
            and remove as many as to reach the minimum tag length,
            "permute" to randomly permute the tags and then remove as
            many as required to reach the minimum tag length or an integer.
            denoting how many tags should be left. If a sample has less
            than this number, it is not touched. Defaults to "generate".
  --savrem|t        Use this option if you wish to save the tags that
            are removed in a different file. The file will be named
            "filename_removed_tags.sav".
  --sort|a      Use this option if you wish to sort the input files
            first. The algorithmic results of the method will not 
            change. Implemented for testing and more options reasons.
  --silent|s        Use this option if you want to turn informative
            messages off.
  --help|h      Display this help text.

The main output of the program is one output BED file for each input BED file
but with random tags removed so as all the output files have the same length.

END
    print $usagetext;
    exit;
}
