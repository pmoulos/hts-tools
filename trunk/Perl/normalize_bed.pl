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
    
    Pod::Usage::pod2usage( -verbose => 1, -exitstatus => 0 ) if ($help);
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ($man);
}

__END__

=pod

=head1 NAME

normalize_bed.pl - Normalize bed files by read downsampling.

=head1 SYNOPSIS

normalize_bed.pl --input file1 [file2, file3, ..., filen] [OPTIONS]

Examples:

=over 4

=item perl normalize_bed.pl --input repli_1.bed repli_2.bed

=item perl normalize_bed.pl --input repli_1.bed repli_2.bed --rand 10000000 --savrem

=back

=head1 DESCRIPTION

A perl program to normalize reads across several samples from compensate for differences 
in sequencing depths by random downsampling. The downsampling can be performed based on
the sample with the lowest reads or down to a given number of reads (e.g. 10000000).Can 
run on Linux or Windows.

=head1 ARGUMENTS

=over 4

=item input B<(required)>

--input or -i

Input bed file(s). Please be careful as there is little checking whether the input
file(s) are indeed bed files. They should contain 3 or 6 columns. Input files can
be sorted or unsorted. Sorted are prefered to maintain uniformity.

=item rand B<(optional)>

--rand or -r

The method of random removal of sequences. Use "generate" to generate random numbers 
between 1..#tags and remove as many as to reach the minimum tag length, "permute" to 
randomly permute the tags and then remove as many as required to reach the minimum 
tag length or an integer. denoting how many tags should be left. If a sample has less
than this number, it is not touched. Defaults to "generate".

=item savrem B<(optional)>

--savrem or -t

Use this option if you wish to save the tags that are removed in a different file. 
The file will be named "filename_removed_tags.sav".

=item sort B<(optional)>

--sort or -s

Use this option if you wish to sort the input files first. The algorithmic results 
of the method will not change.

=item silent B<optional>

--silent or -s

Do not display verbose messages.

=item help B<(optional)>

--help or -h

Display this help text.

=item man B<(optional)>

--man or -m

Display the full manual of the script.

=back

=head1 AUTHOR

Panagiotis Moulos (L<moulos@fleming.gr)>

=cut
