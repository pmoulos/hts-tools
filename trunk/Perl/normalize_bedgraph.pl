#!/usr/bin/perl

# All documentation moves to HTS::Tools::Normalize and HTS::Tools::Normalize::Bedgraph.
# This now a wrapper for HTS::Tools::Normalize and HTS::Tools::Normalize::Bedgraph.
# Use perldoc HTS::Tools::Normalize or HTS::Tools::Normalize::Bedgraph to get instructions.
#
# Author : Panagiotis Moulos (pmoulos@fleming.gr)

use strict;
use Getopt::Long;
use Pod::Usage;

use HTS::Tools;
             
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "normalize_bedgraph.pl";
our @bglist;
our @extnorm;
our $sumto;
our $exportfacs;
our $perlonly;
our @output;
our $prerun;
our $prerunlog;
our $man;
our $log;
our $silent;
our $help;

&check_inputs;

my $tool = HTS::Tools->new({
        "tool" => "normalize",
        "log" => $log,
        "silent" => $silent,
        "params" => {
            "input" => \@bglist,
            "type" => "bedgraph",
            "output" => \@output,
            "extnorm" => \@extnorm,
            "exportfactors" => $exportfacs,
            "perlonly" => $perlonly,
            "prerun" => $prerun,
            "prerunlog" => $prerunlog
        }
});
$tool->run;

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions(
        "input|i=s{,}" => \@bglist,
        "output|o=s{,}" => \@output,
        "extnorm|r=f{,}" => \@extnorm,
        "sumto|s=i" => \$sumto,
        "exportfactors|f=s" => \$exportfacs,
        "perlonly|p" => \$perlonly,
        "prerun|r" => \$prerun,
        "prerunlog|l=s" => \$prerunlog,
        "man|m" => \$man,
        "log|g" => \$log,
        "silent|s" => \$silent,
        "help|h" => \$help
    );

    Pod::Usage::pod2usage( -verbose => 1, -exitstatus => 0 ) if ($help);
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ($man);
}

__END__

=pod

=head1 NAME

normalize_bedgraph.pl - Normalize bedgraph files to a constant total signal

=head1 SYNOPSIS

normalize_bedgraph.pl --input file1 [file2, file3, ..., filen] [OPTIONS]

Examples:

=over 4

=item perl normalize_bedgraph.pl --input repli_1.bg repli_2.bg --output repli_1_norm.bg repli_2_norm.bg --sumto 5000000000 --exportfactors

=item perl normalize_bedgraph.pl --input repli_1.bg repli_2.bg --perlonly --extnorm 0.92 1.28

=back

=head1 DESCRIPTION

Normalize bedgraph signal using to a total signal value (similar to RSeqC
normalize_bigwig.py but faster!) or a set of external normalization factors
(e.g. calculated from DESeq or edgeR).

=head1 ARGUMENTS

=over 4

=item input B<(required)>

--input or -i

Input bedgraph file(s). Please be careful as there is checking whether the input
file(s) are indeed bedgraph files. It's ok if they contain more than 4 columns
but the first four must be bedgraph (chromosome, start, end, signal separated by
tabs). Input files need not to be sorted.

=item output B<(optional)>

--output or -o

Output file names. It can be "stdout" for exporting to STDOUT, a set of file
names equal to the number of input files or nothing for the files to be
auto-generated.

=item sumto B<(optional)>

--sumto or -s

Normalize to --sumto total wig signal. Defaults to 1000000000. It is mutually
exclusive with --extnorm with --extnorm in precedence.

=item extnorm B<(optional)>

--extnorm or -e

A set of external normalization factors (e.g. calculated from DESeq or edgeR).
It is mutually exclusive with --sumto with --extnorm in precedence.

=item exportfactors B<(optional)>

--exportfactors or -f

Export normalization factors and signal sums to a file specified by
--exportfactors.

=item perlonly B<(optional)>

--perlonly or -p

Use pure Perl to run the script, otherwise, uses Linux awk. Useful for e.g.
Windows systems but slower.

=item prerun B<(optional)>

--prerun or -r

If this switch is turned on, the script just counts the total wiggle signal in
the input files and prints it on the screen. This is useful in order for example
to determine the total normalization signal (--sumto).

=item prerunlog B<(optional)>

--prerunlog or -l

Writes the output of --prerun to a file specified by --prerunlog. If only the
--prerunlog is specified, --prerun is executed automatically.

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
