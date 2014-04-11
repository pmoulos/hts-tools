=head1 NAME

HTS::Tools::Track::Signal - A massive motif scanner (not finder!) for short genomic regions

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is a wrapper and significance threshold optimizer for several sequence motif scanners. At the
present point, the scanners supported are pwmscan from the GimmeMotifs suite (van Heeringen and Veenstra, 2010)
and MotifScanner (Thijs et al., 2001). Thus, the installation of these 3rd party tools in the system that
uses HTS::Tools::Track::Signal is necessary (or at least of one of them). Then, the module runs a workflow which
which massively scan the input fasta files (e.g. a set of sequences corresponding to ChIP-Seq peaks) for
the motifs given in the motif file in the form of positional weight matrices (e.g. motifs that have been
characterized as significant during a de novo motif search in the same set of peaks). This workflow determines
an optimal significance threshold based on a set of random background sequences selected from the background
fasta file provided and on the false postive rate (fpr) provided. The module provides three basic outputs:
a simple file with the number of hits of each motif in each input file, gff files with the hits and bed
files with the hits, including also a significance score of the scan hit in the 5th column of the bed file.

    use HTS::Tools::Track::Signal;
    my %params = (
        'input' => ['normal_nfkb_peaks.fa','cancer_nfkb_peaks.fa','mock_peaks.fa']
        'motif' => 'my_motif_matrices.pwm',
        'scanner' => 'pwmscan',
        'background' => 'background_sequences.tab',
        'range' => 0.1:0.1:1,
        'fpr' => 0.05,
        'times' => 10,
        'length' => 400,
        'output' => ['gff','bed','stats']
    )
    my $motifscanner = HTS::Tools::Track::Signal->new(\%params);
    $motifscanner->run;

The acceptable parameters are as follows:

=over 4

=item I<input> B<(required)>

A set of input FASTA file(s) that will be used for matching of the motifs contained in I<motif> file. Each 
of these files will be scanned using the scanner defined by I<scanner> and the hits will be reported either
as a simple tab-delimited stats file, a gff file, or a bed file under conditions (see below).

=item I<motif> B<(required)>

A file containing several motifs in the Position Weight Matrix (PWM) format. The file will be decomposed
in several motifs and each motif will be matched against each input FASTA file.

=item I<background> B<(optional)>

A FASTA file containing background sequences (must be at least the same length as the input sequences)
that will be used for random sampling to define a score cutoff for motif matching significance based on
I<fpr>. It is required if the option I<justscan> is not activated.

=item I<scanner> B<(optional)>

A scanning algorithm to use. Currently, two algorithms are supported: the pwmscan algorithm (van Heeringen
and Veenstra, 2010) and MotifScanner (Thijs et al., 2001).

=item I<range> B<(optional)>

A range of cutoffs that will be used to determine the final matching score cutoff corresponding to I<fpr>.
The range can be given in the form a:b or a:x:b where x is an increment step. If the format a:b is chosen,
the default increment is 1. However, in the latest versions of pwmscan, the matching score is normalized
to one, so this notation will be eventually deprecated.

=item I<fpr> B<(optional)>

The desired False Positive Rate (FPR), that is the percentage of motif matches found in the background
sequences. Defaults to 0.05.

=item  I<times> B<(optional)>

How many times should the background set of sequences that will be used for the determination of the cutoff
score from I<range>, be larger than the input set? It default to 10, which means that if an input FASTA file
contains 100 sequences, the background will contain 1000 random sequences.

=item I<length> B<(optional)>

The length of the background sequences. Generally, it should be larger than the length of input sequences.
It defaults to 400.

=item I<output> B<(optional)>

The output types that the user wished to get. It can be one or more of "stats" for a simple delimited
file containing the hits for each motif and input file, "gff" for the gff output of pwmscan or a related
file from MotifScanner, containing the actual hits and positions in the input sequences, or "bed" for an
output BED file containing the motif matches locations and a score to be used for coloring. The "bed"
output is available only if a set of peak files is given so as to determine the relative location of the
match inside the peak and construct proper bed lines.

=item  I<besthit> B<(optional)>

The number of best hits to be retrieved when I<scanner> is "pwmscan". Defaults to 1.

=item I<uniquestats> B<(optional)>

If the number of besthits is greater than 1, the specifying I<uniquestats> to 1 will cause the output
"stats" file to contain unique hits. Like this you can avoid paradoxes like having more hits than input
FASTA sequences. However, in some situations you might actually want to retrieve multiple hits.

=item I<justscan> B<(optional)>

Set this to 1, to just scan the input sequences for the motifs without defining an FPR. Useful if you
have a very limited number of input sequences (e.g. just one).

=item I<center> B<(optional)>

A set of genomic regions (e.g. peaks) with the SAME IDs as the input FASTA sequences (or a superset of
these) which will be used to assign peak regions to FASTA files in order to determine the proper coordinates 
for the generation of BED output. It is optional, however, if BED output is requested, these files are 
not given, and the FASTA sequence length does not correspond to the length that can be extracted by the 
FASTA ID, the coordinates will be inaccurate.

=item I<colext> B<(optional)>

A vector of length 3, containing the column numbers of peak ID and peak summit, and the length of the 
(possible) extension upstream and downstream of the peak summit. For example I<colext> 4 5 75.

=item I<silent> B<(optional)>

Use this parameter if you want to turn informative messages off.

=head1 OUTPUT

The output of the module is a set of GFF files containing the hits for each motif and each input set of
sequences, a statistics file with the number of hits for each motif and input set of sequences and a set
of BED files with the motif matches that can be used for display in a genome browser.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Track::Signal;

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use Cwd;
use File::Basename;
use File::Temp;
use File::Spec;

use HTS::Tools::Constants;
use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper $const);

our $MODNAME = "HTS::Tools::Track::Signal";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Create an NGS signal track using 3rd party tools.";

BEGIN {
    $helper = HTS::Tools::Utils->new();
    $const = HTS::Tools::Constants->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

The HTS::Tools::Track::Signal object constructor. It accepts a set of parameters that are required to run the
motifscanner and get the output.

    my $motifscanner = HTS::Tools::Track::Signal->new({'input' => ['peaks_1.fa','peaks_2.fa'],'motif' => 'my_motifs.pwm',
        'scanner' => 'pwmscan','justscan' => 1});

=cut

sub new
{
    my ($class,$params) = @_;
    my $self = {};

    # Pass global variables to the helper
    (defined($params->{"silent"})) ? ($helper->set("silent",$params->{"silent"})) :
        ($helper->set("silent",0));
    (defined($params->{"tmpdir"})) ? ($helper->set("tmpdir",$params->{"tmpdir"})) :
        ($helper->set("tmpdir",File::Temp->newdir()));
    $helper->set_logger($params->{"log"}) if (defined($params->{"log"}));
    $helper->advertise($MODNAME,$VERSION,$AUTHOR,$EMAIL,$DESC);

    # Validate the input parameters
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "track", "params" => $params});
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init

HTS::Tools::Track::Signal object initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
    my ($self,$params) = @_;

    # Basic
    foreach my $p (keys(%$params)) { $self->set($p,$params->{$p}); }

    # Global
    $self->set("silent",$helper->get("silent")) unless defined($self->{"silent"});
    $self->set("tmpdir",$helper->get("tmpdir")) unless defined($self->{"tmpdir"});
    
    return($self);
}

=head2 run

The HTS::Tools::Track::Signal run subroutine. It runs the motifscanner with the given parameters in the 
constructor.

    $signaler->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    # Copy some memory-less variables to avoid rewriting the whole thing...
}

=head2 fun2

Parse a file of PWMs and construct different files, specific for the scanner in use. Internal use.

    $signaler->;

=cut

sub fun2 {}

=head2 get

HTS::Tools::Track::Signal object getter.

    my $param_value = $motifscanner->get('param_name')

=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Track::Signal object setter.

    $motifscanner->set('param_name','param_value');
    
=cut

sub set
{
    my ($self,$name,$value) = @_;
    $self->{$name} = $value;
    return($self);
}

=head1 AUTHOR

Panagiotis Moulos, C<< <moulos at fleming.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc HTS::Tools::Track::Signal

You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=HTS-Tools>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/HTS-Tools>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/HTS-Tools>

=item * Search CPAN

L<http://search.cpan.org/dist/HTS-Tools/>

=back

=head1 ACKNOWLEDGEMENTS

=head1 LICENSE AND COPYRIGHT

Copyright 2013 Panagiotis Moulos.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

1; # End of HTS::Tools::Track::Signal
