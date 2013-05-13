=head1 NAME

HTS::Tools::Motifscan - A massive motif scanner (not finder!) for short genomic regions

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is a wrapper and significance threshold optimizer for several sequence motif scanners. At the
present point, the scanners supported are pwmscan from the GimmeMotifs suite (van Heeringen and Veenstra, 2010)
and MotifScanner (Thijs et al., 2001). Thus, the installation of these 3rd party tools in the system that
uses HTS::Tools::Motifscan is necessary (or at least of one of them). Then, the module runs a workflow which
which massively scan the input fasta files (e.g. a set of sequences corresponding to ChIP-Seq peaks) for
the motifs given in the motif file in the form of positional weight matrices (e.g. motifs that have been
characterized as significant during a de novo motif search in the same set of peaks). This workflow determines
an optimal significance threshold based on a set of random background sequences selected from the background
fasta file provided and on the false postive rate (fpr) provided. The module provides three basic outputs:
a simple file with the number of hits of each motif in each input file, gff files with the hits and bed
files with the hits, including also a significance score of the scan hit in the 5th column of the bed file.

	use HTS::Tools::Motifscan;
	my %params = (
		'input' => ['normal_nfkb_peaks.fa','cancer_nfkb_peaks.fa','mock_peaks.fa']
		'motif' => 'my_motif_matrices.pwm',
		'scanner' => 'pwmscan',
		'background' => 'background_sequences.tab',
		'range' => 5:0.1:15,
		'fpr' => 0.05,
		'times' => 10,
		'length' => 400,
		'output' => ['gff','bed','stats']
	)
    my $motifscanner = HTS::Tools::Motifscan->new(\%params);
    $motifscanner->run;

The acceptable parameters are as follows:

=over 4

=item I<input> B<(required)>

A set of input BED-like file(s) to be used as query regions. Each should containing a column with a UNIQUE 
region ID and a column with the region mode (e.g. peak summit, the point with the highest tag pile-up) or 
a location that the user thinks as the point of the region from which the distance to the genes will be 
calculated. If there is no such point, the center of the region may be used (see I<idmode>) parameter below

=item I<motif> B<(required)>

=item I<background> B<(optional)>

=item I<scanner> B<(optional)>

=item I<range> B<(optional)>

=item I<fpr> B<(optional)>

=item  I<times> B<(optional)>

=item I<length> B<(optional)>

=item I<output> B<(optional)>

=item  I<besthit> B<(optional)>

=item I<uniquestats> B<(optional)>

=item I<justscan> B<(optional)>

=item I<center> B<(optional)>

=item I<colext> B<(optional)>

=item  I<silent> B<(optional)>

=item I<silent> B<(optional)>

Use this parameter if you want to turn informative messages off.

=head1 OUTPUT

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Motifscan;

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use Cwd;
use File::Basename;
use File::Temp;
use File::Spec;

use lib '/media/HD4/Fleming/hts-tools/HTS-Tools/lib';
use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper);

our $MODNAME = "HTS::Tools::Count";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Short sequence read counting in genomic regions.";

BEGIN {
	$helper = HTS::Tools::Utils->new();
	select(STDOUT);
	$|=1;
	$SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

=cut

sub new
{
}

=head2 init

HTS::Tools::Motifscan object initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
}

=head2 get_random_seq

Get a random sequence from a set of background sequences.

	$motifscanner->get_random_seq($thefile,$itsindex,$howmany,$length,$num);

=cut

sub get_random_seq
{
	my ($self,$tabfasta,$itsindex,$itslen,$length,$num) = @_;
	my ($line,$id,$seq,$currlen,$start,$end);
	my $count = my $safeswitch = 0;
	my $BIG = 1e+6;
	srand;
	my $outfile = &createOutputFile(" ","sequence");
	open(RANDSEQ,">$outfile");
	open(FASTA,"<$tabfasta");
	open(INDEX,"$itsindex");
	while ($count < $num && $safeswitch < $BIG)
	{
		# Safeswitch in case too many sequences have small lengths, process shouldn't
		# take forever to complete...
		$safeswitch++;
		$line = getIndexedLine(*FASTA,*INDEX,int(rand($itslen)));
		($id,$seq) = split(/\t/,$line);
		$currlen = length($seq);
		next if ($currlen < $length);
		if ($currlen == $length)
		{
			($start,$end) = (1,$length);
			$id .= "_".$start."-".$end;
			&writeSeq(*RANDSEQ,$id,$seq);
			$count++;
		}
		if ($currlen > $length)
		{
			# Restrict the random index generation so that we don't go beyond the sequence end
			$start = int(rand($currlen - $length));
			$end = $start + $length;
			$id .= "_".$start."-".$end;
			&writeSeq(*RANDSEQ,$id,substr($seq,$start,$length));
			$count++;
		}
	}
	close(RANDSEQ);
	close(FASTA);
	close(INDEX);
	if ($safeswitch >= $BIG)
	{
		$helper->disp("Sequence fetching discontinued... $count sequences fetched in total in $outfile...");
		$helper->disp("Probably the FASTA file you supplied to get random sequences from contains too many short sequences.");
		$helper->disp("Try again with larger sequences or smaller length.");
	}
	return($outfile);
}

=head2 convert2bed

Converts the gff output from pwmscan to bed format. The first column of the gff (that is the peak/region ID) 
MUST contain coordinates information in the form chr:start-end (track2fasta) or chr:start:end.
WARNING! If the fasta files used for scanning have been generated with a program like track2fasta from
the GimmeMotifs suite, then the bed co-ordinates for each occurence can be correctly generated. If the 
sequence ids in the fasta files correspond to peak ids rather than exact sequence locations, another file 
with peak ids and peak centers must be provided. The function converts to 6 column bed files. It also
converts the motif score in gff file to the 0-1000 scale of UCSC genome browser so that motif strength can 
be visualized by color. This is done by linear conversion of the form new_value = a*old_value + b and 
by solving the special case of a 2x2 linear system (since we know several of the values):
min(score)*a + b = 0
max(score)*a + b = 1000

	$motifscanner->convert2bed($file_to_convert);

=cut

sub convert2bed
{
	my $f2c = shift @_;
	my %ch = @_;
	my ($base,$dir) = fileparse($f2c,'\..*?');
	my ($line,@lines,@scores,@content,@locs,@newcoord);
	my $bedout = $dir.$base.".bed";
	# In order to determine the coefficients of linear conversion we have to suck in all
	# gff file, the hard way...
	open(F2C,$f2c);
	while ($line = <F2C>)
	{
		$line =~ s/\r|\n$//g;
		push(@lines,$line);
		@content = split(/\t/,$line);
		push(@scores,$content[5]);
	}
	close(F2C);
	# If the scanner was MotifScanner/MotifLocator, there is one xtra line...
	my $tort = shift @scores if (!$scores[0]);
	# Get min, max score
	my ($min,$max) = &minmax(@scores);
	# Get coefficients
	my ($a,$b) = &naiveSolveTwo($min,$max);
	open(BED,">$bedout");
	if (%ch)
	{
		foreach $line (@lines)
		{
			@content = split(/\t/,$line);
			@locs = split(":",$content[0]);
			if ($#locs == 1) # track2fasta format
			{
				my $joined = pop(@locs);
				my @splitted = split("-",$joined);
				push(@locs,@splitted);
			}
			# Shift coordinates to motif occurence
			my $strand = "+";
			$strand = "-" if ($content[6] == -1 || $content[6] eq "R");
			@newcoord = ($locs[0],
						 $ch{$content[0]} - $cntcol[2] + $content[3],
						 $ch{$content[0]} - $cntcol[2] + $content[4],
						 $content[0],$a*$content[5] + $b,$strand);
			print BED join("\t",@newcoord),"\n";
		}
	}
	else
	{
		foreach $line (@lines)
		{
			@content = split(/\t/,$line);
			@locs = split(":",$content[0]);
			if ($#locs == 1) # track2fasta format
			{
				my $joined = pop(@locs);
				my @splitted = split("-",$joined);
				push(@locs,@splitted);
			}
			# Shift coordinates to motif occurence
			my $strand = "+";
			$strand = "-" if ($content[6] == -1 || $content[6] eq "R");
			@newcoord = ($locs[0],$locs[1] + $content[3],$locs[1] + $content[4],
						 $content[0],$a*$content[5] + $b,$strand);
			print BED join("\t",@newcoord),"\n";
		}
	}
	close(BED);
}

=head1 AUTHOR

Panagiotis Moulos, C<< <moulos at fleming.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc HTS::Tools::Motifscan

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

1; # End of HTS::Tools::Motifscan
