=head1 NAME

HTS::Tools::Fetch - The great new HTS::Tools::Fetch!

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

MySQL queries for HTS::Tools::Fetch

    use HTS::Tools::Queries;

    my $q = HTS::Tools::Queries->new();
    my $gene_query = $q->get_query("ucsc_gene_canonical");
    my $utr5_query = $q->get_query("ucsc_utr_canonical",5);
    print "\n$gene_query\n$utr5_query\n";

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Queries;

our $MODNAME = "HTS::Tools::Queries";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "MySQL queries for module HTS::Tools::Fetch.";

=head2 new

Constructor for HTS::Tools::Queries

=cut

sub new
{
	my $class = shift @_;
	my $self = {};
	bless($self,$class);
	return($self);
}

=head2 get_query

Returns the query for the respective elements. See the SYNOPSIS of HTS::Tools::Count for a list of
values for $what
	
	$q->get_query($what);

=cut

sub get_query
{
	my ($self,$what) = @_;
	$self->check($what);
	my $q;

	use v5.14;
	given($what)
	{
		when(/ucsc_canonical_genes/i)
		{
			$q = "SELECT knownCanonical.chrom AS `chrom`, `chromStart`, `chromEnd`, `transcript`, knownGene.exonCount AS `exonCount`, knownGene.strand AS `strand`, `geneName` ".
				"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name ".
				"INNER JOIN `knownToRefSeq` ON knownCanonical.transcript=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name GROUP BY `transcript` ORDER BY `chrom`, `chromStart`";
		}
		when(/ucsc_alternative_genes/i)
		{
			$q = "SELECT knownGene.chrom AS chrom, knownGene.txStart, knownGene.txEnd, knownGene.name, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownGene` INNER JOIN `knownToRefSeq` ON knownGene.name=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.name NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `txStart`";
		}
		when(/refseq_canonical_genes/i)
		{
			$q = "SELECT  refFlat.chrom, refFlat.txStart, refFlat.txEnd, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"GROUP BY `name` ORDER BY `chrom`, `txStart`";
		}
		when(/refseq_alternative_genes/i)
		{
			$q = "SELECT  refFlat.chrom, refFlat.txStart, refFlat.txEnd, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE knownToRefSeq.value NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `txStart`";
		}
		when(/ucsc_canonical_exons/i)
		{
			$q = "SELECT knownGene.chrom, knownGene.exonStarts, knownGene.exonEnds, `name`, `exonCount`, knownGene.strand, `transcript` ".
				"FROM `knownGene` INNER JOIN `knownCanonical` ON knownGene.name=knownCanonical.transcript GROUP BY `name` ORDER BY `chrom`, `exonStarts`";
		}
		when(/ucsc_alternative_exons/i)
		{
			$q = "SELECT knownGene.chrom, knownGene.exonStarts, knownGene.exonEnds, knownGene.name, knownGene.exonCount, knownGene.strand, `name` AS `transcript` ".
				"FROM `knownGene` WHERE `name` NOT IN (SELECT `transcript` FROM knownCanonical)".
				"GROUP BY `name` ORDER BY `chrom`, `exonStarts`";
		}
		when(/refseq_canonical_exons/i)
		{
			$q = "SELECT refFlat.chrom, refFlat.exonStarts, refFlat.exonEnds, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` AS `transcript` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"GROUP BY `name` ORDER BY `chrom`, `txStart`";
		}
		when(/refseq_alternative_exons/i)
		{
			$q = "SELECT  refFlat.chrom, refFlat.exonStarts, refFlat.exonEnds, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` AS `transcript` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE knownToRefSeq.value NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `txStart`";
		}
		when(/ucsc_canonical_5utr/i)
		{
			$q = "SELECT knownGene.chrom AS `chrom`, knownGene.txStart AS `start`, knownGene.cdsStart AS `end`, `transcript`, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name ".
				"INNER JOIN `knownToRefSeq` ON knownCanonical.transcript=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"+\"  AND knownGene.txStart <> knownGene.cdsStart UNION ".
				"SELECT knownGene.chrom AS `chrom`, knownGene.cdsEnd AS `start`, knownGene.txEnd AS `end`, `transcript`, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name ".
				"INNER JOIN `knownToRefSeq` ON knownCanonical.transcript=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"-\" AND knownGene.txEnd <> knownGene.cdsEnd ".
				"GROUP BY `transcript` ORDER BY `chrom`, `start`";
		}
		when(/ucsc_alternative_5utr/i)
		{
			$q = "SELECT knownGene.chrom AS `chrom`, knownGene.txStart AS `start`, knownGene.cdsStart AS `end`, knownGene.name, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownGene` INNER JOIN `knownToRefSeq` ON knownGene.name=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"+\" AND knownGene.txStart <> knownGene.cdsStart ".
				"AND knownGene.name NOT IN (SELECT `transcript` FROM knownCanonical) UNION ".
				"SELECT knownGene.chrom AS `chrom`, knownGene.cdsEnd AS `start`, knownGene.txEnd AS `end`, knownGene.name, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownGene` INNER JOIN `knownToRefSeq` ON knownGene.name=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"-\" AND knownGene.txEnd <> knownGene.cdsEnd ".
				"AND knownGene.name NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `start`";
		}
		when(/ucsc_canonical_3utr/i)
		{
			$q = "SELECT knownGene.chrom AS `chrom`, knownGene.cdsEnd AS `start`, knownGene.txEnd AS `end`, `transcript`, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name ".
				"INNER JOIN `knownToRefSeq` ON knownCanonical.transcript=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"+\" AND knownGene.cdsEnd <> knownGene.txEnd UNION ".
				"SELECT knownGene.chrom AS `chrom`, knownGene.txStart AS `start`, knownGene.cdsStart AS `end`, `transcript`, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name ".
				"INNER JOIN `knownToRefSeq` ON knownCanonical.transcript=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"-\" AND knownGene.txStart <> knownGene.cdsStart ".
				"GROUP BY `transcript` ORDER BY `chrom`, `start`";
		}
		when(/ucsc_alternative_3utr/i)
		{
			$q = "SELECT knownGene.chrom AS `chrom`, knownGene.cdsEnd AS `start`, knownGene.txEnd AS `end`, knownGene.name, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownGene` INNER JOIN `knownToRefSeq` ON knownGene.name=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"+\" AND knownGene.cdsEnd <> knownGene.txEnd ".
				"AND knownGene.name NOT IN (SELECT `transcript` FROM knownCanonical) UNION ".
				"SELECT knownGene.chrom AS `chrom`, knownGene.txStart AS `start`, knownGene.cdsStart AS `end`, knownGene.name, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownGene` INNER JOIN `knownToRefSeq` ON knownGene.name=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.strand = \"-\"  knownGene.txStart <> knownGene.cdsStart ".
				"AND knownGene.name NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `start`";
		}
		when(/refseq_canonical_5utr/i)
		{
			$q = "SELECT refFlat.chrom AS `chrom`, refFlat.txStart AS `start`, refFlat.cdsStart AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"WHERE refFlat.strand = \"+\" AND refFlat.txStart <> refFlat.cdsStart UNION ".
				"SELECT refFlat.chrom AS `chrom`, refFlat.cdsEnd AS `start`, refFlat.txEnd AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"WHERE refFlat.strand = \"-\" AND refFlat.txEnd <> refFlat.cdsEnd ".
				"GROUP BY `name` ORDER BY `chrom`, `start`";
		}
		when(/refseq_alternative_5utr/i)
		{
			$q = "SELECT refFlat.chrom AS `chrom`, refFlat.txStart AS `start`, refFlat.cdsStart AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE refFlat.strand = \"+\" AND refFlat.txStart <> refFlat.cdsStart ".
				"AND knownToRefSeq.value NOT IN (SELECT `transcript` FROM knownCanonical) UNION ".
				"SELECT refFlat.chrom AS `chrom`, refFlat.cdsEnd AS `start`, refFlat.txEnd AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE refFlat.strand = \"-\" AND refFlat.txEnd <> refFlat.cdsEnd ".
				"AND knownToRefSeq.value NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `start`";
		}
		when(/refseq_canonical_3utr/i)
		{
			$q = "SELECT refFlat.chrom AS `chrom`, refFlat.cdsEnd AS `start`, refFlat.txEnd AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"WHERE refFlat.strand = \"+\" AND refFlat.cdsEnd <> refFlat.txEnd UNION ".
				"SELECT refFlat.chrom AS `chrom`, refFlat.txStart AS `start`, refFlat.cdsStart AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"WHERE refFlat.strand = \"-\" AND refFlat.txStart <> refFlat.cdsStart ".
				"GROUP BY `name` ORDER BY `chrom`, `start`";
		}
		when(/refseq_alternative_3utr/i)
		{
			$q = "SELECT refFlat.chrom AS `chrom`, refFlat.cdsEnd AS `start`, refFlat.txEnd AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE refFlat.strand = \"+\" AND refFlat.cdsEnd <> refFlat.txEnd ".
				"AND knownToRefSeq.value NOT IN (SELECT `transcript` FROM knownCanonical) UNION ".
				"SELECT refFlat.chrom AS `chrom`, refFlat.txStart AS `start`, refFlat.cdsStart AS `end`, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE refFlat.strand = \"-\" AND refFlat.txStart <> refFlat.cdsStart ".
				"AND knownToRefSeq.value NOT IN (SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `start`";
		}
		when(/ucsc_canonical_cds/i)
		{
			$q = "SELECT knownCanonical.chrom AS chrom, knownGene.cdsStart, knownGene.cdsEnd, `transcript`, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name ".
				"INNER JOIN `knownToRefSeq` ON knownCanonical.transcript=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.cdsStart <> knownGene.cdsEnd ".
				"GROUP BY `transcript` ORDER BY `chrom`, `cdsStart`";
		}
		when(/ucsc_alternative_cds/i)
		{
			$q = "SELECT knownGene.chrom AS chrom, knownGene.cdsStart, knownGene.cdsEnd, knownGene.name, knownGene.exonCount, knownGene.strand, `geneName` ".
				"FROM `knownGene` INNER JOIN `knownToRefSeq` ON knownGene.name=knownToRefSeq.name ".
				"INNER JOIN `refFlat` ON knownToRefSeq.value=refFlat.name ".
				"WHERE knownGene.cdsStart <> knownGene.cdsEnd AND knownGene.name NOT IN ".
				"(SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `cdsStart`";
		}
		when(/refseq_canonical_cds/i)
		{
			$q = "SELECT  refFlat.chrom, refFlat.cdsStart, refFlat.cdsEnd, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"INNER JOIN knownCanonical ON knownToRefSeq.name=knownCanonical.transcript ".
				"WHERE refFlat.cdsStart <> refFlat.cdsEnd ".
				"GROUP BY `name` ORDER BY `chrom`, `cdsStart`";
		}
		when(/refseq_alternative_cds/i)
		{
			$q = "SELECT  refFlat.chrom, refFlat.txStart, refFlat.txEnd, refFlat.name, refFlat.exonCount, refFlat.strand, `geneName` ".
				"FROM `refFlat` INNER JOIN knownToRefSeq ON refFlat.name=knownToRefSeq.value ".
				"WHERE refFlat.cdsStart <> refFlat.cdsEnd AND knownToRefSeq.value NOT IN ".
				"(SELECT `transcript` FROM knownCanonical) ".
				"GROUP BY `name` ORDER BY `chrom`, `cdsStart`";
		}
	}
	
	return($q);
}

sub check
{
	my ($self,$qc) = @_;
	my $checker = HTS::Tools::Paramcheck->new({"tool" => "queries","params" => {"query" => $qc}});
	$checker->validate;
}

=head1 AUTHOR

Panagiotis Moulos, C<< <moulos at fleming.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc HTS::Tools::Utils


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

1; # End of HTS::Tools::Queries
