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

=head2 new

Constructor for HTS::Tools::Queries

=cut

sub get_query
{
	my ($self,$what,$utr) = @_;
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
