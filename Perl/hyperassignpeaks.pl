#!/usr/bin/perl

# All documentation moves to HTS::Tools::Assign. This now a wrapper for HTS::Tools::Assign. Use perldoc
# HTS::Tools::Assign to get instructions.
#
# Author : Panagiotis Moulos (pmoulos@fleming.gr)

use strict;
use Getopt::Long;

use HTS::Tools;
			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "hyperassignpeaks.pl";
our @peakfile; # The peak file(s) (BED format)
our $sigfile; # Significant regions to be assigned to peaks (BED format)
our $backfile; # Background regions to be assigned to peaks (BED format)
our @span; # Upstream and downstream  (default +/-100k from TSS)
our @sbcols; # Columns containing unique sig/back region ID and strand (default (4,5))
our @pcols; # Columns containing unique peak region ID and mode (default (4,5))
our @expcols; # Columns containing gene expression
our @out; # Output filetype
our $source; # Source for automatic templates
our $splicing; # Splicing for automatic templates
our $test; # Default hypergeometric test
our $pval; # Hypergeometric p-value cutoff (default 0.05)
our $log; # Keep log?
our $silent; # Display verbose messages
our $help; # Help?

# Check inputs
&check_inputs;

my $tool = HTS::Tools->new({
		"tool" => "assign",
		"log" => $log,
		"silent" => $silent,
		"params" => {
			"input" => \@peakfile,
			"region" => $sigfile,
			"background" => $backfile,
			"span" => \@span,
			"idstrand" => \@sbcols,
			"idmode" => \@pcols,
			"expression" => \@expcols,
			"outformat" => \@out,
			"test" => $test,
			"pvalue" => $pval,
			"source" => $source,
			"splicing" => $splicing
		}
});
$tool->run;

# Just parse parameters, checking is now performed by the module
sub check_inputs
{
    GetOptions("input|i=s{,}" => \@peakfile,
    		   "region|r=s" => \$sigfile,
    		   "background|b=s" => \$backfile,
    		   "span|n=i{,}" => \@span,
    		   "idstrand|c=i{,}" => \@sbcols,
    		   "idmode|m=i{,}" => \@pcols,
    		   "expression|e=i{,}" => \@expcols,
    		   "test|t=s" => \$test,
    		   "pvalue|p=f" => \$pval,
    		   "outformat|o=s{,}" => \@out,
    		   "source|u=s" => \$source,
    		   "splicing|x" => \$splicing,
    		   "log|l" => \$log,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    if (!@ARGV || $help)
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
A perl program to assign ChIP-Seq peaks to a set of regulated genes based on
the gene-peak distances between the genes of the regulated set and the gene-
peak distances between the genes of a background set (e.g. the whole genome).
In order to find peaks associatevely enriched in the set of the regulated
genes as compared to the background set, the program uses the hypergeometric
test on the population of peaks at a p-value threshold. There can be multiple
output format, all of them containing the gene-peak associations in different
formats. The hypergeometric method is NOT verified and there are other methods
out there that may perform better. Use at your own risk. The tools works very
nicely to calculate peak-gene distances with a very nice output.

Author : Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname --input peakfile(s) --region regfile --background backfile [OPTIONS]

--- Required ---
  --input|i  file(s)	Peak BED file(s) containing a column with a
  			UNIQUE peak ID and a column with the peak mode (the
  			point with the highest tag pile-up) or a location that
  			the user thinks as the point of the peak from which 
  			the distance to the genes will be calculated.
  --region|r  file	A BED file with the set of regulated genes,
  			containing a column with a UNIQUE gene (or whatever
  			region) ID and a column with the gene strand. Instead
  			of a local file, it can be one of the following, which 
  			will automatically download and use genomic annotations 
  			from the latest version of --source:
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
  --background|b		A BED file with the set of background 
  			genes, containing a column with a UNIQUE gene (or 
  			whatever region) ID and a column with the gene strand.
  			Required when running a statistical test.
  --span|d		Use this option to set the genomic span (distance
  			upstream and downstream from TSS) into which the program
  			will look for any peaks. It should be two values (e.g.
  			--span -50000 50000) and defaults to (-100000,100000).
  --idstrand|t		The columns in BOTH the gene files where their
			unique IDs and strands are. You should provide two values
			(e.g. --idstrand 4 5) where the first denotes the unique
			ID column and the second the strand column. It defaults
			to (4,5).
  --idmode|m		The columns in the peak files where their unique
			IDs and modes are. You should provide two values (e.g
			--idmode 4 5) where the first denotes the unique ID
			column and the second the mode column. It defaults to (4,5).
			Optionally, you can provide three values, where the 3rd 
			represents a peak score if available. This will be reported
			when using the "matrix-peak" output. The values must be 
			provided strictly with the following order: id column, 
			mode column, score column.
  --test|t		What over-representation statistical test to perform.
			Can be one of hypgeom for hypergeometric test, chi2 for
			chi-square test, auto for automatic selection and none for
			no testing. Defaults to hypgeom.
  --pvalue|p		The hypergeometric test p-value threshold. It
			defaults to 0.05.
  --outformat|o		Use this option to determine which output format
			filetype(s) you wish to retrieve.	Possible choices are:
   				"stats" for retrieving the significantly associated 
				peaks with their p-values, Bonferroni corrected p-values
				and enrichment ratios.
   				"gff-peak" for retrieving a peak-based gff file which
   				contains additional columns with peak ids, distances
   				and enrichment ratios. The score column is the p-value.
				"gff-gene" for similar to "gff-peak" but gene-based.
				"gff-peak-db" for same as "gff-peak" but with a header,
				suitable for incorporating to a database.
				"gff-gene-db" for same as "gff-gene" but with a header,
				suitable for incorporating to a database.
				"peak" for a simple file which contains the significantly
				associated peak IDs in the first column and a list of 
				associated genes in the second column.
				"gene" for similar to "peak" but gene based.
				"all-peak" for a simple file which contains ALL (based on
				distance) associated peak IDs in the first column and a 
				list of associated genes in the second column.
				"all-gene" for similar to "all-peak" but gene based.
				"pretty-peak" for retrieving a more human-readable format
				quite self-explicating (please see output).
				"pretty-gene" similar to "pretty-peak" but gene-based
				(please see output).
				"peakdata" for retrieving only the assigned peaks from
				the original peak file.
				"matrix-number" to retrieve a spreadsheet-like file where 
				rows correspond to genes (or the --region file) and columns
				correspond to peak files. The cell (i,j) contains the number
				of peaks in peak file j assigned to gene i.
				"matrix-presence" to retrieve a spreadsheet-like file where 
				rows correspond to genes (or the --region file) and columns
				correspond to peak files. The cell (i,j) contains "+" if peak
				in peak file j assigned to gene i, "-" otherwise.
				"matrix-peaks" to retrieve a spreadsheet-like file where 
				rows correspond to genes (or the --region file) and columns
				correspond to peak files. The cell (i,j) contains the peaks
				in peak file j assigned to gene i, "NP" otherwise.
			Example: --outformat stats gff-peak pretty-gene matrix
  --source|u		Use this option to set the online data source in
			the case of selecting one of the prefefined region templates
			with --region. Can be one of "ucsc", "refseq" or "ensembl".
			Defaults to "ensembl".
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


## Old chunks
#sub printMatrix
#{
#	my %inhash = @_;
#	my ($row,$colhash);
#	my @v;
#	my $outfilename = &createOutputFile($peakfile[0],"matrix");
#	disp("Writing output in $outfilename...");
#	my @rows = keys(%inhash);
#	print Dumper(\@rows);
#	open(OUTPUT,">$outfilename");
#	my $headhash = $inhash{$rows[0]};
#	my @headers = keys(%$headhash);
#	print OUTPUT "GeneID\t",join("\t",@headers),"\n";
#	foreach $row (@rows)
#	{
#		$colhash = $inhash{$row};
#		@v = values(%$colhash);
#		print OUTPUT "$row\t",join("\t",@v),"\n";
#	}
#	#close(OUTPUT);
#}
#
############################### PROTOTYPE SECTION ############################## 
#sub intersect (\@\@)
#{
#	my ($arr1,$arr2) = @_;
#	my %h1 = map{$_ => 1} @{$arr1};
#	my %h2 = map{$_ => 1} @{$arr2};
#	my @inter = grep($h1{$_},@{$arr2});
#	my $len1 = @{$arr1};
#	my $len2 = @{$arr2};
#	my $len3 = @inter;
#	print "\n$len1\t$len2\t$len3\n";
#	
#	return(@inter);
#}
############################## END PROTOTYPE SECTION ############################
