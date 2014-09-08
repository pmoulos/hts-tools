=head1 NAME

HTS::Tools::Assign - Genomic region distance-based assignment, with optional statistical significance.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module assigns one or more sets of genomic regions (queries) to another set of genomic regions (subject) 
using the distance from a fixed point in the queries to another fixed point in the subject regions. One example
of the above is the assignment of a set of ChIP-Seq peaks to a set of regulated genes based on the gene-peak 
distances between the genes of the regulated set. The peak mode (or summit) can be used as the queries fixed
point, while the gene TSS can be used as the fixed point of the subject geneomic regions. The module may
optionally assign a p-value to each query, using a total genomic background subject (e.g. the whole genome
annotation instead of a set of differentially expressed genes) based on an experimental procedure described
below. This procedure has been used in Rao et al., 2011 (http://www.ncbi.nlm.nih.gov/pubmed/18981474) and is
highly experimental. It is optional and you may use it at own risk.
In order to find peaks associatevely enriched in a set of the regulated genes as compared to a background set, 
the module uses the hypergeometric test on the population of peaks at a p-value threshold. This p-value expresses
the probability of a peak being significantly assigned to a set of close regulated genes as compared to the 
total number of genes present in the area specified by the span parameter. Please see the article above for
further details.
There can be multiple outputs, all of them containing the gene-peak associations in different formats which
are described in the parameters below. The hypergeometric method is NOT verified and there are other methods
out there that may perform better. Again, use at your own risk. The tools works very nicely to calculate 
peak-gene distances with a set of very nice and informative outputs.

    use HTS::Tools::Assign;
    my %params = (
        'input' => ['normal_nfkb_peaks.txt','cancer_nfkb_peaks.txt']
        'region' => 'my_experiment_de_genes.txt',
        'span' => [-10000,1000],
        'idstrand' => [4 6],
        'idmode' => [4 5],
        'outformat' => ['pretty-peak','gff-gene','gene-peak-presence','gene-peak-number']
    )
    my $assigner = HTS::Tools::Assign->new(\%params);
    $assigner->run;

The acceptable parameters are as follows:

=over 4

=item I<input> B<(required)>

A set of input BED-like file(s) to be used as query regions. Each should containing a column with a UNIQUE 
region ID and a column with the region mode (e.g. peak summit, the point with the highest tag pile-up) or 
a location that the user thinks as the point of the region from which the distance to the genes will be 
calculated. If there is no such point, the center of the region may be used (see I<idmode>) parameter below

=item I<region> B<(required)>

A BED-like file with the set of subject regions (e.g. a set of regulated genes), containing a column with a 
UNIQUE gene (or whatever region) ID and a column with the gene strand. Alternatively, the I<region> parameter
can be one of the following keywords for a set of predefined region templates, which will be automatically 
downloaded from the source defined by the I<source> parameter:

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

=item I<background> B<(optional)>

A BED file with the set of background regions (e.g. a set of regulated genes), containing a column with a 
UNIQUE gene (or whatever region) ID and a column with the gene strand. Alternatively, the I<region> parameter
can be one of the following keywords for a set of predefined region templates (this file is used only for
the statistical test described above). For these templates, see the I<region> parameter above.

=item I<source> B<(optional)>

Use this option to set the online data source in the case of selecting one of the predefined region 
templates with I<region>. It can be one of "ucsc", "refseq" or "ensembl" and Default to "ensembl".

=item I<splicing> B<(optional)>

Use this option with I<source> to determine whether the canonical or alternatively spliced transcripts will
be used for counting. It can be "canonical" or "alternative" and defaults to "canonical"

=item  I<span> B<(optional)>

Use this parameter to set the genomic span (distance upstream and downstream from the fixed point of the subject
regions) into which the module will look for any queries. It should consist of an array of two values (e.g. 
[-10000,1000]) and defaults to [-10000,10000].

=item I<idstrand> B<(optional)>

The columns in both the subject and the background (if used) file(s) where their unique IDs and strands are. 
You should provide an array of two values (e.g. [4,5]) where the first denotes the unique ID column and the 
second the strand column. It defaults to [4,5].

=item I<expression> B<(optional)>

The columns in the subject file where there are optional expression values, for example if the subject file is
a set of expressed genes. You should provide an array of values (e.g. [7,8,9]) that denote the column number where
the expression values are stored.

=item I<idmode> B<(optional)>

The columns in the query files where their unique IDs and possibly modes are. You should provide at least two values 
(e.g I<idmode> [4,5]) where the first denotes the unique ID column and the second the mode column. If the second
value is not provided, the module assumes that the center of each query region is its mode. It defaults to [4].
Optionally, you can provide three values, where the 3rd represents a peak score if available. This will be reported
when using the "matrix-peak" output. The values must be provided strictly with the following order: id column, mode
column, score column.

=item I<test> B<(optional)>

What over-representation statistical test to perform. Can be one of hypgeom for hypergeometric test, chi2 for
chi-square test, auto for automatic selection and none for no testing. Defaults to none.

=item I<pvalue> B<(optional)>

The statistical test p-value threshold. It defaults to 0.05.

=item I<outformat>  B<(optional)>

Use this parameter to determine which output format filetype(s) you wish to retrieve. Possible choices are:
"stats" for retrieving the significantly associated peaks with their p-values, Bonferroni corrected p-values
and enrichment ratios.
"gff-peak" for retrieving a peak-based gff file which contains additional columns with peak ids, distances
and enrichment ratios. The score column is the p-value.
"gff-gene" for similar to "gff-peak" but gene-based.
"gff-peak-db" for same as "gff-peak" but with a header, suitable for incorporating to a database.
"gff-gene-db" for same as "gff-gene" but with a header, suitable for incorporating to a database.
"peak" for a simple file which contains the significantly associated peak IDs in the first column and a 
list of associated genes in the second column.
"gene" for similar to "peak" but gene based.
"all-peak" for a simple file which contains ALL (based on distance) associated peak IDs in the first 
column and a list of associated genes in the second column.
"all-gene" for similar to "all-peak" but gene based.
"pretty-peak" for retrieving a more human-readable format quite self-explicating (please see output).
"pretty-gene" similar to "pretty-peak" but gene-based (please see output).
"peakdata" for retrieving only the assigned peaks from the original peak file.
"matrix-number" to retrieve a spreadsheet-like file where rows correspond to the subject region file and columns
correspond to query files. The cell (i,j) contains the number of regions in query file j assigned to subject
region i.
"matrix-presence" to retrieve a spreadsheet-like file where rows correspond to the subject region file and columns
correspond to peak files. The cell (i,j) contains "+" if region in query file j assigned to subjecy region i, or"-" 
otherwise.
"matrix-peaks" to retrieve a spreadsheet-like file where rows correspond to subject region file and columns
correspond to query region files. The cell (i,j) contains the regions in query file j assigned to subject region i, 
or "NP" otherwise.

=item I<silent> B<(optional)>

Use this parameter if you want to turn informative messages off.

=head1 OUTPUT

The main output of the program is up to twelve files with information on gene-peak (or more generally, 
region-region) association.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Assign;

use v5.10;
use strict;
use warnings FATAL => 'all';

our $MODNAME = "HTS::Tools::Assign";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Assign query genomic regions to subject genomic regions, possibly based on over-representation.";

use Carp;
use File::Basename;
use File::Temp;
use File::Spec;

use HTS::Tools::Count;
use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper);

BEGIN {
    $helper = HTS::Tools::Utils->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

The HTS::Tools::Assign object constructor. It accepts a set of parameters that are required to run the
asigner and get the output.

    my $assigner = HTS::Tools::Assign->new({'input' => 'my_peaks.txt','region' => 'my_genome.txt'});

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
    my $checker = HTS::Tools::Paramcheck->new();
    $checker->set("tool","assign");
    $checker->set("params",$params);
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init

HTS::Tools::Assign object initialization method. NEVER use this directly, use new instead.

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

The HTS::Tools::Assign run subroutine. It runs the assigner with the given parameters in the constructor.

    $assigner->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    # Copy some memory-less variables to avoid rewriting the whole thing...
    my @input = @{$self->get("input")};
    my $region = $self->get("region");
    my $background = $self->get("background");
    my $test = $self->get("test");
    my $pval = $self->get("pvalue");
    my $source = $self->get("source");
    my $splicing = $self->get("splicing");
    my @pcols = @{$self->get("idmode")};
    my @sbcols = @{$self->get("idstrand")};
    my @expcols = @{$self->get("expression")};
    my @span = @{$self->get("span")};
    my @out = @{$self->get("outformat")};
    
    # Bavard
    my $date = $helper->now;
    $helper->disp("$date - Started...\n");
    $helper->disp("Genomic span from TSS : ($span[0],$span[1])");
    $helper->disp("Over-representation test : $test");
    $helper->disp("p-value threshold : $pval") if ($test ne "none");
    $helper->disp("Chosen output(s) : ",join("\, ",@out),"\n");

    # General indices and some info about gff output
    my ($i,$j,$k);
    my $gffreq = my $pdataout = 0;
    foreach my $o (@out)
    {
        $gffreq = 1 if ($o =~ /gff/);
        $pdataout = 1 if ($o =~ /peakdata|bed/);
    }
    
    # Some intialization
    my (%sigID,%sigStart,%sigEnd);
    my (%backID,%backStart,%backEnd);
    my (@all,$chr,$start,$end,$id,$strand,$expr,$line);
    my $lensig = my $lenback = 0;
    my @lenpeak;

    # Variable for matrix generation
    my (%hasPeak,%hasExpression);
    tie %hasPeak, "Tie::IxHash::Easy" if ($helper->smatch("matrix",@out));
    
    # Initiate a counter in case we have to fetch files
    my $counter = HTS::Tools::Count->new({"tmpdir" => $self->get("tmpdir"),"silent" => 1, "input" => "foo", "region" => $region});

    # Suck in significant region file
    if (! -f $region)
    {
        $region = $counter->fetch_regions($region,$self->get("source"),$self->get("splicing"));
        @sbcols = (4,6);
    }
    open (REG,$region) or croak "\nThe file $region does not exist!\n";;
    $helper->disp("Reading region file $region...");
    $line = <REG>;
    my $reghead = $helper->decide_header($line);
    if (!$reghead)
    {
        seek(REG,0,0);
    }
    else
    {
        if (@expcols)
        {
            my @ht = split("\t",$reghead);
            $hasExpression{"header"} = join("\t",@ht[@expcols]);
        }
    }
    while ($line = <REG>)
    {
        next if ($line =~/^chrM/);
        next if ($line =~/rand/);
        $line =~ s/\r|\n$//g;
        @all = split(/\t/,$line);
        $chr = $all[0];
        $start = $all[1];
        $end = $all[2];
        $id = $all[$sbcols[0]];
        $strand = $all[$sbcols[1]];
        $expr = join("\t",@all[@expcols]) if (@expcols);
        #if ($strand == 1 || $strand eq "+" || $strand eq "F")
        if ($strand =~ m/^[+1F]$/)
        {
            push(@{$sigStart{$chr}},$start);
            push(@{$sigEnd{$chr}},$end);
        }
        #elsif ($strand == -1 || $strand eq "-" || $strand eq "R")
        elsif ($strand =~ m/^(-1)|[-R]$/)
        {
            push(@{$sigStart{$chr}},$end);
            push(@{$sigEnd{$chr}},$start);
        }
        else # Some self-defense
        {
            $helper->disp("Improper strand format... Skipping line $. from $region");
            next;
        }
        push(@{$sigID{$chr}},$id);
        $lensig++;

        # Initiate the hash to keep record of peaks for a peak matrix file generation
        if ($helper->smatch("matrix",@out))
        {
            for ($i=0; $i<@input; $i++)
            {
                $hasPeak{$id}{basename($input[$i])} = ();
            }
            $hasExpression{"data"}{$id} = $expr;
        }
    }
    close(REG);
    
    # Suck in background region file if test is to be performed
    if ($test ne "none")
    {
        if (! -f $background)
        {
            $background = $counter->fetch_regions($background,$self->get("source"),$self->get("splicing"));
            @sbcols = (4,6);
        }
        open (BACK,$background) or croak "\nThe file $background does not exist!\n";
        $helper->disp("Reading background file $background...");
        $line = <BACK>;
        my $backhead = $helper->decide_header($line);
        seek(BACK,0,0) if (!$backhead);
        while ($line = <BACK>)
        {
            next if ($line =~/^chrM/);
            next if ($line =~/rand/);
            $line =~ s/\r|\n$//g;
            @all = split(/\t/,$line);
            $chr = $all[0];
            $start = $all[1];
            $end = $all[2];
            $id = $all[$sbcols[0]];
            $strand = $all[$sbcols[1]];
            if ($strand == 1 || $strand eq "+" || $strand eq "F")
            {
                push(@{$backStart{$chr}},$start);
                push(@{$backEnd{$chr}},$end);
            }
            elsif ($strand == -1 || $strand eq "-" || $strand eq "R")
            {
                push(@{$backStart{$chr}},$end);
                push(@{$backEnd{$chr}},$start);
            }
            else # Some self-defense
            {
                disp("Improper strand format... Skipping line $. from $background");
                next;
            }
            push(@{$backID{$chr}},$id);
            $lenback++;
        }
        close(BACK);
    }
    
    # Some self-defense... Check uniqueness of gene IDs...
    croak "\nThe region IDs in regions file $region are not unique! Exiting...\n" if (!$self->check_unique(%sigID));
    croak "\nThe region IDs in background file $background are not unique! Exiting...\n" if (!$self->check_unique(%backID) && $test ne "none");

    # Read and process peak files
    for ($i=0; $i<@input; $i++)
    {
        $helper->disp("\nReading file $input[$i]...");
        
        my (%peakID,%peakMode,%peakScore,%countSig,%countBack,%allPeakData,%countGenesSig,%countGenesBack);
        my ($pstart,$pend,%peakStarts,%peakEnds) if ($gffreq);
        my (@pall,$pchr,$pmode,$pid,$phead,$pscore);
        my (@starts,@ends,@modes);
        my (@peakchrs,@peakids,@peakscores,@geneids);
        my (@sigallpeakids,@backallpeakids,@sigassgenes,@backassgenes);
        my ($currchr,$currdist,$currpeak,$currgene,$currout,@currgenes);
        my (%genePeaksSig,%peaksGenesSig,%distPeakBased,%distGeneBased,%peakIndex,%geneIndex);
        my ($p,$cp,@elems);
        my %finalPeaks;
        
        open(INPUT,$input[$i]) or croak "\nThe file $input[$i] does not exist!\n";
        $line = <INPUT>;
        $phead = $helper->decide_header($line);
        seek(INPUT,0,0) if (!$phead);
        while ($line = <INPUT>)
        {
            $line =~ s/\r|\n$//g;
            @pall = split(/\t/,$line);
            $pchr = $pall[0];
            $pstart = $pall[1];
            $pend = $pall[2];
            $pid = $pall[$pcols[0]];
            ($pcols[1]) ? ($pmode = $pall[$pcols[1]]) :
            ($pmode = $pstart + $helper->round(($pend - $pstart)/2));
            $pscore = $pall[$pcols[2]] if ($pcols[2]);
            push(@{$peakID{$pchr}},$pid);
            push(@{$peakMode{$pchr}},$pmode);
            push(@{$peakScore{$pchr}},$pscore) if ($pscore);
            if ($gffreq) # Starts and ends required in this case for GFF files
            {
                push(@{$peakStarts{$pchr}},$pstart);
                push(@{$peakEnds{$pchr}},$pend);
            }
            $allPeakData{$pid} = join("\t",@pall) if ($pdataout);
            $lenpeak[$i]++;
        }
        close(INPUT);
        
        $helper->disp("Processing file $input[$i]...");
        croak "\nThe region IDs in file $input[$i] are not unique! Exiting...\n" if (!$self->check_unique(%peakID));
        
        @peakchrs = keys(%peakID);
        
        # Do stuff with siginificant file
        if ($test ne "none")
        {
            $helper->disp("Associating query regions with subject regions in foreground and background region files :");
            $helper->disp("Subject region file : $region");
            $helper->disp("Background region file : $background");
        }
        else
        {
            $helper->disp("Associating query regions with subject regions in region file : $region");
        }
        foreach $currchr (@peakchrs)
        {
            $helper->disp("Queries at $currchr...");
            
            @modes = @{$peakMode{$currchr}};
            @peakids = @{$peakID{$currchr}};
            @peakscores = @{$peakScore{$currchr}} if $pcols[2];
            
            if ($sigID{$currchr}) # Could not have subjects at a specific chromosome
            {
                @geneids = @{$sigID{$currchr}};
                @starts = @{$sigStart{$currchr}};
                @ends = @{$sigEnd{$currchr}};
            
                for ($j=0; $j<@starts; $j++)
                {
                    for ($k=0; $k<@modes; $k++)
                    {
                        $currdist = $self->dist($starts[$j],$ends[$j],$modes[$k]);
                        if ($currdist > $span[0] && $currdist < $span[1])
                        {
                            push(@{$genePeaksSig{$currchr}{$geneids[$j]}},$peakids[$k]);
                            push(@{$peaksGenesSig{$currchr}{$peakids[$k]}},$geneids[$j]);
                            push(@{$distGeneBased{$currchr}{$geneids[$j]}},$currdist);
                            push(@{$distPeakBased{$currchr}{$peakids[$k]}},$currdist);
                            push(@{$peakIndex{$currchr}{$peakids[$k]}},$k);
                            push(@{$geneIndex{$currchr}{$geneids[$j]}},$j);
                            push(@sigallpeakids,$peakids[$k]);
                            push(@sigassgenes,$geneids[$j]);
                            (!$pcols[2]) ? (push(@{$hasPeak{$geneids[$j]}{basename($input[$i])}},$peakids[$k]."_".$currdist)) :
                            (push(@{$hasPeak{$geneids[$j]}{basename($input[$i])}},$peakids[$k]."_".$currdist."_".$peakscores[$k]));
                        }
                    }
                }
            }

            if ($test ne "none") # If no test performed, no need to do anything with background file
            {
                if ($backID{$currchr}) # Could not have genes at chromosome (unlikely for background...)
                {
                    @geneids = @{$backID{$currchr}};
                    @starts = @{$backStart{$currchr}};
                    @ends = @{$backEnd{$currchr}};
                
                    for ($j=0; $j<@starts; $j++)
                    {
                        for ($k=0; $k<@modes; $k++)
                        {
                            $currdist = $self->dist($starts[$j],$ends[$j],$modes[$k]);
                            if ($currdist > $span[0] && $currdist < $span[1])
                            {
                                push(@backallpeakids,$peakids[$k]);
                                push(@backassgenes,$geneids[$j]);
                            }
                        }
                    }
                }
            }
        }
        
        # Get peak counts in significant and background file
        %countSig = $helper->unique(@sigallpeakids);
        %countBack = $helper->unique(@backallpeakids) if ($test ne "none");
        %countGenesSig = $helper->unique(@sigassgenes);
        %countGenesBack = $helper->unique(@backassgenes) if ($test ne "none");
        
        # Run hypergeometric test
        if ($test eq "hypgeom")
        {
            $helper->disp("Running hypergeometric test for each assigned query region..."); 
            @elems = keys(%countSig);
            for ($j=0; $j<@elems; $j++)
            {
                $p = abs(1 - $self->hypergeom_cdf($lensig,$lenback - $lensig,$countBack{$elems[$j]},$countSig{$elems[$j]}));
                ($p*@elems > 1) ? ($cp = 1) : ($cp = $p*@elems); # Bonferroni type correction
                $finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}/$countBack{$elems[$j]}" if ($p < $pval);
            }
        }
        elsif ($test eq "chi2")
        {
            $helper->disp("Running chi-square test for each assigned peak..."); 
            @elems = keys(%countSig);
            for ($j=0; $j<@elems; $j++)
            {
                $p = $self->chisquarecont($countSig{$elems[$j]},$lensig - $countSig{$elems[$j]},$countBack{$elems[$j]},$lenback - $countBack{$elems[$j]});
                ($p*@elems > 1) ? ($cp = 1) : ($cp = $p*@elems); # Bonferroni type correction
                $finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}/$countBack{$elems[$j]}" if ($p < $pval);
            }
        }
        elsif ($test eq "none")
        {
            $helper->disp("No statistical testing performed...");   
            @elems = keys(%countSig);
            for ($j=0; $j<@elems; $j++)
            {
                $p = "NA";
                $cp = "NA";
                $finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}";
            }
        }
        
        my $sp = keys(%finalPeaks);
        my $ap = @elems;
        my $bp = keys(%countBack) if ($test ne "none");
        my $sag = keys(%countGenesSig);
        my $sbg = keys(%countGenesBack) if ($test ne "none");
        if ($test ne "none")
        {
            $helper->disp("\nAssigned peaks in significant list : $ap out of $lenpeak[$i] peaks in $sag out of $lensig genes");
            $helper->disp("Assigned peaks in background : $bp out of $lenpeak[$i] peaks in $sbg out of $lenback genes");
            $helper->disp("Over-represented at p-value<$pval : $sp out of $ap\n");
        }
        else
        {
            $helper->disp("\nAssigned peaks in gene list : $ap out of $lenpeak[$i] peaks in $sag out of $lensig genes");
        }
        
        # Free some memory...
        ($sp,$ap,$bp,$sag,$sbg,@sigassgenes,@backassgenes,%countGenesSig,%countGenesBack) = 
        (undef,undef,undef,undef,undef,undef,undef,undef,undef);
        
        # Construct output
        foreach my $opt (@out)
        {
            if ($opt eq "stats")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                my $co;
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                foreach $co (sort(keys(%finalPeaks)))
                {
                    print OUTPUT "$co\t$finalPeaks{$co}\n";
                }
            }
            if ($opt =~ /gff-peak/)
            { 
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                print OUTPUT "Chromosome\tProgram\tFeature\tStart\tEnd\tp-value\tStrand\tFrame\t".
                             "GeneID\tPeakID\tCorrected p-value\tEnrichment\tDistance\n" if ($opt eq "gff-peak-db");
                foreach $currchr (@peakchrs)
                {
                    foreach $currpeak (keys(%finalPeaks))
                    {   
                        if ($peaksGenesSig{$currchr}{$currpeak})
                        {
                            @currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};
                            my ($q,$cq,$r) = split(/\t/,$finalPeaks{$currpeak});
                            my $pos = 0;
                            foreach $currgene (@currgenes)
                            {   
                                $currout = "$currchr\tHyperAssignPeaks\tPeak\t".
                                           "${$peakStarts{$currchr}}[${$peakIndex{$currchr}{$currpeak}}[$pos]]\t".
                                           "${$peakEnds{$currchr}}[${$peakIndex{$currchr}{$currpeak}}[$pos]]\t".
                                           "$q\t+\t\.\t$currgene\t$currpeak\t$cq\t$r\t".
                                           "${$distPeakBased{$currchr}{$currpeak}}[$pos]\n";
                                print OUTPUT "$currout";
                                $pos++;
                            }
                        }
                    }
                }
                close(OUTPUT);
            }
            if ($opt =~ /gff-gene/)
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                print OUTPUT "Chromosome\tProgram\tFeature\tStart\tEnd\tp-value\tStrand\tFrame\t".
                             "PeakID\tGeneID\tCorrected p-value\tEnrichment\tDistance\n" if ($opt eq "gff-gene-db");
                foreach $currchr (@peakchrs)
                {
                    foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
                    {   
                        my $pos = 0;
                        foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
                        {
                            if ($finalPeaks{$currpeak})
                            {
                                my ($q,$cq,$r) = split(/\t/,$finalPeaks{$currpeak});
                                my $s = ${$sigStart{$currchr}}[${$geneIndex{$currchr}{$currgene}}[$pos]];
                                my $e = ${$sigEnd{$currchr}}[${$geneIndex{$currchr}{$currgene}}[$pos]];
                                my $st = "+";
                                if ($s > $e)
                                {
                                    ($s,$e) = ($e,$s); # Swap them
                                    $st = "-";
                                }
                                $currout = "$currchr\tHyperAssignPeaks\tGene\t$s\t$e\t$q\t$st\t\.\t$currpeak\t".
                                           "$currgene\t$cq\t$r\t${$distGeneBased{$currchr}{$currgene}}[$pos]\n";
                                print OUTPUT "$currout";
                            }
                            $pos++;
                        }
                    }
                }
                close(OUTPUT);
            }
            if ($opt eq "peak")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                foreach $currchr (@peakchrs)
                {
                    foreach $currpeak (keys(%finalPeaks))
                    {   
                        if ($peaksGenesSig{$currchr}{$currpeak})
                        {
                            @currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};            
                            print OUTPUT "$currpeak\t",join("\ ",@currgenes),"\n";
                        }
                    }
                }
                close(OUTPUT);
            }
            if ($opt eq "gene")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                foreach $currchr (@peakchrs)
                {
                    foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
                    {   
                        my @outpeaks;
                        foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
                        {
                            push(@outpeaks,$currpeak) if ($finalPeaks{$currpeak});
                        }
                        print OUTPUT "$currgene\t",join("\ ",@outpeaks),"\n"  if (@outpeaks);
                    }
                }
                close(OUTPUT);
            }
            if ($opt eq "pretty-peak")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                print OUTPUT "PeakID/GeneID\tp-value/Distance\tBonferroni p-value/Strand\tEnrichment\n\n";
                foreach $currchr (@peakchrs)
                {
                    foreach $currpeak (keys(%finalPeaks))
                    {   
                        if ($peaksGenesSig{$currchr}{$currpeak})
                        {
                            print OUTPUT "$currpeak\t$finalPeaks{$currpeak}\n";
                            @currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};
                            my $pos = 0;
                            foreach $currgene (@currgenes)
                            {   
                                my $ct = "+";
                                my $cd = ${$distPeakBased{$currchr}{$currpeak}}[$pos];
                                $ct = "-" if ($cd > 0);
                                print OUTPUT "$currgene\t$cd\t$ct\t\n";
                                $pos++;
                            }
                            print OUTPUT "\n";
                        }
                    }
                }
                close(OUTPUT);
            }
            if ($opt eq "pretty-gene")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                print OUTPUT "GeneID/PeakID\tStrand/p-value\tBonferroni p-value\tEnrichment\tDistance\n\n";
                foreach $currchr (@peakchrs)
                {
                    foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
                    {   
                        my $s = ${$sigStart{$currchr}}[${$geneIndex{$currchr}{$currgene}}[0]];
                        my $e = ${$sigEnd{$currchr}}[${$geneIndex{$currchr}{$currgene}}[0]];
                        my $st = "+";
                        if ($s > $e)
                        {
                            ($s,$e) = ($e,$s); # Swap them
                            $st = "-";
                        }
                        my (@outpeaks,@dis);
                        my $pos = 0;
                        foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
                        {
                            if ($finalPeaks{$currpeak})
                            {
                                push(@outpeaks,$currpeak);
                                push(@dis,${$distGeneBased{$currchr}{$currgene}}[$pos]);        
                            }
                            $pos++;
                        }
                        if (@outpeaks)
                        {
                            print OUTPUT "$currgene\t$st\t\t\t\n";
                            for (my $cpo=0; $cpo<@outpeaks; $cpo++)
                            {
                                print OUTPUT "$outpeaks[$cpo]\t$finalPeaks{$outpeaks[$cpo]}\t$dis[$cpo]\n";
                            }
                            print OUTPUT "\n";
                        }
                    }
                }
                close(OUTPUT);
            }
            if ($opt eq "peakdata")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                print OUTPUT "$phead\n" if ($phead);
                my @fkeys = keys(%finalPeaks);
                foreach $currpeak (@fkeys)
                {
                    print OUTPUT "$allPeakData{$currpeak}\n";
                }
                close(OUTPUT);
            }
            if ($opt eq "bed")
            {
                my $outfile = $self->create_output_file($input[$i],$opt);
                $helper->disp("Writing output in $outfile...");
                open(OUTPUT,">$outfile");
                my @fkeys = keys(%finalPeaks);
                foreach $currpeak (@fkeys)
                {
                    my @spl = split(/\t/,$allPeakData{$currpeak});
                    print OUTPUT $spl[0]."\t".$spl[1]."\t".$spl[2]."\t".$spl[$pcols[0]]."\t".$spl[$pcols[1]]."\t.\n";
                }
                close(OUTPUT);
            }
            $self->print_gene_or_peak($input[$i],"all-peak",%peaksGenesSig) if ($opt eq "all-peak");
            $self->print_gene_or_peak($input[$i],"all-gene",%genePeaksSig) if ($opt eq "all-gene");
        }
    }

    $self->print_matrix(\%hasPeak,"matrix-number",\%hasExpression) if ($helper->smatch("matrix-number",@out));
    $self->print_matrix(\%hasPeak,"matrix-presence",\%hasExpression) if ($helper->smatch("matrix-presence",@out));
    $self->print_matrix(\%hasPeak,"matrix-peaks",\%hasExpression) if ($helper->smatch("matrix-peaks",@out));

    $date = $helper->now;
    $helper->disp("$date - Finished!\n\n");
}

=head2 hypergeom_pdf

Hypergeometric probabilty density function (pdf). There are m "bad" and n "good" balls in an urn.
Pick N of them. The probability of i or more successful selections is: 
(m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
This function is used to perform the hypergeometric test to decide if a genomic region is significantly
"close" to another functional region, depending on the functional regions in the background. This is very
experimental and should not be used until explicitly said in a module update.

    $assigner->hypergeom_pdf($n,$m,$N,$i)

=cut

sub hypergeom_pdf
{   
    my ($self,$n,$m,$N,$i) = @_;
    my $loghyp1 = logfact($m) + logfact($n) + logfact($N) + logfact($m+$n-$N);
    my $loghyp2 = logfact($i) + logfact($n-$i) + logfact($m+$i-$N) + logfact($N-$i) + logfact($m+$n);
    return(exp($loghyp1 - $loghyp2));
}

=head2 hypergeom_pdf

Hypergeometric cumulative distribution function (cdf). This function is used to perform the hypergeometric 
test to decide if a genomic region is significantly "close" to another functional region, depending on 
the functional regions in the background. This is very experimental and should not be used until explicitly 
said in a module update.

    $assigner->hypergeom_cdf($n,$m,$N,$i)

=cut

sub hypergeom_cdf
{  
    my ($self,$n,$m,$N,$i) = @_; 
    my @lessthan = (0..$i);
    my $cum = 0; #;-)
    foreach my $j (@lessthan)
    {
        $cum += hypergeom_pdf($n,$m,$N,$j);
    }
    return($cum);
}

=head2 logfact

Helper function to calculate a hypergeometric test p-value. Efficiently calculates factorial. Internal
use.

=cut

sub logfact 
{
    my ($self,$a);
    return gammaln($a + 1.0);
}

=head2 gammaln

Helper function to calculate a hypergeometric test p-value. Efficiently calculates factorial. Internal
use.

=cut

sub gammaln 
{
    my ($self,$xx) = @_;
    my $eps = 1; 
    my @cof = (76.18009172947146, -86.50532032941677,
               24.01409824083091, -1.231739572450155,
               0.12086509738661e-2, -0.5395239384953e-5);
    my $y = my $x = $xx;
    my $tmp = $x + 5.5;
    $tmp -= ($x + .5) * log($tmp);
    my $ser = 1.000000000190015;
    for my $j (0..5) 
    {
        $ser += $cof[$j]/++$y;
    }
    ($x = $eps) if ($x == 0);
    return(-$tmp + log(2.5066282746310005*$ser/$x));
}

=head2 dist

Distance calculation function among regions, given starting and ending points. Internal use.

    $assigner->dist($start,$end,$anchor)

=cut

sub dist
{
    my ($self,$s,$e,$m) = @_;
    my $d;
    if ($s < $e) # + strand
    {
        #(($m > $s) && ($m < $e)) ? ($d = 0) : ($d = $m - $s);
        $d = $m - $s;
    }
    else # - strand, reversed while reading file
    {
        #(($m < $s) && ($m > $e)) ? ($d = 0) : ($d = $s - $m);
        $d = $s - $m;
    }
    return($d);
}

=head2 chisquarecont

Chi-square contingency table and test.Internal use but may be used also from outside to perform a
chi-square test, given the contingency table (a,b,c,d). The package Math::Cephes is required but this
is checked during module initialization.
    
    $assigner->chisquarecont($a,$b,$c,$d)

=cut

sub chisquarecont
{
    my ($self,$a,$b,$c,$d) = @_;
    my $tot = $a + $b + $c + $d;
    my $A = (($a + $b)*($a + $c))/$tot;
    my $B = (($a + $b)*($b + $d))/$tot;
    my $C = (($c + $d)*($a + $c))/$tot;
    my $D = (($c + $d)*($b + $d))/$tot;
    
    # Chi square statistic
    my $x2 = (($a - $A)**2)/$A + (($b - $B)**2)/$B + (($c - $C)**2)/$C +
             (($d - $D)**2)/$D;

    return(1-chdtr(1,$x2));
}

=head2 print_gene_or_peak

Internal output printing function

    $assigner->print_gene_or_peak($input,$outformat,%resulthash)

=cut

sub print_gene_or_peak
{
    my ($self,$infile,$otype,%inhash) = @_;
    my ($outchr,$ind,$outhash);
    my (@k,@v);
    my $outfilename = $self->create_output_file($infile,$otype);
    $helper->disp("Writing output in $outfilename...");
    my @chrs = keys(%inhash);
    open(OUTPUT,">$outfilename");
    foreach $outchr (sort(@chrs))
    {
        $outhash = $inhash{$outchr};
        @k = keys(%$outhash);
        @v = values(%$outhash);
        for ($ind=0; $ind<@k; $ind++)
        {
            print OUTPUT "$k[$ind]\t";
            print OUTPUT join("\ ",@{$v[$ind]});
            print OUTPUT "\n";
        }
    }
    close(OUTPUT);
}

=head2 print_matrix

Internal output printing function

    $assigner->print_matrix($resulthash,$outformat)

=cut

sub print_matrix
{
    my ($self,$inhash,$type,$exprhash) = @_;
    my ($row,$column,$colhash);
    my $outfilename = $self->create_output_file(${$self->get("input")}[0],$type);
    $helper->disp("Writing output in $outfilename...");
    my @rows = keys(%$inhash);
    open(OUTPUT,">$outfilename");
    my $headhash = $inhash->{$rows[0]};
    my @headers = keys(%$headhash);
    (!defined($exprhash->{"header"})) ? (print OUTPUT "GeneID\t",join("\t",@headers),"\n") :
    (print OUTPUT "GeneID\t",join("\t",@headers),"\t",$exprhash->{"header"},"\n");
    if ($type =~ m/number/i)
    {
        foreach $row (@rows)
        {
            $colhash = $inhash->{$row};
            my @tmp = keys(%$colhash);
            my @v;
            foreach $column (@tmp)
            {
                (defined($colhash->{$column})) ?
                (push(@v,scalar(@{$colhash->{$column}}))) :
                (push(@v,0));
            }
            (!$exprhash->{"data"}->{$row}) ? (print OUTPUT "$row\t",join("\t",@v),"\n") :
            (print OUTPUT "$row\t",join("\t",@v),"\t",$exprhash->{"data"}->{$row},"\n");
        }
    }
    elsif ($type =~ m/presence/i)
    {
        foreach $row (@rows)
        {
            $colhash = $inhash->{$row};
            my @tmp = keys(%$colhash);
            my @v;
            foreach $column (@tmp)
            {
                (defined($colhash->{$column})) ? (push(@v,"+")) : (push(@v,"-"));
            }
            (!$exprhash->{"data"}->{$row}) ? (print OUTPUT "$row\t",join("\t",@v),"\n") :
            (print OUTPUT "$row\t",join("\t",@v),"\t",$exprhash->{"data"}->{$row},"\n");
        }
    }
    elsif ($type =~ m/peaks/i)
    {
        foreach $row (@rows)
        {
            $colhash = $inhash->{$row};
            my @tmp = keys(%$colhash);
            my @v;
            foreach $column (@tmp)
            {
                (defined($colhash->{$column})) ?
                (push(@v,join("; ",@{$colhash->{$column}}))) :
                (push(@v,"NP"));
            }
            (!$exprhash->{"data"}->{$row}) ? (print OUTPUT "$row\t",join("\t",@v),"\n") :
            (print OUTPUT "$row\t",join("\t",@v),"\t",$exprhash->{"data"}->{$row},"\n");
        }
    }
    close(OUTPUT);
}

=head2 create_output_file

Automatic output filename creation for the module outputs. Internal use.

    $assigner->create_output_file($input,$outformat)

=cut

sub create_output_file
{
    my ($self,$in,$type) = @_;
    my $ext;
    my ($base,$dir) = fileparse($in,'\.[^.]*');
    if ($type =~/gff/)
    {
        ($type =~/db/) ? ($ext = ".txt") : ($ext = ".gff");
    }
    elsif ($type =~/bed/) { $ext = ".bed" }
    else { $ext = ".txt" }
    if ($type =~ /matrix/)
    {
        return($dir."gene-peak-$type.txt");
    }
    if ($type =~ /bed/)
    {
        return($dir.$base."_ASSIGNED".$ext);
    }
    else
    {
        return($dir.$base."_ASSIGNED_".$type.$ext);
    }
}

=head2 check_unique

Check uniqueness of hash values. Required for internal control. Internal use.

    $assigner->check_unique(%hash)

=cut

sub check_unique
{
    my ($self,%h) = @_;
    my @vals = values(%h);
    my %ch = $helper->unique(@vals);
    (scalar @vals == scalar keys(%ch)) ? (return(1)) : (return(0));
}

=head2 change_params

Massively change the parameters of an HTS::Tools::Assign object.

    $assigner->change_params({'input' => 'another_file','region' => 'mouse-exon'})
    $assigner->run;
    
=cut

sub change_params
{
    my ($self,$params) = @_;
    
    # Validate the new parameters 
    my $checker = HTS::Tools::Paramcheck->new();
    $checker->set("tool","assign");
    $checker->set("params",$params);
    $params = $checker->validate;
    
    # If validator does not complain, change the parameters
    while (my ($name,$value) = each(%$params))
    {
        $self->set($name,$value);
    }
    return($self);
}

=head2 get

HTS::Tools::Assign object getter

    my $param_value = $assigner->get("param_name")
=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Assign object setter

    $assigner->set("param_name","param_value")
    
=cut

sub set
{
    my ($self,$name,$value) = @_;
    $self->{$name} = $value;
    return($self);
}

=head1 DEPENDENCIES

Tie::IxHash::Easy (optional)
Math::Cephes (optional)

= head1 TODO

=head1 AUTHOR

Panagiotis Moulos, C<< <moulos at fleming.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc HTS::Tools::Assign


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

1; # End of HTS::Tools::Assign
