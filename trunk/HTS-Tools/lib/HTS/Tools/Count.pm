=head1 NAME

HTS::Tools::Count - Count short sequence reads in genomic regions

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Program to count reads in bed file(s) over specific genomic regions in a genome and count. It also 
addresses the problem of partial overlaping between a read and a genomic regions by 3 possible ways:

=over 4

=item *

a. By providing an overlaping coefficient between 0 and 1, the program adds a read into a genomic region 
if the part of the read within the boundaries (start or end) of the genomic region is >= coefficient*read length. 
For example, if a read has length 200 and the overlaping coefficient is 0.9, then 0.9*200 = 180 bases 
should be within the region boundaries to be included in that region. This should provide more accuracy 
when having large reads. It can be used in the case of small reads too, but lucks automation. It can
also split each genomic region in bins of constant length and count reads in individual bins and calculate
several statistics for each genomic region based on the bin stats.

=item *

b. By using a linear probabilistic score. In this case, if a read is located half inside a region, it
is added up to that region, else a score is calculated for the read based on the formula: score = #bps 
inside region/(length of read/2). This score lies between 0 and 1 and it is compared to a uniform random 
number p between 0 and 1. If p<score then this tag is added up to the genomic region else discarded. 
This scoring scheme should be used with reads of small length and should be avoided in the case of 
large reads because of causing possible under-representation in promoter regions, in case one wishes 
to summarize reads in gene regions. However it could perform well when studying promoter or binding regions.

=item *

c. By using an exponential probabilistic score. This case, works as case (b) but the scoring function 
is score = exp(-c^2/#bps inside region) where c is constant. The constant c determines the steepness 
of the score. The lower the value is, the higher is the probability to include in the region tags whose
larger part lies outside the region. This scoring scheme should be used with reads of small length and
for reasons similar to those of (b) but would maybe perform better when used in promoter or binding regions.

=back

Apart from the bed files, the program requires also a bed file (3 first columns should be chromosome, 
start, end and the 4th a UNIQUE region identification, e.g. Ensembl ID, could also contain further 
information in extra columns) which contains the genomic regions of interest, in which number of reads 
should be summarized. The program returns bed file(s) that contain the provided genomic regions with 
any other information together with read counts per region. The region file can be automatically downloaded
among a variety of predefined regions from Ensembl, UCSC or RefSeq. More documentation to come...

    use HTS::Tools::Count;
    my %params = (
        'input' => ['normal_pol2.bed','disease_pol2.bed'],
        'region' => 'human-gene',
        'source' => 'ucsc',
        'split' => 1000
    )
    my $counter = HTS::Tools::Count->new(\%params);
    $counter->run;
    
The acceptable parameters are as follows:

=over 4

=item I<input> B<(required)>

Input file(s). The input bed files can be 3, 6 or n column bedfiles, only the first 3 are used. More
formats will be added soon.

=item I<region> B<(required)>

Genomic regions file. The first 3 columns should be of the same structure as in bed files. The 4th 
column should contain a UNIQUE identifier for each region (e.g. Ensembl IDs). The rest columns can 
contain additional data about the regionsn (e.g. annotation elements, descriptions etc.). Instead of 
a local file, it can be one of the following, which will automatically download and use genomic regions 
from the source defined by the I<source> parameter:

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

=item I<source> B<(optional)>

Use this option to set the online data source in the case of selecting one of the predefined region 
templates with I<region>. It can be one of "ucsc", "refseq" or "ensembl" and Default to "ensembl".

=item I<gversion> B<(optional)>

Use this option to set the version of the genome from which the predefined region templates templates (with
I<region>) will be fetched. For human, it can be "hg18", "hg19". For mouse, it can be "mm9", "mm10". For rat,
it can be "rn5". For fly, it can be "dm3". For zebrafish, it can be "danRer7".

=item I<splicing> B<(optional)>

Use this option with I<source> to determine whether the canonical or alternatively spliced transcripts will
be used for counting. It can be "canonical" or "alternative" and defaults to "canonical"

=item I<sort> B<(optional)>

Use this option if you wish to sort the input files first. This is not obligatory as the new implementation
using the interval trees structure does not require sorting, however, it will run faster if the input files
are already sorted. If you wish to sort the files outside HTS::Tools::Count, in Linux systems, the command
should look like this: 

    sort -k1,1 -k2g,2 inputfile > outputfile 

In Windows systems, it would be better to use a spreadsheet program like Excel to perform sorting by 
columns. Both the input bed files and the region file should be sorted. If the structure of the region 
file is as described (1st column: chromosome, 2nd column: start 3rd column: end, 4th column: unique ID 
etc.) the sorting command is exactly the same as in common bed files (3 or 6 columns). Keep in mind that
in future versions, the I<sort> parameter will be deprecated.

=item I<percent> B<(optional)>

Use this option to provide the overlaping coefficient according to which tags that partially fall outside 
provided genomic regions will be assigned to genomic regions. Defaults to 0.95 and is the default algorithm 
for assigning tags that fall partially inside provided genomic regions.

=item I<lscore> B<(optional)>

Use this option to use the linear probabilistic scoring scheme to address the problem of partial tag
overlap with the provided genomic regions. Please see header of the script file for description. Defaults 
to 0 (not use).

=item I<escore> B<(optional)>

Use this option to use the exponential probabilistic scoring scheme to address the problem of partial tag
overlap with the provided genomic regions. Please see the SYNOPSIS section for a description. Defaults to 0
(not use).

=item I<constant> B<(optional)>

Use this option to provide the constant for the exponential scoring scheme (see description of I<escore>
option). Defaults to 3.

=item I<small> B<(optional)>

Use this option if you wish to take into consideration genomic regions in your genomic regions file which 
are smaller than the tag length (rare but may occur, especially with customized genomic regions). Defaults 
to 0 (not use).

=item I<split> B<(optional)>

Use this option if you wish to further split your genomic regions in smaller areas. In this case, tags 
will be counted per area and the distribution of counts per area will be returned as a column in the 
output file. The argument is the length of sub-areas. The default is 1000 and splits the regions per 
1kb, unless the regions are smaller. In this case, the area will consist of only one sub-area of length 
equal to that of the area. Note also that when using this option, tags that are found inside sub-areas 
are not assigned to those sub-areas based on scoring schemes (options I<percent>, I<lscore> and I<escore>)
but tags are assigned based on the location of their center.

=nbins I<nbins> B<(optional)>

Use this option if you wish to further split your genomic regions in a predefined number (nbins) of smaller
areas. The same things as in I<split> apply.

=item I<stats> B<(optional)>

Use this option to also return basic statistics of counts in the windows used returned by using I<split>.
It should be set to 1 (defaults to 0).

=item I<ncore> B<(optional)>

If the machine has multicore processor(s) and the package Parallel::ForkManager is installed, you can use
parallel processing. Default is 1 and can go up to 12.

=item I<keeporder> B<(optional)>

Use this parameter if you want to force the lines of the output counts table to be in the same order (e.g.
sorted per chromosome or gene name) as the region file. This is accomplished through the use of the module
Tie::IxHash::Easy which must be present in your machine. If the module is not present, the I<keeporder>
option is deactivated. Keep in mind that maintaining the order requires slighlty more memory during runtime.

=item I<output> B<(optional)>

A file to write the output to. If "auto", then it generates an automatic filename in the folder where the
input files are. If not provided, output is written to STDOUT.

=item I<silent> B<(optional)>

Set this to 1 if you want to turn informative messages off.

=back

=head1 OUTPUT

A table-like test file containing the reads inside each genomic region for each input reads file.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Count;

our $MODNAME = "HTS::Tools::Count";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Short sequence read counting in genomic regions.";

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use File::Basename;
use File::Temp;
use File::Spec;
use File::Path qw(make_path remove_tree);
use IntervalTree;

use HTS::Tools::Fetch;
use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper);

use constant MAXCORES => 12;

BEGIN {
    $helper = HTS::Tools::Utils->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

The HTS::Tools::Count object constructor. It accepts a set of parameters that are required to run the
counter and get the output.

    my $counter = HTS::Tools::Count->new({'input' => 'myfile.bed','region' => 'ensembl_genes.txt'});

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
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "count", "params" => $params});
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init($params)

HTS::Tools::Count object initialization method. NEVER use this directly, use new instead.

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

The HTS::Tools::Count run subroutine. It runs the counter with the given parameters in the constructor.

    $counter->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    my $regionfile = $self->get("region");
    my @infile = @{$self->get("input")};
    my @originfile = @infile; # Keep original filenames (verbose purposes)
    
    $regionfile = $self->fetch_regions($regionfile) if (! -f $regionfile);
    $self->set("region",$regionfile) if ($regionfile ne $self->get("region"));
    
    my ($chromosome,$gencounts,$splitcounts,$theHeader) = $self->read_region_file($regionfile,\@originfile);
    ($gencounts,$splitcounts) = $self->count_all_reads($chromosome,$gencounts,$splitcounts,\@infile,\@originfile);
    $self->write_reads($chromosome,$gencounts,$splitcounts,\@originfile,$theHeader);

    $helper->disp("Finished!\n\n");
}

sub fetch_regions
{
    my ($self,$regionfile,$ver,$source,$splicing) = @_;
    my $fetcher = HTS::Tools::Fetch->new({"tmpdir" => $self->get("tmpdir"), "silent" => $self->get("silent")});
    $source = $self->get("source") unless($source);
    $splicing = $self->get("splicing") unless($splicing);
    $ver = $self->get("gversion") unless($ver);
    
    use v5.14;
    given($regionfile)
    {
        when(/human-gene|mouse-gene|rat-gene|fly-gene|zebrafish-gene/i)
        {
            given($source)
            {
                when(/ensembl/i)
                {
                    $regionfile = $fetcher->fetch_ensembl_genes($regionfile);
                    $regionfile = $fetcher->sort_ensembl_genes($regionfile);
                }
                when(/ucsc|refseq/i)
                {
                    $regionfile = $fetcher->fetch_ucsc_genes($regionfile,$ver,$splicing,$source);
                }
            }
        }
        when(/human-exon|mouse-exon|rat-exon|fly-exon|zebrafish-exon/i)
        {
            given($source)
            {
                when(/ensembl/i)
                {
                    $regionfile = $fetcher->fetch_ensembl_exons($regionfile);
                    $regionfile = $fetcher->sort_ensembl_exons($regionfile);
                }
                when(/ucsc|refseq/i)
                {
                    $regionfile = $fetcher->fetch_ucsc_exons($regionfile,$ver,$splicing,$source);
                }
            }
        }
        when(/human-(5|3)utr|mouse-(5|3)utr|rat-(5|3)utr|fly-(5|3)utr|zebrafish-(5|3)utr/i)
        {
            given($source)
            {
                when(/ensembl/i)
                {
                    ($regionfile =~ m/5utr/i) ?
                    ($regionfile = $fetcher->fetch_ensembl_utr($regionfile,5)) :
                    ($regionfile = $fetcher->fetch_ensembl_utr($regionfile,3));
                    $regionfile = $fetcher->sort_ensembl_exons($regionfile);
                }
                when(/ucsc|refseq/i)
                {
                    ($regionfile =~ m/5utr/i) ?
                    ($regionfile = $fetcher->fetch_ucsc_utr($regionfile,$ver,$splicing,$source,5)) :
                    ($regionfile = $fetcher->fetch_ucsc_utr($regionfile,$ver,$splicing,$source,3));
                }
            }
        }
        when(/human-cds|mouse-cds|rat-cds|fly-cds|zebrafish-cds/i)
        {
            given($source)
            {
                when(/ensembl/i)
                {
                    $regionfile = $fetcher->fetch_ensembl_cds($regionfile);
                    $regionfile = $fetcher->sort_ensembl_exons($regionfile);
                }
                when(/ucsc|refseq/i)
                {
                    $regionfile = $fetcher->fetch_ucsc_cds($regionfile,$ver,$splicing,$source);
                }
            }
        }
        default
        {
            ($regionfile,@{$self->get("input")}) = $self->sort_inputs($regionfile,@{$self->get("input")})
                if ($self->get("sort"));
        }
    }

    return($regionfile);
}

=head2 split_area

Creates the bins in the region file genomic regions. Internal use.

    $counter->split_area($area_start,$area_end,$bin_size);
    
=cut

sub split_area
{
    my ($self,$start,$end,$splitlen) = @_;
    
    my $regstart;
    my $regend = $start - 1;
    my $subtree = IntervalTree->new();
    my $regcount = 0;
    
    while ($regend < $end)
    {
        $regstart = $regend + 1;
        $regend = $regstart + $splitlen - 1;
        $subtree->insert($regstart,$regend,{"binid" => $regcount});
        $regcount++;
    }

    # Note: there is no remove function for this CPAN implementation of IntervalTree. It is not such
    # a harm as the overlap of any read will be found even if the subtree goes beyond the splitted
    # genomic region
    
    return $subtree;
}

=head2 read_region_file

Region file parser and constructor of the internal hash representation. Requires Tie::IxHash::Easy.
Internal use.

    $counter->read_region_file($regionfile,[@input_files]);
    
=cut

sub read_region_file
{
    my $self = shift @_;
    my $regionfile = shift @_;
    my @originfile = @{shift @_};
    @originfile = @{$self->get("input")} unless(@originfile);

    my (%chromosome,%gencounts,%splitcounts);
    if ($self->get("keeporder"))
    {
        tie %chromosome, "Tie::IxHash::Easy";
        tie %gencounts, "Tie::IxHash::Easy";
    }
    
    my ($theHeader,$chr,$start,$end,$uid);
    my ($f,$regstart,$regend,$regconts,$regcount,$regline,$step);
    my @rest;
    my %strands = ("+" => 1,"-" => -1,"1" => 1,"-1" => -1,"F" => 1,"R" => -1);

    $helper->disp("Reading genomic regions file...");
    $helper->disp("...also splitting genomic regions in areas of ".$self->get("split")." base pairs...") if ($self->get("split"));
    
    open(REGIONS,$regionfile) or croak "\nThe file $regionfile does not exist!\n";
    $regline = <REGIONS>;
    $theHeader = $helper->decide_header($regline);
    seek(REGIONS,0,0) if (!$theHeader);
    while ($regline = <REGIONS>)
    {
        next if ($regline =~ m/chrM|rand|chrU|hap/i);
        $regline =~ s/\r|\n$//g; # Make sure to remove carriage returns
        ($chr,$start,$end,$uid,@rest) = split(/\t/,$regline);
        
        $chromosome{$chr} = IntervalTree->new() if (!$chromosome{$chr});

        # Fill the IntervalTree
        $chromosome{$chr}->insert_interval(
            IntervalTree::Interval->new(
                $start,$end,{
                    "id" => $uid,
                    "rest" => join("\t",@rest)
                },$chr,
                ($rest[1] && grep {$_ eq $rest[1]} keys(%strands)) ? ($strands{$rest[1]}) : (undef)));

        # Initialize the main count hash
        foreach $f (@originfile)
        {
            ($self->get("ncore") > 1) ? ($gencounts{basename($f)}{$chr}{$uid} = 0) :
                ($gencounts{$chr}{$uid}{basename($f)} = 0);
        }

        # Initialize the bin count hash if requested
        $step = $self->get("split") if ($self->get("split"));
        $step = ceil(($end-$start)/$self->get("nbins")) if ($self->get("nbins"));
        if ($step)
        {
            foreach $f (@originfile)
            {       
                $regcount = 0;
                $regend = $start - 1;
                while ($regend < $end)
                {
                    $regstart = $regend + 1;
                    $regend = $regstart + $step - 1;
                    ($self->get("ncore") > 1) ? ($splitcounts{basename($f)}{$chr}{$uid}{$regcount} = 0) :
                        ($splitcounts{$chr}{$uid}{basename($f)}{$regcount} = 0);
                    $regcount++;
                }
            }
        }
    }
    close(REGIONS);

    return(\%chromosome,\%gencounts,\%splitcounts,$theHeader);
}

=head2 count_all_reads

Wrapper for the main counting subroutine. Internal use.

    $counter->count_all_reads($region_structure,$region_counts,$split_counts,@input_files,[@input_files]);
    
=cut

sub count_all_reads
{   
    my $self = shift @_;
    my $chromosome = shift @_;
    my $gencounts = shift @_;
    my $splitcounts = shift @_;
    my @infile = @{shift @_};
    my @originfile = @{shift @_};
    @originfile = @infile unless(@originfile);
    
    my @base;
    foreach my $of (@originfile)
    {
        my $bb = basename($of);
        push(@base,$bb);
    }

    if ($self->get("ncore") == 1)
    {
        for (my $i=0; $i<@infile; $i++)
        {
            ($gencounts,$splitcounts) = $self->count_reads($chromosome,$gencounts,$splitcounts,$infile[$i],$originfile[$i]);
        }
    }
    else
    {
        my (%vhash1,%vhash2);
        foreach (my $i=0; $i<@infile; $i++)
        {
            $vhash1{$base[$i]} = $infile[$i];
            $vhash2{$base[$i]} = $originfile[$i];
        }
        my $pm = new Parallel::ForkManager($self->get("ncore"));
        $pm -> run_on_finish(sub {
            my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data_structure_reference) = @_;
            if (defined($data_structure_reference)) {
                my @hr = @{$data_structure_reference};
                while(my($chr,$genehash) = each(%{$hr[0]}))
                {       
                    while(my($genename,$count) = each(%{$genehash}))
                    {
                        $gencounts->{$ident}->{$chr}->{$genename} = $count;
                    }
                }
                if ($self->get("split"))
                {
                    while(my($chr,$genehash) = each(%{$hr[1]}))
                    {       
                        while(my($genename,$idhash) = each(%{$genehash}))
                        {
                            while(my($id,$scount) = each(%{$idhash}))
                            {
                                $splitcounts->{$ident}->{$chr}->{$genename}->{$id} = $scount;
                            }
                        }
                    }
                }
            } else {
                print STDERR qq|No message received from child process $pid!\n|;
            }
        });
        foreach my $ff (@base)
        {
            my $pid = $pm->start($ff) and next;
            my @arr = $self->count_reads_multi($chromosome,$gencounts->{$ff},$splitcounts->{$ff},$vhash1{$ff},$vhash2{$ff});
            $pm->finish(0,\@arr);
        }
        $pm->wait_all_children;
    }
    
    return($gencounts,$splitcounts);
}

=head2 count_reads

Single core main read counter. Internal use.

    $counter->count_reads($region_structure,$region_counts,$split_counts,$input_file,[$input_file]);
    
=cut

sub count_reads
{
    my ($self,$chromosome,$gencounts,$splitcounts,$infile,$originfile) = @_;
    $originfile = $infile unless($originfile);

    my ($bedchr,$bedstart,$bedend,$bedrest);
    my ($i,$j,$currstart,$currend,$gene);
    my ($overstru,$sostru,$subtree,$bin);
    my ($diff,$score,$p,$filename);

    $bin = $self->get("split") if ($self->get("split"));
    $bin = $self->get("nbins") if ($self->get("nbins"));
    my $c = $self->get("constant");
    
    $helper->disp("Reading and processing file $originfile...");

    $filename = basename($originfile);
    open(BEDFILE,$infile);

    while (<BEDFILE>)
    {
        chomp $_;
        ($bedchr,$bedstart,$bedend,$bedrest) = split(/\t/,$_);

        $overstru = $chromosome->{$bedchr}->find($bedstart,$bedend);

        # If overlap found, check conditions that satisfy our parameters
        if ($overstru)
        {
            for ($i=0; $i< scalar @$overstru; $i++)
            {
                $currstart = $overstru->[$i]->{"start"};
                $currend = $overstru->[$i]->{"end"};
                $gene = $overstru->[$i]->{"value"}->{"id"};
            
                # Security, too small region, tag does not fit? Mark that there is something there if user wishes and proceed
                if ($currstart >= $bedstart && $currend <= $bedend)
                {
                    $gencounts->{$bedchr}->{$gene}->{$filename}++ if ($self->get("small"));
                    if ($bin && $self->get("small")) # If splitting of regions
                    {
                        $subtree = $self->split_area($currstart,$currend,$bin);
                        $sostru = $subtree->find($bedstart,$bedend);
                        for ($j=0; $j< scalar @$sostru; $j++)
                        {
                            $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                        }
                    }
                }
                # Ideal case... whole tag inside gene... no problem
                elsif ($currstart <= $bedstart && $currend >= $bedend)
                {
                    $gencounts->{$bedchr}->{$gene}->{$filename}++;
                    if ($bin) # If splitting of regions
                    {
                        $subtree = $self->split_area($currstart,$currend,$bin);
                        $sostru = $subtree->find($bedstart,$bedend);
                        for ($j=0; $j< scalar @$sostru; $j++)
                        {
                            $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                        }
                    }
                }
                # Part of tag outside after region
                elsif ($currstart < $bedstart && $currend < $bedend && $bedstart < $currend)
                {
                    $diff = $currend - $bedstart;
                    if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                    {
                        if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside no need for further check
                        {
                            $gencounts->{$bedchr}->{$gene}->{$filename}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                        else # Use probabilistic scoring scheme
                        {
                            $score = ($currend - $bedstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                            $score = exp(-$c**2/($currend - $bedstart)) if ($self->get("escore")); # Exponential
                            $p = rand();
                            if ($p < $score)
                            {
                                $gencounts->{$bedchr}->{$gene}->{$filename}++; 
                                if ($bin) # If splitting of regions
                                {
                                    $subtree = $self->split_area($currstart,$currend,$bin);
                                    $sostru = $subtree->find($bedstart,$bedend);
                                    for ($j=0; $j< scalar @$sostru; $j++)
                                    {
                                        $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                    }
                                }
                            }
                        }
                    }
                    else # Simple overlap
                    {
                        if ($diff >= $self->get("percent")*($bedend - $bedstart))
                        {
                            $gencounts->{$bedchr}->{$gene}->{$filename}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                    }
                }
                # Part of tag outside before region
                elsif ($currstart > $bedstart && $currend > $bedend && $bedend > $currstart)
                {
                    $diff = $bedend - $currstart;
                    if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                    {
                        if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside
                        {
                            $gencounts->{$bedchr}->{$gene}->{$filename}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                        else # Use probabilistic scoring scheme
                        {
                            $score = ($bedend - $currstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                            $score = exp(-$c**2/($bedend - $currstart)) if ($self->get("escore")); # Exponential
                            $p = rand();
                            if ($p < $score)
                            {
                                $gencounts->{$bedchr}->{$gene}->{$filename}++;
                                if ($bin) # If splitting of regions
                                {
                                    $subtree = $self->split_area($currstart,$currend,$bin);
                                    $sostru = $subtree->find($bedstart,$bedend);
                                    for ($j=0; $j< scalar @$sostru; $j++)
                                    {
                                        $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                    }
                                }
                            }
                        }
                    }
                    else # Simple overlap
                    {
                        if ($diff >= $self->get("percent")*($bedend - $bedstart))
                        {
                            $gencounts->{$bedchr}->{$gene}->{$filename}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$bedchr}->{$gene}->{$filename}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    close(BEDFILE);

    return($gencounts,$splitcounts);
}

=head2 count_reads_multi

Multiple core main read counter. Internal use.

    $counter->count_reads_multi($region_structure,$region_counts,$split_counts,$input_file,[$input_file]);
    
=cut

sub count_reads_multi
{
    my ($self,$chromosome,$gencounts,$splitcounts,$infile,$originfile) = @_;
    $originfile = $infile unless($originfile);

    my ($bedchr,$bedstart,$bedend,$bedrest);
    my ($i,$j,$currstart,$currend,$gene);
    my ($overstru,$sostru,$subtree,$bin);
    my ($diff,$score,$p,$filename);

    $bin = $self->get("split") if ($self->get("split"));
    $bin = $self->get("nbins") if ($self->get("nbins"));
    my $c = $self->get("constant");
    
    $helper->disp("Reading and processing file $originfile...");

    $filename = basename($originfile);
    open(BEDFILE,$infile) or croak "Can't open $infile! Does it exist?\n";

    while (<BEDFILE>)
    {
        chomp $_;
        ($bedchr,$bedstart,$bedend,$bedrest) = split(/\t/,$_);

        ($chromosome->{$bedchr}) ?
        ($overstru = $chromosome->{$bedchr}->find($bedstart,$bedend)) :
        (next);

        # If overlap found, check conditions that satisfy our parameters
        if ($overstru)
        {
            for ($i=0; $i< scalar @$overstru; $i++)
            {
                $currstart = $overstru->[$i]->{"start"};
                $currend = $overstru->[$i]->{"end"};
                $gene = $overstru->[$i]->{"value"}->{"id"};
            
                # Security, too small region, tag does not fit? Mark that there is something there if user wishes and proceed
                if ($currstart >= $bedstart && $currend <= $bedend)
                {
                    $gencounts->{$bedchr}->{$gene}++ if ($self->get("small"));
                    if ($bin && $self->get("small")) # If splitting of regions
                    {
                        $subtree = $self->split_area($currstart,$currend,$bin);
                        $sostru = $subtree->find($bedstart,$bedend);
                        for ($j=0; $j< scalar @$sostru; $j++)
                        {
                            $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                        }
                    }
                }
                # Ideal case... whole tag inside gene... no problem
                elsif ($currstart <= $bedstart && $currend >= $bedend)
                {
                    $gencounts->{$bedchr}->{$gene}++;
                    if ($bin) # If splitting of regions
                    {
                        $subtree = $self->split_area($currstart,$currend,$bin);
                        $sostru = $subtree->find($bedstart,$bedend);
                        for ($j=0; $j< scalar @$sostru; $j++)
                        {
                            $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                        }
                    }
                }
                # Part of tag outside after region
                elsif ($currstart < $bedstart && $currend < $bedend && $bedstart < $currend)
                {
                    $diff = $currend - $bedstart;
                    if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                    {
                        if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside no need for further check
                        {
                            $gencounts->{$bedchr}->{$gene}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                        else # Use probabilistic scoring scheme
                        {
                            $score = ($currend - $bedstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                            $score = exp(-$c**2/($currend - $bedstart)) if ($self->get("escore")); # Exponential
                            $p = rand();
                            if ($p < $score)
                            {
                                $gencounts->{$bedchr}->{$gene}++; 
                                if ($bin) # If splitting of regions
                                {
                                    $subtree = $self->split_area($currstart,$currend,$bin);
                                    $sostru = $subtree->find($bedstart,$bedend);
                                    for ($j=0; $j< scalar @$sostru; $j++)
                                    {
                                        $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                    }
                                }
                            }
                        }
                    }
                    else # Simple overlap
                    {
                        if ($diff >= $self->get("percent")*($bedend - $bedstart))
                        {
                            $gencounts->{$bedchr}->{$gene}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                    }
                }
                # Part of tag outside before region
                elsif ($currstart > $bedstart && $currend > $bedend && $bedend > $currstart)
                {
                    $diff = $bedend - $currstart;
                    if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                    {
                        if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside
                        {
                            $gencounts->{$bedchr}->{$gene}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                        else # Use probabilistic scoring scheme
                        {
                            $score = ($bedend - $currstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                            $score = exp(-$c**2/($bedend - $currstart)) if ($self->get("escore")); # Exponential
                            $p = rand();
                            if ($p < $score)
                            {
                                $gencounts->{$bedchr}->{$gene}++;
                                if ($bin) # If splitting of regions
                                {
                                    $subtree = $self->split_area($currstart,$currend,$bin);
                                    $sostru = $subtree->find($bedstart,$bedend);
                                    for ($j=0; $j< scalar @$sostru; $j++)
                                    {
                                        $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                    }
                                }
                            }
                        }
                    }
                    else # Simple overlap
                    {
                        if ($diff >= $self->get("percent")*($bedend - $bedstart))
                        {
                            $gencounts->{$bedchr}->{$gene}++;
                            if ($bin) # If splitting of regions
                            {
                                $subtree = $self->split_area($currstart,$currend,$bin);
                                $sostru = $subtree->find($bedstart,$bedend);
                                for ($j=0; $j< scalar @$sostru; $j++)
                                {
                                    $splitcounts->{$filename}->{$bedchr}->{$gene}->{$sostru->[$j]->{"value"}->{"binid"}}++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    close(BEDFILE);

    return($gencounts,$splitcounts);
}

=head2 write_reads

Main output writer. Internal use.

    $counter->write_reads($region_structure,$region_counts,$split_counts,$input_file,[$input_file],$header);
    
=cut

sub write_reads
{
    my $self = shift @_;
    my $chromosome = shift @_;
    my $gencounts = shift @_;
    my $splitcounts = shift @_;
    my @originfile = @{shift @_};
    my $theHeader = shift @_;
    @originfile = @{$self->get("input")} unless(@originfile);
    
    my ($filename,$retrhash,$outhead,$outsubhead,$chr,$tree,$finalout);
    my ($mean,$median,$stdev,$mad,$jdistr);
    my (@base,@distr);
    
    foreach my $of (@originfile)
    {
        my $bb = fileparse($of,'\.[^.]*');
        push(@base,$bb);
    }

    my $out;
    if ($self->get("output"))
    {
        open(OUTPUT,">".$self->get("output"));
        $out = *OUTPUT;
    }
    else { $out = *STDOUT; }

    if ($self->get("output"))
    {
        ($self->get("stats")) ? ($helper->disp("Writing reads and calculating statistics per genomic regions in ".$self->get("output")."...")) :
            $helper->disp("Writing reads per genomic regions for file in ".$self->get("output")."...");
    }
    else
    {
        ($self->get("stats")) ? ($helper->disp("Writing reads and calculating statistics per genomic regions to standard output...\n")) :
            $helper->disp("Writing reads per genomic regions for file to standard output...\n");
    }

    if ($theHeader) # Add a header based on what we did plus additional data given
    {
        $outhead = $theHeader;
        if (scalar @base == 1)
        {
            $outhead .= "\t".$base[0];
            $outhead .= "\tSub-area Counts/$self->get(\"split\") bp" if ($self->get("split")); # Return also sub-area distributions
            $outhead .= "\tMean Counts/$self->get(\"split\") bp\tMedian Counts/$self->get(\"split\") bp\tStDev Counts\tMAD Counts" if ($self->get("stats"));
            print $out "$outhead\n";
        }
        else
        {
            foreach my $b (@base)
            {
                $outhead .= "\t".$b;
                $outhead .= "\t" if ($self->get("split"));
                $outhead .= "\t"x4 if ($self->get("stats"));
            }
            # We must also print a subheader in case of reporting stats
            $outsubhead = "Counts";
            $outsubhead .= "\tSub-area Counts/$self->get(\"split\") bp" if ($self->get("split"));
            $outsubhead .= "\tMean Counts/$self->get(\"split\") bp\tMedian Counts/$self->get(\"split\") bp\tStDev Counts\tMAD Counts" if ($self->get("stats"));
            $outsubhead = "$outsubhead"x(scalar @base);
            $outsubhead = ("\t"x(scalar split(/\t/,$theHeader))).$outsubhead;
            print $out "$outhead\n";
            print $out "$outsubhead\n" if ($self->get("split"));
        }
    }

    if ($self->get("ncore") == 1)
    {
        while(($chr,$tree) = each(%$chromosome))
        {
            $tree->traverse(sub {
                if ($_[0]->{"interval"}->{"value"}->{"rest"}) {
                    $finalout = $chr."\t".
                        $_[0]->{"interval"}->{"start"}."\t".
                        $_[0]->{"interval"}->{"end"}."\t".
                        $_[0]->{"interval"}->{"value"}->{"id"}."\t".
                        $_[0]->{"interval"}->{"value"}->{"rest"};
                }
                else {
                    $finalout = $chr."\t". 
                        $_[0]->{"interval"}->{"start"}."\t".
                        $_[0]->{"interval"}->{"end"}."\t".
                        $_[0]->{"interval"}->{"value"}->{"id"};
                }
                foreach $filename (@originfile) {
                    if ($self->get("split") || $self->get("nbins")) { # Retrieve calculated distributions per region 
                        $retrhash = $splitcounts->{$chr}->{$_[0]->{"interval"}->{"value"}->{"id"}}->{basename($filename)};
                        @distr = values(%$retrhash);
                        $jdistr = join(" ",@distr);
                        if ($self->get("stats")) { # Calculate stats
                            $mean = $helper->mean(@distr);
                            $median = $helper->median(@distr);
                            $stdev = $helper->stdev(@distr);
                            $mad = $helper->mad(@distr);
                        }
                    }
                    $finalout .= "\t".$gencounts->{$chr}->{$_[0]->{"interval"}->{"value"}->{"id"}}->{basename($filename)};
                    $finalout .= "\t".$jdistr if ($self->get("split") || $self->get("nbins"));
                    $finalout .= "\t".$mean."\t".$median."\t".$stdev."\t".$mad if ($self->get("stats"));
                }
                print $out $finalout,"\n";
            });
        }
    }
    else
    {
        while(($chr,$tree) = each(%$chromosome))
        {
            $tree->traverse(sub {
                if ($_[0]->{"interval"}->{"value"}->{"rest"}) {
                    $finalout = $chr."\t".
                        $_[0]->{"interval"}->{"start"}."\t".
                        $_[0]->{"interval"}->{"end"}."\t".
                        $_[0]->{"interval"}->{"value"}->{"id"}."\t".
                        $_[0]->{"interval"}->{"value"}->{"rest"};
                }
                else {
                    $finalout = $chr."\t".
                        $_[0]->{"interval"}->{"start"}."\t".
                        $_[0]->{"interval"}->{"end"}."\t".
                        $_[0]->{"interval"}->{"value"}->{"id"};
                }
                foreach $filename (@originfile) {
                    if ($self->get("split") || $self->get("nbins")) { # Retrieve calculated distributions per region 
                        $retrhash = $splitcounts->{basename($filename)}->{$chr}->{$_[0]->{"interval"}->{"value"}->{"id"}};
                        @distr = values(%$retrhash);
                        $jdistr = join(" ",@distr);
                        if ($self->get("stats")) { # Calculate stats
                            $mean = $helper->mean(@distr);
                            $median = $helper->median(@distr);
                            $stdev = $helper->stdev(@distr);
                            $mad = $helper->mad(@distr);
                        }
                    }
                    $finalout .= "\t".$gencounts->{basename($filename)}->{$chr}->{$_[0]->{"interval"}->{"value"}->{"id"}};
                    $finalout .= "\t".$jdistr if ($self->get("split") || $self->get("nbins"));
                    $finalout .= "\t".$mean."\t".$median."\t".$stdev."\t".$mad if ($self->get("stats"));
                }
                print $out $finalout,"\n";
            });
        }
    }
    close(OUTPUT) if ($self->get("output"));
}

=head2 sort_inputs

Helper sorting function for input files and region. Internal use.

    $counter->sort_inputs($region_file,@input_files);
    
=cut

sub sort_inputs
{
    my ($self,$regionfile,@infile) = @_;
    my $tmpdir = $self->get("tmpdir");
    my $tmpfile;
    
    if ($^O !~ /MSWin/) # Case of linux, easy sorting
    {
        for (my $i=0; $i<@infile; $i++)
        {
            $helper->disp("Sorting bed file $infile[$i]...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp"."$i".".in$$");
            `sort -k1,1 -k2g,2 $infile[$i] > $tmpfile `;
            $infile[$i] = $tmpfile;
        }
        $helper->disp("Sorting region file $regionfile...");
        $tmpfile = File::Spec->catfile($tmpdir,"tempreg.in$$");
        `sort -k1,1 -k2g,2 $regionfile > $tmpfile `;
        $regionfile = $tmpfile;
    }
    else # We are in Windows... package required
    {
        $helper->try_module("File::Sort","sort_file");
        eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
        for (my $i=0; $i<@infile; $i++)
        {
            $helper->disp("Sorting file $infile[$i]...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp"."$i".".tmp");
            sort_file(
            {
                I => $infile[$i],
                o => $tmpfile,
                k => ['1,1','2n,2'],
                t => "\t"
            });
            $infile[$i] = $tmpfile;
        }
        $helper->disp("Sorting region file $regionfile...");
        $tmpfile = File::Spec->catfile($tmpdir,"tempreg.tmp");
        sort_file(
        {
            I => $regionfile,
            o => $tmpfile,
            k => ['1,1','2n,2'],
            t => "\t"
        });
        $regionfile = $tmpfile;
    }

    return($regionfile,@infile);
}

=head2 sort_one

Helper sorting function. Internal use.

    $counter->sort_one($input_file);
    
=cut

sub sort_one
{
    my ($self,$file) = @_;
    my $tmpdir = $self->get("tmpdir");
    my $tmpfile;
    
    if ($^O !~ /MSWin/) # Case of linux, easy sorting
    {
        $helper->disp("Sorting file $file...");
        $tmpfile = File::Spec->catfile($tmpdir,"temptss.in$$");
        `sort -k1,1 -k2g,2 $file > $tmpfile `;
        $file = $tmpfile;
    }
    else # We are in Windows... package required
    {
        $helper->try_module("File::Sort","sort_file");
        eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
        $helper->disp("Sorting file $file...");
        $tmpfile = File::Spec->catfile($tmpdir,"temptss.tmp");
        sort_file(
        {
            I => $file,
            o => $tmpfile,
            k => ['1,1','2n,2'],
            t => "\t"
        });
        $file = $tmpfile;
    }

    return($file);
}

=head2 change_params

Massively change the parameters of an HTS::Tools::Count object.

    $counter->change_params({'input' => 'another_file','region' => 'mouse-exon'})
    
=cut

sub change_params
{
    my ($self,$params) = @_;
    
    # Validate the new parameters 
    my $checker = HTS::Tools::Paramcheck->new();
    $checker->set("tool","count");
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

HTS::Tools::Count object getter

    my $param_value = $counter->get("param_name")
=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Count object setter

    $counter->set("param_name","param_value")
    
=cut

sub set
{
    my ($self,$name,$value) = @_;
    $self->{$name} = $value;
    return($self);
}

=head1 DEPENDENCIES

Tie::IxHash::Easy (optional)
File::Sort (optional)

= head1 TODO

Check how we can use gtf2tree.pl from the ngsplot package (https://code.google.com/p/ngsplot/)

=head1 AUTHOR

Panagiotis Moulos, C<< <moulos at fleming.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc HTS::Tools::Count


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

1; # End of HTS::Tools::Count

############################### LEGACY SUBROUTINES, USING BINARY SEARCH ############################

#sub read_region_file
#{
    #my $self = shift @_;
    #my $regionfile = shift @_;
    #my @originfile = @{shift @_};
    #@originfile = @{$self->get("input")} unless(@originfile);
    
    #my ($theHeader,$chr,$start,$end,$uid,$score,$strand,$rest);
    #my @arest;

    #my (%chromosome,%gencounts,%splitcounts);
    #tie %chromosome, "Tie::IxHash::Easy";
    #tie %gencounts, "Tie::IxHash::Easy";

    #my ($f,$regstart,$regend,$regconts,$regcount,$regline);
    #my @tempconts;

    #$helper->disp("Reading genomic regions file...");
    #$helper->disp("...also splitting genomic regions in areas of ".$self->get("split")." base pairs...") if ($self->get("split"));
    
    #open(REGIONS,$regionfile) or croak "\nThe file $regionfile does not exist!\n";
    #$regline = <REGIONS>;
    #$theHeader = $helper->decide_header($regline);
    #seek(REGIONS,0,0) if (!$theHeader);
    #while ($regline = <REGIONS>)
    #{
        #$regline =~ s/\r|\n$//g; # Make sure to remove carriage returns
        
        #($chr,$start,$end,$uid,@arest) = split(/\t/,$regline);
        #$rest = join("\t",@arest);
        #($rest eq "") ? ($chromosome{$chr}{$uid} = $start."\t".$end) :
            #($chromosome{$chr}{$uid} = $start."\t".$end."\t".$rest);
        
        #foreach $f (@originfile)
        #{
            #($self->get("ncore") > 1) ? ($gencounts{basename($f)}{$chr}{$uid} = 0) :
                #($gencounts{$chr}{$uid}{basename($f)} = 0);

            #if ($self->get("split"))
            #{
                #$regcount = 0;
                #$regend = $start - 1;
                #while ($regend < $end)
                #{
                    #$regstart = $regend + 1;
                    #$regend = $regstart + $self->get("split") - 1;
                    #push(@tempconts,$regstart." ".$regend."\t");
                    #($self->get("ncore") > 1) ? ($splitcounts{basename($f)}{$chr}{$uid}{$regcount} = 0) :
                        #($splitcounts{$chr}{$uid}{basename($f)}{$regcount} = 0);
                    #$regcount++;
                #}
            #}
        #}
    #}
    #close(REGIONS);

    #return(\%chromosome,\%gencounts,\%splitcounts,$theHeader);
#}

#sub split_area
#{
    #my ($self,$start,$end,$splitlen) = @_;
    
    #my $regstart;
    #my $regend = $start - 1;
    #my @subareas;
    
    #while ($regend < $end)
    #{
        #$regstart = $regend + 1;
        #$regend = $regstart + $splitlen - 1;
        #push(@subareas,$regstart."\t".$regend);
    #}
    #if ($regend > $end)
    #{
        #pop(@subareas);
        #$regend = $end;
        #push(@subareas,$regstart."\t".$regend);
    #}
    
    #return @subareas;
#}
#
#sub count_reads
#{
    #my ($self,$chromosome,$gencounts,$splitcounts,$infile,$originfile) = @_;
    #$originfile = $infile unless($originfile);

    #my (@k,@v);
    #my ($gene,$coords,$diff,$score,$p,$filename,$bprog,$numberlines);
    #my ($tmp1,$tmp2);
    #my ($bedchr,$bedstart,$bedend,$bedrest);
    #my ($currhash,$currchr,$currstart,$currend,$currrest);
    #my (@regarr,@bsr);
    #my $c = $self->get("constant");
    
    #$helper->disp("Reading and processing file $originfile...");
    
    #$filename = basename($originfile);
    #open(BEDFILE,$infile);
    #my $nextchr = "";

    #while (<BEDFILE>)
    #{
        #chomp $_;
        #($bedchr,$bedstart,$bedend,$bedrest) = split(/\t/,$_);

        #if ($bedchr ne $nextchr)
        #{
            #$currhash = $chromosome->{$bedchr};
            #@k = keys(%$currhash);
            #@v = values(%$currhash);
            #$nextchr = $bedchr;
        #}

        ## Start binary search
        #my $i;
        #my ($l,$u) = (0,$#k);
        #while ($l <= $u)
        #{
            #$i = int(($l + $u)/2);
            #my $gene = $k[$i];
            #my $coords = $v[$i];
            #($currstart,$currend,$currrest) = split(/\t/,$coords);

            ## Conditions that satisfy our search
            ## Security, too small region, tag does not fit?
            ## Mark that there is something there if user wishes and proceed
            #if ($currstart >= $bedstart && $currend <= $bedend)
            #{
                #$gencounts->{$bedchr}->{$gene}->{$filename}++ if ($self->get("small"));
                #if ($self->get("split") && $self->get("small")) # If splitting of regions
                #{
                    #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                    #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                    #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                #}
                #last;
            #}
            ## Ideal case... whole tag inside gene... no problem
            #elsif ($currstart <= $bedstart && $currend >= $bedend)
            #{
                #$gencounts->{$bedchr}->{$gene}->{$filename}++;
                #if ($self->get("split")) # If splitting of regions
                #{
                    #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                    #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                    #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                #}
                #last;
            #}
            ## Part of tag outside after region
            #elsif ($currstart < $bedstart && $currend < $bedend && $bedstart < $currend)
            #{
                #$diff = $currend - $bedstart;
                #if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                #{
                    #if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside no need for further check
                    #{
                        #$gencounts->{$bedchr}->{$gene}->{$filename}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                    #else # Use probabilistic scoring scheme
                    #{
                        #$score = ($currend - $bedstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                        #$score = exp(-$c**2/($currend - $bedstart)) if ($self->get("escore")); # Exponential
                        #$p = rand();
                        #if ($p < $score)
                        #{
                            #$gencounts->{$bedchr}->{$gene}->{$filename}++; 
                            #if ($self->get("split")) # If splitting of regions
                            #{
                                #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                                #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                                #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                            #}
                        #}
                    #}
                #}
                #else # Simple overlap
                #{
                    #if ($diff >= $self->get("percent")*($bedend - $bedstart))
                    #{
                        #$gencounts->{$bedchr}->{$gene}->{$filename}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                #}
                #last;
            #}
            ## Part of tag outside before region
            #elsif ($currstart > $bedstart && $currend > $bedend && $bedend > $currstart)
            #{
                #$diff = $bedend - $currstart;
                #if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                #{
                    #if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside
                    #{
                        #$gencounts->{$bedchr}->{$gene}->{$filename}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                    #else # Use probabilistic scoring scheme
                    #{
                        #$score = ($bedend - $currstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                        #$score = exp(-$c**2/($bedend - $currstart)) if ($self->get("escore")); # Exponential
                        #$p = rand();
                        #if ($p < $score)
                        #{
                            #$gencounts->{$bedchr}->{$gene}->{$filename}++;
                            #if ($self->get("split")) # If splitting of regions
                            #{
                                #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                                #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                                #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                            #}
                        #}
                    #}
                #}
                #else # Simple overlap
                #{
                    #if ($diff >= $self->get("percent")*($bedend - $bedstart))
                    #{
                        #$gencounts->{$bedchr}->{$gene}->{$filename}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                #}
                #last;
            #}
            #else # If none of them met, reduce the searching subset
            #{
                #$u = $i - 1 if ($bedend <= $currstart);
                #$l = $i + 1 if ($bedstart >= $currend);
            #}

        #}
    #}
    #close(BEDFILE);

    #return($gencounts,$splitcounts);
#}
#
#sub count_reads_multi
#{
    #my ($self,$chromosome,$gencounts,$splitcounts,$infile,$originfile) = @_;
    #$originfile = $infile unless($originfile);

    #my (@k,@v);
    #my ($gene,$coords,$diff,$score,$p,$filename,$bprog,$numberlines);
    #my ($bedchr,$bedstart,$bedend,$bedrest);
    #my ($currhash,$currchr,$currstart,$currend,$currrest);
    #my (@regarr,@bsr);
    #my $c = $self->get("constant");
    
    #$helper->disp("Reading and processing file $originfile...");

    #$filename = basename($originfile);
    #open(BEDFILE,$infile);
    #my $nextchr = "";

    #while (<BEDFILE>)
    #{
        #chomp $_;
        #($bedchr,$bedstart,$bedend,$bedrest) = split(/\t/,$_);
        
        #if ($bedchr ne $nextchr)
        #{
            #$currhash = $chromosome->{$bedchr};
            #@k = keys(%$currhash);
            #@v = values(%$currhash);
            #$nextchr = $bedchr;
        #}

        ## Start binary search
        #my $i;
        #my ($l,$u) = (0,$#k);
        #while ($l <= $u)
        #{
            #$i = int(($l + $u)/2);
            #my $gene = $k[$i];
            #my $coords = $v[$i];
            #($currstart,$currend,$currrest) = split(/\t/,$coords);

            ## Conditions that satisfy our search
            ## Security, too small region, tag does not fit?
            ## Mark that there is something there if user wishes and proceed
            #if ($currstart >= $bedstart && $currend <= $bedend)
            #{
                #$gencounts->{$bedchr}->{$gene}++ if ($self->get("small"));
                #if ($self->get("split") && $self->get("small")) # If splitting of regions
                #{
                    #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                    #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                    #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                #}
                #last;
            #}
            ## Ideal case... whole tag inside gene... no problem
            #elsif ($currstart <= $bedstart && $currend >= $bedend)
            #{
                #$gencounts->{$bedchr}->{$gene}++;
                #if ($self->get("split")) # If splitting of regions
                #{
                    #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                    #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                    #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                #}
                #last;
            #}
            ## Part of tag outside after region
            #elsif ($currstart < $bedstart && $currend < $bedend && $bedstart < $currend)
            #{
                #$diff = $currend - $bedstart;
                #if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                #{
                    #if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside no need for further check
                    #{
                        #$gencounts->{$bedchr}->{$gene}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                    #else # Use probabilistic scoring scheme
                    #{
                        #$score = ($currend - $bedstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                        #$score = exp(-$c**2/($currend - $bedstart)) if ($self->get("escore")); # Exponential
                        #$p = rand();
                        #if ($p < $score)
                        #{
                            #$gencounts->{$bedchr}->{$gene}++; 
                            #if ($self->get("split")) # If splitting of regions
                            #{
                                #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                                #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                                #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                            #}
                        #}
                    #}
                #}
                #else # Simple overlap
                #{
                    #if ($diff >= $self->get("percent")*($bedend - $bedstart))
                    #{
                        #$gencounts->{$bedchr}->{$gene}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                #}
                #last;
            #}
            ## Part of tag outside before region
            #elsif ($currstart > $bedstart && $currend > $bedend && $bedend > $currstart)
            #{
                #$diff = $bedend - $currstart;
                #if ($self->get("lscore") || $self->get("escore")) # Linear or exponential scoring
                #{
                    #if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside
                    #{
                        #$gencounts->{$bedchr}->{$gene}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                    #else # Use probabilistic scoring scheme
                    #{
                        #$score = ($bedend - $currstart)/(($bedend - $bedstart)/2) if ($self->get("lscore")); # Linear
                        #$score = exp(-$c**2/($bedend - $currstart)) if ($self->get("escore")); # Exponential
                        #$p = rand();
                        #if ($p < $score)
                        #{
                            #$gencounts->{$bedchr}->{$gene}++;
                            #if ($self->get("split")) # If splitting of regions
                            #{
                                #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                                #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                                #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                            #}
                        #}
                    #}
                #}
                #else # Simple overlap
                #{
                    #if ($diff >= $self->get("percent")*($bedend - $bedstart))
                    #{
                        #$gencounts->{$bedchr}->{$gene}++;
                        #if ($self->get("split")) # If splitting of regions
                        #{
                            #@regarr = $self->split_area($currstart,$currend,$self->get("split"));
                            #@bsr = $self->bin_search_loc($bedstart,$bedend,@regarr);
                            #$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
                        #}
                    #}
                #}
                #last;
            #}
            #else # If none of them met, reduce the searching subset
            #{
                #$u = $i - 1 if ($bedend <= $currstart);
                #$l = $i + 1 if ($bedstart >= $currend);
            #}
        #}
    #}
    #close(BEDFILE);

    #return($gencounts,$splitcounts);
#}
#
#sub write_reads
#{
    #my $self = shift @_;
    #my $chromosome = shift @_;
    #my $gencounts = shift @_;
    #my $splitcounts = shift @_;
    #my @originfile = @{shift @_};
    #my $theHeader = shift @_;
    #@originfile = @{$self->get("input")} unless(@originfile);
    
    #my (@outrest,@distr,@filenames,@outchrs,@outgeneids);
    #my ($outgene,$files,$filename,$b,$retrhash,$outhead,$outsubhead);
    #my ($outchr,$outgeneid,$outcounts,$outstart,$outend,$outrest,$finalout);
    #my ($mean,$median,$stdev,$mad,$jdistr);
    
    #my @base;
    #foreach my $of (@originfile)
    #{
        #my $bb = fileparse($of,'\.[^.]*');
        #push(@base,$bb);
    #}

    #my $out;
    #if ($self->get("output"))
    #{
        #open(OUTPUT,">".$self->get("output"));
        #$out = *OUTPUT;
    #}
    #else { $out = *STDOUT; }

    #if ($self->get("output"))
    #{
        #($self->get("stats")) ? ($helper->disp("Writing reads and calculating statistics per genomic regions in ".$self->get("output")."...")) :
            #$helper->disp("Writing reads per genomic regions for file in ".$self->get("output")."...");
    #}
    #else
    #{
        #($self->get("stats")) ? ($helper->disp("Writing reads and calculating statistics per genomic regions to standard output...\n")) :
            #$helper->disp("Writing reads per genomic regions for file to standard output...\n");
    #}

    #if ($theHeader) # Add a header based on what we did plus additional data given
    #{
        #$outhead = $theHeader;
        #if (scalar @base == 1)
        #{
            #$outhead .= "\t".$base[0];
            #$outhead .= "\tSub-area Counts/$self->get(\"split\") bp" if ($self->get("split")); # Return also sub-area distributions
            #$outhead .= "\tMean Counts/$self->get(\"split\") bp\tMedian Counts/$self->get(\"split\") bp\tStDev Counts\tMAD Counts" if ($self->get("stats"));
            #print $out "$outhead\n";
        #}
        #else
        #{
            #foreach my $b (@base)
            #{
                #$outhead .= "\t".$b;
                #$outhead .= "\t" if ($self->get("split"));
                #$outhead .= "\t"x4 if ($self->get("stats"));
            #}
            ## We must also print a subheader in case of reporting stats
            #$outsubhead = "Counts";
            #$outsubhead .= "\tSub-area Counts/$self->get(\"split\") bp" if ($self->get("split"));
            #$outsubhead .= "\tMean Counts/$self->get(\"split\") bp\tMedian Counts/$self->get(\"split\") bp\tStDev Counts\tMAD Counts" if ($self->get("stats"));
            #$outsubhead = "$outsubhead"x(scalar @base);
            #$outsubhead = ("\t"x(scalar split(/\t/,$theHeader))).$outsubhead;
            #print $out "$outhead\n";
            #print $out "$outsubhead\n" if ($self->get("split"));
        #}
    #}

    #if ($self->get("ncore") == 1)
    #{
        #while(($outchr,$outgene) = each(%$gencounts))
        #{
            #while (($outgeneid,$files) = each(%$outgene))
            #{  
                #($outstart,$outend,@outrest) = split(/\t/,$chromosome->{$outchr}->{$outgeneid});
                #if (scalar @outrest != 0)
                #{
                    #$outrest = join("\t",@outrest);
                    #$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid."\t".$outrest;
                #}
                #else
                #{
                    #$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid;
                #}

                #while (($filename,$outcounts) = each(%$files))
                #{
                    #if ($self->get("split")) # Retrieve calculated distributions per region
                    #{
                        #$retrhash = $splitcounts->{$outchr}->{$outgeneid}->{basename($filename)};
                        #@distr = values(%$retrhash);
                        #$jdistr = join(" ",@distr);
                        #if ($self->get("stats")) # Calculate stats
                        #{
                            #$mean = $helper->mean(@distr);
                            #$median = $helper->median(@distr);
                            #$stdev = $helper->stdev(@distr);
                            #$mad = $helper->mad(@distr);
                        #}
                    #}
                    #$finalout .= "\t".$outcounts;
                    #$finalout .= "\t".$jdistr if ($self->get("split"));
                    #$finalout .= "\t".$mean."\t".$median."\t".$stdev."\t".$mad if ($self->get("stats"));
                #}
                #print $out "$finalout\n";
            #}
        #}
    #}
    #else
    #{
        #@outchrs = keys(%$chromosome);
        #@filenames = keys(%$gencounts);
        #foreach my $outchr (@outchrs)
        #{
            #@outgeneids = keys(%{$chromosome->{$outchr}});
            #foreach $outgeneid (@outgeneids)
            #{
                #($outstart,$outend,@outrest) = split(/\t/,$chromosome->{$outchr}->{$outgeneid});
                #if (scalar @outrest != 0)
                #{
                    #$outrest = join("\t",@outrest);
                    #$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid."\t".$outrest;
                #}
                #else
                #{
                    #$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid;
                #}

                #foreach $filename (@filenames)
                #{
                    #$outcounts = $gencounts->{$filename}->{$outchr}->{$outgeneid};
                    #if ($self->get("split")) # Retrieve calculated distributions per region
                    #{
                        #$retrhash = $splitcounts->{basename($filename)}->{$outchr}->{$outgeneid};
                        #@distr = values(%$retrhash);
                        #$jdistr = join(" ",@distr);
                        #if ($self->get("stats")) # Calculate stats
                        #{
                            #$mean = $helper->mean(@distr);
                            #$median = $helper->median(@distr);
                            #$stdev = $helper->stdev(@distr);
                            #$mad = $helper->mad(@distr);
                        #}
                    #}
                    #$finalout .= "\t".$outcounts;
                    #$finalout .= "\t".$jdistr if ($self->get("split"));
                    #$finalout .= "\t".$mean."\t".$median."\t".$stdev."\t".$mad if ($self->get("stats"));
                #}
                #print $out "$finalout\n";
            #}
        #}
    #}
    #close(OUTPUT) if ($self->get("output"));
#}
#
#=head2 bin_search_loc

#Binary search for genomic regions. Internal use.

    #$counter->bin_search_loc($start,$end,@array_of_coordinates);
    
#=cut

#sub bin_search_loc
#{
    #my ($self,$start,$end,@areas) = @_;
    
    #my ($ind,$currstart,$currend,$center);
    #my ($l,$u) = (0,$#areas);
    #while ($l <= $u)
    #{
        #$ind = int(($l + $u)/2);
        #($currstart,$currend) = split(/\t/,$areas[$ind]);
        #$center = $start + $helper->round(($end - $start)/2); # Location of center of the tag
        #if ($currstart < $center && $currend > $center)
        #{
            #return (1,$ind);
        #}
        #elsif ($currstart == $center) # Randomly assign to one window
        #{
            #(int(10*rand(1)) < 5) ? (return(1,$ind-1)) : (return(1,$ind));
        #}
        #elsif ($currend == $center)
        #{
            #(int(10*rand(1)) < 5) ? (return(1,$ind)) : (return(1,$ind+1));
        #}
        #else
        #{
            #$u = $ind - 1 if ($center <= $currstart);
            #$l = $ind + 1 if ($center >= $currend);
        #}
    #}
    #return (0,-1);
#}
