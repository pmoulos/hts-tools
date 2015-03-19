=head1 NAME

HTS::Tools::Motifscan - A massive motif scanner (not finder!) for short genomic 
regions

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is a wrapper and significance threshold optimizer for several 
sequence motif scanners. At the present point, the scanners supported are 
pwmscan from the GimmeMotifs suite (van Heeringen and Veenstra, 2010) and 
MotifScanner (Thijs et al., 2001). Thus, the installation of these 3rd party 
tools in the system that uses HTS::Tools::Motifscan is necessary (or at least of
one of them). Then, the module runs a workflow which which massively scan the 
input fasta files (e.g. a set of sequences corresponding to ChIP-Seq peaks) for
the motifs given in the motif file in the form of positional weight matrices 
(e.g. motifs that have been characterized as significant during a de novo motif 
search in the same set of peaks). This workflow determines an optimal 
significance threshold based on a set of random background sequences selected 
from the background fasta file provided and on the false postive rate (fpr) 
provided. The module provides three basic outputs: a simple file with the number 
of hits of each motif in each input file, gff files with the hits and bed files 
with the hits, including also a significance score of the scan hit in the 5th 
column of the bed file.

    use HTS::Tools::Motifscan;
    my %params = (
        'input' => ['normal_nfkb_peaks.fa','cancer_nfkb_peaks.fa',
            'mock_peaks.fa']
        'motif' => 'my_motif_matrices.pwm',
        'scanner' => 'pwmscan',
        'background' => 'background_sequences.tab',
        'range' => 0.1:0.1:1,
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

A set of input FASTA file(s) that will be used for matching of the motifs 
contained in I<motif> file. Each of these files will be scanned using the 
scanner defined by I<scanner> and the hits will be reported either as a simple 
tab-delimited stats file, a gff file, or a bed file under conditions (see below).

=item I<motif> B<(required)>

A file containing several motifs in the Position Weight Matrix (PWM) format. The
file will be decomposed in several motifs and each motif will be matched against
each input FASTA file.

=item I<background> B<(optional)>

A FASTA file containing background sequences (must be at least the same length
as the input sequences) that will be used for random sampling to define a score
cutoff for motif matching significance based on I<fpr>. It is required if the 
option I<justscan> is not activated.

=item I<scanner> B<(optional)>

A scanning algorithm to use. Currently, two algorithms are supported: the 
pwmscan algorithm (van Heeringen and Veenstra, 2010) and MotifScanner (Thijs et 
al., 2001).

=item I<sigmethod> B<(optional)>

It can be "bootstrap", "converge" or "none". If "bootstrap", the input sequence 
set is scanned using the I<justscan> threshold and the hits are counted. Then, 
using the same  threshold, I<times> background sequences of the same length are 
scanned and the number of hits is counted. The p-value is the number of times 
that the number of hits is larger in the background than the input sequence set 
divided by I<times>. When "converge", the input sequences are scanned with 
I<range> thresholds until less hits remain in I<times> background  sequences 
until I<fpr> is reached. If "none", then just a scan is performed with the 
I<justscan> threshold without any statistical significance. The default is 
"bootstrap".

=item I<ncore> B<(optional)>

If the machine has multicore processor(s) and the package Parallel::Iterator is 
installed, you can use parallel processing when I<sigmethod> is "bootstrap". 
Default is 1 and can go up to 12.

=item I<range> B<(optional)>

A range of cutoffs that will be used to determine the final matching score 
cutoff corresponding to I<fpr>. The range can be given in the form a:b or a:x:b 
where x is an increment step. If the format a:b is chosen, the default increment
is 1. However, in the latest versions of pwmscan, the matching score is 
normalized to one, so this notation will be eventually deprecated.

=item I<fpr> B<(optional)>

The desired False Positive Rate (FPR), that is the percentage of motif matches 
found in the background sequences. Defaults to 0.05.

=item  I<times> B<(optional)>

How many times should the background set of sequences that will be used for the 
determination of the cutoff score from I<range>, be larger than the input set? 
It defaulta to 10, which means that if an input FASTA file contains 100 sequences, 
the background will contain 1000 random sequences.

=item I<length> B<(optional)>

The length of the background sequences. Generally, it should be larger than the 
length of input sequences. It defaults to 400.

=item I<output> B<(optional)>

The output types that the user wished to get. It can be one or more of "stats" 
for a simple delimited file containing the hits for each motif and input file, 
"gff" for the gff output of pwmscan or a related file from MotifScanner, 
containing the actual hits and positions in the input sequences, or "bed" for an
output BED file containing the motif matches locations and a score to be used 
for coloring. The "bed" output is available only if a set of peak files is given
so as to determine the relative location of the match inside the peak and 
construct proper bed lines.

=item  I<besthit> B<(optional)>

The number of best hits to be retrieved when I<scanner> is "pwmscan". Defaults 
to 1.

=item I<uniquestats> B<(optional)>

If the number of besthits is greater than 1, the specifying I<uniquestats> to 1 
will cause the output "stats" file to contain unique hits. Like this you can 
avoid paradoxes like having more hits than input FASTA sequences. However, in 
some situations you might actually want to retrieve multiple hits.

=item I<justscan> B<(optional)>

Set this to a real number between 0 and 1 (preferably close to 0.9), to just 
scan the input sequences using this score cutoff for the motifs without defining 
an FPR. Useful if you have a very limited number of input sequences (e.g. just 
one).

=item I<center> B<(optional)>

A set of genomic regions (e.g. peaks) with the SAME IDs as the input FASTA 
sequences (or a superset of these) which will be used to assign peak regions to 
FASTA files in order to determine the proper coordinates for the generation of 
BED output. It is optional, however, if BED output is requested, these files are 
not given, and the FASTA sequence length does not correspond to the length that 
can be extracted by the FASTA ID, the coordinates will be inaccurate.

=item I<colext> B<(optional)>

A vector of length 3, containing the column numbers of peak ID and peak summit,
and the length of the (possible) extension upstream and downstream of the peak 
summit. For example I<colext> 4 5 75.

=item I<silent> B<(optional)>

Use this parameter if you want to turn informative messages off.

=head1 OUTPUT

The output of the module is a set of GFF files containing the hits for each 
motif and each input set of sequences, a statistics file with the number of hits
for each motif and input set of sequences and a set of BED files with the motif 
matches that can be used for display in a genome browser.

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

use HTS::Tools::Constants;
use HTS::Tools::Convert;
use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper $const);

our $MODNAME = "HTS::Tools::Motifscan";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Motif scan in a set of sequences using 3rd party tools.";

BEGIN {
    $helper = HTS::Tools::Utils->new();
    $const = HTS::Tools::Constants->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

The HTS::Tools::Motifscan object constructor. It accepts a set of parameters 
that are required to run the motifscanner and get the output.

    my $motifscanner = HTS::Tools::Motifscan->new({'input' => ['peaks_1.fa',
        'peaks_2.fa'],'motif' => 'my_motifs.pwm','scanner' => 'pwmscan',
        'justscan' => 1});

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
    $checker->set("tool","motifscan");
    $checker->set("params",$params);
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init

HTS::Tools::Motifscan object initialization method. NEVER use this directly, use 
new instead.

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

The HTS::Tools::Motifscan run subroutine. It runs the motifscanner with the 
given parameters in the constructor.

    $motifscanner->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    # Copy some memory-less variables to avoid rewriting the whole thing...
    our @seqfile = @{$self->get("input")};
    our @output = @{$self->get("output")};
    our @cntfile = @{$self->get("center")};
    our @cntcol = @{$self->get("colext")};
    our @range = @{$self->get("range")};
    our $motfile = $self->get("motif"); 
    our $backfile = $self->get("background");
    our $scanner = $self->get("scanner");
    our $sigmethod = $self->get("sigmethod");
    our $ncore = $self->get("ncore");
    our $fpr = $self->get("fpr");
    our $times = $self->get("times");
    our $length = $self->get("length");
    our $besthit = $self->get("besthit");
    our $unistats = $self->get("uniquestats");
    our $justscan = $self->get("justscan");
    
    my ($ogff,$obed,$ostat);
    our $olog;
    foreach my $o (@output)
    {
        $ogff = 1 if ($o eq "gff");
        $obed = 1 if ($o eq "bed");
        $ostat = 1 if ($o eq "stats");
        $olog = 1 if ($o eq "log");
    }

    # Start logging if chosen
    our @log;
    push(@log,"##########     HTS::Tools::Motifscan logfile     ##########\n") if ($olog);

    # Bavard
    my $date = $helper->now;
    $helper->disp("$date - Started...\n");

    # More bavard... I like information...
    $helper->disp("Selected scanner : $scanner");
    $helper->disp("Enrichment significance assessment : $sigmethod");
    if ($sigmethod eq "converge")
    {
        $helper->disp("Background sequences files larger by factor of : $times");
    }
    elsif ($sigmethod eq "bootstrap")
    {
        $helper->disp("Bootstrap iterations with background : $times");
    }
    ($justscan) ? ($helper->disp("Cutoff range : $justscan")) : 
    ($helper->disp("Cutoff range : ",$range[0],"to",$range[$#range]));
    $helper->disp("False Positive Rate : $fpr") if ($fpr);
    $helper->disp("Chosen output(s) : ",join("\, ",@output),"\n");

    # If bed output chosen and cntfile given, read them and merge them to database
    my %cnthash;
    if ($obed && @cntfile)
    {
        $helper->disp("Merging files containing peak centers...");
        foreach my $c (@cntfile)
        {
            my ($line,@columns);
            open(CNT,"$c") or croak "\nCannot open peak center file $c.\n";
            $line = <CNT>;
            my $header = $helper->decide_header($line);
            seek(CNT,0,0) if (!$header);
            while ($line = <CNT>)
            {
                @columns = split(/\t/,$line);
                $cnthash{$columns[$cntcol[0]]} = $columns[$cntcol[1]];
            }
            close(CNT);
        }
    }
    elsif ($obed && !(@cntfile || @cntcol))
    {
        $helper->disp("WARNING! BED output requested but no files containing peak centers are given...");
        $helper->disp("Proceeding only with co-ordinates as will can be derived from sequence IDs...");
        $helper->disp("Assuming sequence ID format : chr:start-end or chr:start:end...\n");
    }

    my ($i,$j,$k); # General indices
    my (@stats,@pvalues,@means,@sds); # Stats arrays in case

    # Parse the motifs file, split it into different files
    my @mfiles = $self->parse_motifs($motfile,$scanner);

    # If running the full process, not just scanning given only one threshold
    if (!$justscan && $sigmethod eq "converge")
    {
        # Quick and dirty checking about tabular file
        my $chk = $helper->check_tabseq($backfile);
        my $backbackup = $backfile;
        if ($chk)
        {
            $helper->disp("Background file $backfile does not appear to be a tabular sequence file...");
            $helper->disp("Checking if $backfile is in FASTA format...");
            my $chkchk = $helper->check_fasta($backfile);
            croak "\nBackground file $backfile does not appear to be in FASTA format either! Exiting...\n" if ($chkchk);
            $helper->disp("$backfile is in FASTA format. Converting to tabular format...");
            $backbackup = $backfile;
            my $converter = HTS::Tools::Convert->new();
            $backfile = $converter->fasta2tab($backfile);
        }

        # Count number of background sequences
        $helper->disp("Counting sequences in background file $backfile...");
        my $nback = $helper->count_lines($backfile); 
        $helper->disp("Background file $backfile contains $nback sequences.");

        # Construct background index file for quick reference
        my ($bbi,$bbd) = fileparse($backfile,'\.[^.]*');
        $bbi = $bbd.$bbi.".idx";
        if (! -e $bbi)
        {
            $helper->disp("Constructing index file for background file $backfile...");
            open(BACK,"< $backfile");
            open(INDEX,"+>$bbi") or croak "Can't open $bbi for read/write: $!\n";
            $self->build_index(*BACK,*INDEX);
            close(BACK);
            close(INDEX);
            $helper->disp("Index file for background constructed in $bbi.");
        }
        else
        {
            $helper->disp("Index file $bbi for background file $backfile already exists. Proceeding...\n");
        }
        
        if ($scanner =~ /MotifScanner/i)
        {
            $helper->disp("Constructing Markov background model to be used with MotifScanner...");
            $self->construct_MS_background($backbackup);
        }

        # Do job for peak-region-whatever files
        for ($i=0; $i<@seqfile; $i++)
        {
            $helper->disp("Processing $seqfile[$i]...");
            my $chkat1 = 0;
            my $chk = $helper->check_fasta($seqfile[$i]);
            if ($chk)
            {
                $helper->disp("File $seqfile[$i] does not appear to be a FASTA file. Proceeding to next...");
                next;
            }
            else { $chkat1 = 1; }
            
            # Open, count lines, call background etc.
            $helper->disp("Counting sequences in file $seqfile[$i]...");
            my $nseqs = $helper->count_fasta($seqfile[$i]);
            $helper->disp("File $seqfile[$i] has $nseqs sequences.");
            
            # Generate adjusted background for each input file
            my $n = $times*$nseqs;
            $helper->disp("Generating $n background sequences to be used for scanning of $seqfile[$i]...");
            my $bfile = $self->get_random_seq($backfile,$bbi,$nback,$length,$n);
            $helper->disp("$n background sequences written in file $bfile...\n");
            
            # For each motif, calculate threshold and scan according to selected scanner
            use v5.10;
            if ($scanner =~ m/pwmscan/i)
            {
                my ($sl,$cutoff);
                for ($j=0; $j<@mfiles; $j++)
                {
                    $helper->disp("Calculating score cutoff based on sequences in $bfile for motif $mfiles[$j] using pwmscan...");
                    $helper->disp("Output will be written in bgcut.gff.\n");
                    # Use range vector
                    for  ($k=0; $k<@range; $k++)
                    {
                        $helper->disp("Now testing threshold $range[$k]...");
                        `gimme scan $bfile $mfiles[$j] -c $range[$k] -n $besthit > bgcut.gff `;
                        if ($besthit > 1) 
                        {
                            ($unistats) ? ($sl = $helper->count_unique_lines("bgcut.gff")) : ($sl = $helper->count_lines("bgcut.gff"));
                        }
                        else { $sl = $helper->count_lines("bgcut.gff"); }
                        if ($sl <= $fpr*$n)
                        {
                            $cutoff = $range[$k];
                            last;
                        }
                        $helper->disp("FPR $fpr not reached... $sl matches out of $n sequences remain in background.");
                    }
                    if ($k == @range)
                    {
                        $helper->disp("Last number of sequences remaining in background : $sl");
                        $helper->disp("No cutoff within given range satisfies FPR criteria... Skipping to next.\n");
                        $cutoff = 0;
                        $stats[$i][$j] = 0 if ($ostat);
                        next;
                    } 
                    else 
                    { 
                        $helper->disp("Cutoff for FPR $fpr determined at $cutoff.");
                        $helper->disp("$sl matches out of $n sequences remain in background.\n");
                    }
                    $helper->disp("Scanning $seqfile[$i] for motif $mfiles[$j] using pwmscan... Cutoff: $cutoff.");
                    my $mbase = fileparse($mfiles[$j],'\.[^.]*');
                    my $currout = $self->create_output_file($seqfile[$i],"output",$mbase);
                    #`python gimme scan -i $seqfile[$i] -c $cutoff -s $spacer[$j] -p $mfiles[$j] -n $besthit > $currout `;
                    `gimme scan $seqfile[$i] $mfiles[$j] -c $cutoff -n $besthit > $currout `;
                    $self->convert2bed($currout,$cntcol[2],%cnthash) if ($obed);
                    my $nmatch;
                    if ($besthit > 1) 
                    {
                        ($unistats) ? ($nmatch = $helper->count_unique_lines($currout)) : ($nmatch = $helper->count_lines($currout));
                    }
                    else { $nmatch = $helper->count_lines($currout); }
                    $stats[$i][$j] = $nmatch if ($ostat);
                    $helper->disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
                    unlink($currout) if (!$ogff);
                }
            }                
            elsif ($scanner =~ m/MotifScanner/i)
            {
                my ($sl,$cutoff);
                for ($j=0; $j<@mfiles; $j++)
                {
                    $helper->disp("Calculating score cutoff based on sequences in $bfile for motif $mfiles[$j] using MotifScanner...");
                    $helper->disp("Output will be written in bgcut.gff.\n");
                    # Use range vector
                    for  ($k=0; $k<@range; $k++)
                    {
                        $helper->disp("Now testing p-value threshold $range[$k]...\n");
                        ($^O !~ /MSWin/) ? (`./MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $range[$k] -s 0 -o bgcut.gff `) :
                        (`MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $range[$k] -s 0 -o bgcut.gff `);
                        ($unistats) ? ($sl = $helper->count_unique_lines("bgcut.gff")) : ($sl = $helper->count_lines("bgcut.gff"));
                        $sl--; # One line header of MotifScanner output
                        if ($sl <= $fpr*$n) 
                        {
                            $cutoff = $range[$k];
                            last;
                        }
                        $helper->disp("FPR $fpr not reached... $sl matches out of $n sequences remain in background.");
                    }
                    if ($k == @range)
                    {
                        $helper->disp ("Last number of sequences remaining in background : $sl");
                        $helper->disp ("No cutoff within given range satisfies FPR criteria... Skipping to next.\n");
                        $cutoff = 0;
                        $stats[$i][$j] = 0 if ($ostat);
                        next;
                    } 
                    else 
                    { 
                        $helper->disp("Cutoff for FPR $fpr determined at $cutoff.");
                        $helper->disp("$sl matches out of $n sequences remain in background.\n");
                    }
                    $helper->disp("Scanning $seqfile[$i] for motif $mfiles[$j]... Cutoff: $cutoff.");
                    my $mbase = fileparse($mfiles[$j],'\.[^.]*');
                    my $currout = $self->create_output_file($seqfile[$i],"output",$mbase);
                    ($^O !~ /MSWin/) ? (`./MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $cutoff -s 0 -o $currout `) :
                    (`MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $cutoff -s 0 -o $currout `);
                    $self->convert2bed($currout,$cntcol[2],%cnthash) if ($obed);
                    my $nmatch;
                    ($unistats) ? ($nmatch = $helper->count_unique_lines($currout)) : ($nmatch = $helper->count_lines($currout));
                    $nmatch--; # One line header of MotifScanner output
                    $stats[$i][$j] = $nmatch if ($ostat);
                    $helper->disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
                    unlink($currout) if (!$ogff);
                }
            }
            $helper->disp(" ");
        }
    }
    else # Just scan the sequence files using one defined threshold or bootstrap
    {
        my ($bbi,$bbd,$nback,$backbackup,$chk);
        if ($sigmethod eq "bootstrap")
        {
            # Quick and dirty checking about tabular file
            $chk = $helper->check_tabseq($backfile);
            $backbackup  = $backfile;
            if ($chk)
            {
                $helper->disp("Background file $backfile does not appear to be a tabular sequence file...");
                $helper->disp("Checking if $backfile is in FASTA format...");
                my $chkchk = $helper->check_fasta($backfile);
                croak "\nBackground file $backfile does not appear to be in FASTA format either! Exiting...\n" if ($chkchk);
                $helper->disp("$backfile is in FASTA format. Converting to tabular format...");
                $backbackup = $backfile;
                my $converter = HTS::Tools::Convert->new();
                $backfile = $converter->fasta2tab($backfile);
            }
            # Construct background index file for quick reference
            ($bbi,$bbd) = fileparse($backfile,'\.[^.]*');
            $bbi = $bbd.$bbi.".idx";
            if (! -e $bbi)
            {
                $helper->disp("Constructing index file for background file $backfile...");
                open(BACK,"< $backfile");
                open(INDEX,"+>$bbi") or croak "Can't open $bbi for read/write: $!\n";
                $self->build_index(*BACK,*INDEX);
                close(BACK);
                close(INDEX);
                $helper->disp("Index file for background constructed in $bbi.");
            }
            else
            {
                $helper->disp("Index file $bbi for background file $backfile already exists. Proceeding...\n");
            }
            $helper->disp("Counting sequences in background file $backfile...");
            $nback = $helper->count_lines($backfile); 
            $helper->disp("Background file $backfile contains $nback sequences.");
        }
        
        if ($scanner =~ /MotifScanner/i)
        {
            $helper->disp("Constructing Markov background model to be used with MotifScanner...");
            $self->construct_MS_background($backbackup);
        }
        
        # Do job for peak-region-whatever files
        for ($i=0; $i<@seqfile; $i++)
        {
            $helper->disp("Processing $seqfile[$i]...");
            my $chkat1 = 0;
            my $chk = $helper->check_fasta($seqfile[$i]);
            if ($chk)
            {
                $helper->disp("File $seqfile[$i] does not appear to be a FASTA file. Proceeding to next...");
                next;
            }
            else { $chkat1 = 1; }

            # For each motif, calculate threshold and scan according to selected scanner
            if ($scanner =~ m/pwmscan/i)
            {
                for ($j=0; $j<@mfiles; $j++)
                {
                    $helper->disp("Scanning $seqfile[$i] for motif $mfiles[$j] using pwmscan... Cutoff: $justscan.");
                    my $mbase = fileparse($mfiles[$j],'\.[^.]*');
                    my $currout = $self->create_output_file($seqfile[$i],"output",$mbase);
                    `gimme scan $seqfile[$i] $mfiles[$j] -c $justscan -n $besthit > $currout `;
                    $self->convert2bed($currout,$cntcol[2],%cnthash) if ($obed);
                    my $nmatch;
                    if ($besthit > 1) 
                    {
                        ($unistats) ? ($nmatch = $helper->count_unique_lines($currout)) : ($nmatch = $helper->count_lines($currout));
                    }
                    else { $nmatch = $helper->count_lines($currout); }
                    $stats[$i][$j] = $nmatch;
                    $helper->disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
                    unlink($currout) if (!$ogff);
                }
                if ($sigmethod eq "bootstrap")
                {
                    $helper->disp("Counting sequences in file $seqfile[$i]...");
                    my $nseqs = $helper->count_fasta($seqfile[$i]);
                    $helper->disp("File $seqfile[$i] has $nseqs sequences.");
                    $helper->disp("\nAssessing statistical significance for tne input motifs for $seqfile[$i] using resampling... It might take some time...\n");
                    for ($j=0; $j<@mfiles; $j++)
                    {
                        my (@bhits,@nums);
                        my $sl;
                        if ($ncore == 1)
                        {                           
                            for ($k=1; $k<=$times; $k++)
                            {
                                $helper->disp("  Motif $mfiles[$j] -- Iteration $k");
                                $helper->disp("    Generating $nseqs background sequences to be used for motif scanning...");
                                my $bfile = $self->get_random_seq($backfile,$bbi,$nback,$length,$nseqs);
                                $helper->disp("    $nseqs background sequences written in file $bfile...");
                                $helper->disp("    Now scanning file $bfile...");
                                `gimme scan $bfile $mfiles[$j] -c $justscan -n $besthit > bgcut.gff `;
                                if ($besthit > 1) 
                                {
                                    ($unistats) ? ($sl = $helper->count_unique_lines("bgcut.gff")) : ($sl = $helper->count_lines("bgcut.gff"));
                                }
                                else { $sl = $helper->count_lines("bgcut.gff"); }
                                $helper->disp("    Found $sl hits");
                                push(@bhits,$sl);
                                unlink($bfile);
                            }
                        }
                        else
                        {
                            @nums = (1..$times);
                            $helper->disp("Parallel boostraping enrichment for motif $mfiles[$j]...");
                            @bhits = Parallel::Iterator->iterate_as_array({workers => $ncore},sub{
                                my ($id,$job) = @_;
                                my $sl;
                                my $bfile = $self->get_random_seq($backfile,$bbi,$nback,$length,$nseqs,$job);
                                `gimme scan $bfile $mfiles[$j] -c $justscan -n $besthit > bgcut_$job.gff `;
                                if ($besthit > 1) 
                                {
                                    ($unistats) ? ($sl = $helper->count_unique_lines("bgcut_$job.gff")) : 
                                    ($sl = $helper->count_lines("bgcut_$job.gff"));
                                }
                                else { $sl = $helper->count_lines("bgcut_$job.gff"); }
                                unlink($bfile);
                                unlink("bgcut_$job.gff");
                                return($sl);
                            },\@nums);
                        }
                        
                        my $ng = scalar grep { $_ > $stats[$i][$j] } @bhits;
                        $means[$i][$j] = $helper->mean(@bhits);
                        $sds[$i][$j] = $helper->stdev(@bhits);
                        $pvalues[$i][$j] = $ng/$times;
                        $helper->disp("  Average background hits for $mfiles[$j]: $means[$i][$j] +/- $sds[$i][$j]");
                        $helper->disp("  p-value for $mfiles[$j]: $pvalues[$i][$j]\n");
                    }
                }
            }
            if ($scanner =~ m/MotifScanner/i)
            {
                for ($j=0; $j<@mfiles; $j++)
                {
                    $helper->disp("Scanning $seqfile[$i] for motif $mfiles[$j]... Cutoff: $justscan.\n");
                    my $mbase = fileparse($mfiles[$j],'\.[^.]*');
                    my $currout = $self->create_output_file($seqfile[$i],"output",$mbase);
                    ($^O !~ /MSWin/) ? (`./MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o $currout `) :
                    (`MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o $currout `);
                    $self->convert2bed($currout,$cntcol[2],%cnthash) if ($obed);
                    my $nmatch;
                    ($unistats) ? ($nmatch = $helper->count_unique_lines($currout)) : ($nmatch = $helper->count_lines($currout));
                    $nmatch--; # One line header of MotifScanner output
                    $stats[$i][$j] = $nmatch if ($ostat);
                    $helper->disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
                    unlink($currout) if (!$ogff);
                }
                if ($sigmethod eq "bootstrap")
                {
                    $helper->disp("Counting sequences in file $seqfile[$i]...");
                    my $nseqs = $helper->count_fasta($seqfile[$i]);
                    $helper->disp("File $seqfile[$i] has $nseqs sequences.");
                    $helper->disp("\nAssessing statistical significance for tne input motifs for $seqfile[$i] using resampling... It might take some time...\n");
                    
                    for ($j=0; $j<@mfiles; $j++)
                    {
                        my (@bhits,@nums);
                        my $sl;
                        if ($ncore == 1)
                        {
                            for ($k=1; $k<=$times; $k++)
                            {
                                $helper->disp("  Motif $mfiles[$j] -- Iteration $k");
                                $helper->disp("    Generating $nseqs background sequences to be used for motif scanning...");
                                my $bfile = $self->get_random_seq($backfile,$bbi,$nback,$length,$nseqs);
                                $helper->disp("    $nseqs background sequences written in file $bfile...");
                                $helper->disp("    Now scanning file $bfile...");
                                ($^O !~ /MSWin/) ? (`./MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o bgcut.gff `) :
                                (`MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o bgcut.gff `);
                                ($unistats) ? ($sl = $helper->count_unique_lines("bgcut.gff")) : ($sl = $helper->count_lines("bgcut.gff"));
                                $sl--; # One line header of MotifScanner output
                                $helper->disp("    Found $sl hits");
                                push(@bhits,$sl);
                                unlink($bfile);
                            }
                        }
                        else
                        {
                            @nums = (1..$times);
                            $helper->disp("Parallel boostraping enrichment for motif $mfiles[$j]...");
                            @bhits = Parallel::Iterator->iterate_as_array(sub{
                                my ($id,$job) = @_;
                                my $sl;
                                my $bfile = $self->get_random_seq($backfile,$bbi,$nback,$length,$nseqs,$job);
                                ($^O !~ /MSWin/) ? (`./MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o bgcut_$job.gff `) :
                                (`MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o bgcut.gff `);
                                ($unistats) ? ($sl = $helper->count_unique_lines("bgcut_$job.gff")) : ($sl = $helper->count_lines("bgcut_$job.gff"));
                                $sl--; # One line header of MotifScanner output
                                unlink($bfile);
                                unlink("bgcut_$job.gff");
                                return($sl);
                            },\@nums);
                        }
                        
                        my $ng = scalar grep { $_ > $stats[$i][$j] } @bhits;
                        $means[$i][$j] = $helper->mean(@bhits);
                        $sds[$i][$j] = $helper->stdev(@bhits);
                        $pvalues[$i][$j] = $ng/$times;
                        $helper->disp("  Average background hits for $mfiles[$j]: $means[$i][$j] +/- $sds[$i][$j]\n");
                        $helper->disp("  p-value for $mfiles[$j]: $pvalues[$i][$j]\n");
                    }
                }
            }
            $helper->disp(" ");
        }
    }

    # Print stats if requested
    if ($ostat)
    {
        my $outstat = $self->create_output_file(" ","stats");
        open(STATS,">$outstat");
        print STATS "\t",join("\t",@mfiles),"\n";
        for ($i=0; $i<@seqfile; $i++)
        {
            print STATS "$seqfile[$i]";
            for ($j=0; $j<@mfiles; $j++)
            {
                print STATS "\t$stats[$i][$j]";
            }
            print STATS "\n\t";
            for ($j=0; $j<@mfiles; $j++)
            {
                print STATS "\t$pvalues[$i][$j]";
            }
            print STATS "\n\t";
            for ($j=0; $j<@mfiles; $j++)
            {
                print STATS "\t$means[$i][$j]";
            }
            print STATS "\n\t";
            for ($j=0; $j<@mfiles; $j++)
            {
                print STATS "\t$sds[$i][$j]";
            }
            print STATS "\n";
        }
        close(STATS);
    }

    $date = $helper->now;
    $helper->disp("$date - Finished!\n\n");

    if ($olog)
    {
        my $logfile = $self->create_output_file(" ","log");
        open(LOG,">$logfile");
        print LOG join("\n",@log),"\n";
        close(LOG);
    }
}

=head2 parse_motifs

Parse a file of PWMs and construct different files, specific for the scanner in 
use. Internal use.

    $motifscanner->parse_motifs($motifsfile,$scanner);

=cut

sub parse_motifs
{
    my ($self,$infile,$scn) = @_;
    my ($base,$dir) = fileparse($infile,'\.[^.]*');
    my ($line,$om,$f,@fhs,@outnames);
    my $c = 0;
    open(MOTIF,"$infile");
    
    if ($scn =~ m/pwmscan/i)
    {
        $line = <MOTIF>;
        die "\nMotif file (pwmscan) does not appear to have the correct syntax!\n" if ($line !~ /^>/);
        seek(MOTIF,0,0);
        while ($line = <MOTIF>)
        {
            if ($line =~ /^>/) # Signal to open new file
            {
                $line =~ s/\r|\n$//g;
                my $bline = $line;
                $line =~ s/^>//g;
                $om = $self->create_output_file($infile,$line);
                push(@outnames,$om);
                local *PWM;
                open(PWM,">$om");
                $fhs[$c] = *PWM;
                $c++;
                print PWM $bline."\n";
            }
            else 
            { 
                $f = $fhs[$c-1];
                print $f $line;
            }
        }
        # Close the split motifs
        foreach my $fh (@fhs) { close($fh); }
    }
    elsif ($scn =~ m/MotifScanner/i)
    {
        # Some checking
        my $fline = <MOTIF>;
        $fline =~ s/\r|\n$//g;
        die "\nMotif file (MotifScanner) does not appear to have the correct syntax (line 1)!\n" if ($fline !~ /^#INCLUSive/i);
        $c = 0; # Reset counter
        my %motifs;
        while ($line = <MOTIF>)
        {
            $line =~ s/\r|\n$//g;
            if ($line =~ /^#ID/) # Signal to open new file
            { 
                my @sep = split("=",$line);
                $sep[1] =~ s/^\s+|\s+$//g;
                my $name = $sep[1];
                $om = $self->create_output_file($infile,$name);
                push(@outnames,$om);
                local *PWM;
                open(PWM,">$om");
                $fhs[$c] = *PWM;
                push(@{$motifs{$c}},$fline);
                $c++;
            }
            next if ($line =~ /^#W/); # Skip, we will create afterwards...
            push(@{$motifs{$c-1}},$line);       
        }
        foreach my $k (sort { $a <=> $b } keys(%motifs))
        {
            my @fileconts = @{$motifs{$k}};
            my $mw = $#fileconts - 1;
            splice(@fileconts,2,0,"#W = $mw");
            $f = $fhs[$k];
            print $f join("\n",@fileconts);
            print $f "\n";
            close($f)
        }
    }
    
    return(@outnames);
}

=head2 construct_MS_background

Construct a set of background sequences, suitable to use with MotifSampler. 
Internal use.

    $motifscanner->construct_MS_background($backgroundfile);

=cut

sub construct_MS_background
{
    my ($self,$infile) = @_;
    if (-e "MSmodel.bkg")
    {
        $helper->disp("Background model MSmodel.bkg already exists. Proceeding...");
        return 0;
    }
    my ($base,$dir,$ext) = fileparse($infile,'\.[^.]*');
    my $chkfa = $helper->check_fasta($infile);
    if ($chkfa && ! -e $dir.$base.".fa")
    {
        $helper->disp("Background file $infile does not appear to be in FASTA format...");
        $helper->disp("Checking if $infile is in tabular format...");
        my $chktab = $helper->check_tabseq($infile);
        croak "\nBackground file $infile does not appear to be in tabular format either! Exiting...\n" if ($chktab);
        $helper->disp("$infile is in tabular format. Converting to FASTA format...");
        my $converter = HTS::Tools::Convert->new();
        $converter->tab2fasta($infile);
    }
    else { $infile = File::Spec->catfile($dir,$base.".fa"); }
    $helper->disp("Creating Markov model background from $infile... Output will be written in MSmodel.bkg\n");
    `./CreateBackgroundModel -f \"$infile\" -b MSmodel.bkg > temp.std `;
    unlink("temp.std");
    if (! -e "MSmodel.bkg")
    {
        $helper->disp("Markov background model not created!");
        $helper->disp("Probably program CreateBackgroundModel cannot be located on the system! Exiting...");
        croak;
    }
}

=head2 get_random_seq

Get a random sequence from a set of background sequences.

    $motifscanner->get_random_seq($thefile,$itsindex,$howmany,$length,$num);

=cut

sub get_random_seq
{
    my ($self,$tabfasta,$itsindex,$itslen,$length,$num,$job) = @_;
    my ($line,$id,$seq,$currlen,$start,$end);
    my $count = my $safeswitch = 0;
    my $BIG = 1e+6;
    srand;
    my $outfile = $self->create_output_file(" ","sequence");
    $outfile = $outfile."_$job"  if ($job);
    open(RANDSEQ,">$outfile");
    open(FASTA,"<$tabfasta");
    open(INDEX,"$itsindex");
    while ($count < $num && $safeswitch < $BIG)
    {
        # Safeswitch in case too many sequences have small lengths, process shouldn't
        # take forever to complete...
        $safeswitch++;
        $line = $self->get_indexed_line(*FASTA,*INDEX,int(rand($itslen))+1);
        ($id,$seq) = split(/\t/,$line);
        $currlen = length($seq);
        next if ($currlen < $length);
        if ($currlen == $length)
        {
            ($start,$end) = (1,$length);
            $id .= "_".$start."-".$end;
            $self->write_seq(*RANDSEQ,$id,$seq);
            $count++;
        }
        if ($currlen > $length)
        {
            # Restrict the random index generation so that we don't go beyond the sequence end
            $start = int(rand($currlen - $length));
            $end = $start + $length;
            $id .= "_".$start."-".$end;
            $self->write_seq(*RANDSEQ,$id,substr($seq,$start,$length));
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

=head2 write_seq

Write a random sequence, randomly from a set of background sequences. Internal 
use.

    $motifscanner->write_seq($thefile,$theseqid,$theseq);

=cut

sub write_seq
{
    my ($self,$file,$id,$seq) = @_;
    $id = ">".$id if ($id !~ /^>/);
    print $file "$id\n";
    while ($seq) # Write sequences of length 100
    {
        my $wseq = substr($seq,0,100,"");
        print $file "$wseq\n";
    }
}

=head2 convert2bed

Converts the gff output from pwmscan to bed format. The first column of the gff 
(that is the peak/region ID) MUST contain coordinates information in the form 
chr:start-end (track2fasta) or chr:start:end. WARNING! If the fasta files used 
for scanning have been generated with a program like track2fasta from the 
GimmeMotifs suite, then the bed co-ordinates for each occurence can be correctly 
generated. If the sequence ids in the fasta files correspond to peak ids rather 
than exact sequence locations, another file with peak ids and peak centers must 
be provided. The function converts to 6 column bed files. It also converts the 
motif score in gff file to the 0-1000 scale of UCSC genome browser so that motif 
strength can be visualized by color. This is done by linear conversion of the 
form new_value = a*old_value + b and by solving the special case of a 2x2 linear
system (since we know several of the values): 
min(score)*a + b = 0 max(score)*a + b = 1000

    $motifscanner->convert2bed($file_to_convert);

=cut

sub convert2bed
{
    my ($self,$f2c,$cntcol,%ch) = @_;
    my ($base,$dir) = fileparse($f2c,'\.[^.]*');
    my ($line,@lines,@scores,@content,@prelocs,@locs,@newcoord);
    my $bedout = $dir.$base.".bed";
    # In order to determine the coefficients of linear conversion we have to 
    # suck in all gff file, the hard way...
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
    my ($min,$max) = $helper->minmax(@scores);
    # Get coefficients
    my ($a,$b) = $self->naive_solve_two($min,$max);
    open(BED,">$bedout");
    if (%ch)
    {
        foreach $line (@lines)
        {
            @content = split(/\t/,$line);
            @prelocs = split(/\s/,$content[0]);
            @locs = split(":",$prelocs[0]);
            if ($#locs == 1) # track2fasta format
            {
                my $joined = pop(@locs);
                my @splitted = split("-",$joined);
                push(@locs,@splitted);
            }
            # Shift coordinates to motif occurence
            my $strand = "+";
            $strand = "-" if ($content[6] eq "-1" || $content[6] eq "R" || $content[6] eq "-");
            @newcoord = ($locs[0],
                        $ch{$prelocs[1]} - $cntcol + $content[3],
                        $ch{$prelocs[1]} - $cntcol + $content[4],
                        $prelocs[1],$a*$content[5] + $b,$strand);
            print BED join("\t",@newcoord),"\n";
        }
    }
    else
    {
        foreach $line (@lines)
        {
            @content = split(/\t/,$line);
            @prelocs = split(/\s/,$content[0]);
            @locs = split(":",$prelocs[0]);
            if ($#locs == 1) # track2fasta format
            {
                my $joined = pop(@locs);
                my @splitted = split("-",$joined);
                push(@locs,@splitted);
            }
            # Shift coordinates to motif occurence
            my $strand = "+";
            $strand = "-" if ($content[6] eq "-1" || $content[6] eq "R" || $content[6] eq "-");
            @newcoord = ($locs[0],$locs[1] + $content[3],$locs[1] + $content[4],
                         $content[0],$a*$content[5] + $b,$strand);
            print BED join("\t",@newcoord),"\n";
        }
    }
    close(BED);
}

=head2 naive_solve_two

Naive solution of a 2x2 system for proportionally scaling motif scores to BED 
scores. Internal use.

    my ($x,$y) = $motifscanner->naive_solve_two($minscore,$maxscore);

=cut

sub naive_solve_two
{
    my ($self,$min,$max) = @_;
    my $eps = 0.000001;
    my $div;
    (abs($min-$max) < $eps) ? ($div = $eps) : ($div = $max/$min);
    $min = $eps if (!$min);
    my $y = 1000/(1-($div));
    my $x = -$y/$min;
    return($x,$y);
}

=head2 build_index

Index a text file for quick access. Used internally to index the file of 
background sequences.

    my $index = $motifscanner->build_index(*DATAHANDLE,*INDEXHANDLE);

=cut

sub build_index 
{
    my ($self,$datafile,$indexfile) = @_;
    my $offset = 0;
    while (<$datafile>) 
    {
        print $indexfile pack("N",$offset);
        $offset = tell($datafile);
    }
}

=head2 get_indexed_line

Get the line of an indexed file. Returns line or undef if the requested line is 
out of range. Internal use.

    my $iline = $motifscanner->get_indexed_line(*DATAHANDLE,*INDEXHANDLE,
        $theline);

=cut

sub get_indexed_line 
{
    my ($self,$datafile,$indexfile,$linenumber) = @_;
    my $size;    # Size of an index entry
    my $ioffset; # Offset into the index of the entry
    my $entry;   # Index entry
    my $doffset; # Offset into the data file

    $size = length(pack("N",0));
    $ioffset = $size * ($linenumber - 1);
    seek($indexfile,$ioffset,0) or return;
    read($indexfile,$entry,$size);
    $doffset = unpack("N",$entry);
    seek($datafile,$doffset,0);
    return scalar(<$datafile>);
}

=head2 create_output_file

Automatic creation of the output file name depending on the output type. 
Internal use.

    my $name = $motifscanner->create_output_file($thefile,$itstype,$itssubtype);

=cut

sub create_output_file
{
    my ($self,$in,$type,$subtype) = @_;
    my $cdir = getcwd;
    my $date = $helper->now("machine");
    
    if ($type =~ m/sequence/i)
    {
        return(File::Spec->catfile($cdir,"randseq"."$date".".fa"));
    }
    elsif ($type =~ m/stats/i)
    {
        return(File::Spec->catfile($cdir,"stats"."$date".".txt"));
    }
    elsif ($type =~ m/log/i)
    {
        return(File::Spec->catfile($cdir,"log"."$date".".txt"));
    }
    elsif ($type =~ m/output/i)
    {
        my $base = fileparse($in,'\.[^.]*');
        return(File::Spec->catfile($cdir,$base."_".$subtype.".gff"));
    }
    else # Motif file
    {
        my ($base,$dir) = fileparse($in,'\.[^.]*');
        return(File::Spec->catfile($dir,$type."_motif.pwm"));
    }
}

=head2 get

HTS::Tools::Motifscan object getter.

    my $param_value = $motifscanner->get('param_name')

=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Motifscan object setter.

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

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, 
or through the web interface at 
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, 
and then you'll automatically be notified of progress on your bug as I make 
changes.

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
