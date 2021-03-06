=head1 NAME

HTS::Tools::Track::Signal - Conver among several NGS track formats with visualization purposes.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is 

    use HTS::Tools::Track::Signal;
    my %params = (
        'input' => 'wt.bam'
        'source' => 'bam',
        'destination' => 'bigwig',
        'options' => (
            'color' => '0,160,0',
            'visbility' => 'dense'
        )
    )
    my $track = HTS::Tools::Track::Signal->new(\%params);
    $track->run;

The acceptable parameters are as follows:

=over 4

=item I<input> B<(required)>

Ipse deus

=item I<source> B<(required)>

Velox discurrere

=item I<destination> B<(required)>

Gaudet

=item I<dir> B<(optional)>

In altis montibus

=item I<urlbase> B<(required when destination is bigbed, bigwig, bam)>

Et subitas

=item I<org> B<(required when destination is bigbed, bigwig, wig, bedgraph)>

Concipit

=item I<gversion> B<(required when destination is bigbed, bigwig, wig, bedgraph)>

Ipse fugas

=item I<cleanlevel> B<(optional)>

The cleanlevel parameter controls what filtering will be applied to the raw reads so as to produce
the signal track. It can have three values: 0 for not cleaning anything (reporting reads as they
are, no unique and no removal of unlocalized regions and mitochondrial DNA reads), 1 for removing
unlocalized regions (chrU, hap, random etc.), 2 for removing reads of level 1 plus mitochondrial
reads (chrM) and 3 for removing reads of level 2 plus returning unique reads only. The default is
level 1.

=item  I<sort> B<(optional)>

Sort the BAM or BED files according to co-ordinates. This process is required for some conversions
so if you do not supply this parameter, make sure that the source tracks are sorted.

=item  I<options> B<(optional)>

Ipse Deus

=item I<silent> B<(optional)>

Use this parameter if you want to turn informative messages off.

=head1 OUTPUT

The output of the module is usually a UCSC Genome Browser track.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Track::Signal;

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use File::Basename;
use File::Temp;
use File::Spec;

use HTS::Tools::Constants;
use HTS::Tools::Fetch;
use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper $const $fetcher);

our $MODNAME = "HTS::Tools::Track::Signal";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Create an NGS signal track using 3rd party tools.";

BEGIN {
    $helper = HTS::Tools::Utils->new();
    $const = HTS::Tools::Constants->new();
    $fetcher = HTS::Tools::Fetch->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

The HTS::Tools::Track::Signal object constructor. It accepts a set of parameters that are required to run the
motifscanner and get the output.

    my $track = HTS::Tools::Track::Signal->new({'input' => 'wt.bam','output' => 'wt.wig',
        options => {'color' => '0,160,0'}});
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
    $helper->advertise($MODNAME,$VERSION,$AUTHOR,$EMAIL,$DESC);

    # Validate the input parameters
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "track_signal",
        "params" => $params});
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

The HTS::Tools::Track::Signal run subroutine. It runs the track converter with the given parameters in the 
constructor.

    $track->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    # Copy some memory-less variables to avoid rewriting the whole thing...
    my $input = $self->get("input");
    my $source = lc($self->get("source"));
    my $destination = lc($self->get("destination"));
    my $dir = $self->get("dir");
    my $org = $self->get("org");
    my $ver = $self->get("gversion");
    my $clevel = $self->get("cleanlevel");
    my $sort = $self->get("sort");
    my $options = $self->get("options");
    
    #print "\n\n$source\n$destination\n";
    #print $source."2".$destination."\n";
    #if ($source."2".$destination eq "bedgraph2bigwig") {
    #    print "Fucking matches!\n\n";
    #}

    # Check the existence of proper Constants
    $self->check_constants($source,$destination);

    # Actual job must be performed, organized in subroutines...
    my ($track,$header);
    if ($source."2".$destination eq "bam2bedgraph")
    {
        ($track,$header) = $self->bam2bedgraph($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bam2bigbed")
    {
        ($track,$header) = $self->bam2bigbed($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bam2bigwig")
    {
        ($track,$header) = $self->bam2bigwig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bam2wig")
    {
        ($track,$header) = $self->bam2wig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bed2bam")
    {
        ($track,$header) = $self->bed2bam($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bed2bedgraph")
    {
        ($track,$header) = $self->bed2bedgraph($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bed2bigbed")
    {
         ($track,$header) = $self->bed2bigbed($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bed2bigwig")
    {
        ($track,$header) = $self->bed2bigwig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bed2wig")
    {
        ($track,$header) = $self->bed2wig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bedgraph2bigwig")
    {
        ($track,$header) = $self->bedgraph2bigwig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bedgraph2wig")
    {
        ($track,$header) = $self->bedgraph2wig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bigbed2bam")
    {
        ($track,$header) = $self->bigbed2bam($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bigbed2bed")
    {
        ($track,$header) = $self->bigbed2bed($input,$dir,$options);
    }
    elsif ($source."2".$destination eq "bigbed2bedgraph")
    {
        ($track,$header) = $self->bigbed2bedgraph($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bigbed2bigwig")
    {
        ($track,$header) = $self->bigbed2bigwig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "bigwig2bedgraph")
    {
        ($track,$header) = $self->bigwig2bedgraph($input,$dir,$options);
    }
    elsif ($source."2".$destination eq "bigwig2wig")
    {
        ($track,$header) = $self->bigwig2wig($input,$dir,$options);
    }
    elsif ($source."2".$destination eq "wig2bigwig")
    {
        ($track,$header) = $self->wig2bigwig($input,$dir,$org,$ver,$options);
    }
    elsif ($source."2".$destination eq "wig2bedgraph")
    {
        ($track,$header) = $self->bigbed2bigwig($input,$dir,$org,$ver,$options);
    }
    elsif ($source."2".$destination eq "sam2bedgraph")
    {
        ($track,$header) = $self->sam2bedgraph($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "sam2bigbed")
    {
        ($track,$header) = $self->sam2bigbed($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "sam2bigwig")
    {
        ($track,$header) = $self->sam2bigwig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }
    elsif ($source."2".$destination eq "sam2wig")
    {
        ($track,$header) = $self->sam2wig($input,$dir,$org,$ver,$clevel,$sort,$options);
    }

    # Then, header creation and append or write to file according to destination type
    return($track,$header);
    
}

=head2 bam2bedgraph

bam2bedgraph converter using 3rd party tools.

    $track->bam2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bed
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $bam2bed = File::Spec->catfile($bedtoolshome,"bedtools bamtobed -split -i");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".tmpbed");
    $helper->disp("Converting to bed...");
    my $fail = system($bam2bed." ".$input." > ".$output);
    return(0) if ($fail);

    $header = "track type=bed";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".bed");
    open(OUTPUT,$output);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output);
    $track = $outfinal;
    
    return($track,$header);
}

=head2 bam2bedgraph

bam2bedgraph converter using 3rd party tools.

    $track->bam2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bedgraph
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $bam2bed = File::Spec->catfile($bedtoolshome,"bedtools bamtobed -split -i");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bam2bed." ".$input." > ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".tmp");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bedgraph...\n");
    my $fail2 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output2);
    ($fail2) ? (return(0)) : ($track = $output2);

    $header = "track type=bedGraph";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".bedGraph");
    open(OUTPUT,$output2);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output1);
    unlink($output2);
    $track = $outfinal;
    
    return($track,$header);
}

=head2 bam2bigbed

bam2bigbed converter using 3rd party tools.

    $track->bam2bigbed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bigbed
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bam2bed = File::Spec->catfile($bedtoolshome,"bedtools bamtobed -split -i");
    my $bed2bigbed = File::Spec->catfile($kenthome,"bedToBigBed");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bam2bed." ".$input." > ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bigBed");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bigbed...");
    my $fail2 = system($bed2bigbed." ".$sorted." ".$chromsize." ".$output2);
    ($fail2) ? (return(0)) : ($track = $output2);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigBed";
    $header = "track type=bigBed";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bbh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    unlink($output1);
    
    return($track,$header);
}

=head2 bam2bigwig

bam2bigwig converter using 3rd party tools.

    $track->bam2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bigwig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bam2bed = File::Spec->catfile($bedtoolshome,"bedtools bamtobed -split -i");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bam2bed." ".$input." > ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bedGraph");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bedgraph...");
    my $fail2 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output2);
    return(0) if ($fail2);

    my $output3 = File::Spec->catfile($dir,$basename.".bigWig");
    $helper->disp("Converting to bigwig...\n");
    my $fail3 = system($bedgraph2bigwig." ".$output2." ".$chromsize." ".$output3);
    ($fail3) ? (return(0)) : ($track = $output3);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigWig";
    $header = "track type=bigWig";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bwh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    unlink($output1);
    unlink($output2);
    
    return($track,$header);
}

=head2 bam2wig

bam2wig converter using 3rd party tools.

    $track->bam2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2wig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bam2bed = File::Spec->catfile($bedtoolshome,"bedtools bamtobed -split -i");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");
    my $bigwig2wig = File::Spec->catfile($kenthome,"bigWigToWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bam2bed." ".$input." > ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bedGraph");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bedgraph...");
    my $fail2 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output2);
    return(0) if ($fail2);

    my $output3 = File::Spec->catfile($dir,$basename.".bigWig");
    $helper->disp("Converting to bigwig...");
    my $fail3 = system($bedgraph2bigwig." ".$output2." ".$chromsize." ".$output3);
    return(0) if ($fail3);

    my $output4 = File::Spec->catfile($dir,$basename.".tmp");
    $helper->disp("Converting to wig...\n");
    my $fail4 = system($bigwig2wig." ".$output3." ".$output4);
    return(0) if ($fail4);

    $header = "track type=wiggle_0";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".wig");
    open(OUTPUT,$output4);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output1);
    unlink($output2);
    unlink($output3);
    unlink($output4);
    $track = $outfinal;

    return($track,$header);
}

=head2 bed2bam

bed2bam converter using 3rd party tools.

    $track->bed2bam($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bam
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $samtoolshome = $const->get("SAMTOOLS_HOME");
    my $bed2bam = File::Spec->catfile($bedtoolshome,"bedtools bedtobam -i");
    my $samtoolsindex = File::Spec->catfile($samtoolshome,"samtools index");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".bam");
    my $sorted = $self->clean_bedstar($input,$clevel,$sort);
    $helper->disp("Converting to bam...");
    my $fail = system($bed2bam." ".$sorted." -g ".$chromsize." > ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    $helper->disp("Indexing the bam...\n");
    my $fi = system($samtoolsindex." ".$track);
    return(0) if ($fi);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bam";
    $header = "track type=bam";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bmh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    
    return($track,$header);
}

=head2 bed2bedgraph

bed2bedgraph converter using 3rd party tools.

    $track->bed2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bedgraph
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".tmp");
    my $sorted = $self->clean_bedstar($input,$clevel,$sort);
    $helper->disp("Converting to bedgraph...\n");
    my $fail = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    $header = "track type=bedGraph";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".bedGraph");
    open(OUTPUT,$output);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output);
    $track = $outfinal;
    
    return($track,$header);
}

=head2 bed2bigbed

bed2bigbed converter using 3rd party tools.

    $track->bed2bigbed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bigbed
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bed2bigbed = File::Spec->catfile($kenthome,"bedToBigBed");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".bigBed");
    my $sorted = $self->clean_bedstar($input,$clevel,$sort);

    $helper->disp("Converting to bigbed...\n");
    my $fail = system($bed2bigbed." ".$sorted." ".$chromsize." ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    # Construct track header
    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigBed";
    $header = "track type=bigBed";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bbh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    
    return($track,$header);
}

=head2 bed2bigwig

bed2bigwig converter using 3rd party tools.

    $track->bed2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bigwig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $kenthome = $const->get("KENTBIN_HOME");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bedGraph");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bedgraph...");
    my $fail1 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bigWig");
    $helper->disp("Converting to bigwig...\n");
    my $fail2 = system($bedgraph2bigwig." ".$output1." ".$chromsize." ".$output2);
    ($fail2) ? (return(0)) : ($track = $output2);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigWig";
    $header = "track type=bigWig";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bwh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    unlink($output1);
    
    return($track,$header);
}

=head2 bed2wig

bed2wig converter using 3rd party tools.

    $track->bed2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2wig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $kenthome = $const->get("KENTBIN_HOME");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");
    my $bigwig2wig = File::Spec->catfile($kenthome,"bigWigToWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bedGraph");
    my $sorted = $self->clean_bedstar($input,$clevel,$sort);
    $helper->disp("Converting to bedgraph...");
    my $fail1 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bigWig");
    $helper->disp("Converting to bigwig...");
    my $fail2 = system($bedgraph2bigwig." ".$output1." ".$chromsize." ".$output2);
    return(0) if ($fail2);

    my $output3 = File::Spec->catfile($dir,$basename.".tmp");
    $helper->disp("Converting to wig...\n");
    my $fail3 = system($bigwig2wig." ".$output2." ".$output3);
    return(0) if ($fail3);

    $header = "track type=wiggle_0";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".wig");
    open(OUTPUT,$output3);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output1);
    unlink($output2);
    unlink($output3);
    $track = $outfinal;

    return($track,$header);
}

=head2 bedgraph2bigwig

bedgraph2bigwig converter using 3rd party tools.

    $track->bedgraph2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bedgraph2bigwig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);
    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".bigWig");
    my $sorted = $self->clean_bedstar($input,$clevel,$sort);
    $helper->disp("Converting to bigWig...\n");
    my $fail = system($bedgraph2bigwig." ".$sorted." ".$chromsize." ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigWig";
    $header = "track type=bigWig";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bwh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    
    return($track,$header);
}

=head2 bedgraph2wig

bedgraph2wig converter using 3rd party tools.

    $track->bedgraph2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bedgraph2wig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");
    my $bigwig2wig = File::Spec->catfile($kenthome,"bigWigToWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bigWig");
    my $sorted = $self->clean_bedstar($input,$clevel,$sort);
    $helper->disp("Converting to bigwig...");
    my $fail1 = system($bedgraph2bigwig." ".$sorted." ".$chromsize." ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".tmp");
    $helper->disp("Converting to wig...\n");
    my $fail2 = system($bigwig2wig." ".$output1." ".$output2);
    return(0) if ($fail2);

    $header = "track type=wiggle_0";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".wig");
    open(OUTPUT,$output2);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output1);
    unlink($output2);
    $track = $outfinal;

    return($track,$header);
}

=head2 bigbed2bam

bigbed2bam converter using 3rd party tools.

    $track->bigbed2bam($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bam
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $samtoolshome = $const->get("SAMTOOLS_HOME");
    my $bigbed2bed = File::Spec->catfile($kenthome,"bigBedToBed");
    my $bed2bam = File::Spec->catfile($bedtoolshome,"bedtools bedtobam -i");
    my $samtoolsindex = File::Spec->catfile($samtoolshome,"samtools index");

    print "\n\n$input\n\n";
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bigbed2bed." ".$input." ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bam");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bam...");
    my $fail2 = system($bed2bam." ".$sorted." -g ".$chromsize." > ".$output2);
    ($fail2) ? (return(0)) : ($track = $output2);

    $helper->disp("Indexing the bam...\n");
    my $fi = system($samtoolsindex." ".$track);
    return(0) if ($fi);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bam";
    $header = "track type=bam";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bmh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    unlink($output1);
    
    return($track,$header);
}

=head2 bigbed2bed

bigbed2bed converter using 3rd party tools.

    $track->bigbed2bed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bed
{
    my ($self,$input,$dir,$options) = @_;
    my ($track,$header,$tracktmp);

    my $kenthome = $const->get("KENTBIN_HOME");
    my $bigbed2bed = File::Spec->catfile($kenthome,"bigBedToBed");
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".tmpbed");
    
    $helper->disp("Converting to bed...\n");
    my $fail = system($bigbed2bed." ".$input." ".$output);
    ($fail) ? (return(0)) : ($tracktmp = $output);

    # Construct track header
    $header = "track type=bed";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outfinal = File::Spec->catfile($dir,$basename.".bed");
    open(OUTPUT,$output);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output);
    $track = $outfinal;
    
    return($track,$header);
}

=head2 bigbed2bedgraph

bigbed2bedgraph converter using 3rd party tools.

    $track->bigbed2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bedgraph
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $bigbed2bed = File::Spec->catfile($kenthome,"bigBedToBed");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bigbed2bed." ".$input." ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".tmp");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bedgraph...\n");
    my $fail2 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output2);
    return(0) if ($fail2);

    $header = "track type=bedGraph";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }

    my $outfinal = File::Spec->catfile($dir,$basename.".bedGraph");
    open(OUTPUT,$output2);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output1);
    unlink($output2);
    $track = $outfinal;

    return($track,$header);
}

=head2 bigbed2bigwig

bigbed2bigwig converter using 3rd party tools.

    $track->bigbed2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bigwig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bedtoolshome = $const->get("BEDTOOLS_HOME");
    my $bigbed2bed = File::Spec->catfile($kenthome,"bigBedToBed");
    my $genomecov = File::Spec->catfile($bedtoolshome,"bedtools genomecov -bg -i");
    my $bedgraph2bigwig = File::Spec->catfile($kenthome,"bedGraphToBigWig");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output1 = File::Spec->catfile($dir,$basename.".bed");
    $helper->disp("Converting to bed...");
    my $fail1 = system($bigbed2bed." ".$input." ".$output1);
    return(0) if ($fail1);

    my $output2 = File::Spec->catfile($dir,$basename.".bedGraph");
    my $sorted = $self->clean_bedstar($output1,$clevel,$sort);
    $helper->disp("Converting to bedgraph...");
    my $fail2 = system($genomecov." ".$sorted." -g ".$chromsize." > ".$output2);
    return(0) if ($fail2);

    my $output3 = File::Spec->catfile($dir,$basename.".bigWig");
    $helper->disp("Converting to bigWig...\n");
    my $fail3 = system($bedgraph2bigwig." ".$output2." ".$chromsize." ".$output3);
    ($fail3) ? (return(0)) : ($track = $output3);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigWig";
    $header = "track type=bigWig";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bwh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    unlink($output1);
    unlink($output2);
    
    return($track,$header);
}

=head2 bigwig2bedgraph

bigwig2bedgraph converter using 3rd party tools.

    $track->bigwig2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigwig2bedgraph
{
    my ($self,$input,$dir,$options) = @_;
    my ($track,$header);

    my $kenthome = $const->get("KENTBIN_HOME");
    my $bigwig2bedgraph = File::Spec->catfile($kenthome,"bigWigToBedGraph");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".tmpbedgraph");
    $helper->disp("Converting to bedgraph...\n");
    my $fail = system($bigwig2bedgraph." ".$input." ".$output);
    return(0) if ($fail);

    $header = "track type=bedGraph";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outfinal = File::Spec->catfile($dir,$basename.".bedGraph");
    open(OUTPUT,$output);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output);
    $track = $outfinal;
    
    return($track,$header);
}

=head2 bigwig2wig

bigwig2wig converter using 3rd party tools.

    $track->bigwig2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigwig2wig
{
    my ($self,$input,$dir,$options) = @_;
    my ($track,$header);

    my $kenthome = $const->get("KENTBIN_HOME");
    my $bigwig2wig = File::Spec->catfile($kenthome,"bigWigToWig");
    
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".tmpwig");
    $helper->disp("Converting to wig...\n");
    my $fail = system($bigwig2wig." ".$input." ".$output);
    return(0) if ($fail);

    $header = "track type=wiggle_0";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outfinal = File::Spec->catfile($dir,$basename.".wig");
    open(OUTPUT,$output);
    open(OUTFINAL,">$outfinal");
    print OUTFINAL $header,"\n";
    while (<OUTPUT>)
    {
        print OUTFINAL $_;
    }
    close(OUTPUT);
    close(OUTFINAL);
    unlink($output);
    $track = $outfinal;
    
    return($track,$header);
}

=head2 sam2bam

sam2bam converter using 3rd party tools.

    $track->sam2bam($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bam
{
    my ($self,$input,$dir,$options) = @_;
    my ($track,$header);

    my $samtoolshome = $const->get("SAMTOOLS_HOME");
    my $sam2bam = File::Spec->catfile($samtoolshome,"samtools view -bh");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".bam");
    $helper->disp("Converting to bam...\n");
    my $fail = system($sam2bam." ".$input." -o ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bam";
    $header = "track type=bam";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bmh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    
    return($track,$header);
}

=head2 sam2bedgraph

sam2bedgraph converter using 3rd party tools.

    $track->sam2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bedgraph
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my $tmptrack = $self->sam2bam($input,$dir,$options);
    return($self->bam2bedgraph($tmptrack,$dir,$org,$ver,$clevel,$sort,$options));
}

=head2 sam2bigbed

sam2bigbed converter using 3rd party tools.

    $track->sam2bigbed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bigbed
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my $tmptrack = $self->sam2bam($input,$dir,$options);
    return($self->bam2bigbed($tmptrack,$dir,$org,$ver,$clevel,$sort,$options));
}

=head2 sam2bigwig

sam2bigwig converter using 3rd party tools.

    $track->sam2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bigwig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my $tmptrack = $self->sam2bam($input,$dir,$options);
    return($self->bam2bigwig($tmptrack,$dir,$org,$ver,$clevel,$sort,$options));
}

=head2 sam2wig

sam2wig converter using 3rd party tools.

    $track->sam2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2wig
{
    my ($self,$input,$dir,$org,$ver,$clevel,$sort,$options) = @_;
    my $tmptrack = $self->sam2bam($input,$dir,$options);
    return($self->bam2wig($tmptrack,$dir,$org,$ver,$clevel,$sort,$options));
}

=head2 wig2bigwig

wig2bigwig converter using 3rd party tools.

    $track->wig2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub wig2bigwig
{
    my ($self,$input,$dir,$org,$ver,$options) = @_;
    my ($track,$header);

    my $chromsize = $self->get_chrom_size($org,$ver);
    my $kenthome = $const->get("KENTBIN_HOME");
    my $wig2bigwig = File::Spec->catfile($kenthome,"wigToBigWig");

    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".bigWig");
    $helper->disp("Converting to bigWig...\n");
    my $fail = system($wig2bigwig." ".$input." ".$chromsize." ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigWig";
    $header = "track type=bigWig";
    while (my ($key,$value) = each (%{$options}))
    {
        $header .= " ".$key."=".$value;
    }
    my $outheader = File::Spec->catfile($dir,$basename.".bwh");
    open(HEADER,">$outheader");
    print HEADER $header,"\n";
    close(HEADER);
    
    return($track,$header);
}

=head2 wig2bedgraph

wig2bedgraph converter using 3rd party tools.

    $track->wig2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub wig2bedgraph
{
    my ($self,$input,$dir,$org,$ver,$options) = @_;
    my $tmptrack = $self->wig2bigwig($input,$dir,$org,$ver,$options);
    return($self->bigwig2bedgraph($tmptrack,$dir,$options));
}

=head2 check_constants

Check if required 3rd party tools exist for the required conversion. Internal use.

    $signaler->check_constants("sam","bigwig");

=cut

sub check_constants
{
    my $self = shift @_;
    my ($source,$dest) = @_;
    
    if ($source."2".$dest =~ m/bam2bedgraph/i)
    {
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bam2bigbed/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bam2bigwig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/bam2wig/i)
    {
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/bed2bam/i)
    {
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bed2bedgraph/i)
    {
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bed2bigbed/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/bed2bigwig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bed2wig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bedgraph2bigwig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/bedgraph2wig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/bigbed2bam/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bigbed2bed/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
    }
    elsif ($source."2".$dest =~ m/bigbed2bedgraph/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/bigbed2bigwig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/bigwig2wig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
    }
    elsif ($source."2".$dest =~ m/bigwig2bedgraph/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
    }
    elsif ($source."2".$dest =~ m/sam2bedgraph/i)
    {
        croak "SAM tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/sam2bigbed/i)
    {
        croak "SAM tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
    }
    elsif ($source."2".$dest =~ m/sam2bigwig/i)
    {
        croak "SAM tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
         croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/sam2wig/i)
    {
        croak "SAM tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
        croak "BED tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
    elsif ($source."2".$dest =~ m/wig2bigwig/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
    }
    elsif ($source."2".$dest =~ m/wig2bedgraph/i)
    {
        croak "Kent tools not found in 3rd party tools! Please install them and retry!"
            if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
            if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
    }
}

=head2 clean_bedstar

Helper cleaning/sorting function for BED-like files (bed, bedgraph). Internal use.

    $track->sort_bedstar($file,$cleanlevel,$sort);
    
=cut

sub clean_bedstar
{
    my ($self,$infile,$clevel,$sort) = @_;
    my $tmpdir = $self->get("tmpdir");
    my $tmpfile;
    my ($lev1,%seen); # Vars for pure Perl cleaning on Windows

    if ($sort)
    {
        if ($^O !~ /MSWin/) # Case of linux, easy sorting
        {
            $helper->disp("Cleaning and sorting bed-like file $infile...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp.in$$");
            if ($clevel == 0)
            {
                `sort -k1,1 -k2g,2 -k3g,3 $infile > $tmpfile `;
            }
            elsif ($clevel == 1)
            {
                `grep -vP 'chrU|chrG|rand|hap|loc|cox' $infile | sort -k1,1 -k2g,2 -k3g,3  > $tmpfile `;
            }
            elsif ($clevel == 2)
            {
                `grep -vP 'chrM|chrU|chrG|rand|hap|loc|cox' $infile | sort -k1,1 -k2g,2 -k3g,3  > $tmpfile `;
            }
            elsif ($clevel == 3)
            {
                `grep -vP 'chrM|chrU|chrG|rand|hap|loc|cox' $infile | sort -k1,1 -k2g,2 -k3g,3 -u  > $tmpfile `;
            }
            $infile = $tmpfile;
        }
        else # We are in Windows... package required
        {
            $helper->try_module("File::Sort","sort_file");
            eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
            $helper->disp("Cleaning and sorting file $infile...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp.tmp");
            if ($clevel != 0)
            {
                $lev1 = File::Spec->catfile($tmpdir,"lev1.tmp");
                open(LEV1INPUT,$infile);
                open(LEV1OUTPUT,">$lev1");
                if ($clevel != 3)
                {
                    while (<LEV1INPUT>)
                    {
                        next if ($_ =~ m/chrU|chrG|rand|hap|loc|cox/i && $clevel == 1);
                        next if ($_ =~ m/chrM|chrU|chrG|rand|hap|loc|cox/i && $clevel == 2);
                        print LEV1OUTPUT $_;
                    }
                }
                else
                {
                    while (my $line = <LEV1INPUT>)
                    {
                        my @cols = split(/\t/,$line);
                        next if ($line =~ m/chrM|chrU|chrG|rand|hap|loc|cox/i && $seen{join("\t",@cols[0..2])});
                        $seen{join("\t",@cols[0..2])}++;
                        print LEV1OUTPUT $line;
                    }
                }
                close(LEV1INPUT);
                close(LEV1OUTPUT);
                $infile = $lev1;
            }
            sort_file({
                I => $infile,
                o => $tmpfile,
                k => ['1,1','2n,2','3n,3'],
                t => "\t"
            }); 
            $infile = $tmpfile;
        }
    }
    else
    {
        if ($^O !~ /MSWin/)
        {
            $helper->disp("Cleaning bed-like file $infile...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp.in$$");
            if ($clevel == 0)
            {
                $tmpfile = $infile;
            }
            elsif ($clevel == 1)
            {
                `grep -vP 'chrU|chrG|rand|hap|loc|cox' $infile > $tmpfile `;
            }
            elsif ($clevel == 2)
            {
                `grep -vP 'chrM|chrG|rand|hap|loc|cox' $infile > $tmpfile `;
            }
            elsif ($clevel == 3) # In this case, sorting is forced...
            {
                `grep -vP 'chrM|chrG|rand|hap|loc|cox' $infile | sort -k1,1 -k2g,2 -k3g,3 -u  > $tmpfile `;
            }
            $infile = $tmpfile;
        }
        else # We are in Windows... package required
        {
            $helper->disp("Cleaning file $infile...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp.tmp");
            if ($clevel != 0)
            {
                $lev1 = File::Spec->catfile($tmpdir,"lev1.tmp");
                open(LEV1INPUT,$infile);
                open(LEV1OUTPUT,">$lev1");
                if ($clevel != 3)
                {
                    while (<LEV1INPUT>)
                    {
                        next if ($_ =~ m/chrG|rand|hap|loc|cox/i && $clevel == 1);
                        next if ($_ =~ m/chrM|chrG|rand|hap|loc|cox/i && $clevel == 2);
                        print LEV1OUTPUT $_;
                    }
                }
                else
                {
                    while (my $line = <LEV1INPUT>)
                    {
                        my @cols = split(/\t/,$line);
                        next if ($line =~ m/chrM|chrG|rand|hap|loc|cox/i && $seen{join("\t",@cols[0..2])});
                        $seen{join("\t",@cols[0..2])}++;
                        print LEV1OUTPUT $line;
                    }
                }
                close(LEV1INPUT);
                close(LEV1OUTPUT);
                $infile = $lev1;
            }
            else
            {
                $tmpfile = $infile;
            }
            $infile = $tmpfile;
        }
    }

    return($infile);
}

=head2 get_chrom_size

Fetch/find in path a file containing the size of each chromosome, according to the requested organism
and genome version.

    $track->get_chrom_size($organism,$version);

=cut

sub get_chrom_size
{
    my ($self,$org,$ver) = @_;
    if ($const->get("IGENOMES_HOME") && $const->get("REMOTE_HOST"))
    {
        return($self->format_igenomes_chrom_size($org,$ver));
    }
    elsif ($const->get("IGENOMES_HOME") && !$const->get("REMOTE_HOST"))
    {
        return($self->format_igenomes_chrom_size($org,$ver));
    }
    elsif (!$const->get("IGENOMES_HOME") && $const->get("REMOTE_HOST"))
    {
        return($fetcher->fetch_chrom_info($org));
    }
}

=head2 format_igenomes

Construct a string representing with the correct nomenclature, species for each supported iGenomes
chromosomal size annotation. Mostly for internal use.

    $track->format_igenomes($organism,$version);

=cut

sub format_igenomes_chrom_size
{
    my ($self,$org,$ver) = @_;
    my $base = $const->get("IGENOMES_HOME");
    if ($org =~ m/human/)
    {
        if ($ver =~ m/hg19/)
        {
            return(File::Spec->catfile($base,"Homo_sapiens","UCSC","hg19","Annotation","Genes","ChromInfo.txt"));
        }
        elsif ($ver =~ m/hg18/)
        {
            return(File::Spec->catfile($base,"Homo_sapiens","UCSC","hg18","Annotation","Genes","ChromInfo.txt"));
        }
    }
    elsif ($org =~ m/mouse/)
    {
        if ($ver =~ m/mm10/)
        {
            return(File::Spec->catfile($base,"Mus_musculus","UCSC","mm10","Annotation","Genes","ChromInfo.txt"));
        }
        elsif ($ver =~ m/mm9/)
        {
            return(File::Spec->catfile($base,"Mus_musculus","UCSC","mm9","Annotation","Genes","ChromInfo.txt"));
        }
    }
    elsif ($org =~ m/rat/)
    {
        if ($ver =~ m/rn5/)
        { 
            return(File::Spec->catfile($base,"Rattus_norvegicus","UCSC","rn5","Annotation","Genes","ChromInfo.txt"));
        }
        elsif ($ver =~ m/rn4/)
        { 
            return(File::Spec->catfile($base,"Rattus_norvegicus","UCSC","rn4","Annotation","Genes","ChromInfo.txt"));
        }
    }
    elsif ($org =~ m/fly/)
    {
        if ($ver =~ m/dm3/)
        {
            return(File::Spec->catfile($base,"Drosophila_melanogaster","UCSC","dm3","Annotation","Genes","ChromInfo.txt"));
        }
    }
    elsif ($org =~ m/zebrafish/)
    {
        if ($ver =~ m/danrer7/i)
        {
            return(File::Spec->catfile($base,"Danio_rerio","UCSC","danRer7","Annotation","Genes","ChromInfo.txt"));
        }
    }
}

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
