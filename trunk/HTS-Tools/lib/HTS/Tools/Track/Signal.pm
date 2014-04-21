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

Ipse deus

=item I<org> B<(required when destination is bigbed, bigwig, wig, bedgraph)>

Nudus nudos

=item I<gversion> B<(required when destination is bigbed, bigwig, wig, bedgraph)>

Iubet ire

=item  I<options> B<(optional)>

Ministrus

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

The HTS::Tools::Track::Signal run subroutine. It runs the track converter with the given parameters in the 
constructor.

    $track->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    # Copy some memory-less variables to avoid rewriting the whole thing...
    my $input = $self->get("input");
    my $source = $self->get("source");
    my $destination = $self->get("destination");
    my $dir = $self->get("dir");
    my $org = $self->get("org");
    my $ver = $self->get("gversion");
    my $options = $self->get("options");

    # Check the existence of propoer Constants
    $self->check_constants($source,$destination);

    # Actual job must be performed, organized in subroutines...
    my ($track,$header);
    use v5.14;
    given($source."2".$destination)
    {
        when(/bam2bedgraph/i)
        {
        }
        when(/bam2bigbed/i)
        {
        }
        when(/bam2bigwig/i)
        {
        }
        when(/bam2wig/i)
        {
        }
        when (/bed2bam/i)
        {
        }
        when (/bed2bedgraph/i)
        {
        }
        when (/bed2bigbed/i)
        {
             ($track,$header) = $self->bed2bigbed($input,$dir,$org,$ver,$options);
        }
        when (/bed2bigwig/i)
        {
        }
        when (/bed2wig/i)
        {
        }
        when(/bedgraph2bigwig/i)
        {
        }
        when(/bedgraph2wig/i)
        {
        }
        when(/bigbed2bam/i)
        {
        }
        when(/bigbed2bed/i)
        {
            ($track,$header) = $self->bigbed2bed($input,$dir,$org,$options);
        }
        when(/bigbed2bedgraph/i)
        {
        }
        when(/bigbed2bigwig/i)
        {
        }
        when(/bigwig2wig/i)
        {
        }
        when(/bigwig2bedgraph/i)
        {
        }
        when(/wig2bigwig/i)
        {
        }
        when(/wig2bedgraph/i)
        {
        }
        when(/sam2bedgraph/i)
        {
        }
        when(/sam2bigbed/i)
        {
        }
        when(/sam2bigwig/i)
        {
        }
        when(/sam2wig/i)
        {
        }
    }

    # Then, header creation and append or write to file according to destination type
    
}

=head2 bam2bedgraph

bam2bedgraph converter using 3rd party tools.

    $track->bam2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bedgraph
{
}

=head2 bam2bigbed

bam2bigbed converter using 3rd party tools.

    $track->bam2bigbed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bigbed
{
}

=head2 bam2bigwig

bam2bigwig converter using 3rd party tools.

    $track->bam2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2bigwig
{
}

=head2 bam2wig

bam2wig converter using 3rd party tools.

    $track->bam2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bam2wig
{
}

=head2 bed2bam

bed2bam converter using 3rd party tools.

    $track->bed2bam($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bam
{
}

=head2 bed2bedgraph

bed2bedgraph converter using 3rd party tools.

    $track->bed2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bedgraph
{
}

=head2 bed2bigbed

bed2bigbed converter using 3rd party tools.

    $track->bed2bigbed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bigbed
{
    my ($self,$input,$dir,$org,$ver,$options) = @_;
    my ($track,$header);

    # Construct file
    my $chromsize;
    if ($const->get("IGENOMES_HOME") && $const->get("REMOTE_HOST"))
    {
        $chromsize = $self->format_igenomes_chrom_size($org,$ver);
    }
    elsif ($const->get("IGENOMES_HOME") && !$const->get("REMOTE_HOST"))
    {
        $chromsize = $self->format_igenomes_chrom_size($org,$ver);
    }
    elsif (!$const->get("IGENOMES_HOME") && $const->get("REMOTE_HOST"))
    {
        $chromsize = $fetcher->fetch_chrom_info($org);
    }
    my $kenthome = $const->get("KENTBIN_HOME");
    my $bed2bigbed = File::Spec->catfile($kenthome,"bedToBigBed");
    my ($basename,$dirname,$ext) = fileparse($input,'\.[^.]*');
    my $output = File::Spec->catfile($dir,$basename.".bigBed");
    my $sorted = $self->sort_bedstar($input);
    my $fail = system($bed2bigbed." ".$sorted." ".$chromsize." ".$output);
    ($fail) ? (return(0)) : ($track = $output);

    # Construct track header
    $options->{"bigDataUrl"} = $options->{"bigDataUrl"}."/".$basename.".bigBed";
    $header = "track type=bigBed";
	while ($key,$value) = each (%{$options})
	{
		$header .= " ".$key."=".$value;
    }
    
    return($track,$header);
}

=head2 bed2bigwig

bed2bigwig converter using 3rd party tools.

    $track->bed2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2bigwig
{
}

=head2 bed2wig

bed2wig converter using 3rd party tools.

    $track->bed2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bed2wig
{
}

=head2 bedgraph2bigwig

bedgraph2bigwig converter using 3rd party tools.

    $track->bedgraph2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bedgraph2bigwig
{
}

=head2 bedgraph2wig

bedgraph2wig converter using 3rd party tools.

    $track->bedgraph2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bedgraph2wig
{
}

=head2 bigbed2bam

bigbed2bam converter using 3rd party tools.

    $track->bigbed2bam($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bam
{
}

=head2 bigbed2bed

bigbed2bed converter using 3rd party tools.

    $track->bigbed2bed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bed
{
}

=head2 bigbed2bedgraph

bigbed2bedgraph converter using 3rd party tools.

    $track->bigbed2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bedgraph
{
}

=head2 bigbed2bigwig

bigbed2bigwig converter using 3rd party tools.

    $track->bigbed2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigbed2bigwig
{
}

=head2 bigwig2bedgraph

bigwig2bedgraph converter using 3rd party tools.

    $track->bigwig2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigwig2bedgraph
{
}

=head2 bigwig2wig

bigwig2wig converter using 3rd party tools.

    $track->bigwig2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub bigwig2wig
{
}

=head2 sam2bedgraph

sam2bedgraph converter using 3rd party tools.

    $track->sam2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bedgraph
{
}

=head2 sam2bigbed

sam2bigbed converter using 3rd party tools.

    $track->sam2bigbed($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bigbed
{
}

=head2 sam2bigwig

sam2bigwig converter using 3rd party tools.

    $track->sam2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2bigwig
{
}

=head2 sam2wig

sam2wig converter using 3rd party tools.

    $track->sam2wig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub sam2wig
{
}

=head2 wig2bigwig

wig2bigwig converter using 3rd party tools.

    $track->wig2bigwig($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub wig2bigwig
{
}

=head2 wig2bedgraph

wig2bedgraph converter using 3rd party tools.

    $track->wig2bedgraph($input,$dir,$org,\%options);

The subroutine outputs the filename of the new track and its header when needed.

=cut

sub wig2bedgraph
{
}

=head2 check_constants

Check if required 3rd party tools exist for the required conversion. Internal use.

    $signaler->check_constants("sam","bigwig");

=cut

sub check_constants
{
    my $self = shift @_;
    my ($source,$dest) = @_;
    
    use v5.14;
    given($source."2".$dest)
    {
        when(/bam2bedgraph/i)
        {
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/bam2bigbed/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/bam2bigwig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when(/bam2wig/i)
        {
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when (/bed2bam/i)
        {
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when 
        (/bed2bedgraph/i)
        {
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when (/bed2bigbed/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when (/bed2bigwig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when (/bed2wig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/bedgraph2bigwig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when(/bedgraph2wig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when(/bigbed2bam/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/bigbed2bed/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        }
        when(/bigbed2bedgraph/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/bigbed2bigwig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when(/bigwig2wig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        }
        when(/bigwig2bedgraph/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        }
        when(/sam2bedgraph/i)
        {
            croak "SAM tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/sam2bigbed/i)
        {
            croak "SAM tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
        }
        when(/sam2bigwig/i)
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
        when(/sam2wig/i)
        {
            croak "SAM tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("SAMTOOLS_HOME") || (! -d $const->get("SAMTOOLS_HOME")));
            croak "BED tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("BEDTOOLS_HOME") || (! -d $const->get("BEDTOOLS_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
        when(/wig2bigwig/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
        }
        when(/wig2bedgraph/i)
        {
            croak "Kent tools not found in 3rd party tools! Please install them and retry!"
                if (!$const->get("KENTBIN_HOME") || (! -d $const->get("KENTBIN_HOME")));
            croak "At least one of Illumina iGenomes directory or HOST-USERNAME pair for accessing UCSC databases must be provided! Please check and retry!"
                if ((!$const->get("IGENOMES_HOME") || (! -d $const->get("IGENOMES_HOME"))) &&
                    (!$const->get("REMOTE_HOST") || !$const->get("REMOTE_USER")));
        }
    }
}

=head2 sort_bedstar

Helper sorting function for BED-like files (bed, bedgraph). Internal use.

    $track->sort_bedstar($file);
    
=cut

sub sort_bedstar
{
    my ($self,$infile) = @_;
    my $tmpdir = $self->get("tmpdir");
    my $tmpfile;
    
    if ($^O !~ /MSWin/) # Case of linux, easy sorting
    {
        $helper->disp("Sorting bed-like file $infile...");
        $tmpfile = File::Spec->catfile($tmpdir,"temp.in$$");
        `sort -k1,1 -k2g,2 $infile > $tmpfile `;
        $infile = $tmpfile;
    }
    else # We are in Windows... package required
    {
        $helper->try_module("File::Sort","sort_file");
        eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
        $helper->disp("Sorting file $infile...");
        $tmpfile = File::Spec->catfile($tmpdir,"temp.tmp");
        sort_file(
        {
            I => $infile,
            o => $tmpfile,
            k => ['1,1','2n,2'],
            t => "\t"
        });
        $infile = $tmpfile;
    }

    return($infile);
}

=head2 format_igenomes

Construct a string representing with the correct nomenclature, species for each supported iGenomes
chromosomal size annotation. Mostly for internal use.

    $track->format_igenomes($organism,$version);

=cut

sub format_igenomes_chrom_size
{
    use v5.14;
    my ($self,$org,$ver) = @_;
    my $base = $const->get("IGENOMES_HOME");
    given($org)
    {
        when(/human/)
        {
            given($ver)
            {
                when(/hg19/)
                {
                    return(File::Spec->catfile($base,"Homo_sapiens","UCSC","hg19","Annotation","Genes","ChromInfo.txt"));
                }
                when(/hg18/)
                {
                    return(File::Spec->catfile($base,"Homo_sapiens","UCSC","hg18","Annotation","Genes","ChromInfo.txt"));
                }
            }
        }
        when(/mouse/)
        {
            given($ver)
            {
                when(/mm10/)
                {
                    return(File::Spec->catfile($base,"Mus_musculus","UCSC","mm10","Annotation","Genes","ChromInfo.txt"));
                }
                when(/mm9/)
                {
                    return(File::Spec->catfile($base,"Mus_musculus","UCSC","mm9","Annotation","Genes","ChromInfo.txt"));
                }
            }
        }
        when(/rat/)
        {
            given($ver)
            {
                when(/rn5/)
                { 
                    return(File::Spec->catfile($base,"Rattus_norvegicus","UCSC","rn5","Annotation","Genes","ChromInfo.txt"));
                }
                when(/rn4/)
                { 
                    return(File::Spec->catfile($base,"Rattus_norvegicus","UCSC","rn4","Annotation","Genes","ChromInfo.txt"));
                }
            }
        }
        when(/fly/)
        {
            given($ver)
            {
                when(/dm3/)
                {
                    return(File::Spec->catfile($base,"Drosophila_melanogaster","UCSC","dm3","Annotation","Genes","ChromInfo.txt"));
                }
            }
        }
        when(/zebrafish/)
        {
            given($ver)
            {
                when(/danrer7/i)
                {
                    return(File::Spec->catfile($base,"Danio_rerio","UCSC","danRer7","Annotation","Genes","ChromInfo.txt"));
                }
            }
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
