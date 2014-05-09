=head1 NAME

HTS::Tools::Fetch - The great new HTS::Tools::Paramcheck!

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Parameter checker and validator for HTS::Tools

    use HTS::Tools::Paramcheck;

    my $checker = HTS::Tools::Check->new($tool,\%params);

    # Correct
    my $checker = HTS::Tools::Check->new("count",{"input" => "myfile.bed", "region" => "human-gene"});
    $checker->validate;
    # Blows up
    my $checker = HTS::Tools::Check->new("count",{"input" => "myfile.bed", "regions" => "human-gene"});
    $checker->validate;

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Paramcheck;

our $MODNAME = "HTS::Tools::Paramcheck";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Parameter checking module for module HTS::Tools.";


use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use File::Spec;
use File::Basename;

use HTS::Tools::Utils;

use vars qw($helper);

BEGIN {
    $helper = HTS::Tools::Utils->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

use constant MAXCORES => 12;

=head2 new($tool,$args)

The HTS::Tools::Paramcheck object constructor. See the SYNOPSIS for usage examples.

=cut

sub new
{
    my ($class,$args) = @_;
    my $self = {};

    # Pass global variables to the helper
    (defined($args->{"params"}->{"silent"})) ? ($helper->set("silent",$args->{"params"}->{"silent"})) :
        ($helper->set("silent",0));
    $helper->set_logger($args->{"params"}->{"log"}) if (defined($args->{"params"}->{"log"}));

    bless($self,$class);
    $self->init($args);
    return($self);
}

=head2 init($params)

HTS::Tools::Count initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
    my ($self,$args) = @_;
    my @accept = ("tool","params");
    foreach my $p (keys(%$args))
    {
        croak "Unknown parameter for initialization of HTS::Tools::Paramcheck: $p\n" if (!($p ~~ @accept));
        $self->set($p,$args->{$p});
    }
    return($self);
}

=head2 validate

The main parameter validator function of the module

=cut

sub validate
{
    use v5.14;
    my $self = shift @_;
    given($self->{"tool"})
    {
        when(/assign/i)
        {
            $self->validate_assign;
        }
        when(/constants/i)
        {
            $self->validate_constants;
        }
        when(/convert/i)
        {
            $self->validate_convert;
        }
        when(/count/i)
        {
            $self->validate_count;
        }
        when(/fetch/i)
        {
            $self->validate_fetch;
        }
        when(/intersect/i)
        {
            $self->validate_intersect;
        }
        when(/motifscan/i)
        {
            $self->validate_motifscan;
        }
        when(/multisect/i)
        {
            $self->validate_multisect;
        }
        when(/normalize/i)
        {
            $self->validate_normalize;
        }
        when(/profile/i)
        {
            $self->validate_profile;
        }
        when(/qc/i)
        {
            $self->validate_qc;
        }
        when(/queries/i)
        {
            $self->validate_queries;
        }
        when(/track/i)
        {
            $self->validate_track;
        }
    }
}

=head2 validate_assign

The parameter validator function of the HTS::Assign module. Do not use this directly, use the validate
function instead

=cut

sub validate_assign
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Assign"; 
    my $status;
    
    my @accept = ("input","region","background","span","idstrand","idmode","test","pvalue","outformat",
        "source","splicing","expression","gversion","log","silent","tmpdir");
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input query region file(s) ---\n" if (!@{$self->{"params"}->{"input"}});
    $stop .= "--- Please specify significant region file ---\n" if (!$self->{"params"}->{"region"});
    $stop .= "--- Please specify background region file ---\n" if (!$self->{"params"}->{"background"} && $self->{"params"}->{"test"} && $self->{"params"}->{"test"} ne "none");
    $stop .= "--- The supported genomes for the region parameter are organism-type, where organism is human, mouse, rat, fly or zebrafish and type is gene, exon, 5utr, 3utr or cds! Alternatively, it must be a file ---\n"
        if ( ! -f $self->{"params"}->{"region"}
            && $self->{"params"}->{"region"} !~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i);
    $stop .= "--- The supported genomes for the background parameter are organism-type, where organism is human, mouse, rat, fly or zebrafish and type is gene, exon, 5utr, 3utr or cds! Alternatively, it must be a file ---\n"
        if (($self->{"params"}->{"background"} && ! -f $self->{"params"}->{"background"}) && ($self->{"params"}->{"test"} && $self->{"params"}->{"test"} ne "none")
            && $self->{"params"}->{"background"} !~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i);
    if ($stop)
    {
        $helper->disp("$stop\n");
        croak "Type perldoc $modname for help in usage.\n\n";
        exit;
    }
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    # Check required packages
    if ($self->{"params"}->{"test"} && $self->{"params"}->{"test"} eq "chi2")
    {
        $status = eval { $helper->try_module("Math::Cephes") };
        if ($status)
        {
            $helper->disp("Module Math::Cephes is required to perform the chi-square test! Using default (none)...");
            $self->{"params"}->{"test"} = "none";
        }
        else { use Math::Cephes; }
    }
    if (defined($self->{"params"}->{"outformat"}) && @{$self->{"params"}->{"outformat"}} && @{$self->{"params"}->{"outformat"}} ~~ /matrix/)
    {
        $status = eval { $helper->try_module("Tie::IxHash::Easy") };
        if ($status)
        {
            $helper->disp("Module Tie::IxHash::Easy is required for one or more of the selected outputs! Using default (gff-peak)...");
            @{$self->{"params"}->{"outformat"}} = ("gff-peak");
        }
        else { use Tie::IxHash::Easy; }
    }
    
    # Check the rest
    # Check statistical test
    if ($self->{"params"}->{"test"} && $self->{"params"}->{"test"} ne "hypgeom" && $self->{"params"}->{"test"} ne "chi2" && 
        $self->{"params"}->{"test"} ne "none" && $self->{"params"}->{"test"} ne "auto")
    {
        $helper->disp("test parameter should be one of \"hypgeom\", \"chi2\", \"none\" or \"auto\"! Using default (none)...");
        $self->{"params"}->{"test"} = "none";
    }
    else # A test must be set
    {
        $self->{"params"}->{"test"} = "none";
    }
    # Check if span given
    if ((defined($self->{"params"}->{"span"}) && !@{$self->{"params"}->{"span"}}) || !defined($self->{"params"}->{"span"}))
    {
        $helper->disp("Search range from region start points (e.g. TSS) not given! Using defaults (-10kbp,10kbp)");
        @{$self->{"params"}->{"span"}} = (-10000,10000);
    }
    # Check if id and strand columns given for sig/back files
    if ((defined($self->{"params"}->{"idstrand"}) && !@{$self->{"params"}->{"idstrand"}}) || !defined($self->{"params"}->{"idstrand"}))
    {
        $helper->disp("Unique ID and strand columns for region and background files not given! Using defaults as from BED format (4,6)...");
        @{$self->{"params"}->{"idstrand"}} = (3,5);
    }
    else # Proper perl indexing
    {
        ${$self->{"params"}->{"idstrand"}}[0]--;
        ${$self->{"params"}->{"idstrand"}}[1]--;
    }
    # Check if id and mode columns given for peak files
    if ((defined($self->{"params"}->{"idmode"}) && !@{$self->{"params"}->{"idmode"}}) || !defined($self->{"params"}->{"idmode"}))
    {
        $helper->disp("Unique ID and mode columns for query region files not given! Using default ID as from BED format (4) and");
        $helper->disp("query modes will be assumed to be the centers of the input regions...");
        @{$self->{"params"}->{"idmode"}} = (3);
    }
    else # Proper perl indexing
    {
        ${$self->{"params"}->{"idmode"}}[0]--;
        if (!${$self->{"params"}->{"idmode"}}[1])
        {
            $helper->disp("Input query regions modes not given! The region centers will be used instead...");
        }
        else
        {
            ${$self->{"params"}->{"idmode"}}[1]--;
        }
        if (!${$self->{"params"}->{"idmode"}}[2])
        {
            $helper->disp("Additional query regions scores not given! They will not be reported...");
        }
        else
        {
            ${$self->{"params"}->{"idmode"}}[2]--;
        }
    }
    # Check if expression column is given in genes file and set proper Perl indexing
    if (defined($self->{"params"}->{"expression"}) && @{$self->{"params"}->{"expression"}})
    {
        my $l = scalar @{$self->{"params"}->{"expression"}};
        for (my $i=0; $i<$l; $i++)
        {
             ${$self->{"params"}->{"expression"}}[$i]--;
        }
    }
    # Check proper output format
    if (@{$self->{"params"}->{"outformat"}})
    {
        foreach my $c (@{$self->{"params"}->{"outformat"}})
        {
            if ($c ne "bed" && $c ne "stats" && $c ne "gff-peak" && $c ne "gff-gene" && $c ne "peak" &&  
                $c ne "gene" && $c ne "all-peak" && $c ne "all-gene" && $c ne "pretty-peak" && 
                $c ne "pretty-gene" && $c ne "gff-peak-db" && $c ne "gff-gene-db" && $c ne "peakdata"
                && $c ne "matrix-number" && $c ne "matrix-presence" && $c ne "matrix-peaks")
            {
                my $msg = "WARNING! outformat parameter options should be one or more of \"bed\", \"gff-peak\", \"gff-gene\",\n".
                          "\"peak\", \"gene\", \"all-peak\", \"all-gene\", \"pretty-peak\", \"pretty-gene\",\n".
                          "\"gff-peak-db\", \"gff-gene-db\", \"peakdata\", \"stats\", \"matrix-number\", \"matrix-presence\"".
                          "or \"matrix-peaks\"! \nUsing default (\"gff-peak\")...";
                $helper->disp($msg);
                @{$self->{"params"}->{"outformat"}} = ("gff-peak");
            }
        }
    }
    else
    {
        $helper->disp("Output format not given! Using default (gff-peak)");
        @{$self->{"params"}->{"outformat"}} = ("gff-peak");
    }
    if ($self->{"params"}->{"region"} =~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i)
    {
        my %sources = ("ucsc" => "UCSC","refseq" => "RefSeq","ensembl" => "Ensembl");
        my $source = $self->{"params"}->{"source"};
        my $splicing = $self->{"params"}->{"splicing"};
        my $gversion = $self->{"params"}->{"gversion"};
        if ($source)
        {
            $source = lc($source);
            if (grep {$_ eq $source} keys(%sources))
            {
                $helper->disp("Selected template regions source: ",$sources{$source});
            }
            else
            {
                $helper->disp("Source for template region files is not well-defined! Using default (ensembl)...");
                $self->{"params"}->{"source"} = "ensembl";
            }
        }
        else
        {
            $helper->disp("Source for template region files not given! Using default (ensembl)...");
            $self->{"params"}->{"source"} = "ensembl";
        }
        if ($splicing)
        {
            $splicing = lc($splicing);          
            if (grep {$_ eq $splicing} ("canonical","alternative"))
            {
                $helper->disp("Selected splicing for template regions source: ",$splicing);
            }
            else
            {
                $helper->disp("Splicing for template region files is not well-defined! Using default (canonical)...");
                $self->{"params"}->{"splicing"} = "canonical";
            }       
        }
        else
        {
            if ($source eq "ucsc" || $source eq "refseq")
            {
                $helper->disp("Splicing for template region files required but not defined! Using default (canonical)...");
                $self->{"params"}->{"splicing"} = "canonical";
            }
            else
            {
                $helper->disp("Splicing is not supported for Ensembl!");
                delete $self->{"params"}->{"splicing"};
            }
        }
        if ($gversion)
        {
            $gversion = lc($gversion);
            if ($self->{"params"}->{"region"} =~ m/human/i)
            {
                if ($gversion ne "hg18" && $gversion ne "hg19")
                {
                    $helper->disp("Human genome version wrong or not supported! Using default for human (hg19)...");
                    $self->{"params"}->{"gversion"} = "hg19";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/mouse/i)
            {
                if ($gversion ne "mm10" && $gversion ne "mm9")
                {
                    $helper->disp("Mouse genome version wrong or not supported! Using default for mouse (mm10)...");
                    $self->{"params"}->{"gversion"} = "mm10";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/rat/i)
            {
                if ($gversion ne "rn5")
                {
                    $helper->disp("Rat genome version wrong or not supported! Using default for rat (rn5)...");
                    $self->{"params"}->{"gversion"} = "rn5";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/fly/i)
            {
                if ($gversion ne "dm3")
                {
                    $helper->disp("Fruitfly genome version wrong or not supported! Using default for fruitfly (dm3)...");
                    $self->{"params"}->{"gversion"} = "dm3";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/zebrafish/i)
            {
                if ($gversion ne "danrer7")
                {
                    $helper->disp("Zebrafish genome version wrong or not supported! Using default for zebrafish (danRer7)...");
                    $self->{"params"}->{"gversion"} = "danRer7";
                }
            }
        }
        else
        {
            if ($self->{"params"}->{"region"} =~ m/human/i)
            {
                $helper->disp("Human genome version not given! Using default for human (hg19)...");
                $self->{"params"}->{"gversion"} = "hg19";
            }
            if ($self->{"params"}->{"region"} =~ m/mouse/i)
            {
                $helper->disp("Mouse genome version not given! Using default for mouse (mm10)...");
                $self->{"params"}->{"gversion"} = "mm10";
            }
            if ($self->{"params"}->{"region"} =~ m/rat/i)
            {
                $helper->disp("Rat genome version not given! Using default for rat (rn5)...");
                $self->{"params"}->{"gversion"} = "rn5";
            }
            if ($self->{"params"}->{"region"} =~ m/fly/i)
            {
                $helper->disp("Fruitfly genome version not given! Using default for fruitfly (dm3)...");
                $self->{"params"}->{"gversion"} = "dm3";
            }
            if ($self->{"params"}->{"region"} =~ m/zebrafish/i)
            {
                $helper->disp("Zebrafish genome version not given! Using default for zebrafish (danRer7)...");
                $self->{"params"}->{"gversion"} = "danRer7";
            }
        }
    }
    
    return($self->{"params"});
}

=head2 validate_constants

The parameter validator function of the HTS::Constants module. Do not use this directly, use the validate
function instead

=cut

sub validate_constants
{
    my ($self,$params) = @_;
    my $modname = "HTS::Tools::Constants";
    my $wrongs = 0; 

    my @accept = ();

    foreach my $c (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized constant : $c   --- Ignoring...") if (!($c ~~ @accept));
        $wrongs++;
    }
    
    if ($wrongs == scalar keys(%{$self->{"params"}}))
    {
        $helper->disp("No valid constant was given for $modname\n");
        croak "Type perldoc $modname for help in usage.\n\n";
    }
}

=head2 validate_convert

The parameter validator function of the HTS::Convert module. Do not use this directly, use the validate
function instead

=cut

sub validate_convert
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Convert";
    
    my @accept = ();
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    return($self->{"params"});
}

=head2 validate_count

The parameter validator function of the HTS::Count module. Do not use this directly, use the validate
function instead

=cut

sub validate_count
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Count";
    my $status;

    my @accept = ("input","region","sort","percent","lscore","escore","constant","small","split","nbins",
        "stats","output","ncore","source","gversion","splicing","keeporder","log","silent","tmpdir");

    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!$self->{"params"}->{"input"});
    $stop .= "--- Please specify region file ---\n" if (!$self->{"params"}->{"region"});
    $stop .= "--- The supported genomes are organism-type, where organism is human, mouse, rat, fly or zebrafish and type is gene, exon, 5utr, 3utr or cds! Alternatively, it must be a file ---\n"
        if (($self->{"params"}->{"region"} && ! -f $self->{"params"}->{"region"})
            && $self->{"params"}->{"region"} !~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i);
    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }

    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    # Check required packages
    $helper->try_module("IntervalTree");
    if ($self->{"params"}->{"keeporder"})
    {
        $status = eval { $helper->try_module("Tie::IxHash::Easy") };
        if ($status)
        {
            $helper->disp("Module Tie::IxHash::Easy is required for the keeporder parameter! Deactivating...");
            $self->{"params"}->{"keeporder"} = 0;
        }
        else { use Tie::IxHash::Easy; }
    }
    # Check the rest
    if (!$self->{"params"}->{"percent"} && !$self->{"params"}->{"lscore"} && !$self->{"params"}->{"escore"}) # If both given, use exponential scoring
    {
        $helper->disp("You did not define a partial overlap scheme! Using default (95% percent overlap)...");
        $self->{"params"}->{"percent"} = 0.95;
    }
    if ($self->{"params"}->{"lscore"} && $self->{"params"}->{"escore"}) # If both given, use exponential scoring
    {
        $helper->disp("You chose both linear and exponential scoring. Only exponential will be used...");
        $self->{"params"}->{"lscore"} = 0;
    }
    if ($self->{"params"}->{"stats"} && !($self->{"params"}->{"split"} || $self->{"params"}->{"nbins"}))
    {
        $helper->disp("You can't calculate area statistics without splitting to sub-areas or defining a number of bins! Option deactivated...");
        $self->{"params"}->{"stats"} = 0;
    }
    if ($self->{"params"}->{"ncore"})
    {
        $status = eval { $helper->try_module("Parallel::ForkManager") };
        if ($status)
        {
            $helper->disp("Module Parallel::ForkManager not found, proceeding with one core...");
            $self->{"params"}->{"ncore"} = 1;
        }
        else { use Parallel::ForkManager; }
        if ($self->{"params"}->{"ncore"} > MAXCORES)
        {
            my $c = MAXCORES;
            $helper->disp("The maximum number of cores allowed is $c...");
            $self->{"params"}->{"ncore"} = MAXCORES;
        }
    }
    else { $self->{"params"}->{"ncore"} = 1; }
    if ($self->{"params"}->{"output"})
    {
        if ($self->{"params"}->{"output"} eq "auto")
        {
            my ($base,$dir,$ext) = fileparse(${$self->{"params"}->{"input"}}[0],'\.[^.]*');
            $self->{"params"}->{"output"} = File::Spec->catfile($dir,$base."_REGIONCOUNTSTATS".$ext);
        }
    }
    if ($self->{"params"}->{"region"} =~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i)
    {
        my %sources = ("ucsc" => "UCSC","refseq" => "RefSeq","ensembl" => "Ensembl");
        my $source = $self->{"params"}->{"source"};
        my $splicing = $self->{"params"}->{"splicing"};
        my $gversion = $self->{"params"}->{"gversion"};
        if ($source)
        {
            $source = lc($source);
            if (grep {$_ eq $source} keys(%sources))
            {
                $helper->disp("Selected template regions source: ",$sources{$source});
            }
            else
            {
                $helper->disp("Source for template region files is not well-defined! Using default (ensembl)...");
                $self->{"params"}->{"source"} = "ensembl";
            }
        }
        else
        {
            $helper->disp("Source for template region files not given! Using default (ensembl)...");
            $self->{"params"}->{"source"} = "ensembl";
        }
        if ($splicing)
        {
            $splicing = lc($splicing);
            if ($self->{"params"}->{"source"} eq "ensembl")
            {
                $helper->disp("Splicing is not supported for Ensembl!");
                delete $self->{"params"}->{"splicing"};
            }
            else
            {
                if (grep {$_ eq $splicing} ("canonical","alternative"))
                {
                    $helper->disp("Selected splicing for template regions source: ",$splicing);
                }
                else
                {
                    $helper->disp("Splicing for template region files is not well-defined! Using default (canonical)...");
                    $self->{"params"}->{"splicing"} = "canonical";
                }
            }
        }
        else
        {
            if ($self->{"params"}->{"source"} eq "ucsc" || $self->{"params"}->{"source"} eq "refseq")
            {
                $helper->disp("Splicing for template region files required but not defined! Using default (canonical)...");
                $self->{"params"}->{"splicing"} = "canonical";
            }
        }
        if ($gversion)
        {
            $gversion = lc($gversion);
            if ($self->{"params"}->{"region"} =~ m/human/i)
            {
                if ($gversion ne "hg18" && $gversion ne "hg19")
                {
                    $helper->disp("Human genome version wrong or not supported! Using default for human (hg19)...");
                    $self->{"params"}->{"gversion"} = "hg19";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/mouse/i)
            {
                if ($gversion ne "mm10" && $gversion ne "mm9")
                {
                    $helper->disp("Mouse genome version wrong or not supported! Using default for mouse (mm10)...");
                    $self->{"params"}->{"gversion"} = "mm10";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/rat/i)
            {
                if ($gversion ne "rn5")
                {
                    $helper->disp("Rat genome version wrong or not supported! Using default for rat (rn5)...");
                    $self->{"params"}->{"gversion"} = "rn5";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/fly/i)
            {
                if ($gversion ne "dm3")
                {
                    $helper->disp("Fruitfly genome version wrong or not supported! Using default for fruitfly (dm3)...");
                    $self->{"params"}->{"gversion"} = "dm3";
                }
            }
            if ($self->{"params"}->{"region"} =~ m/zebrafish/i)
            {
                if ($gversion ne "danrer7")
                {
                    $helper->disp("Zebrafish genome version wrong or not supported! Using default for zebrafish (danRer7)...");
                    $self->{"params"}->{"gversion"} = "danRer7";
                }
            }
        }
        else
        {
            if ($self->{"params"}->{"region"} =~ m/human/i)
            {
                $helper->disp("Human genome version not given! Using default for human (hg19)...");
                $self->{"params"}->{"gversion"} = "hg19";
            }
            if ($self->{"params"}->{"region"} =~ m/mouse/i)
            {
                $helper->disp("Mouse genome version not given! Using default for mouse (mm10)...");
                $self->{"params"}->{"gversion"} = "mm10";
            }
            if ($self->{"params"}->{"region"} =~ m/rat/i)
            {
                $helper->disp("Rat genome version not given! Using default for rat (rn5)...");
                $self->{"params"}->{"gversion"} = "rn5";
            }
            if ($self->{"params"}->{"region"} =~ m/fly/i)
            {
                $helper->disp("Fruitfly genome version not given! Using default for fruitfly (dm3)...");
                $self->{"params"}->{"gversion"} = "dm3";
            }
            if ($self->{"params"}->{"region"} =~ m/zebrafish/i)
            {
                $helper->disp("Zebrafish genome version not given! Using default for zebrafish (danRer7)...");
                $self->{"params"}->{"gversion"} = "danRer7";
            }
        }
    }
    if ($self->{"params"}->{"nbins"})
    {
        $helper->disp("Number of genomic bins must be a positive integer! Using default (20)...")
            if ($self->{"params"}->{"nbins"} < 0 || $self->{"params"}->{"nbins"} !~ m/\d+/);
        $self->{"params"}->{"nbins"} = 20;
    }
    if ($self->{"params"}->{"nbins"} && $self->{"params"}->{"split"})
    {
        $helper->disp("nbins and split parameters are mutually exclusive! Using nbins...");
        delete $self->{"params"}->{"split"};
    }

    return($self->{"params"});
}

=head2 validate_fetch

The parameter validator function of the HTS::Fetch module. Do not use this directly, use the validate
function instead

=cut

sub validate_fetch
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Fetch";
    
    # Check required packages
    $helper->try_module("Tie::IxHash::Easy");

    my @accept = ("log","silent","tmpdir","output");
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }

    return($self->{"params"});
}

=head2 validate_intersect

The parameter validator function of the HTS::Intersect module. Do not use this directly, use the validate
function instead

=cut

sub validate_intersect
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Intersect";
    
    my @accept = ("inputA","inputB","sort","percent","any","extend","mode","autoextend","both","exact","keeporder","maxud",
            "reportonce","gap","output","multi","dryrun","waitbar","log","silent","tmpdir");
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!$self->{"params"}->{"inputA"} || !$self->{"params"}->{"inputB"});

    if (defined($self->{"params"}->{"percent"}) && @{$self->{"params"}->{"percent"}})
    {
        if (${$self->{"params"}->{"percent"}}[0] =~ /\d\:\d+/) 
        {
            my ($s,$e) = split(":",${$self->{"params"}->{"percent"}}[0]);
            @{$self->{"params"}->{"percent"}} = ($s..$e);
        }
        foreach my $cpp (@{$self->{"params"}->{"percent"}})
        {
            if ($cpp < 0 || $cpp > 100)
            {
                $stop .= "--- Overlap percentage should be a value between 0 and 100 ---\n";
                last;
            }
        }
    }
    if ($stop)
    {
        $helper->disp("$stop\n");
        croak "Type perldoc $modname for help in usage.\n\n";
        exit;
    }

    # At least the "any" parameter must be explicitly specified...
    if (!defined($self->{"params"}->{"any"}) || !$self->{"params"}->{"any"})
    {
        if (!defined($self->{"params"}->{"percent"}) || !@{$self->{"params"}->{"percent"}})
        {
            $self->{"params"}->{"any"} = 1;
        }
    }
    # Check gap
    disp("The gap parameter should be >=0! Using default (10000)...") if ($self->{"params"}->{"gap"} && $self->{"params"}->{"gap"} < 0);
    # Mode
    if ($self->{"params"}->{"mode"})
    {
        if ($self->{"params"}->{"mode"} < 3) # Can't be 0,1,2 as we are talking about bed-like file
        {
            $helper->disp("Invalid peak mode column: ".$self->{"params"}->{"mode"}." It will not be used...");
            delete $self->{"params"}->{"mode"};
        }
        else
        {
            $self->{"params"}->{"mode"} -= 4;
        }
    }
    # Check what is given on extend
    if ($self->{"params"}->{"extend"} && @{$self->{"params"}->{"extend"}})
    {
        if (!$self->{"params"}->{"mode"})
        {
            $helper->disp("Region mode/summit column must be provided with extend lengths! Assuming the middle of the region...");
            $self->{"params"}->{"mode"} = -1;
        }
        if (!${$self->{"params"}->{"extend"}}[1])
        {
            $helper->disp("The extend parameter has one argument... Assuming same extension in both sides.");
            ${$self->{"params"}->{"extend"}}[1] = ${$self->{"params"}->{"extend"}}[0];
        }
        if ($self->{"params"}->{"autoextend"})
        {
            $helper->disp("extend and autoextend options cannot be given together! Ignoring autoextend...");
            $self->{"params"}->{"autoextend"} = 0;
        }
    }
    else
    {
        if ($self->{"params"}->{"autoextend"} && !$self->{"params"}->{"mode"})
        {
            $helper->disp("Region mode/summit column must be provided with autoextend option! Assuming the middle of the region...");
            $self->{"params"}->{"mode"} = -1;
        }
    }
    if ($self->{"params"}->{"maxud"} && $self->{"params"}->{"maxud"} < 0)
    {
        $helper->disp("The maxud parameter must be greater than 0! Using default (5)...");
        $self->{"params"}->{"maxud"} = 5;
    }
    # Check proper output format
    if ($self->{"params"}->{"output"} && @{$self->{"params"}->{"output"}})
    {
        foreach my $c (@{$self->{"params"}->{"output"}})
        {
            if ($c ne "overlapA" && $c ne "overlapB" && $c ne "onlyA" && $c ne "onlyB" && $c ne "overpairs" &&
                $c ne "nonpairs")
            {
                my $msg = "WARNING! Output options should be one or more of \"overlapA\", \"overlapB\",".
                          " \"onlyA\", \"onlyB\", \"overpairs\" or \"nonpairs\",\n".
                          "Using default (\"overlapA\")...";
                $helper->disp($msg);
                @{$self->{"params"}->{"output"}} = ("overlapA");
            }
            if ($c eq "nonpairs" || $c eq "overpairs")
            {
                if (!$self->{"params"}->{"mode"})
                {
                    $helper->disp("Region mode/summit column must be provided with nonpairs and overpairs output options! Assuming the middle of the region...");
                    $self->{"params"}->{"mode"} = -1;
                }
            }
            if ($c eq "nonpairs")
            {
                if (!$self->{"params"}->{"gap"})
                {
                    $helper->disp("Acceptable gap for non-overlapping regions not specified! Using default (1000)...");
                    $self->{"params"}->{"gap"} = 1000;
                }
                if (!$self->{"params"}->{"maxud"})
                {
                    $helper->disp("Maximum number of non-overlapping regions to report not specified! Using default (5)...");
                    $self->{"params"}->{"maxud"} = 5;
                }
            }
        }
    }
    else
    {
        if (!${$self->{"params"}->{"percent"}}[1])
        {
            $helper->disp("Output file type not given! Using default (overlapA)...") ;
            @{$self->{"params"}->{"output"}} = ("overlapA");
        }
    }
    if ($self->{"params"}->{"keeporder"})
    {
        my $status = eval { $helper->try_module("Tie::IxHash::Easy") };
        if ($status)
        {
            $helper->disp("Module Tie::IxHash::Easy is required for the keeporder parameter! Deactivating...");
            $self->{"params"}->{"keeporder"} = 0;
        }
        else { require Tie::IxHash::Easy; }
    }
    
    return($self->{"params"});
}

=head2 validate_motifscan

The parameter validator function of the HTS::Motifscan module. Do not use this directly, use the validate
function instead

=cut

sub validate_motifscan
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Motifscan";
    
    my @accept = ("input","motif","background","scanner","center","colext","range","fpr","times","length",
            "besthit","output","uniquestats","justscan","log","silent","tmpdir");
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input sequence file(s) ---\n" if (!$self->{"params"}->{"input"} && !@{$self->{"params"}->{"input"}});
    $stop .= "--- Please specify motifs file ---\n" if (!$self->{"params"}->{"motif"});
    $stop .= "--- Please specify background sequences file ---\n" if (!$self->{"params"}->{"background"} && !$self->{"params"}->{"justscan"});

    if ($stop)
    {
        $helper->disp("$stop\n");
        croak "Type perldoc $modname for help in usage.\n\n";
        exit;
    }
    
    # Check the scanner
    if (defined($self->{"params"}->{"scanner"}) && $self->{"params"}->{"scanner"} !~ /pwmscan/i && $self->{"params"}->{"scanner"} !~ /MotifScanner/i)
    {
        $helper->disp("The scanner should be one of \"pwmscan\" or \"MotifScanner\". Using default (pwmscan)...");
        $self->{"params"}->{"scanner"} = "pwmscan";
    }
    else
    {
        $helper->disp("Motif canner not set! Using default (pwmscan)...");
        $self->{"params"}->{"scanner"} = "pwmscan";
    }
    # Check range - increase per real number, very simple expression, use with caution
    if (defined($self->{"params"}->{"range"}) && @{$self->{"params"}->{"range"}})
    {
        my @range = @{$self->{"params"}->{"range"}};
        if (@range == 1 && $range[0] =~ /\d+\:\d+(\.\d*)?\:\d+/)
        {
            my ($s,$inc,$e) = split(":",$range[0]);
            if (($s > 1 || $inc > 1 || $e > 1) && $self->{"params"}->{"scanner"} =~ /MotifScanner/i)
            {
                disp("Range for MotifScanner should be <1. Using default (0.2:0.01:0.7)...");
                @range = $helper->range_vector(0.2,0.7,0.01);
            }
            else 
            { 
                @range = $helper->range_vector($s,$e,$inc); 
            }
            $self->{"params"}->{"range"} = \@range;
        }
        elsif (@range == 1 && $range[0] =~ /\d+\:\d+/) # Increase per 1 if pwmscan or per 0.1 if MotifScanner
        {
            my ($s,$e) = split(":",$range[0]);
            use v5.10;
            given($self->{"params"}->{"scanner"})
            {
                when(/pwmscan/i)
                {
                    @range = ($s..$e);
                }
                when(/MotifScanner/i)
                {
                    if ($s > 1 || $e > 1)
                    {
                        $helper->disp("Range for MotifScanner should be <1. Using default (0.2:0.1:0.7)...");
                        @range = $helper->range_vector(0.5,0.01,0.01);
                    }
                    else 
                    { 
                        @range = $helper->range_vector($s,$e,0.1); 
                    }
                }
            }
            $self->{"params"}->{"range"} = \@range;
        }
    }
    else
    {
        my @range;
        if (!defined($self->{"params"}->{"justscan"}) || !$self->{"params"}->{"justscan"})
        {
            $helper->disp("Cutoff range not given. Using default (0..1)...");
            @range = $helper->range_vector(0,1,0.1);
        } 
        else 
        { 
            $range[0] = $self->{"params"}->{"justscan"}; 
        }
        $self->{"params"}->{"range"} = \@range;
    }
    # Check fpr
    if ((!defined($self->{"params"}->{"fpr"}) || !$self->{"params"}->{"fpr"}) && !$self->{"params"}->{"justscan"})
    {
        $helper->disp("FPR not defined! Using default (0.05)...");
        $self->{"params"}->{"fpr"} = 0.05;
    }
    elsif (defined($self->{"params"}->{"fpr"}) && ($self->{"params"}->{"fpr"} < 0 || $self->{"params"}->{"fpr"} > 1)
        && !$self->{"params"}->{"justscan"})
    {
        $helper->disp("FPR should be a number between 0 and 1! Using default (0.05)...");
        $self->{"params"}->{"fpr"} = 0.05;
    }
    if (defined($self->{"params"}->{"output"}) && @{$self->{"params"}->{"output"}})
    {
        foreach my $c (@{$self->{"params"}->{"output"}})
        {
            if ($c ne "gff" && $c ne "bed" && $c ne "stats" && $c ne "log")
            {
                my $msg = "WARNING! --output options should be one or more of \"gff\", \"bed\", \"stats\" or \"log\"\n".
                          "Using default (\"gff\")...";
                $helper->disp($msg);
                $self->{"params"}->{"output"} = ["gff"];
            }
        }
    }
    else
    {
        $helper->disp("Output file type(s) not defined! Using default (\"gff\")...");
        $self->{"params"}->{"output"} = ["gff"];
    }
    if (defined($self->{"params"}->{"colext"}) && @{$self->{"params"}->{"colext"}})
    {
        my @cntcol = @{$self->{"params"}->{"colext"}};
        my $l = @cntcol;
        if ($l != 3)
        {
            $helper->disp("ID, center columns and downstream extension not given properly... Using defaults (1,2,200)...");
            @cntcol = (0,1,200);
        }
        else
        {
            $cntcol[0]--;
            $cntcol[1]--;
        }
        $self->{"params"}->{"colext"} = \@cntcol;
    }
    if (defined($self->{"params"}->{"besthit"}) && $self->{"params"}->{"besthit"} < 0)
    {
        $helper->disp("The number of best hits should be a positive integer! Using default (1)...");
        $self->{"params"}->{"besthit"} = 1;
    }
    else
    {
        $helper->disp("Number of best hits not defined! Using default (1)...");
        $self->{"params"}->{"besthit"} = 1;
    }
    if (defined($self->{"params"}->{"times"}) && $self->{"params"}->{"times"} < 0)
    {
        $helper->disp("The times parameter must be a positive integer! Using default (10)...");
        $self->{"params"}->{"times"} = 10;
    }
    else
    {
        $helper->disp("times parameter not defined! Using default (10)...");
        $self->{"params"}->{"times"} = 10;
    }
    if (!defined($self->{"params"}->{"uniquestats"}))
    {
        $self->{"params"}->{"uniquestats"} = 0;
    }
    if (!defined($self->{"params"}->{"justscan"}))
    {
        $self->{"params"}->{"justscan"} = 0;
    }
    
    return($self->{"params"});
}

=head2 validate_normalize

The parameter validator function of the HTS::Normalize module. Do not use this directly, use the validate
function instead

=cut

sub validate_normalize
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Normalize";
    
    my @accept = ("type","input");
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!@{$self->{"params"}->{"input"}});
    $stop .= "--- Please specify file type ---\n" if (!$self->{"params"}->{"type"});
    $stop .= "--- File type must be one of \"bed\", \"bedgraph\" ---\n" # Other will be added at some point
        if ($self->{"params"}->{"type"} ne "bed" && $self->{"params"}->{"type"} ne "bedgraph");
    
    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }
    
    if ($self->{"params"}->{"type"} eq "bed")
    {
        $self->validate_normalize_bed;
    }
    elsif ($self->{"params"}->{"type"} eq "bedgraph")
    {
        $self->validate_normalize_bedgraph;
    }
}

=head2 validate_normalize_bed

The parameter validator function of the HTS::Normalize::Bed module. Do not use this directly, use the validate
function instead

=cut
sub validate_normalize_bed
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Normalize::Bed";
    
    my @accept = ("type","input","sumto","sort","savrem","sort","log","silent","tmpdir");
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!$self->{"params"}->{"input"});
    
    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }

    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }

    # Check normalization option
    if (!$self->{"params"}->{"sumto"})
    {
        $helper->disp("The normalization mode not given! Using default(\"generate\")...");
        $self->{"params"}->{"sumto"} = "generate";
    }
    if ($self->{"params"}->{"sumto"} && $self->{"params"}->{"sumto"} ne "generate" && $self->{"params"}->{"sumto"} ne "permute" &&
        $self->{"params"}->{"sumto"} !~ m/^[1-9]\d*$/)
    {
        $helper->disp("The normalization mode must be one of \"generate\", \"permute\" or a positive integer for randomly removing tags! Using \"generate\"");
        $self->{"params"}->{"sumto"} = "generate";
    }

    # Check required packages
    if ($self->{"params"}->{"sumto"} eq "generate" || $self->{"params"}->{"sumto"} eq "permute")
    {
        my $status = eval { $helper->try_module("Math::Random") };
        if ($status)
        {
            $helper->disp("Module Math::Random is required for the \"generate\" or \"permute\" normalization types! Deactivating and normalizing to 10M tags...");
            $self->{"params"}->{"sumto"} = 10000000;
        }
        else { require Math::Random; }
    }
    
    return($self->{"params"});
}

=head2 validate_normalize_bedgraph

The parameter validator function of the HTS::Normalize::Bedgraph module. Do not use this directly, use the validate
function instead

=cut
sub validate_normalize_bedgraph
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Normalize::Bedgraph";
    
    my @accept = ("input","type","output","extnorm","sumto","exportfactors","perlonly","prerun",
            "prerunlog","log","silent","tmpdir");
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!$self->{"params"}->{"input"});
    
    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }

    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    # Check signal summarization
    if (!$self->{"params"}->{"sumto"})
    {
        $helper->disp("The normalization total signal not given! Using default(1000000000)...");
        $self->{"params"}->{"sumto"} = 1000000000;
    }
    if ($self->{"params"}->{"sumto"} !~ m/^[1-9]\d*$/ || $self->{"params"}->{"sumto"} <= 0)
    {
        $helper->disp("The signal sum parameter must be a positive integer! Using default (1000000000)...");
        $self->{"params"}->{"sumto"} = 1000000000;
    }
    # Check presence of external normalization factors
    if (@{$self->{"params"}->{"extnorm"}})
    {
        my $e = @{$self->{"params"}->{"extnorm"}};
        my $b = @{$self->{"params"}->{"input"}};
        if ($e != $b)
        {
            $helper->disp("The number of external normalization factors given must be equal to the number of input files! Ignoring...");
            @{$self->{"params"}->{"extnorm"}} = ();
        }
        if (@{$self->{"params"}->{"extnorm"}} && $self->{"params"}->{"sumto"})
        {
            $helper->disp("Normalizing signal sum (sumto) and external normalization factors (extnorm) are mutually exclusive! Ignoring sumto...");
            $self->{"params"}->{"sumto"} = 0;
        }
    }
    # Check dry run and output files
    if (!$self->{"params"}->{"prerun"} && !$self->{"params"}->{"prerunlog"})
    {
        if (@{$self->{"params"}->{"output"}} && ${$self->{"params"}->{"output"}}->[0] ne "stdout")
        {
            my $o = @{$self->{"params"}->{"output"}};
            my $b = @{$self->{"params"}->{"input"}};
            if ($o != $b)
            {
                $helper->disp("The number of output files must be equal to the number of input files! Ignoring and autogenerating...");
                @{$self->{"params"}->{"output"}} = ();
            }
        }
        elsif (!@{$self->{"params"}->{"output"}})
        {
            $helper->disp("Output filenames will be autogenerated...");
            @{$self->{"params"}->{"output"}} = ();
        }
    }

    # Check if we are on Linux for the usage of awk
    if (!$self->{"params"}->{"perlonly"})
    {
        if ($^O =~ /MSWin/) # Windows... bad news...
        {
            $helper->disp("Windows OS detected! Switching to pure Perl for file streaming...");
            $self->{"params"}->{"perlonly"} = 1;
        }
    }

    return($self->{"params"});
}

=head2 validate_profile

The parameter validator function of the HTS::Profile module. Do not use this directly, use the validate
function instead

=cut

sub validate_profile
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Profile";
    
    my @accept = ();
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    return($self->{"params"});
}

=head2 validate_qc

The parameter validator function of the HTS::QC module. Do not use this directly, use the validate
function instead

=cut

sub validate_qc
{
    my $self = shift @_;
    my $modname = "HTS::Tools::QC";
    
    my @accept = ();
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }
    
    return($self->{"params"});
}

=head2 validate_query

The parameter validator function of the HTS::Queries module. Do not use this directly, use the validate
function instead

=cut

sub validate_queries
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Queries";
    my $query = $self->{"params"}->{"query"};

    my @accept = (
        "ucsc_canonical_genes",
        "ucsc_alternative_genes",
        "refseq_canonical_genes",
        "refseq_alternative_genes",
        "ucsc_canonical_exons",
        "ucsc_alternative_exons",
        "refseq_canonical_exons",
        "refseq_alternative_exons",
        "ucsc_canonical_5utr",
        "ucsc_alternative_5utr",
        "ucsc_canonical_3utr",
        "ucsc_alternative_3utr",
        "refseq_canonical_5utr",
        "refseq_alternative_5utr",
        "refseq_canonical_3utr",
        "refseq_alternative_3utr",
        "ucsc_canonical_cds",
        "ucsc_alternative_cds",
        "refseq_canonical_cds",
        "refseq_alternative_cds",
        "chrom_info"
    );

    if (!($query ~~ @accept))
    {
        $helper->disp("Unkown query for $modname: $query\n");
        croak "Type perldoc $modname for help in usage.\n\n";
    }
}

=head2 validate_track

The parameter validator function of the HTS::Track module. Do not use this directly, use the validate
function instead

=cut

sub validate_track
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Track";
    
    my @accept = ("type","input");
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!$self->{"params"}->{"input"});
    $stop .= "--- Please specify track set type ---\n" if (!$self->{"params"}->{"type"});
    $stop .= "--- File type must be one of \"signal\", \"hub\" ---\n"
        if (!$self->{"params"}->{"type"} || ($self->{"params"}->{"type"} ne "signal" && $self->{"params"}->{"type"} ne "hub"));
    
    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }
    
    if ($self->{"params"}->{"type"} eq "signal")
    {
        $self->validate_track_signal;
    }
    elsif ($self->{"params"}->{"type"} eq "hub")
    {
        $self->validate_track_hub;
    }
}

=head2 validate_track_signal

The parameter validator function of the HTS::Tools::Track::Signal module. Do not use this directly, use the validate
function instead

=cut

sub validate_track_signal
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Track::Signal";
    
    my @accept = ("type","input","source","destination","dir","urlbase","org","gversion","cleanlevel","sort","options");
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }

    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file ---\n" if (!$self->{"params"}->{"input"});
    $stop .= "--- Please specify source track type ---\n" if (!$self->{"params"}->{"source"});
    $stop .= "--- Please specify destination track type ---\n" if (!$self->{"params"}->{"destination"});
    $stop .= "--- Please specify big data url base ---\n"
        if (($self->{"params"}->{"destination"} eq "bigbed" || $self->{"params"}->{"destination"} eq "bigwig" ||
            $self->{"params"}->{"destination"} eq "bam" || $self->{"params"}->{"destination"} eq "sam") &&
            !$self->{"params"}->{"urlbase"});
    $stop .= "--- Please specify organism ---\n"
        if (($self->{"params"}->{"destination"} eq "bigbed" || $self->{"params"}->{"destination"} eq "bigwig" ||
            $self->{"params"}->{"destination"} eq "bedgraph" || $self->{"params"}->{"destination"} eq "wig") &&
            !$self->{"params"}->{"org"});
    $stop .= "--- Please specify genome version ---\n"
        if (($self->{"params"}->{"destination"} eq "bigbed" || $self->{"params"}->{"destination"} eq "bigwig" ||
            $self->{"params"}->{"destination"} eq "bedgraph" || $self->{"params"}->{"destination"} eq "wig") &&
            !$self->{"params"}->{"gversion"});
    $stop .= "--- Source track type must be one of \"bigbed\", \"bigwig\", \"wig\", \"bedgraph\", \"bed\", \"sam\" or \"bam\" ---\n"
        if ($self->{"params"}->{"source"} ne "bigbed" && $self->{"params"}->{"source"} ne "bigwig" &&
            $self->{"params"}->{"source"} ne "wig" && $self->{"params"}->{"source"} ne "bedgraph" &&
            $self->{"params"}->{"source"} ne "bed" && $self->{"params"}->{"source"} ne "bam" &&
            $self->{"params"}->{"source"} ne "sam");
    $stop .= "--- Destination track type must be one of \"bigbed\", \"bigwig\", \"wig\", \"bedgraph\", \"bed\", \"sam\" or \"bam\" ---\n"
        if ($self->{"params"}->{"destination"} ne "bigbed" && $self->{"params"}->{"destination"} ne "bigwig" &&
            $self->{"params"}->{"destination"} ne "wig" && $self->{"params"}->{"destination"} ne "bedgraph" &&
            $self->{"params"}->{"destination"} ne "bed" && $self->{"params"}->{"destination"} ne "bam" &&
            $self->{"params"}->{"source"} ne "sam");

    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }

    # We must have some defaults in this case...
    my %options;
    if ($self->{"params"}->{"destination"} eq "bigbed")
    {
        %options = (
            "name" => $self->{"params"}->{"input"},
            "description" => $self->{"params"}->{"input"},
            "color" => "0,0,160",
            "maxHeightPixels" => "128:64:16",
            "visibility" => "pack",
            "colorByStrand" => "\"255,0,0 0,0,255\"",
            "bigDataUrl" => $self->{"params"}->{"urlbase"}
        )
    }
    if ($self->{"params"}->{"destination"} eq "bigwig")
    {
        %options = (
            "name" => $self->{"params"}->{"input"},
            "description" => $self->{"params"}->{"input"},
            "color" => "0,0,160",
            "maxHeightPixels" => "128:64:16",
            "visibility" => "full",
            "bigDataUrl" => $self->{"params"}->{"urlbase"}
        )
    }
    if ($self->{"params"}->{"destination"} eq "bedgraph")
    {
        %options = (
            "name" => $self->{"params"}->{"input"},
            "description" => $self->{"params"}->{"input"},
            "color" => "0,0,160",
            "maxHeightPixels" => "128:64:16",
            "visibility" => "full"
        )
    }
    if ($self->{"params"}->{"destination"} eq "wig")
    {
        %options = (
            "name" => $self->{"params"}->{"input"},
            "description" => $self->{"params"}->{"input"},
            "color" => "0,0,160",
            "maxHeightPixels" => "128:64:16",
            "visibility" => "dense"
        )
    }
    if ($self->{"params"}->{"destination"} eq "bed")
    {
        %options = (
            "name" => $self->{"params"}->{"input"},
            "description" => $self->{"params"}->{"input"},
            "color" => "0,0,160",
            "colorByStrand" => "\"255,0,0 0,0,255\"",
            "maxHeightPixels" => "128:64:16",
            "visibility" => "pack"
        )
    }
    if ($self->{"params"}->{"destination"} eq "bam" || $self->{"params"}->{"destination"} eq "sam")
    {
        %options = (
            "name" => $self->{"params"}->{"input"},
            "description" => $self->{"params"}->{"input"},
            "color" => "0,0,160",
            "maxHeightPixels" => "128:64:16",
            "visibility" => "pack",
            "bigDataUrl" => $self->{"params"}->{"urlbase"},
            "bamColorMode" => "strand",
            "bamGrayMode" => "aliQual"
        )
    }
    #$self->{"params"}->{"options"} = \%options;

    # Define possible pairs
    my %pairs = (
        "bam2bed" => 1,
        "bam2bedgraph" => 1,
        "bam2bigbed" => 1,
        "bam2bigwig" => 1,
        "bam2wig" => 1,
        "bed2bam" => 1,
        "bed2bedgraph" => 1,
        "bed2bigbed" => 1,
        "bed2bigwig" => 1,
        "bed2wig" => 1,
        "bedgraph2bam" => 0,
        "bedgraph2bigbed" => 0,
        "bedgraph2bigwig" => 1,
        "bedgraph2wig" => 1,
        "bigbed2bam" => 1,
        "bigbed2bed" => 1,
        "bigbed2bedgraph" => 1,
        "bigbed2bigwig" => 1,
        "bigwig2bam" => 0,
        "bigwig2bedgraph" => 1,
        "bigwig2bigbed" => 0,
        "bigwig2wig" => 1,
        "wig2bigbed" => 0,
        "wig2bigwig" => 1,
        "wig2bedgraph" => 1,
        "wig2bam" => 0,
        "sam2bam" => 1,
        "sam2bedgraph" => 1,
        "sam2bigbed" => 1,
        "sam2bigwig" => 1,
        "sam2wig" => 1
    );

    if (!$self->{"params"}->{"dir"})
    {
        $helper->disp("No output directory specified! The input directory will be used...");
        my ($base,$dir,$ext) = fileparse($self->{"params"}->{"input"},'\.[^.]*');
        $self->{"params"}->{"dir"} = $dir;
    }
    if (!$self->{"params"}->{"cleanlevel"})
    {
        $helper->disp("No clean level specified! Using default (2: Removal of unlocalized/random chromosome regions and mitochondrial DNA)...");
        $self->{"params"}->{"cleanlevel"} = 2;
    }
    if (!$self->{"params"}->{"sort"})
    {
        $helper->disp("No co-rdinate sorting specified! If source track is not sorted, some conversions might crash...");
        $self->{"params"}->{"sort"} = 0;
    }
    else
    {
        $self->{"params"}->{"sort"} = 1;
    }
    if (!$self->{"params"}->{"options"})
    {
        $helper->disp("No additional options specified! The defaults will be used according to track types...");
        $self->{"params"}->{"options"} = \%options;
    }
    else
    {
        if (!$self->{"params"}->{"options"}->{"name"})
        {
            $helper->disp("No track name specified! The file name will be used...");
            $self->{"params"}->{"options"}->{"name"} = $self->{"params"}->{"input"};
        }
        if (!$self->{"params"}->{"options"}->{"description"})
        {
            $helper->disp("No track description specified! The track name will be used...");
            ($self->{"params"}->{"options"}->{"name"}) ?
            ($self->{"params"}->{"options"}->{"description"} = $self->{"params"}->{"options"}->{"name"}) :
            ($self->{"params"}->{"options"}->{"description"} = $self->{"params"}->{"input"});
        }
        if (!$self->{"params"}->{"options"}->{"color"})
        {
            $helper->disp("No track color specified! The default will be used (0,0,160)...");
            $self->{"params"}->{"options"}->{"color"} = "0,0,160";
        }
        if (!$self->{"params"}->{"options"}->{"maxHeightPixels"})
        {
            $helper->disp("No default track height area specified! The default used (128:64:16)...");
            $self->{"params"}->{"options"}->{"maxHeightPixels"} = "128:64:16";
        }
        if (!$self->{"params"}->{"options"}->{"bigDataUrl"} && ($self->{"params"}->{"destination"} eq "bigbed" ||
            $self->{"params"}->{"destination"} eq "bigwig" || $self->{"params"}->{"destination"} eq "bam" ||
            $self->{"params"}->{"destination"} eq "sam"))
        {
            $helper->disp("No different big data url specified! Using default from input...");
            $self->{"params"}->{"options"}->{"bigDataUrl"} = $self->{"params"}->{"urlbase"};
        }
    }

    # In every case, we add quotes to name and description
    $self->{"params"}->{"options"}->{"name"} = "\"".$self->{"params"}->{"options"}->{"name"}."\"";
    $self->{"params"}->{"options"}->{"description"} = "\"".$self->{"params"}->{"options"}->{"description"}."\"";

    my $gversion = $self->{"params"}->{"gversion"};
    if ($gversion)
    {
        $gversion = lc($gversion);
        if ($self->{"params"}->{"org"} =~ m/human/i)
        {
            if ($gversion ne "hg18" && $gversion ne "hg19")
            {
                $helper->disp("Human genome version wrong or not supported! Using default for human (hg19)...");
                $self->{"params"}->{"gversion"} = "hg19";
            }
        }
        if ($self->{"params"}->{"org"} =~ m/mouse/i)
        {
            if ($gversion ne "mm10" && $gversion ne "mm9")
            {
                $helper->disp("Mouse genome version wrong or not supported! Using default for mouse (mm10)...");
                $self->{"params"}->{"gversion"} = "mm10";
            }
        }
        if ($self->{"params"}->{"org"} =~ m/rat/i)
        {
            if ($gversion ne "rn5")
            {
                $helper->disp("Rat genome version wrong or not supported! Using default for rat (rn5)...");
                $self->{"params"}->{"gversion"} = "rn5";
            }
        }
        if ($self->{"params"}->{"org"} =~ m/fly/i)
        {
            if ($gversion ne "dm3")
            {
                $helper->disp("Fruitfly genome version wrong or not supported! Using default for fruitfly (dm3)...");
                $self->{"params"}->{"gversion"} = "dm3";
            }
        }
        if ($self->{"params"}->{"org"} =~ m/zebrafish/i)
        {
            if ($gversion ne "danrer7")
            {
                $helper->disp("Zebrafish genome version wrong or not supported! Using default for zebrafish (danRer7)...");
                $self->{"params"}->{"gversion"} = "danRer7";
            }
        }
    }
    else
    {
        if ($self->{"params"}->{"org"} && $self->{"params"}->{"org"} =~ m/human/i)
        {
            $helper->disp("Human genome version not given! Using default for human (hg19)...");
            $self->{"params"}->{"gversion"} = "hg19";
        }
        if ($self->{"params"}->{"org"} && $self->{"params"}->{"org"} =~ m/mouse/i)
        {
            $helper->disp("Mouse genome version not given! Using default for mouse (mm10)...");
            $self->{"params"}->{"gversion"} = "mm10";
        }
        if ($self->{"params"}->{"org"} && $self->{"params"}->{"org"} =~ m/rat/i)
        {
            $helper->disp("Rat genome version not given! Using default for rat (rn5)...");
            $self->{"params"}->{"gversion"} = "rn5";
        }
        if ($self->{"params"}->{"org"} && $self->{"params"}->{"org"} =~ m/fly/i)
        {
            $helper->disp("Fruitfly genome version not given! Using default for fruitfly (dm3)...");
            $self->{"params"}->{"gversion"} = "dm3";
        }
        if ($self->{"params"}->{"org"} && $self->{"params"}->{"org"} =~ m/zebrafish/i)
        {
            $helper->disp("Zebrafish genome version not given! Using default for zebrafish (danRer7)...");
            $self->{"params"}->{"gversion"} = "danRer7";
        }
    }

    # Check combination
    my $source = $self->{"params"}->{"source"};
    my $dest = $self->{"params"}->{"destination"};
    my $conv = $source."2".$dest;
    croak "The requested conversion from $source to $dest is not possible! Please select another conversion..."
        if (!$pairs{$conv});

    return($self->{"params"});
}

=head2 validate_track_hub

The parameter validator function of the HTS::Tools::Track::Hub module. Do not use this directly, use the validate
function instead

=cut

sub validate_track_hub
{
    my $self = shift @_;
    my $modname = "HTS::Tools::Track::Hub";
    
    my @accept = ("type","config","hubid","hubname","hubdesc","hubgenomes","hubmail","hubbase","tracks");
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$self->{"params"}}))
    {
        $helper->disp("Unrecognized parameter : $p   --- Ignoring...") if (!($p ~~ @accept));
    }

    # Check fatal
    my $stop;
    $stop .= "--- Please specify configuration YAML or the rest of required parameters ---\n"
        if (!$self->{"params"}->{"config"} && !$self->{"params"}->{"hubid"} && !$self->{"params"}->{"hubmail"}
            && !$self->{"params"}->{"hubbase"} && !$self->{"params"}->{"tracks"});
    $stop .= "--- Please specify track hub unique identification ---\n"
        if (!$self->{"params"}->{"config"} && !$self->{"params"}->{"hubid"});
    $stop .= "--- Please specify track hub administrator e-mail ---\n"
        if (!$self->{"params"}->{"config"} && !$self->{"params"}->{"hubmail"});
    $stop .= "--- Please specify track hub public html base directory ---\n"
        if (!$self->{"params"}->{"config"} && !$self->{"params"}->{"hubbase"});
    $stop .= "--- Please specify track hub tracks file ---\n"
        if (!$self->{"params"}->{"config"} && !$self->{"params"}->{"tracks"});

    if ($stop)
    {
        $helper->disp("$stop\n");
        $helper->disp("Type perldoc $modname for help in usage.\n\n");
        exit;
    }

    if ($self->{"params"}->{"config"} && ($self->{"params"}->{"hubid"} || !$self->{"params"}->{"hubmail"}
        || !$self->{"params"}->{"hubbase"} || !$self->{"params"}->{"tracks"}))
    {
        $helper->disp("Both YAML configuration file and manual parameters given! Will use only YAML file...");
        delete($self->{"params"}->{"hubid"}) if ($self->{"params"}->{"hubid"});
        delete($self->{"params"}->{"hubname"}) if ($self->{"params"}->{"hubname"});
        delete($self->{"params"}->{"hubdesc"}) if ($self->{"params"}->{"hubdesc"});
        delete($self->{"params"}->{"hubbase"}) if ($self->{"params"}->{"hubbase"});
        delete($self->{"params"}->{"hubmail"}) if ($self->{"params"}->{"hubmail"});
    }

    my $status;
    if ($self->{"params"}->{"config"})
    {
        $status = eval { $helper->try_module("YAML") };
        if ($status)
        {
            croak "Module YAML is required to continue with the track hub creation! Otherwise, specify the run parameters individually...";
        }
        else { use YAML; }
    }

    my $testval;
    if (!$self->{"params"}->{"config"})
    {
        if ($self->{"params"}->{"hubid"})
        {
            $testval = $self->{"params"}->{"hubid"};
            if ($testval =~ m/\s/g)
            {
                $helper->disp("The hub unique identifier cannot contain spaces! Collapsing...");
                $testval =~ s/\s/_/g;
                $self->{"params"}->{"hubid"} = $testval;
            }
        }
        $self->{"params"}->{"hubname"} = $self->{"params"}->{"hubid"} if (!$self->{"params"}->{"hubid"});
        $self->{"params"}->{"hubdesc"} = $self->{"params"}->{"hubid"} if (!$self->{"params"}->{"hubdesc"});
    }
    
    return($self->{"params"});
}

=head2 get

HTS::Tools::Paramcheck object getter

    my $param_value = $helper->get("param_name")
=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Paramcheck object setter

    $helper->set("param_name","param_value")
    
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

    perldoc HTS::Tools::Fetch


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

1; # End of HTS::Tools::Paramcheck
