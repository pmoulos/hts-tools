=head1 NAME

HTS::Tools::Fetch - Annotation downloader for several widely used genomes and genomic features.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

HTS::Tools::Fetch is the module which serves the automatic downloading and formating of the genomic
region templates mentioned in the SYNOPSIS of HTS::Tools::Count module. See that part for an explanation
of the module possibilities. This module can be used independently to fetch these regions for additional
applications. It also serves the module HTS::Tools::Profile.

Perhaps a little code snippet.

    use HTS::Tools::Fetch;
    my $fetcher = HTS::Tools::Fetch->new();
    my $fetcher = HTS::Tools::Fetch->new({'tmpdir' => 'my_output_dir'});
    $fetcher->fetch_ensembl_genes('human');
    $fetcher->fetch_ucsc_exons('mouse','canonical','ucsc');
    $fetcher->fetch_ucsc_utr('rat','alternative','refseq',5);

Note that the 'tmpdir' parameter must be set if this module is called independently so as to store the
region templates that will be fetched, otherwise, they are written in a temporary directory which is
destroyed after the module has finished. The 'silent' option is also available. See the SYNOPSIS of
HTS::Tools::Count.

=head1 EXPORT

Functions that can be exported:

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Fetch;

our $MODNAME = "HTS::Tools::Fetch";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Annotation download module for module HTS::Tools.";

#################### SPECIAL SECTION ####################
sub cosort
{
    my @one = split("-",$a);
    my @two = split("-",$b);
    return 1 if ($one[0] > $two[0]);
    return -1 if ($one[0] < $two[0]);
    if ($one[0] == $two[0])
    {
        return 1 if ($one[1] > $two[1]);
        return 0 if ($one[1] == $two[1]);
        return -1 if ($one[1] < $two[1]);
    }
}
################## END SPECIAL SECTION ##################

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use DBI;
use File::Basename;
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp;
use HTTP::Request;
use LWP::UserAgent;

use HTS::Tools::Constants;
use HTS::Tools::Queries;
use HTS::Tools::Utils;

use vars qw($helper $qobj $const);

#use constant BIOMART_PATH => "http://www.biomart.org/biomart/martservice?";
#use constant REMOTE_HOST => "genome-mysql.cse.ucsc.edu";
#use constant REMOTE_USER => "genome";

BEGIN {
    $helper = HTS::Tools::Utils->new();
    $qobj = HTS::Tools::Queries->new();
    $const = HTS::Tools::Constants->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

use constant BIOMART_PATH => $const->get("BIOMART_PATH");
use constant REMOTE_HOST => $const->get("REMOTE_HOST");
use constant REMOTE_USER => $const->get("REMOTE_USER");

# Make sure output is unbuffered
select(STDOUT);
$|=1;

=head2 new

Constructor of HTS::Tools::Fetch module

    use HTS::Tools::Fetch;
    my $fetcher = HTS::Tools::Fetch->new($params);

See the SYNOPSIS for possible parameters.

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
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "fetch","params" => $params});
    $params = $checker->validate;
    
    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init($params)

HTS::Tools::Fetch object initialization method. NEVER use this directly, use new instead.

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

=head2 fetch_ensembl_genes

Fetch a BED6+1 headed file with the latest Ensembl genes annotation for one of the suported organisms.
Mostly for internal use but it can easily be used for fetching annotation too.

    $fetcher->fetch_ensembl_genes($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count.

=cut

sub fetch_ensembl_genes
{
    my ($self,$org) = @_;
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $sp = $self->format_species($org,"","ensembl");
    my $xml = $self->get_xml_genes_query($sp);
    my $path = BIOMART_PATH;
    my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
    my $ua = LWP::UserAgent->new;
    my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my $response;

    # We need to put data in a temporary file because it's scrambled by asynchronicity
    $helper->disp("Querying Biomart...");
    $ua->request($request,sub {
        my ($data,$response) = @_;
        if ($response->is_success) {
            print $tmpfh "$data";
        }
        else {
            warn ("Problems with the web server: ".$response->status_line);
        }
    },1000);

    seek($tmpfh,0,SEEK_SET);
    my %strand = ( 1 => "+", -1 => "-" );
    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\tgene_id\tgc_content\tstrand\tgene_name\tbiotype\n";
    while (my $line = <$tmpfh>)
    {
        next if ($line !~ m/^[0-9XY]/);
        $line =~ s/\r|\n$//g;
        my @cols = split(/\t/,$line);
        print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\t$cols[7]\n";
    }
    close(REGS);
    close($tmpfh);

    return($regs);
}

=head2 fetch_ensembl_exons

Fetch a BED6+1 headed file with the latest Ensembl exons annotation for one of the suported organisms.
Mostly for internal use but it can easily be used for fetching annotation too.

    $fetcher->fetch_ensembl_exons($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count.

=cut

sub fetch_ensembl_exons
{
    my ($self,$org) = @_;
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $sp = $self->format_species($org,"","ensembl");
    my $xml = $self->get_xml_exons_query($sp);
    my $path = BIOMART_PATH;
    my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
    my $ua = LWP::UserAgent->new;
    my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my $response;

    # We need to put data in a temporary file because it's scrambled by asynchronicity
    $helper->disp("Querying Biomart...");
    $ua->request($request,sub {   
        my ($data,$response) = @_;
        if ($response->is_success) {
            print $tmpfh "$data";
        }
        else {
            warn ("Problems with the web server: ".$response->status_line);
        }
    },1000);

    seek($tmpfh,0,SEEK_SET);
    my %strand = ( 1 => "+", -1 => "-" );
    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\texon_id\tgene_id\tstrand\tgene_name\tbiotype\n";
    while (my $line = <$tmpfh>)
    {
        next if ($line !~ m/^[0-9XY]/);
        $line =~ s/\r|\n$//g;
        my @cols = split(/\t/,$line);
        print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\t$cols[7]\n";
    }
    close(REGS);
    close($tmpfh);

    return($regs);
}

=head2 fetch_ensembl_utr

Fetch a BED6+1 headed file with the latest Ensembl UTR annotation for one of the suported organisms.
Mostly for internal use but it can easily be used for fetching annotation too.

    $fetcher->fetch_ensembl_utr($organism,$utr);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $utr can be 5 or 3.

=cut

sub fetch_ensembl_utr
{
    my ($self,$org,$utr) = @_;
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $sp = $self->format_species($org,"","ensembl");
    my $xml = $self->get_xml_utr_query($sp,$utr);
    my $path = BIOMART_PATH;
    my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
    my $ua = LWP::UserAgent->new;
    my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my $response;

    # We need to put data in a temporary file because it's scrambled by asynchronicity
    $helper->disp("Querying Biomart...");
    $ua->request($request,sub {
        my ($data,$response) = @_;
        if ($response->is_success) {
            print $tmpfh "$data";
        }
        else {
            warn ("Problems with the web server: ".$response->status_line);
        }
    },1000);

    seek($tmpfh,0,SEEK_SET);
    my %strand = ( 1 => "+", -1 => "-" );
    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\tgene_id\ttranscript_count\tstrand\tgene_name\tbiotype\n";
    while (my $line = <$tmpfh>)
    {
        next if ($line !~ m/^[0-9XY]/);
        $line =~ s/\r|\n$//g;
        my @cols = split(/\t/,$line);
        next if (!$cols[1]);
        print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\t$cols[7]\n";
    }
    close(REGS);
    close($tmpfh);

    return($regs);
}

=head2 fetch_ensembl_cds

Fetch a BED6+1 headed file with the latest Ensembl CDS annotation for one of the suported organisms.
Mostly for internal use but it can easily be used for fetching annotation too.

    $fetcher->fetch_ensembl_cds($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count.

=cut

sub fetch_ensembl_cds
{
    my ($self,$org) = @_;
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $sp = $self->format_species($org,"","ensembl");
    my $xml = $self->get_xml_cds_query($sp);
    my $path = BIOMART_PATH;
    my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
    my $ua = LWP::UserAgent->new;
    my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my $response;

    # We need to put data in a temporary file because it's scrambled by asynchronicity
    $helper->disp("Querying Biomart...");
    $ua->request($request,sub {   
        my ($data,$response) = @_;
        if ($response->is_success) {
            print $tmpfh "$data";
        }
        else {
            warn ("Problems with the web server: ".$response->status_line);
        }
    },1000);

    seek($tmpfh,0,SEEK_SET);
    my %strand = ( 1 => "+", -1 => "-" );
    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\texon_id\tgene_id\tstrand\tgene_name\tbiotype\n";
    while (my $line = <$tmpfh>)
    {
        next if ($line !~ m/^[0-9XY]/);
        $line =~ s/\r|\n$//g;
        my @cols = split(/\t/,$line);
        next if (!$cols[1]);
        print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\t$cols[7]\n";
    }
    close(REGS);
    close($tmpfh);

    return($regs);
}

=head2 fetch_ucsc_genes

Fetch a BED6+1 headed file with the UCSC genes annotation for one of the suported organisms, by connecting
either to UCSC or a local version of the UCSC Genome Browser database. Mostly for internal use but it 
can easily be used for fetching annotation too.

    $fetcher->fetch_ucsc_genes($organism,$splicing_type,$db_source);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $splicing_type can be 'canonical'
or 'alternative', $db_source can be 'ucsc' or 'refseq'.

=cut

sub fetch_ucsc_genes
{
    my ($self,$org,$ver,$type,$source) = @_;
    my $sp = $self->format_species($org,$ver,$source);
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my ($q,$data);

    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\tucsc_id\texon_count\tstrand\tgene_name\n";

    if ($source eq "ucsc")
    {
        if ($type eq "canonical")
        {
            $helper->disp("Querying UCSC database for UCSC canonical genes...");
            $q = $qobj->get_query("ucsc_canonical_genes");
        }
        elsif ($type eq "alternative")
        {
            $helper->disp("Querying UCSC database for UCSC alternative genes...");
            $q = $qobj->get_query("ucsc_alternative_genes");
        }
    }
    elsif ($source eq "refseq")
    {
        if ($type eq "canonical")
        {
            $helper->disp("Querying UCSC database for RefSeq canonical genes...");
            $q = $qobj->get_query("refseq_canonical_genes");
        }
        elsif ($type eq "alternative")
        {
            $helper->disp("Querying UCSC database for RefSeq alternative genes...");
            $q = $qobj->get_query("refseq_alternative_genes");
        }
    }
    my $conn = $self->open_connection($sp);
    my $sth = $conn->prepare($q);
    $sth->execute();
    while ($data = $sth->fetchrow_hashref())
    {
        print REGS $data->{"chrom"}."\t".$data->{"chromStart"}."\t".$data->{"chromEnd"}."\t".
            $data->{"transcript"}."\t".$data->{"exonCount"}."\t".$data->{"strand"}."\t".$data->{"geneName"}."\n";
    }
    $self->close_connection($conn);
    close(REGS);
    return($regs);
}

=head2 fetch_ucsc_exons

Fetch a BED6+1 headed file with the UCSC exons annotation for one of the suported organisms, by connecting
either to UCSC or a local version of the UCSC Genome Browser database. Mostly for internal use but it 
can easily be used for fetching annotation too.

    $fetcher->fetch_ucsc_exons($organism,$splicing_type,$db_source);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $splicing_type can be 'canonical'
or 'alternative', $db_source can be 'ucsc' or 'refseq'.

=cut

sub fetch_ucsc_exons
{
    my ($self,$org,$ver,$type,$source) = @_;
    my $sp = $self->format_species($org,$ver,$source);
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my ($q,$data,$i);
    my (@starts,@ends,@coords,@ses);

    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\texon_id\texon_count\tstrand\ttranscript_name\n";

    if ($source eq "ucsc")
    {
        if ($type eq "canonical")
        {
            $helper->disp("Querying UCSC database for UCSC canonical exons...");
            $q = $qobj->get_query("ucsc_canonical_exons");
        }
        elsif ($type eq "alternative")
        {
            $helper->disp("Querying UCSC database for UCSC alternative exnos...");
            $q = $qobj->get_query("ucsc_alternative_exons");
        }
    }
    elsif ($source eq "refseq")
    {
        if ($type eq "canonical")
        {
            $helper->disp("Querying UCSC database for RefSeq canonical exons...");
            $q = $qobj->get_query("refseq_canonical_exons");
        }
        elsif ($type eq "alternative")
        {
            $helper->disp("Querying UCSC database for RefSeq alternative exons...");
            $q = $qobj->get_query("refseq_alternative_exons");
        }
    }
    my $conn = $self->open_connection($sp);
    my $sth = $conn->prepare($q);
    $sth->execute();
    while ($data = $sth->fetchrow_hashref())
    {
        @starts = split(",",$data->{"exonStarts"});
        @ends = split(",",$data->{"exonEnds"});
        for ($i=0; $i<@starts; $i++)
        {
            push(@coords,$starts[$i]."-".$ends[$i]);
        }
        my %u = $helper->unique(@coords);
        @coords = sort cosort keys(%u);
        foreach my $k (@coords)
        {
            @ses = split("-",$k);
            print REGS $data->{"chrom"}."\t".$ses[0]."\t".$ses[1]."\t".$data->{"name"}."\t".
            $data->{"exonCount"}."\t".$data->{"strand"}."\t".$data->{"transcript"}."\n";
        }
    }
    $self->close_connection($conn);
    close(REGS);
    return($regs);
}

=head2 fetch_ucsc_utr

Fetch a BED6+1 headed file with the UCSC UTR annotation for one of the suported organisms, by connecting
either to UCSC or a local version of the UCSC Genome Browser database. Mostly for internal use but it 
can easily be used for fetching annotation too.

    $fetcher->fetch_ucsc_utr($organism,$splicing_type,$db_source,$utr);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $splicing_type can be 'canonical'
or 'alternative', $db_source can be 'ucsc' or 'refseq'. $utr can be 5 or 3.

=cut

sub fetch_ucsc_utr
{
    my ($self,$org,$ver,$type,$source,$utr) = @_;
    my $sp = $self->format_species($org,$ver,$source);
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my ($q,$data);

    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\tucsc_id\texon_count\tstrand\tgene_name\n";

    if ($source eq "ucsc")
    {
        if ($type eq "canonical")
        {
            if ($utr == 5)
            {
                $helper->disp("Querying UCSC database for UCSC canonical 5'UTRs...");
                $q = $qobj->get_query("ucsc_canonical_5utr");
            }
            elsif ($utr == 3)
            {
                $helper->disp("Querying UCSC database for UCSC canonical 3'UTRs...");
                $q = $qobj->get_query("ucsc_canonical_3utr");
            }
        }
        elsif ($type == "alternative")
        {
            if ($utr == 5)
            {
                $helper->disp("Querying UCSC database for UCSC alternative 5'UTRs...");
                $q = $qobj->get_query("ucsc_alternative_5utr");
            }
            elsif ($utr == 3)
            {
                $helper->disp("Querying UCSC database for UCSC alternative 3'UTRs...");
                $q = $qobj->get_query("ucsc_alternative_3utr");
            }
        }
    }
    elsif ($source eq "refseq")
    {
        if ($type eq "canonical")
        {
            if ($utr == 5)
            {
                $helper->disp("Querying UCSC database for RefSeq canonical 5'UTRs...");
                $q = $qobj->get_query("refseq_canonical_5utr");
            }
            elsif ($utr == 3)
            {
                $helper->disp("Querying UCSC database for RefSeq canonical 3'UTRs...");
                $q = $qobj->get_query("refseq_canonical_3utr");
            }
        }
        elsif ($type == "alternative")
        {
            if ($utr == 5)
            {
                $helper->disp("Querying UCSC database for RefSeq alternative 5'UTRs...");
                $q = $qobj->get_query("refseq_alternative_5utr");
            }
            elsif ($utr == 3)
            {
                $helper->disp("Querying UCSC database for RefSeq alternative 3'UTRs...");
                $q = $qobj->get_query("refseq_alternative_3utr");
            }
        }
    }
    
    my $conn = $self->open_connection($sp);
    my $sth = $conn->prepare($q);
    $sth->execute();
    while ($data = $sth->fetchrow_hashref())
    {
        print REGS $data->{"chrom"}."\t".$data->{"start"}."\t".$data->{"end"}."\t".
            $data->{"name"}."\t".$data->{"exonCount"}."\t".$data->{"strand"}."\t".$data->{"geneName"}."\n";
    }
    $self->close_connection($conn);
    close(REGS);
    return($regs);
}

=head2 fetch_ucsc_cds

Fetch a BED6+1 headed file with the UCSC CDS annotation for one of the suported organisms, by connecting
either to UCSC or a local version of the UCSC Genome Browser database. Mostly for internal use but it 
can easily be used for fetching annotation too.

    $fetcher->fetch_ucsc_genes($organism,$splicing_type,$db_source);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $splicing_type can be 'canonical'
or 'alternative', $db_source can be 'ucsc' or 'refseq'.

=cut

sub fetch_ucsc_cds
{
    my ($self,$org,$ver,$type,$source) = @_;
    my $sp = $self->format_species($org,$ver,$source);
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $regs = File::Spec->catfile($tmpdir,"regions.txt");
    my ($q,$data);

    open(REGS,">$regs");
    print REGS "chromosome\tstart\tend\tucsc_id\texon_count\tstrand\tgene_name\n";

    if ($source eq "ucsc")
    {
        if ($type eq "canonical")
        {
            disp("Querying UCSC database for UCSC canonical CDS...");
            $q = $qobj->get_query("ucsc_canonical_cds");
        }
        elsif ($type eq "alternative")
        {
            disp("Querying UCSC database for UCSC alternative CDS...");
            $q = $qobj->get_query("ucsc_alternative_cds");
        }
    }
    elsif ($source eq "refseq")
    {
        if ($type eq "canonical")
        {
            disp("Querying UCSC database for RefSeq canonical CDS...");
            $q = $qobj->get_query("refseq_canonical_cds");
        }
        elsif ($type eq "alternative")
        {
            disp("Querying UCSC database for RefSeq alternative CDS...");
            $q = $qobj->get_query("refseq_alternative_cds");
        }
    }
    my $conn = $self->open_connection($sp);
    my $sth = $conn->prepare($q);
    $sth->execute();
    while ($data = $sth->fetchrow_hashref())
    {
        print REGS $data->{"chrom"}."\t".$data->{"chromStart"}."\t".$data->{"chromEnd"}."\t".
            $data->{"transcript"}."\t".$data->{"exonCount"}."\t".$data->{"strand"}."\t".$data->{"geneName"}."\n";
    }
    $self->close_connection($conn);
    close(REGS);
    return($regs);
}

=head2 fetch_chrom_info

Fetch a 2-column file with chromosome name and length for one of the suported organisms, by connecting
either to UCSC or a local version of the UCSC Genome Browser database. Mostly for internal use but it 
can easily be used for fetching chromosomal information too.

    $fetcher->fetch_chrom_info($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count.

=cut

sub fetch_chrom_info
{
    my ($self,$org,$ver) = @_;
    my $sp = $self->format_species($org,$ver,"ucsc");
    print "\n\n$org\n\n$ver\n\n$sp\n\n";
    my $tmpdir = (defined($self->get("output"))) ? ($self->get("output")) : ($self->get("tmpdir"));
    my $chrinfo = File::Spec->catfile($tmpdir,"chrom_info_$org.txt");
    my ($q,$data);

    open(CHRINFO,">$chrinfo");

    $helper->disp("Querying UCSC database for chromosomal information...");
    $q = $qobj->get_query("chrom_info");
    my $conn = $self->open_connection($sp);
    my $sth = $conn->prepare($q);
    $sth->execute();
    while ($data = $sth->fetchrow_hashref())
    {
        print CHRINFO $data->{"chrom"}."\t".$data->{"size"}."\n";
    }
    $self->close_connection($conn);
    close(CHRINFO);
    return($chrinfo);
}

=head2 get_xml_genes_query

Construct an XML string used for the genes query in Ensembl with Biomart. Mostly for internal use but
may be used to display the query.

    $fetcher->get_xml_exons_query($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $organism must be properly 
formatted through the format_species function.

=cut

sub get_xml_genes_query
{
    my ($self,$species) = @_;
    my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
              "<!DOCTYPE Query>".
              "<Query virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
              "<Dataset name = \"$species\" interface = \"default\" >\n".
              "<Attribute name = \"chromosome_name\" />".
              "<Attribute name = \"start_position\" />".
              "<Attribute name = \"end_position\" />".
              "<Attribute name = \"ensembl_gene_id\" />".
              "<Attribute name = \"percentage_gc_content\" />".
              "<Attribute name = \"strand\" />".
              "<Attribute name = \"external_gene_id\" />".
              "<Attribute name = \"gene_biotype\" />".
              "</Dataset>\n".
              "</Query>\n";
    return($xml);
}

=head2 get_xml_exons_query

Construct an XML string used for the exons query in Ensembl with Biomart. Mostly for internal use but
may be used to display the query.

    $fetcher->get_xml_exons_query($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $organism must be properly 
formatted through the format_species function.

=cut

sub get_xml_exons_query
{
    my ($self,$species) = @_;
    my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
          "<!DOCTYPE Query>".
          "<Query virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
          "<Dataset name = \"$species\" interface = \"default\" >\n".
          "<Attribute name = \"chromosome_name\" />".
          "<Attribute name = \"exon_chrom_start\" />".
          "<Attribute name = \"exon_chrom_end\" />".
          "<Attribute name = \"ensembl_exon_id\" />".
          "<Attribute name = \"ensembl_gene_id\" />".
          "<Attribute name = \"strand\" />".
          "<Attribute name = \"external_gene_id\" />".
          "<Attribute name = \"gene_biotype\" />".
          "</Dataset>\n".
          "</Query>\n";
    return($xml);
}

=head2 get_xml_utr_query

Construct an XML string used for the UTR query in Ensembl with Biomart. Mostly for internal use but
may be used to display the query.

    $fetcher->get_xml_utr_query($organism,$utr);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $organism must be properly 
formatted through the format_species function. $utr can be 5 or 3.

=cut

sub get_xml_utr_query
{
    my ($self,$species,$utr) = @_;
    my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
        "<!DOCTYPE Query>\n".
        "<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
        "<Dataset name = \"$species\" interface = \"default\" >\n".
        "<Attribute name = \"chromosome_name\" />\n".
        "<Attribute name = \"".$utr."_utr_start\" />\n".
        "<Attribute name = \"".$utr."_utr_end\" />\n".
        "<Attribute name = \"ensembl_gene_id\" />\n".
        "<Attribute name = \"transcript_count\" />\n".
        "<Attribute name = \"strand\" />\n".
        "<Attribute name = \"external_gene_id\" />\n".
        "<Attribute name = \"gene_biotype\" />\n".
        "</Dataset>\n".
        "</Query>\n";
    return($xml);
}

=head2 get_xml_cds_query

Construct an XML string used for the CDS query in Ensembl with Biomart. Mostly for internal use but
may be used to display the query.

    $fetcher->get_xml_cds_query($organism);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $organism must be properly 
formatted through the format_species function.

=cut

sub get_xml_cds_query
{
    my ($self,$species) = @_;
    my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
        "<!DOCTYPE Query>\n".
        "<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
        "<Dataset name = \"$species\" interface = \"default\" >\n".
        "<Attribute name = \"chromosome_name\" />\n".
        "<Attribute name = \"genomic_coding_start\" />\n".
        "<Attribute name = \"genomic_coding_end\" />\n".
        "<Attribute name = \"ensembl_exon_id\" />\n".
        "<Attribute name = \"ensembl_gene_id\" />\n".
        "<Attribute name = \"strand\" />\n".
        "<Attribute name = \"external_gene_id\" />\n".
        "<Attribute name = \"gene_biotype\" />\n".
        "</Dataset>\n".
        "</Query>\n";
    return($xml);
}

=head2 format_species

Construct a string representing with the correct nomenclature, species for each external annotaton
source database. Mostly for internal use.

    $fetcher->format_species($organism,$source);
    
For a list of supported organisms, see the SYNOPSIS of HTS::Tools::Count. $source can be one of 'ucsc',
'refseq' or 'ensembl'.

=cut

sub format_species
{
    my ($self,$org,$ver,$source) = @_;
    if ($source =~ m/ucsc|refseq/)
    {
        if ($org =~ m/human/)
        {
            if ($ver =~ m/hg19/) { return("hg19"); }
            elsif ($ver =~ m/hg18/) { return("hg18"); }
        }
        elsif ($org =~ m/mouse/)
        {
            if ($ver =~ m/mm10/) { return("mm10"); }
            elsif ($ver =~ m/mm9/) { return("mm9"); }
        }
        elsif ($org =~ m/rat/)
        {
            if ($ver =~ m/rn5/) { return("rn5"); }
        }
        elsif ($org =~ m/fly/)
        {
            if ($ver =~ m/dm3/) { return("dm3"); }
        }
        elsif ($org =~ m/zebrafish/)
        {
            if ($ver =~ m/danrer7/i) { return("danRer7"); }
        }
    }
    elsif ($source =~ m/ensembl/)
    {
        if ($org =~ m/human/) { return("hsapiens_gene_ensembl"); }
        elsif ($org =~ m/mouse/) { return("mmusculus_gene_ensembl"); }
        elsif ($org =~ m/rat/) { return("rnorvegicus_gene_ensembl"); }
        elsif ($org =~ m/fly/) { return("dmelanogaster_gene_ensembl"); }
        elsif ($org =~ m/zebrafish/) { return("drerio_gene_ensembl"); }
    }
}

=head2 sort_ensembl_genes

Sort the genes file fetched from Ensembl through Biomart. Internal use.

    $fetcher->sort_ensembl_genes($file);

=cut

sub sort_ensembl_genes
{
    my ($self,$infile) = @_;
    my $tmpdir = $self->get("tmpdir");
    my $tmpfile;
    
    if ($^O !~ /MSWin/) # Case of linux, easy sorting
    {
        $helper->disp("Sorting bed file $infile...");
        $tmpfile = File::Spec->catfile($tmpdir,"temp".".in$$");
        `awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $infile > $tmpfile `;
        $infile = $tmpfile;
    }
    else # We are in Windows... package required not able to sort file with header...
    {
        my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
                    "messing up data. Please sort files outside $MODNAME first (e.g. using\n".
                    "Excel or something similar.";
        die "\n$dmsg\n\n";
    }

    return($infile);
}

=head2 sort_ensembl_exons

Sort the exons file fetched from Ensembl through Biomart. Internal use.

    $fetcher->sort_ensembl_exons($file);

=cut

sub sort_ensembl_exons
{
    my ($self,$infile) = @_;
    my $tmpdir = $self->get("tmpdir");
    my $tmpfile;
    
    if ($^O !~ /MSWin/) # Case of linux, easy sorting
    {
        $helper->disp("Sorting bed file $infile...");
        $tmpfile = File::Spec->catfile($tmpdir,"temp".".in$$");
        `awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2 -k3g,3 -u\"}' $infile > $tmpfile `;
        $infile = $tmpfile;
    }
    else # We are in Windows... package required not able to sort file with header...
    {
        my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
                    "messing up data. Please sort files outside $MODNAME first (e.g. using\n".
                    "Excel or something similar.";
        die "\n$dmsg\n\n";
    }

    return($infile);
}

=head2 open_connection($db,@dbdata)

Open a connection to a local ore remote database given a host and database connection credits. Do not
directly use this function, it serves only internal purposes of retrieving data from UCSC database
in order to annotate, read count and plot.

    my $conn = $helper->open_connection("hg19","gbuser","gbpass");

=cut

sub open_connection
{   
    my ($self,$database,@dbdata) = @_;

    # Generally, dbdata should be read from HTS::Tools::Constants
    #@dbdata = ("user","password");
    
    my ($hostname,$conn);
    #if (@dbdata && $self->check_db_existence($database,@dbdata))
    #{
    #    $hostname = "localhost";
    #    $conn = DBI->connect("dbi:mysql:database=$database;host=$hostname;port=3306",$dbdata[0],$dbdata[1]);
    #}
    #else # Connect to the public MySQL host at UCSC
    {
        $hostname = REMOTE_HOST;
        $conn = DBI->connect("dbi:mysql:database=$database;host=$hostname;port=3306",REMOTE_USER);
    }
    return $conn;
}

=head2 close_connection($db,@dbdata)

Close the connection to a local ore remote database. Do not directly use this function, it serves only
internal purposes of retrieving data from UCSC database in order to annotate, read count and plot.

    $helper->close_connection($conn);

=cut

sub close_connection
{ 
    my ($self,$conn) = @_;
    $conn->disconnect();
}

=head2 check_existence($db,@dbdata)

Check if a local or remote database exists. Do not directly use this function, it serves only internal
purposes of retrieving data from UCSC database in order to annotate, read count and plot.

    $helper->check_db_existence("arbDB","user","pass");

=cut

sub check_db_existence
{
    my ($self,$dbcheck,@dbdata) = @_;
    my $out = 1;
    my $conn = DBI->connect("dbi:mysql:database=information_schema;host=localhost;port=3306",$dbdata[0],$dbdata[1]);
    my $query = "SELECT `SCHEMA_NAME` FROM `SCHEMATA` WHERE `SCHEMA_NAME` = \"$dbcheck\"";
    my $sth = $conn->prepare($query);
    $sth->execute();
    $out = 0 if (!$sth->rows());
    $sth->finish();
    $self->close_connection($conn);
    return($out);
}

=head2 get

HTS::Tools::Fetch object getter

    my $param_value = $fetcher->get("param_name")
    
=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Fetch object setter

    $fetcher->set("param_name","param_value")
    
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

1; # End of HTS::Tools::Fetch
