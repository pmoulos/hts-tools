=head1 NAME

HTS::Tools::Track::Hub - A massive motif scanner (not finder!) for short genomic regions

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module constructs a UCSC track hub using either a set of basic track hub parameters (e.g. the hub
base, the hub name and description etc.) together with a text tab delimited file with the tracks to be
generated and their properties (some are required, see HTS::Tools::Track::Signal), or with the use of
a single YAML file containing all the parameters required to build the hub

    use HTS::Tools::Track::Hub;
    my %params = (
        'config' => 'trackhub.yml'
    )
    my $huber = HTS::Tools::Track::Hub->new(\%params);
    $huber->run;

    # OR

    %params = (
        'hubid' => 'test_hub',
        'hubname' => 'A test hub',
        'hubdesc' => 'A hub demonstrating HTS::Tools::Track::Hub',
        'hubmail' => hub_admin@trackhubs.org,
        'hubbase' => '/the/public/html/directory/you/store/hubs',
        'tracks' => '/the/tab/delimited/file/describing/tracks/and/options'
    )
    my $huber2 = HTS::Tools::Track::Hub->new(\%params);
    $huber2->run;

The acceptable parameters are as follows:

=over 4

=item I<config> B<(required)>

A YAML file with the parameters of hub creation. It can look like the following:

HUB_ID: test_hub
HUB_NAME: The first test hub
HUB_DESCRIPTION: This is the first hub made with HTS::Tools::Track::Hub module
HUB_BASE: /media/raid/tracks/test
HUB_MAIL: moulos@fleming.gr
TRACKS:
 TCF4:
  filename: /media/HD4/Fleming/play/trackhub_play/TCF4.bed
  source: bed
  destination: bigwig
  name: TCF4
  urlbase: http://epigenomics.fleming.gr/tracks/test
  description: TCF4 ChIP-Seq signal
  color: 0,0,160
  maxHeightPixels: 128:64:16
  visibility: full
  boxedCfg: on
  autoScale: on
  group: user
  priority: auto
 TCF4_peaks:
  filename: /media/HD4/Fleming/play/trackhub_play/TCF4_peaks.bed
  source: bed
  destination: bigbed
  name: TCF4 peaks
  urlbase: http://epigenomics.fleming.gr/tracks/test
  description: TCF4 ChIP-Seq peaks
  color: 0,0,160
  maxHeightPixels: 128:64:16
  visibility: dense
  boxedCfg: on
  autoScale: on
  group: user
  priority: auto
 CON_BR1:
  filename: /media/HD4/Fleming/play/trackhub_play/CON_BR1.bam
  source: bam
  destination: bigwig
  name: CON BR1
  urlbase: http://epigenomics.fleming.gr/tracks/test
  description: CON RNA-Seq signal replicate 1
  color: 120,120,120
  maxHeightPixels: 128:64:16
  visibility: full
  boxedCfg: on
  autoScale: on
  group: user
  priority: auto
 DOX_BR1:
  filename: /media/HD4/Fleming/play/trackhub_play/DOX_BR1.bam
  source: bam
  destination: bigwig
  name: DOX BR1
  urlbase: http://epigenomics.fleming.gr/tracks/test
  description: DOX RNA-Seq signal replicate 1
  color: 0,160,0
  maxHeightPixels: 128:64:16
  visibility: full
  boxedCfg: on
  autoScale: on
  group: user
  priority: auto

If config is given, no other parameters are required and if given, will be ignored.

=item I<hubid> B<(required)>

The unique identifier of the hub to be created.

=item I<hubname> B<(optional)>

A short name for the hub to be created.

=item I<hubdesc> B<(optional)>

A longer description for the hub to be created.

=item I<hubmail> B<(required)>

The track hub admin e-mail.

=item  I<hubbase> B<(required)>

The public html directory where the hub will live.

=item I<tracks> B<(required)>

A text tab-delimited file containing information for the tracks that will be included in the track hub.
See also HTS::Tools::Track::Signal for details.

=item I<silent> B<(optional)>

Use this parameter if you want to turn informative messages off.

=head1 OUTPUT

The output of the module is the track hub directory, properly structured and ready to be loaded in an
instance of the UCSC Genome Browser, as long as the hub base is public.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Track::Hub;

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use Cwd;
use File::Basename;
use File::Path qw(make_path remove_tree);
use File::Temp;
use File::Spec;

use HTS::Tools::Paramcheck;
use HTS::Tools::Utils;

use vars qw($helper $const);

our $MODNAME = "HTS::Tools::Track::Hub";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Create a UCSC Genome Browser track hub.";

BEGIN {
    $helper = HTS::Tools::Utils->new();
    select(STDOUT);
    $|=1;
    $SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new

The HTS::Tools::Track::Hub object constructor. It accepts a set of parameters that are required to run the
motifscanner and get the output.

    my $huber = HTS::Tools::Track::Hub->new({'config' => 'myhub.yml'});

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
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "track_hub",
        "params" => $params});
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init

HTS::Tools::Track::Hub object initialization method. NEVER use this directly, use new instead.

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

The HTS::Tools::Track::Hub run subroutine. It runs the motifscanner with the given parameters in the 
constructor.

    $huber->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    # Copy some memory-less variables to avoid rewriting the whole thing...
    my ($hubid,$hubname,$hubdesc,$hubmail,$hubbase,$hubtracks);
    if ($self->get("config"))
    {
        my ($opts,$tracks) = $self->read_yaml($self->get("config"));
        croak "Corrupted or not valid configuration file" if (!$opts);
        $hubid = $opts->{"hubid"};
        $hubname = $opts->{"hubname"};
        $hubdesc = $opts->{"hubdesc"};
        $hubmail = $opts->{"hubmail"};
        $hubbase = $opts->{"hubbase"};
        $hubtracks = $tracks;
    }
    else
    {
        $hubid = $self->get("hubid");
        $hubname = $self->get("hubname");
        $hubdesc = $self->get("hubdesc");
        $hubmail = $self->get("hubmail");
        $hubbase = $self->get("hubbase");
        $hubtracks = $self->get("tracks");
    }
    
    # And implement the workflow
    my $hubdir = $self->create_hub_tree($hubbase,$hubtracks);
    my $hubfile = $self->create_hub_file($hubdir);
    my $genomefile = $self->create_genomes_file($hubdir,$hubtracks);
    my @signals = $self->create_signals;
}

=head2 read_yaml

Read the YAML file with all configurations for the hub creation.

    my $config = $huber->read_yaml($configfile);

=cut

sub read_yaml
{
    my ($self,$file) = @_;
    my ($yfh,$yaml);
    use YAML qw(LoadFile Dump);
    eval
    {
        open($yfh,"<",$file);
        $yaml = LoadFile($yfh);
        close($yfh);
    };
    
    croak "Bad YAML hub parameters file! Commiting suicide..." if ($@);
    return($yaml);
}

=head2 read_tab

Read the tab-delimited file with all details regarding the tracks in the track hub.

    my $tracks = $huber->read_tab($trackfile);

=cut

sub read_tab
{
    my ($self,$file) = @_;
    my $tracks;
    my $line;
    my @cols;
    open(TRACKS,$file) or croak "$file does not exist! Commiting suicide...";
    my $header = <TRACKS>;
    $header =~ s/\r|\n$//g;
    # Header should determine the structure of the file and validation checks
    while ($line = <TRACKS>)
    {
        $line =~ s/\r|\n$//g;
        @cols = split(/\t/,$line);
    }
}

=head2 create_hub_tree

Create the hub tree directory structure.

    my $file = $huber->create_hub_tree($hubbase,$hubid);

=cut

sub create_hub_tree {}

=head2 create_hub_file

Create the hub description file.

    my $file = $huber->create_hub_file($hubbase,$hubid,$hubname,$hubdesc,$hubmail);

=cut

sub create_hub_file {}

=head2 create_genomes_file

Create the hub genomes file.

    my $file = $huber->create_genomes_file($hubbase,$hubgenomes);

=cut

sub create_genomes_file {}

=head2 parse_track_line

Parse a track line and create the corresponding trackDb elements.

    my %trackopts = $huber->parse_track_line($header);

=cut

sub parse_track_line {}

=head2 convert_genome

Convert between genomes (of the same organism) to create tracks for different genome versions. The
genome versions are the ones supported by HTS::Tools::Fetch.

    my $file = $huber->convert_genomes($track,$from,$to);

=cut

sub convert_genome {}

=head2 create_track

Create the hub tracks using HTS::Tools::Track::Signal.

    my ($track,$header) = $huber->create_track($opts);

=cut

sub create_track {}

=head2 get

HTS::Tools::Track::Hub object getter.

    my $param_value = $motifscanner->get('param_name')

=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Track::Hub object setter.

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

    perldoc HTS::Tools::Track::Hub

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

1; # End of HTS::Tools::Track::Hub
