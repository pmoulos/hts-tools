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

use lib '/media/HD4/Fleming/dev/HTS-Tools/lib';
use HTS::Tools::Utils;

use vars qw($helper);

BEGIN {
	$helper = HTS::Tools::Utils->new();
	select(STDOUT);
	$|=1;
	$SIG{INT} = sub { $helper->catch_cleanup; }
}

use constant MAXCORES => 12;

=head2 new($tool,%params)

=cut

sub new
{
	my ($class,$args) = @_;
	my $self = {};

	# Pass global variables to the helper
	(defined($args->{"silent"})) ? ($helper->set("silent",$args->{"silent"})) :
		($helper->set("silent",0));
		
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
	return($self->{"params"});
}

=head2 validate_convert

The parameter validator function of the HTS::Convert module. Do not use this directly, use the validate
function instead

=cut

sub validate_convert
{
	my $self = shift @_;
	my $modname = "HTS::Tools::Convert";
	return($self->{"params"});
}

=head2 validate_count

The parameter validator function of the HTS::Count module. Do not use this directly, use the validate
function instead

=cut

# TODO: Check for unrecognized parameters
sub validate_count
{
	my $self = shift @_;
	my $modname = "HTS::Tools::Count";

	# Check required packages
	$helper->try_module("Tie::IxHash::Easy");

	my @accept = ("input","region","sort","percent","lscore","escore","constant","small","split","stats",
		"output","ncore","source","silent","tmpdir");

	# Check fatal
	my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!$self->{"params"}->{"input"});
    $stop .= "--- Please specify region file ---\n" if (!$self->{"params"}->{"region"});
    $stop .= "--- The supported genomes are organism-type, where organism is human, mouse, rat, fly or zebrafish and type is gene, exon, 5utr, 3utr or cds! ---\n"
		if ( !-f $self->{"params"}->{"region"}
			&& $self->{"params"}->{"region"} !~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i);
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
    if ($self->{"params"}->{"stats"} && !$self->{"params"}->{"split"})
    {
    	$helper->disp("You can't calculate area statistics without splitting to sub-areas. Option deactivated...");
    	$self->{"params"}->{"stats"} = 0;
    }
    if ($self->{"params"}->{"ncore"})
    {
		#$helper->try_module("Parallel::ForkManager");
		my $status = eval { $helper->try_module("Parallel::ForkManager") };
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
		if ($source)
		{
			$source = lc($source);
			if ($source ~~ keys(%sources))
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
	return($self->{"params"});
}

=head2 validate_multisect

The parameter validator function of the HTS::Multisect module. Do not use this directly, use the validate
function instead

=cut

sub validate_multisect
{
	my $self = shift @_;
	my $modname = "HTS::Tools::Multisect";
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
	return($self->{"params"});
}

=head2 validate_query

The parameter validator function of the HTS::Queries module. Do not use this directly, use the validate
function instead

=cut

sub validate_queries
{
	my ($self,$params) = @_;
	my $modname = "HTS::Tools::Queries";
	my $query = $params->{"query"};

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
		"refseq_alternative_cds"
	);

	if (!($query ~~ @accept))
	{
		$helper->disp("Unkown query for $modname: $query\n");
		croak "Type perldoc $modname for help in usage.\n\n";
	}
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

1; # End of HTS::Tools::Fetch
