=head1 NAME

HTS::Tools::Constants - A set of constants (e.g. external tool paths) for the HTS::Tools module.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module initiates a set of standard parameters that can be used from HTS::Tools like paths to external
tools, local and remote database hosts, usernames, passwords, number of cores to use etc. A full list of
possible constants will be provided soon...
		
	# Just load the module so that constants can be used
    use HTS::Tools::Constants;
    my $const->HTS::Tools::Constants->new();
    # Load constants from external YAML file
    my $const->HTS::Tools::Constants->new({'file' => 'my_constans.yml'});
	# Change constant
    $const->set('BEDTOOLS_HOME','/opt/BEDTools');
    # Massively set new constants
    my %new_constants = (
		'BEDTOOLS_HOME' => '/opt/BEDTools',
		'SAMTOOLS_HOME' => '/opt/SAMTools',
		'MAX_CORES' => 8
	)
    $const->change_constants(\%new_constants);
    # Get constant
    $const->get('MAX_CORES');
    # Reload external constants from a YAML file
    $const->load_constants('my_new_constants.yml');

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Constants;

our $MODNAME = "HTS::Tools::Constants";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Several constants for module HTS::Tools.";

#use lib 'D:/Software/hts-tools/HTS-Tools/lib';
use lib '/media/HD4/Fleming/hts-tools/HTS-Tools/lib';
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

Constructor for HTS::Tools::Constants

=cut

sub new
{
	my ($class,$params) = shift @_;
	my $self = {};
	bless($self,$class);
	$self->init($params);
	return($self);
}

=head2 init

HTS::Tools::Constants object initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
	my ($self,$params) = @_;
	my ($constants,$userdef,$loaded,$checker);
	
	# First load the defaults as the external file might not be full
	$constants = $self->load_default_constants;
	# Then check the YAML file and validate
	if (defined($params->{"file"}))
	{
		($userdef,$loaded) = $self->load_user_constants($params->{"file"});
		# In the case that an external file has succesfully loaded, we have to validate, else, no need
		# as the defaults are loaded, which are error-free of course! 
		if ($loaded)
		{
			$checker = HTS::Tools::Paramcheck->new({"tool" => "constants","params" => $external});
			$userdef = $checker->validate;
			$self->change_constants($external,"skip");
		}
	}
	
	return($self);
}

=head2 load_constants

Load constants from an external YAML file.

	$const->load_constants($yaml_file)

=cut

sub load_user_constants
{
	my ($self,$file) = @_;
	my ($pfh,$phref,$hasloaded);
	use YAML qw(LoadFile Dump);
	eval
	{
		open($pfh,"<",$file);
		$phref = LoadFile($pfh);
		close($pfh);
		$hasloaded = 1;
	};
	if ($@)
	{
		$helper->disp("Bad constants file! Ignoring...");
		$phref = $self->load_default_constants;
		$hasloaded = 0;
	}
	return($phref,$hasloaded);
}

=head2 load_default_constants

Load default constants in the absence of an external YAML file.

	$const->load_default_constants({'BEDTOOLS_HOME' => '/opt/BEDTools','SAMTOOLS_HOME' => '/opt/SAMTools'})

=cut

sub load_default_constants
{	
	my $self = shift @_;
	
	my ($pfh,$phref);
	use YAML qw(LoadFile Dump);
	if (-f "../config.yml")
	{
		eval
		{
			open($pfh,"<","../config.yml");
			$phref = LoadFile($pfh);
			close($pfh);
		};
		if ($@)
		{
			$helper->disp("Error loading constants file! Falling back to minimums...");
			$phref = $self->load_hard_constants;
		}
	}
	else
	{
		$helper->disp("Constants configuration file does not exist! Falling back to minimums...");
		$phref = $self->load_hard_constants;
	}
	
	return($phref);
}

=head2 load_hard_constants

Load hard coded constants as a last resort. There must be a default fallback with very basic constants.

	$const->load_hard_constants;

=cut

sub load_hard_constants
{
	my $self = shift @_;
	
	my $constants = {
		"LOCAL_HOST" => "localhost",
		"LOCAL_USER" => "user",
		"LOCAL_PASS" => "password",
		"BIOMART_PATH" => "http://www.biomart.org/biomart/martservice?",
		"REMOTE_HOST" => "genome-mysql.cse.ucsc.edu",
		"REMOTE_USER" => "genome",
		"BEDTOOLS_HOME" => "/opt/NGSTools/BEDTools/bin",
		"SAMTOOLS_HOME" => "/opt/NGSTools/SAMTools",
		"GENOMICTOOLS_HOME" => "/opt/NGSTools/GenomicTools",
		"MAX_CORES" => 12,
		"MACS_HOME" => "/usr/bin/macs14",
		"MACS2_HOME" => "/usr/local/bin/macs2",
		"SICER_HOME" => "/opt/NGSTools/SICER",
		"RSEG_HOME" => "/opt/NGSTools/rseg",
		"GIMMEMOTIFS_HOME" => "/usr/bin",
		"MOTIFSCANNER_HOME" => "",
		"LOCAL_GENOMES" => "/opt/genomes"
	};
	
	return($constants);
}

=head2 change_constants

Massively change the parameters of an HTS::Tools::Constants object.

	$const->change_constants({'BEDTOOLS_HOME' => '/opt/BEDTools','SAMTOOLS_HOME' => '/opt/SAMTools'})

=cut

sub change_constants
{
	my ($self,$params,$check) = @_;
	$check = "do" unless($check);
	
	if ($check eq "do") # Validate the new parameters
	{
		my $checker = HTS::Tools::Paramcheck->new({"tool" => "constants","params" => $params});
		$params = $checker->validate;
	}
	
	# If validator does not complain or if we are loading the defaults, change/set the constants
	while (my ($name,$value) = each(%$params))
	{
		$self->set($name,$value);
	}
	return($self);
}

=head2 get

HTS::Tools::Constants object getter.

	my $param_value = $const->get('param_name')

=cut

sub get
{
	my ($self,$name) = @_;
	return($self->{$name});
}

=head2 set

HTS::Tools::Constants object setter.

	$const->set('param_name','param_value');
	
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

    perldoc HTS::Tools::Utils


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

1; # End of HTS::Tools::Constants

#$status = eval { $helper->try_module("YAML") };
#if ($status)
#{
#	$helper->disp("Module YAML is required to read an external constants file! Will try with the defaults...");
#	$self->load_default_constants;
#}
#else
#{ }
