=head1 NAME

HTS::Tools::Constants - A set of constants (e.g. external tool paths) for the HTS::Tools module.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

MySQL queries for HTS::Tools::Constants
	
	# Just load the module so that constants can be used
    use HTS::Tools::Constants;
    my $const->HTS::Tools::Constants->new();
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
    # Load external constants from a YAML file
    $const->load_constants('my_constants.yml');

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Constants;

our $MODNAME = "HTS::Tools::Constants";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Several constants for module HTS::Tools.";

=head2 new

Constructor for HTS::Tools::Constants

=cut

sub new
{
	my $class = shift @_;
	my $self = {};
	bless($self,$class);
	return($self);
}

=head2 init

HTS::Tools::Constants object initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
	my ($self,$params) = @_;
	
	# Here we must initiate a set of standard parameters like paths for external tools, local and
	# remote database hosts, usernames, passwords, number of cores to use etc.
	
	return($self);
}

=head2 load_constants

Load constants from an external YAML file.

	$const->change_constants({'BEDTOOLS_HOME' => '/opt/BEDTools','SAMTOOLS_HOME' => '/opt/SAMTools'})

=cut

sub load_constants
{
	my ($self,$file) = @_;
	
	# Read the YAML file using the module YAML, create a parameters hash and pass it to change_constants
	# which will do the rest of the work like validation 
	
	return($self);
}

=head2 change_constants

Massively change the parameters of an HTS::Tools::Constants object.

	$const->change_constants({'BEDTOOLS_HOME' => '/opt/BEDTools','SAMTOOLS_HOME' => '/opt/SAMTools'})

=cut

sub change_constants
{
	my ($self,$params) = @_;
	
	# Validate the new parameters 
	my $checker = HTS::Tools::Paramcheck->new();
	$checker->set("tool","constants");
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

HTS::Tools::Constants object getter

	my $param_value = $count->get('param_name')

=cut

sub get
{
	my ($self,$name) = @_;
	return($self->{$name});
}

=head2 set

HTS::Tools::Constants object setter

	$intersecter->set('param_name','param_value');
	
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
