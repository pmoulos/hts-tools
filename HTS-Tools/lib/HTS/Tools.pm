=head1 NAME

HTS::Tools - A Perl collection for processing results from next generation sequencing experiments!

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is just a wrapper for the rest of the HTS::Tools

    use HTS::Tools;

    my $tool = HTS::Tools->new($tool,%params);
    $tool->run;

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools;

our $MODNAME = "HTS::Tools";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Tool collection for the analysis of Next Generation Sequencing experiments.";

use v5.10;
use strict;
use warnings FATAL => 'all';

use File::Temp;

#use lib 'D:/Software/hts-tools/HTS-Tools/lib';
use lib '/media/HD4/Fleming/hts-tools/HTS-Tools/lib';
use HTS::Tools::Assign;
use HTS::Tools::Convert;
use HTS::Tools::Count;
use HTS::Tools::Fetch;
use HTS::Tools::Intersect;
use HTS::Tools::Motifscan;
use HTS::Tools::Normalize;
use HTS::Tools::Profile;
use HTS::Tools::QC;
use HTS::Tools::Utils;

use vars qw($tmpdir $helper);

BEGIN {
	$helper = HTS::Tools::Utils->new();
	# On Ctrl-C or die, do cleanup
	$SIG{INT} = sub { $helper->catch_cleanup; }
}

=head2 new($args)

=cut

sub new
{
	my ($class,$args) = @_;
	my $self = {};
	bless($self,$class);
	$self->init($args);
	return($self);
}

=head2 init($args)

HTS::Tools object initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
	my ($self,$args) = @_;

	# Init self with global variables such as tmpdir and verbosity
	(defined($args->{"tmpdir"})) ? ($self->set("tmpdir",$args->{"tmpdir"})) :
		($self->set("tmpdir",File::Temp->newdir()));
	(defined($args->{"silent"})) ? ($self->set("silent",$args->{"silent"})) :
		($self->set("silent",0));
	(defined($args->{"log"})) ? ($self->set("log",$args->{"log"})) :
		($self->set("log",0));

	# Pass the above global parameters to the parameter structure for each tool. We do it like this
	# because each module is supposed to be used also independently of the wrapper.
	$args->{"params"}->{"tmpdir"} = $self->get("tmpdir");
	$args->{"params"}->{"silent"} = $self->get("silent");
	$args->{"params"}->{"log"} = $self->get("log");
	
	use v5.14;
	given($args->{"tool"})
	{
		when(/assign/i) {
			$self->set("tool",HTS::Tools::Assign->new($args->{"params"}));
		}
		when(/convert/i) {
			$self->set("tool",HTS::Tools::Convert->new($args->{"params"}));
		}
		when(/count/i) {
			$self->set("tool",HTS::Tools::Count->new($args->{"params"}));
		}
		when(/fetch/i) {
			$self->set("tool",HTS::Tools::Fetch->new($args->{"params"}));
		}
		when(/intersect/i) {
			$self->set("tool",HTS::Tools::Intersect->new($args->{"params"}));
		}
		when(/motifscan/i) {
			$self->set("tool",HTS::Tools::Motifscan->new($args->{"params"}));
		}
		when(/normalize/i) {
			$self->set("tool",HTS::Tools::Normalize->new($args->{"params"}));
		}
		when(/profile/i) {
			$self->set("tool",HTS::Tools::Profile->new($args->{"params"}));
		}
		when(/qc/i) {
			$self->set("tool",HTS::Tools::QC->new($args->{"params"}));
		}
	}

	return($self);
}

=head2 run

# The main running function of HTS::Tools

=cut

sub run	
{
	my $self = shift @_;
	my $tool = $self->get("tool");
	$self->get("tool")->run;
	$helper->cleanup;
}

=head2 get

HTS::Tools::Tools object getter

	my $param_value = $helper->get("param_name")
=cut

sub get
{
	my ($self,$name) = @_;
	return($self->{$name});
}

=head2 set

HTS::Tools::Tools object setter

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

    perldoc HTS::Tools


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

1; # End of HTS::Tools
