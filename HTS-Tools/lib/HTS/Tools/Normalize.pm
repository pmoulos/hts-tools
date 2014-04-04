=head1 NAME

This module is just a wrapper for the rest of the HTS::Tools

    use HTS::Tools;

    my $normalizer = HTS::Tools::Normalize->new($what,%params);
    $normalizer->run;

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use HTS::Tools::Normalize;

    my $foo = HTS::Tools::Normalize->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Normalize;

our $MODNAME = "HTS::Tools::Normalize";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Track normalization wrapper for HTS::Tools.";

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use HTS::Tools::Normalize::Bed;
use HTS::Tools::Normalize::Bedgraph;
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

The HTS::Tools::Normalize object constructor. It accepts a set of parameters that are required to run
the normalizer and get the output.

    my $normalizer = HTS::Tools::Normalize->new({'type' => 'bedgraph','input' => 'myfile.bedGraph',
        'sumto' => 1000000000});

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
    $checker->set("tool","normalize");
    $checker->set("params",$params);
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}


=head2 init($args)

HTS::Tools::Normalize object initialization method. NEVER use this directly, use new instead.

=cut

sub init
{
    my ($self,$args) = @_;

    ## Pass the above global parameters to the parameter structure for each tool. We do it like this
    ## because each module is supposed to be used also independently of the wrapper.
    #$args->{"params"}->{"tmpdir"} = $self->get("tmpdir");
    #$args->{"params"}->{"silent"} = $self->get("silent");
    #$args->{"params"}->{"log"} = $self->get("log");

    # Can't use matching as bed and bedgraph would cause problem
    if ($args->{"type"} eq "bed")
    {
        $self->set("tool",HTS::Tools::Normalize::Bed->new($args));
    }
    elsif ($args->{"type"} eq "bedgraph")
    {
        $self->set("tool",HTS::Tools::Normalize::Bedgraph->new($args));
    }
    
    return($self);
}

=head2 run

# The main running function of HTS::Tools::Normalize

=cut

sub run 
{
    my $self = shift @_;
    my $tool = $self->get("tool");
    $tool->run;
    $helper->cleanup;
}

=head2 get

HTS::Tools::Normalize object getter

    my $param_value = $tool->get("param_name")
=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Normalize object setter

    $tool->set("param_name","param_value")
    
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

    perldoc HTS::Tools::Normalize


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

1; # End of HTS::Tools::Normalize
