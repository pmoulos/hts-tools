=head1 NAME

HTS::Tools::Utils - Helper functions for the HTS::Tools module

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module includes mostly helper functions which are commonly used within the HTS::Tools module

    use HTS::Tools::Utils;

    my $helper = HTS::Tools::Utils->new();
    $helper->now; # Current day and time formatted in a human readable format
    $helper->now("machine"); # Same but in a machine friendly format, useful for naming files etc.
    
    my @array_without_duplicates = keys($helper->unique(@array_with_duplicates));

    my $median_absolute_deviation = $helper->mad(@numerical_array);

    my $nlines = $helper->count_lines($filename);
    my $nlines = $helper->count_lines($filehandle);

    # ...and more helper functions, read the documentation...

=head1 EXPORT

All the functions in this module are exportable

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Utils;

our $MODNAME = "HTS::Tools::Utils";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Utilities for module HTS::Tools.";

use v5.14;
use strict;
use warnings FATAL => 'all';

use Carp;
use File::Spec;
use File::Path qw(make_path remove_tree);
use DBI;

use constant REMOTE_HOST => "genome-mysql.cse.ucsc.edu";
use constant REMOTE_USER => "genome";

BEGIN {
	select(STDOUT);
	$|=1;
}

=head2 new

Simple class constructor

	my $helper = HTS::Tools::Utils->new;

=cut

sub new
{
	my $class = $_[0];
	my $self = {};
	bless($self,$class);
	return($self);
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
	my ($hostname,$conn);
	if (&check_existence($database))
	{
		$hostname = "localhost";
		$conn = DBI->connect("dbi:mysql:database=$database;host=$hostname;port=3306",$dbdata[0],$dbdata[1]);
	}
	else # Connect to the public MySQL host at UCSC
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
	&close_connection($conn);
	return($out);
}

=head2 now($format)

Formats the current time in a human readable or machine friendly (just one string) format. $format
can be "human" (default) or "machine"

	my $date_human = $helper->now;
	my $date_machine = $helper->now("machine");

=cut

sub now
{
	my ($self,$format) = shift @_;
	$format = "human" if (!$format);
	my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
	$year += 1900;
	$month++;
	$month = "0".$month if (length($month)==1);
	$day = "0".$day if (length($day)==1);
	$hour = "0".$hour if (length($hour)==1);
	$min = "0".$min if (length($min)==1);
	$sec = "0".$sec if (length($sec)==1);
	($format ne "machine") ? (return($day."/".$month."/".$year." ".$hour.":".$min.":".$sec)) :
	(return($year.$month.$day.$hour.$min.$sec));
}

=head2 count_lines($file)

Counts the lines of a (non-binary) file

	my $number_of_lines = $helper->count_lines($file);

=cut

sub count_lines
{
	open(IN,$_[1]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
}

=head2 round($number)

Scientifically rounds a number to the closest integer

	my $rounded = $helper->round($real);
	
=cut

sub round
{
	my ($self,$number) = @_;
	return int($number + .5*($number <=> 0));
}

=head2 mean(@array)

Calculates the arithmetic mean of a numeric array

	my $mean = $helper->mean(@array);
	
=cut

sub mean 
{
	my $self = shift @_;
	my $result;
	foreach (@_) { $result += $_ ;}
	return $result/@_;
}

=head2 stdev(@array)

Calculates the standard deviation of a numeric array

	my $stdev = $helper->stdev(@array);
	
=cut

sub stdev 
{
	my $self = shift @_;
	my $mean = mean(@_);
	my @elemsquared;
	foreach (@_)
	{
		push (@elemsquared,($_**2));
	}
	return sqrt(mean(@elemsquared) - ($mean**2));
}	

=head2 median(@array)

Calculates the median of a numeric array

	my $med = $helper->median(@array);
	
=cut

sub median
{
	my $self = shift @_;
    my @pole = sort(@_);
    my $ret;
    if((@pole % 2) == 1)
    {
    	$ret = $pole[((@pole+1)/2) - 1];
    }
    else 
    {
        $ret = ($pole[(@pole/2) - 1] + $pole[@pole/2])/2;
    }
	return $ret;
}

=head2 mad(@array)

Calculates the median absolute deviation of a numeric array

	my $MAD = $helper->mad(@array);
	
=cut

sub mad
{
	my $self = shift @_;
	my @absdiff;
	my $med = median(@_);
	foreach (@_)
	{
		push(@absdiff,($_ - $med));
	}
	return abs(median(@absdiff));
}

=head2 try_module($module)

Checks if a required by the package module exists and dies with additional info in the case that the
module does not exist

	my $exists = $helper->try_module($module);
	
=cut

sub try_module
{
	my ($self,$module,@fun) = @_;
	eval "require $module";
	if ($@)
	{
		my $killer = "Module $module is required to continue with the execution. If you are in\n". 
					 "Windows and you have ActiveState Perl installed, use the Package Manager\n".
					 "to get the module. If you are under Linux, log in as a super user (or use\n".
					 "sudo under Ubuntu) and type \"perl -MCPAN -e shell\" (you will possibly have\n".
					 "to answer some questions). After this type \"install $module\" to install\n".
					 "the module. If you don't know how to install the module, contact your\n".
					 "system administrator.";
		die "\n$killer\n\n";
	}
	else
	{
		if (@fun)
		{
			my $funs = join(" ",@fun);
			eval "use $module qw($funs)";
		}
		else { eval "use $module"; }
	}
}

=head2 swap(@array_two_members)

Simply swaps the two first elements of an array, ignores the rest. It is intended to be used ONLY with
arrays with two elements.

	my @swapped_array = $helper->swap(@array);
	
=cut

sub swap
{
	return($_[2],$_[1]);
}

=head2 disp(@array_of_messages)

Simply displays a user-defined message in STDERR if verbosity is requested from a higher level. It also
prints messages to a log file handle, if requested from a higher level.

	$helper->disp("Hello world!");
	
=cut

sub disp
{
	my $self = shift @_;
	print STDERR "\n@_" if (!$self->get("silent"));
	#print $logfilehandle "\n@_" if (!$silent && $log);
}

=head2 disp(@array_of_messages)

Prints some credits for the respective module

	my $mod = "My::Module"
	my $ver = "0.01";
	my $auth = "John Doe";
	my $email = "john.doe\@example.com"
	my $adtext = "My extraordinary module."
	my $helper = HTS::Tools::Utils->new();
	$helper->advertise($mod,$ver,$auth,$email,$adtext);
	
=cut

sub advertise
{
	my ($self,$mod,$ver,$auth,$email,$adtext) = @_;
	$mod = "My::Module" unless($mod);
	$ver = "0.01" unless($ver);
	$auth = "John Doe" unless($auth);
	$email = "john.doe\@example.com" unless($email);
	$adtext = "My extraordinary module" unless($adtext);
	my $credits = $adtext." Copyright: ".$auth." ".$email;
	use Term::ANSIColor;
	print color 'bold yellow on_blue';
	$self->disp($mod." ".$ver);
	print color 'bold green';
	$self->disp($credits."\n");
	print color 'reset';
}

=head2 catch_cleanup

Ctrl-C signal catcher that ensures removing of temporary directories. Never use this directly as it
will throw an error

	SIG{INT} = sub { $helper->catch_cleanup }

=cut

sub catch_cleanup 
{
	my $self = shift @_;
	$self->cleanup;
	croak "\nCatching Ctrl-C, cleaning temporary files!\n";
}

=head2 cleanup

The cleanup function of the signal catcher

	$helper->cleanup

=cut

sub cleanup 
{
	my $self = shift @_;
	remove_tree($self->get("tmpdir"));
}

=head2 get

HTS::Tools::Utils object getter

	my $param_value = $helper->get("param_name")
=cut

sub get
{
	my ($self,$name) = @_;
	return($self->{$name});
}

=head2 set

HTS::Tools::Utils object setter

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

1; # End of HTS::Tools::Utils
