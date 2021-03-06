
=head1 NAME

HTS::Tools::Convert - A set of format conversion utilities for HTS::Tools

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use HTS::Tools::Convert;

    my $foo = HTS::Tools::Convert->new();

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Convert;

use v5.10;
use strict;
use warnings FATAL => 'all';

our $MODNAME = "HTS::Tools::Convert";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Set of format conversion utilities.";

use Carp;
use File::Basename;
use File::Spec;

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
    
    # Validate the input parameters
    my $checker = HTS::Tools::Paramcheck->new();
    $checker->set("tool","convert");
    $checker->set("params",$params);
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init($params)

HTS::Tools::Convert object initialization method. NEVER use this directly, use 
new instead.

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

=head2 fasta2tab($fasta_file)

Convert a FASTA file to sequences in tabular format

    my $tab_converted = $converter->fasta2tab($fasta_file);

=cut

sub fasta2tab
{
    my ($self,$infile) = @_;
    my $currid;
    open(INPUT,$infile) or croak "\nThe file $infile does not exist!\n";
    my ($base,$dir) = fileparse($infile,'\.[^.]*');
    my $outfile = File::Spec->catfile($dir,$base.".tab");
    open(OUTPUT,">$outfile");
    while (my $line = <INPUT>)
    {
        $line =~ s/\r|\n$//g;
        if ($line =~ /^>/)
        {
            print OUTPUT "\n" if ($currid);
            $currid = $line;
            $currid =~ s/^>//g;
            print OUTPUT "$currid\t";
        }
        elsif ($currid)
        {
            print OUTPUT "$line";
        }
    }
    print OUTPUT "\n";
    close(INPUT);
    close(OUTPUT);
    return($outfile);
}

=head2 tab2fasta($tab_file)

Convert a tabular sequences file to FASTA format

    my $fasta_converted = $converter->tab2fasta($tab_file);

=cut

sub tab2fasta
{
    my ($self,$infile) = @_;
    my @conts;
    open(INPUT,$infile) or croak "\nThe file $infile does not exist!\n";
    my ($base,$dir) = fileparse($infile,'\.[^.]*');
    my $outfile = File::Spec->catfile($dir,$base.".fa");
    open(OUTPUT,">$outfile");
    while (my $line = <INPUT>)
    {
        $line =~ s/\r|\n$//g;
        @conts = split(/\t/,$line);
        print OUTPUT ">$conts[0]\n";
        my $seq = $conts[1];
        while ($seq)
        {
            my $olin = substr($seq,0,100,"");
            print OUTPUT "$olin\n";
        }
    }
    close(INPUT);
    close(OUTPUT);
}

=head2 gff2bed6($tab_file)

Convert a GFF file to BED6 (no groups) file format

    my $gff2bed6 = $converter->gff2bed6($gff_file);

=cut

sub gff2bed6
{
    my ($self,$infile) = @_;
    my @conts;
    open(INPUT,$infile) or croak "\nThe file $infile does not exist!\n";
    my ($base,$dir) = fileparse($infile,'\.[^.]*');
    my $outfile = File::Spec->catfile($dir,$base.".bed");
    my $fetcount = 0;
    open(OUTPUT,">$outfile");
    while (my $line = <INPUT>)
    {
        $fetcount++;
        $line =~ s/\r|\n$//g;
        @conts = split(/\t/,$line);
        print OUTPUT "$conts[0]\t$conts[4]\t$conts[4]\t$conts[3].$fetcount\t$conts[5]\t$conts[6]\n";
    }
    close(INPUT);
    close(OUTPUT);
}

=head2 gff2bed6($tab_file)

Convert a GFF file to BED12 file format

    my $gff2bed12 = $converter->gff2bed12($gff_file);

=cut

sub gff2bed12
{
    my ($self,$infile) = @_;
    my @conts;
    open(INPUT,$infile) or croak "\nThe file $infile does not exist!\n";
    my ($base,$dir) = fileparse($infile,'\.[^.]*');
    my $outfile = File::Spec->catfile($dir,$base.".bed");
    my $fetcount = 0;
    open(OUTPUT,">$outfile");
    while (my $line = <INPUT>)
    {
        $fetcount++;
        $line =~ s/\r|\n$//g;
        @conts = split(/\t/,$line);
        # This remains to be filled...
    }
    close(INPUT);
    close(OUTPUT);
}

=head2 get

HTS::Tools::Convert object getter

    my $param_value = $converter->get("param_name")
=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Convert object setter

    $converter->set("param_name","param_value")
    
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

Please report any bugs or feature requests to C<bug-hts-tools at rt.cpan.org>, 
or through the web interface at 
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=HTS-Tools>.  I will be notified, 
and then you'll automatically be notified of progress on your bug as I make 
changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc HTS::Tools::Convert


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

1; # End of HTS::Tools::Convert
