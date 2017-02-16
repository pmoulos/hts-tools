=head1 NAME

HTS::Tools::Normalize::Bed - Normalize UCSC BED files by read downsampling

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

A perl module to normalize reads across several samples from compensate for differences 
in sequencing depths by random downsampling. The downsampling can be performed based on
the sample with the lowest reads or down to a given number of reads (e.g. 10000000).Can 
run on Linux or Windows.

    use HTS::Tools::Normalize::Bed;
    my %params = (
        'input' => ['normal_rnaseq_1.bed','normal_rnaseq_2.bed',
            'disease_rnaseq_1.bed','disease_rnaseq_2.bed'],
        'rand' => '20000000'
    )
    my $bed = HTS::Tools::Normalize::Bed->new(\%params);
    $bed->run;

The acceptable parameters are as follows:

=over 4

=item input B<(required)>

Input bed file(s). Please be careful as there is little checking whether the input
file(s) are indeed bed files. They should contain 3 or 6 columns. Input files can
be sorted or unsorted. Sorted are prefered to maintain uniformity.

=item rand B<(optional)>

The method of random removal of sequences. Use "generate" to generate random numbers 
between 1..#tags and remove as many as to reach the minimum tag length, "permute" to 
randomly permute the tags and then remove as many as required to reach the minimum 
tag length or an integer. denoting how many tags should be left. If a sample has less
than this number, it is not touched. Defaults to "generate".

=item savrem B<(optional)>

Use this option if you wish to save the tags that are removed in a different file. 
The file will be named "filename_removed_tags.sav".

=item sort B<(optional)>

Use this option if you wish to sort the input files first. The algorithmic results 
of the method will not change.

=item silent B<optional>

Do not display verbose messages.

=back

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Normalize::Bed;

our $MODNAME = "HTS::Tools::Normalize::Bed";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Bed track normalization function.";

use v5.10;
use strict;
use warnings FATAL => 'all';

use Carp;
use File::Copy;
use File::Basename;
use File::Spec;
use File::Temp;

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

The HTS::Tools::Normalize::Bed object constructor. It accepts a set of parameters that are required to run the
counter and get the output.

    my $bednormer = HTS::Tools::Normalize::Bed->new({'input' => ('myfile1.bed','myfile2.bed')});
    $bednormer = HTS::Tools::Normalize::Bed->new({'input' => ('myfile1.bed','myfile2.bed'), 'sumto' => 10000000});

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
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "normalize_bed",
        "params" => $params});
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init($params)

HTS::Tools::Normalize::Bed object initialization method. NEVER use this directly, use new instead.

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

The HTS::Tools::Normalize::Bed run subroutine. It runs the normalizer with the given parameters in the constructor.

    $bednorm->run;
    
=cut

sub run
{
    my $self = shift @_;

    use Data::Dumper;
    
    my @infile = @{$self->get("input")};
    my $savrem = $self->get("savrem");
    my $rand = $self->get("sumto");
    my $sort = $self->get("sort");

    my @originfile = @infile; # Keep original filenames (verbose purposes)
    @infile = $self->sort_inputs(@infile) if ($sort);
    my (@base,@dir,@ext);
    for my $i (0..$#originfile)
    {
        ($base[$i],$dir[$i],$ext[$i]) = fileparse($originfile[$i],'\.[^.]*');
    }

    # Some important variables
    my $numberfiles = @infile;

    # Get the number of lines for each file, performing at the same time a check to see if it
    # is a valid bed file
    my @linesfiles = $self->check_and_count(@infile);

    # Sort the array which contains the number of lines so as to get the min number of lines
    my @sortedlines = sort {$a <=> $b} @linesfiles;
    our $minlines = $sortedlines[0];

    my (@outfile,@randnum,@sortidx,@indrem);
    my ($i,$j); # General index
    my $nochange;
    for ($i=0; $i<$numberfiles; $i++)
    {
        $helper->disp("Now processing file $infile[$i]...");
        $nochange = 0;

        if ($rand eq "generate")
        {
            $helper->disp("Generating random numbers...");
            @randnum = Math::Random::random_uniform($linesfiles[$i],0,1);
            $helper->disp("Sorting them...");
            @sortidx = $helper->sort_by_index(@randnum);
            $helper->disp("Getting required and writing output...");
            @indrem = sort {$a <=> $b} @sortidx[0..$linesfiles[$i] - $minlines - 1];
            $nochange = 1 if ($#indrem == -1);
        }
        elsif ($rand eq "permute")
        {
            $helper->disp("Generating random indices...");
            @randnum = Math::Random::random_permuted_index($linesfiles[$i]);
            $helper->disp("Getting required and writing output...");
            @indrem = sort {$a <=> $b} @randnum[0..$linesfiles[$i] - $minlines - 1];
            $nochange = 1 if ($#indrem == -1);
        }
        else # $rand is an integer
        {
            if ($linesfiles[$i] <= $rand)
            {
                $helper->disp("$infile[$i] has $linesfiles[$i] records, less than $rand. Will be reported as is...");
                $nochange = 1;
            }
            else
            {
                $helper->disp("Generating $rand random indices...");
                @randnum = Math::Random::random_permuted_index($linesfiles[$i]);
                $helper->disp("Getting required and writing output...");
                @indrem = sort {$a <=> $b} @randnum[0..$linesfiles[$i] - $rand -1];
            }
        }

        if ($savrem) # Save removed tags?
        {
            $helper->disp("Saving also removed tags...");
            my $savetags = $dir[$i].$base[$i]."_removed_tags.sav";
            open(SAVEREMTAGS,">$savetags");
        }

        $outfile[$i] = $dir[$i].$base[$i]."_norm".$ext[$i];
        open(INPUT,$infile[$i]) or die "\nThe file $infile[$i] does not exist!\n";
        open(OUTPUT,">$outfile[$i]");
        $j = 1;

        if ($nochange)
        {
            copy($infile[$i],$outfile[$i]);
        }
        else
        {
            while (<INPUT>)
            {
                #if (($. - 1) == $indrem[$j] && $j <= $#indrem)
                if ($. == $indrem[$j-1] && $j <= $#indrem)
                {
                    $j++;
                    print SAVEREMTAGS "$indrem[$j]"."\t".$_ if ($savrem);
                    #next;
                }
                else
                {
                    print OUTPUT $_;
                }
            }
        }

        close(INPUT);
        close(OUTPUT);
        close(SAVEREMTAGS) if ($savrem);
    }

    $helper->disp("Finished!\n\n");
}

=head2 sort_inputs

Sort input bed files if requested (not necessary) before bed downsampling

    $bednorm->sort_inputs(@infile);
    
=cut

sub sort_inputs
{
    my ($self,@infile) = @_;
    my $tmpfile;
    my $tmpdir = $helper->get("tmpdir");
    
    if ($^O !~ /MSWin/) # Case of linux, easy sorting
    {
        for (my $i=0; $i<@infile; $i++)
        {
            $helper->disp("Sorting bed file $infile[$i]...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp"."$i".".in$$");
            `sort -k1,1 -k2g,2 $infile[$i] > $tmpfile `;
            $infile[$i] = $tmpfile;
        }
    }
    else # We are in Windows... package required
    {
        $helper->try_module("File::Sort","sort_file");
        eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
        for (my $i=0; $i<@infile; $i++)
        {
            $helper->disp("Sorting file $infile[$i]...");
            $tmpfile = File::Spec->catfile($tmpdir,"temp"."$i".".tmp");
            File::Sort->sort_file(
            {
                I => $infile[$i],
                o => $tmpfile,
                k => ['1,1','2n,2'],
                t => "\t"
            });
            $infile[$i] = $tmpfile;
        }
    }

    return(@infile);
}

=head2 check_and_count

Check the bed file (if bed file) and count its lines

    $bednorm->check_and_count(@infile);
    
=cut

sub check_and_count
{
    my ($self,@infile) = @_;
    my @linesfiles;
    for (my $i=0; $i<@infile; $i++)
    {
        open(INPUT,$infile[$i]) or croak "\nThe file $infile[$i] does not exist!\n";
        my $fline = <INPUT>;
        my @check = split(/\t/,$fline);
        my $len = @check;
        if ($len != 3 && $len !=6 && $len !=12)
        {
            croak "\nThe file $infile[$i] does not appear to be a .bed file.\n";
        }
        sysseek(INPUT,0,0); # Go to the beginning of the file, since 1 line read for checking
        $helper->disp("Counting lines in file $infile[$i]..."); # Get number of lines
        my $totlines=0;
        $totlines += tr/\n/\n/ while sysread(INPUT,$_,2**16);
        seek(INPUT,0,0);
        $linesfiles[$i] = $totlines;
        close(INPUT);
    }
    return(@linesfiles);
}

=head2 get

HTS::Tools::Normalize::Bed object getter

    my $param_value = $bednorm->get("param_name")

=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Normalize::Bed object setter

    $bednorm->set("param_name","param_value")
    
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

1; # End of HTS::Tools::Normalize::Bed
