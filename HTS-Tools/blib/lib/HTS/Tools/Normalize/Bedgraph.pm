=head1 NAME

HTS::Tools::Normalize::Bedgraph - Normalize UCSC BedGraph files to a total signal

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module normalized a UCSC BedGraph signal file either to a total signal summed over the
entire genome or using some externally calculated normalization factors (e.g. from software
like DESeq or edgeR). This is a different normalization approach than read count downsampling
as it also takes into account read lenghts (downsampling read counts is not safe when there
are read of variable length as the removal of longer reads destroyes the uniformity of the 
process). The process is either to normalized the total signal sum over the genome to a given
constant (e.g. normalize signals of 1.2e+10 and 9.8e+8 to 1e+9) or use externally calculated
normalization factors such as the ones returned by the calcNormFactors function of the egdeR
Bioconductor package. The input is a set of BedGraph files and the output is the same set of
BedGraph files with normalized signal.

    use HTS::Tools::Normalize::Bedgraph;
    my %params = (
        'input' => ['normal_rnaseq_1.bedGraph','normal_rnaseq_2.bedGraph',
            'disease_rnaseq_1.bedGraph','disease_rnaseq_2.bedGraph'],
        'sumto' => '10000000000',
        'exportfactors' => 'normfactors.txt'
    )
    my $bg = HTS::Tools::Normalize::Bedgraph->new(\%params);
    $bg->run;
    
The acceptable parameters are as follows:

=over 4

=item input B<(required)>

Input bedgraph file(s). Please be careful as there is checking whether the input
file(s) are indeed bedgraph files. It's ok if they contain more than 4 columns
but the first four must be bedgraph (chromosome, start, end, signal separated by
tabs). Input files need not to be sorted.

=item output B<(optional)>

Output file names. It can be "stdout" for exporting to STDOUT, a set of file
names equal to the number of input files or nothing for the files to be
auto-generated.

=item sumto B<(optional)>

Normalize to --sumto total wig signal. Defaults to 1000000000. It is mutually
exclusive with --extnorm with --extnorm in precedence. It can be negative e.g.
for creating normalized stranded tracks where the minus strand shows negative
signal.

=item extnorm B<(optional)>

A set of external normalization factors (e.g. calculated from DESeq or edgeR).
It is mutually exclusive with --sumto with --extnorm in precedence.

=item exportfactors B<(optional)>

Export normalization factors and signal sums to a file specified by the parameter.

=item perlonly B<(optional)>

Use pure Perl to run the script, otherwise, uses Linux awk. Useful for e.g.
Windows systems but slower.

=item prerun B<(optional)>

If this switch is turned on, the script just counts the total wiggle signal in
the input files and prints it on the screen. This is useful in order for example
to determine the total normalization signal (sumto).

=item prerunlog B<(optional)>

Writes the output of prerun option to a file specified by prerunlog. If only the
prerunlog is specified, prerun is assumed and executed automatically.

=item ncores B<(optional)>

How many cores to use. Requires the presence of Parallel::Loops.

=item silent B<optional>

Do not display verbose messages.

=back

=head1 OUTPUT

BedGraph files with normalized signal.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Normalize::Bedgraph;

our $MODNAME = "HTS::Tools::Normalize::Bedgraph";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Bedgraph track normalization function.";

use v5.10;
use strict;
use warnings FATAL => 'all';

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

The HTS::Tools::Normalize::Bedgraph object constructor. It accepts a set of parameters that are required to run the
counter and get the output.

    my $bgnormer = HTS::Tools::Normalize::Bedgraph->new({'input' => ('myfile1.bedgraph','myfile2.bedgraph'),
        'sumto' => 5000000000});

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
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "normalize_bedgraph",
        "params" => $params});
    $params = $checker->validate;

    # After validating, bless and initialize
    bless($self,$class);
    $self->init($params);
    return($self);
}

=head2 init($params)

HTS::Tools::Normalize::Bedgraph object initialization method. NEVER use this directly, use new instead.

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

The HTS::Tools::Normalize::Bedgraph run subroutine. It runs the normalizer.

    $bgnorm->run;
    
=cut

sub run
{
    my $self = shift @_;
    
    my @bglist = @{$self->get("input")};
    my @extnorm = (defined($self->get("extnorm"))) ? (@{$self->get("extnorm")}) : (());
    my $sumto = $self->get("sumto");
    my $exportfacs = $self->get("exportfactors");
    my $perlonly = $self->get("perlonly");
    my @output = (defined($self->get("output"))) ? (@{$self->get("output")}) : (());
    my $prerun = $self->get("prerun");
    my $prerunlog = $self->get("prerunlog");
    my $ncores = $self->get("ncores");

    if (!@output)
    {
        foreach my $f (@bglist)
        {
            push(@output,$self->create_output_file($f));
        }
    }
    elsif ($helper->smatch("stdout",@output))
    {
        @output = ();
    }

    my ($line,$out,$chr,$start,$end,$signal);
    my (%wigsum,%normfactor);

    if ($perlonly)
    {
        # Calculate normalization factors or use external (e.g. from edgeR)
        if (@extnorm)
        {
            for (my $i=0; $i<@bglist; $i++)
            {
                $normfactor{$bglist[$i]} = $extnorm[$i];
            }
        }
        else
        {
            if ($ncores==1)
            {
                for (my $i=0; $i<@bglist; $i++)
                {
                    $helper->disp("Reading ".basename($bglist[$i])."...");
                    $wigsum{$bglist[$i]} = 0;
                    open(BGIN,$bglist[$i]);
                    while ($line = <BGIN>)
                    {
                        $helper->disp("  Read $. signals from ".basename($bglist[$i])."...") if ($.%1000000 == 0);
                        $line =~ s/\r|\n$//g;
                        ($chr,$start,$end,$signal) = split(/\t/,$line);
                        $wigsum{$bglist[$i]} += ($end - $start)*$signal;
                    }
                    close(BGIN);
                    $normfactor{$bglist[$i]} = $sumto/$wigsum{$bglist[$i]};
                }
            }
            else
            {
                my $pl1 = Parallel::Loops->new($ncores);
                $pl1->share(\%wigsum,\%normfactor);
                $pl1->foreach(\@bglist, sub {
                    $helper->disp("Reading ".basename($_)."...");
                    $wigsum{$_} = 0;
                    open(BGIN,$_);
                    while ($line = <BGIN>)
                    {
                        $line =~ s/\r|\n$//g;
                        ($chr,$start,$end,$signal) = split(/\t/,$line);
                        $wigsum{$_} += ($end - $start)*$signal;
                    }
                    close(BGIN);
                    $normfactor{$_} = $sumto/$wigsum{$_};
                });
            }
        }

        if (!$prerun && !$prerunlog)
        {
            # Reparse the bedgraph, rescale and write
            if ($ncores==1)
            {
                for (my $i=0; $i<@bglist; $i++)
                {
                    if ($output[$i])
                    {
                        $helper->disp("Writing ".basename($output[$i])."...");
                        open(OUTPUT,">$output[$i]");
                        $out = *OUTPUT;
                    }
                    else
                    {
                        $helper->disp("Writing to stdout...");
                        $out = *STDOUT;
                    }
                    open(BGIN,$bglist[$i]);
                    while ($line = <BGIN>)
                    {
                        $helper->disp("  Wrote $. normalized signals for ".basename($bglist[$i])."...") if ($.%1000000 == 0);
                        $line =~ s/\r|\n$//g;
                        ($chr,$start,$end,$signal) = split(/\t/,$line);
                        $signal = $helper->round($signal*$normfactor{$bglist[$i]},2);
                        print $out "$chr\t$start\t$end\t$signal\n";
                    }
                    close(BGIN);
                }
            }
            else
            {
                my $pl2 = Parallel::Loops->new($ncores);
                $pl2->share(\%wigsum,\%normfactor);
                $pl2->foreach(\@bglist, sub {
                    my $outp = $self->create_output_file($_);
                    $helper->disp("Writing ".basename($outp)."...");
                    open(OUTPUT,">$out");
                    $out = *OUTPUT;
                    
                    open(BGIN,$_);
                    while ($line = <BGIN>)
                    {
                        $line =~ s/\r|\n$//g;
                        ($chr,$start,$end,$signal) = split(/\t/,$line);
                        $signal = $helper->round($signal*$normfactor{$_},2);
                        print $out "$chr\t$start\t$end\t$signal\n";
                    }
                    close(BGIN);
                });
            }
        }
        else
        {
            if ($prerunlog)
            {
                open(LOG,">$prerunlog");
                print LOG "filename\ttotal signal\n";
                foreach my $f (keys(%wigsum))
                {
                    print LOG basename($f)."\t$wigsum{$f}\n";
                }
                close(LOG);
            }
            else
            {
                $helper->disp("\nfilename\ttotal signal");
                foreach my $f (keys(%wigsum))
                {
                    $helper->disp(basename($f)."\t$wigsum{$f}");
                }
                $helper->disp("")
            }
        }
    }
    else # Use awk
    { 
        # Calculate normalization factors or use external (e.g. from edgeR)
        if (@extnorm)
        {
            for (my $i=0; $i<@bglist; $i++)
            {
                $normfactor{$bglist[$i]} = $extnorm[$i];
            }
        }
        else
        {
            if ($ncores==1)
            {
                for (my $i=0; $i<@bglist; $i++)
                {
                    $helper->disp("Reading ".basename($bglist[$i])."...");
                    $wigsum{$bglist[$i]} = `awk '{sum += (\$3-\$2)*\$4; if (FNR \% 1000000==0) { printf("\\n  Read %d signals from %s",FNR,FILENAME) | "cat 1>&2"}} END {print sum}' $bglist[$i]`;
                    $normfactor{$bglist[$i]} = $sumto/$wigsum{$bglist[$i]};
                }
            }
            else
            {
                my $pl3 = Parallel::Loops->new($ncores);
                $pl3->share(\%wigsum,\%normfactor);
                $pl3->foreach(\@bglist, sub {
                    $helper->disp("Reading ".basename($_)."...");
                    $wigsum{$_} = `awk '{sum += (\$3-\$2)*\$4; } END {print sum}' $_`;
                    $normfactor{$_} = $sumto/$wigsum{$_};
                });
            }
        }

        if (!$prerun && !$prerunlog)
        {
            # Reparse the bedgraph, rescale and write
            if ($ncores==1)
            {
                for (my $i=0; $i<@bglist; $i++)
                {
                    if ($output[$i])
                    {
                        $helper->disp("Writing ".basename($output[$i])."...");
                        #`awk -v var="$output[$i]" '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4*$normfactor{$bglist[$i]}; if (FNR \% 1000000==0) { printf("\\n  Wrote %d signals to %s",FNR,var) | "cat 1>&2"}}' $bglist[$i] > $output[$i]`;
                        `awk -v var="$output[$i]" '{printf \"%s\\t%s\\t%s\\t%.2f\\n\", \$1,\$2,\$3,\$4*$normfactor{$bglist[$i]}; if (FNR \% 1000000==0) { printf("\\n  Wrote %d signals to %s",FNR,var) | "cat 1>&2"}}' $bglist[$i] > $output[$i]`;
                        
                    }
                    else {
                        $helper->disp("Writing to stdout...");
                        #`awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4*$normfactor{$bglist[$i]}}' $bglist[$i]`;
                        `awk '{printf \"%s\\t%s\\t%s\\t%.2f\\n\", \$1,\$2,\$3,\$4*$normfactor{$bglist[$i]}}' $bglist[$i]`;
                    }
                }
            }
            else
            {
                my $pl4 = Parallel::Loops->new($ncores);
                $pl4->share(\%wigsum,\%normfactor);
                $pl4->foreach(\@bglist, sub {
                    my $outp = $self->create_output_file($_);
                    $helper->disp("Writing ".basename($outp)."...");
                    `awk -v var="$outp" '{printf \"%s\\t%s\\t%s\\t%.2f\\n\", \$1,\$2,\$3,\$4*$normfactor{$_};}' $_ > $outp`;
                });
            }
        }
        else
        {
            if ($prerunlog)
            {
                open(LOG,">$prerunlog");
                print LOG "filename\ttotal signal\n";
                foreach my $f (keys(%wigsum))
                {
                    print LOG basename($f)."\t$wigsum{$f}";
                }
                close(LOG);
            }
            else
            {
                print STDERR "\nfilename\ttotal signal\n";
                foreach my $f (keys(%wigsum))
                {
                    print STDERR basename($f)."\t$wigsum{$f}";
                }
            }
        }
    }

    if ($exportfacs)
    {
        $helper->disp("Exporting normalization factors...");
        open(FACS,">$exportfacs");
        print FACS "file\tnormalization factor\ttotal signal\n";
        foreach my $f (keys(%normfactor))
        {
            print FACS basename($f)."\t".$helper->round($normfactor{$f},6)."\t$wigsum{$f}";
            print FACS "\n" if ($perlonly); # For some reason, awk maintains the EOL
        }
        close(FACS);
    }

    $helper->disp("Finished!\n\n");
    
    return(@output);
}

=head2 create_output_file

The HTS::Tools::Normalize::Bedgraph output filename creation subroutine. It constructs output filenames for
normalized bedgraphs if they are not provided.

    $bgnorm->create_output_file($infile);
    
=cut

sub create_output_file
{
    my ($self,$in) = @_;
    my ($base,$dir,$ext) = fileparse($in,'\.[^.]*');
    return(File::Spec->catfile($dir,$base."_norm".$ext));
}

=head2 get

HTS::Tools::Normalize::Bedgraph object getter

    my $param_value = $bgnorm->get("param_name")

=cut

sub get
{
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set

HTS::Tools::Normalize::Bedgraph object setter

    $bgnorm->set("param_name","param_value")
    
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

1; # End of HTS::Tools::Normalize::Bedgraph
