=head1 NAME

HTS::Tools::Normalize::Bedgraph - The great new HTS::Tools::Normalize::Bedgraph!

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use HTS::Tools::Normalize::Bedgraph;

    my $foo = HTS::Tools::Normalize::Bedgraph->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

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
    my $checker = HTS::Tools::Paramcheck->new({"tool" => "normalize", "params" => $params});
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

    if (!@output)
    {
        foreach my $f (@bglist)
        {
            push(@output,$self->create_output_file($f));
        }
    }
    elsif ("stdout" ~~ (@output))
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

        if (!$prerun && !$prerunlog)
        {
            # Reparse the bedgraph, rescale and write
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
                    $signal = $helper->round($signal*$normfactor{$bglist[$i]},6);
                    print $out "$chr\t$start\t$end\t$signal\n";
                }
                close(BGIN);
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
            for (my $i=0; $i<@bglist; $i++)
            {
                $helper->disp("Reading ".basename($bglist[$i])."...");
                $wigsum{$bglist[$i]} = `awk '{sum += (\$3-\$2)*\$4; if (FNR \% 1000000==0) { printf("\\n  Read %d signals from %s",FNR,FILENAME) | "cat 1>&2"}} END {print sum}' $bglist[$i]`;
                $normfactor{$bglist[$i]} = $sumto/$wigsum{$bglist[$i]};
            }
        }

        if (!$prerun && !$prerunlog)
        {
            # Reparse the bedgraph, rescale and write
            for (my $i=0; $i<@bglist; $i++)
            {
                if ($output[$i])
                {
                    $helper->disp("Writing ".basename($output[$i])."...");
                    `awk -v var="$output[$i]" '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4*$normfactor{$bglist[$i]}; if (FNR \% 1000000==0) { printf("\\n  Wrote %d signals to %s",FNR,var) | "cat 1>&2"}}' $bglist[$i] > $output[$i]`;
                }
                else {
                    $helper->disp("Writing to stdout...");
                    `awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4*$normfactor{$bglist[$i]}}' $bglist[$i]`;
                }
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
