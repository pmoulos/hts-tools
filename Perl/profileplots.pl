#!/usr/bin/perl

# Program to create various profile plots for NGS experiments. More documentation later
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 08 - 03 - 2013 (dd - mm - yyyy)
# Last Update : XX - XX - XXXX (dd - mm - yyyy)
# Version     : 1.0
#
# CHANGELOG
#
# TODO: Use the --split property to create TSS profile heatmaps and plots. Continue implementation of this...
# TODO: 5'UTR = cdsStart - txStart or chr:txStart-cdsStart, 3'UTR = txEnd - cdsEnd or chr:cdsEnd-txEnd

# The 5' UTR is the region between the txStart and cdsStart and the 3' UTR is the region between the cdsEnd and the txEnd

use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
use File::Spec;
use File::Path qw(make_path remove_tree);
use Parallel::Loops;
use LWP::UserAgent;

use constant MAXCORES => 12;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# On Ctrl-C or die, do cleanup
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "profileplots.pl";
our @infile;
our $regionfile;
our @profile;
our @proftype;
our @width;
our $sort = 0;
our $binsize;
our $binumber;
our $ncore = 1;
our $silent = 0;
our $help = 0;
# readcounter.pl specific, see it for help
our $percentage = 0.95;
our $lscore = 0;
our $escore = 0;
our $c = 3;
our $countsmall;
our $output;

# Advertise
&advertise;

# Check inputs
&check_inputs;

# Global variables
our $tmpdir = File::Temp->newdir();

# MAIN IDEA
# Use the --split property of readcounter.pl to create profile (tss, genebody etc.) heatmaps
# and plots. Use ensembl's interface to get annotation or even at some point in the future
# a connection to UCSC or local for RefSeq (see chipseqPQC or covsat.pl for this).

$regionfile = &prepare_tss($regionfile) if (@profile);

sub check_inputs
{
	my $stop;
    GetOptions("input|i=s{,}" => \@infile,
    		   "region|r=s" => \$regionfile,
    		   "profile|x=s{,}" => \@profile,
    		   "proftype|y=s{,}" => \@proftype,
    		   "sort|a" => \$sort,
    		   "sbin|b=i" => \$binsize,
    		   "nbin|n=i" => \$binumber,
    		   "ncore|t=i" => \$ncore,
    		   "percent|p=f" => \$percentage,
    		   "lscore|l" => \$lscore,
    		   "escore|e" => \$escore,
    		   "constant|c=f" => \$c,
    		   "small|m" => \$countsmall,
               #"output|o=s" => \$output,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&program_usage;
    	exit;
    }
    my @availableprofiles = ("tss","genebody","exon","5utr","3utr"); # More to come...
    $stop .= "--- Please specify input file(s) ---\n" if (!@infile);
    $stop .= "--- Please specify region file ---\n" if (!$regionfile);
    $stop .= "--- The supported genomes so far are only human-type, mouse-type, rat-type or fly-type, zebrafish-type where type is gene or exon ---\n"
		if ( !-f $regionfile && $regionfile !~ m/human-(gene|exon)|mouse-(gene|exon)|rat-(gene|exon)|fly-(gene|exon)|zebrafish-(gene|exon)/i);
    foreach my $pro (@profile)
	{
		#if ($pro !~ m/tss/i && $pro !~ m/genebody/i && $pro !~ m/exon/i)
		if (!($pro ~~ @availableprofiles))
		{
			$stop .= "--- --profile can be ".join(", ",@availableprofiles)." or more than one ---";
			last;
		}
	}
    if ($stop)
    {
            print STDERR "\n$stop\n";
            print STDERR "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }	
	if (@proftype)
	{
		foreach my $prt (@proftype)
		{
			if ($prt !~ m/heatmap/i && $prt !~ m/coverage/i)
			{
				disp("--profile can be \"heatmap\", \"coverage\" or both. Assuming default (coverage)...");
				@proftype = ("coverage");
				last;
			}
		}
	}
	if (@profile && !@width)
	{
		disp("Width for TSS (or start of a region) profile not given. Assuming default (1000 upstream and downstream)...");
		@width = (1000,1000);
	}
	# Silently assume the same upstream and downstream distance if only one number given
	@width = ($width[0],$width[0]) if (@width && scalar @width<2);
}

sub prepare_tss
{
	my $regionfile = shift @_;
	my $tssfile = File::Spec->catfile($tmpdir,"tss_$$.txt");
	my ($theHeader,$line,$chr,$start,$end,$uid,$score,$strand,$rest);
	my @arest;

	disp("Preparing profile file...");
	open(REGIONS,$regionfile) or die "\nThe file $regionfile does not exist!\n";
	open(FORTSS,">$tssfile");
	$line = <REGIONS>;
	$theHeader = &decide_header($line);
	($theHeader) ? (print FORTSS $theHeader."\n") : (seek(REGIONS,0,0));
	while ($line = <REGIONS>)
	{
		$line =~ s/\r|\n$//g; # Make sure to remove carriage returns
		($chr,$start,$end,$uid,$score,$strand,@arest) = split(/\t/,$line);
		$rest = join("\t",@arest);
		($start,$end) = &swap($start,$end) if ($strand eq "-" || $strand == -1 || $strand eq "R");
		$start = $start - $width[0];
		$end = $start + $width[0] + $width[1];
		($rest eq "") ? (print FORTSS "$chr\t$start\t$end\t$uid\t$score\t$strand\n") :
			(print FORTSS "$chr\t$start\t$end\t$uid\t$score\t$strand\t$rest\n");
	}
	close(REGIONS);
	close(FORTSS);

	($theHeader) ? ($tssfile = &sort_ensembl_genes($tssfile)) :
		($tssfile = &sort_one($tssfile));

	return($tssfile);
}

sub construct_optargs
{
	my $str;
	$str.= " --sort" if ($sort);
	$str.= " --percent $percentage" if ($percentage);
	$str.= " --lscore $lscore" if ($lscore);
	$str.= " --escore $escore" if ($escore);
	$str.= " --constant $c" if ($c);
	$str.= " --small" if ($countsmall);
	$str.= " --output $output" if ($output);
	return($str);
}

sub swap
{
	return($_[1],$_[0]);
}

sub disp
{
	print STDERR "\n@_" if (!$silent);
}

sub catch_cleanup 
{
	print STDERR "\nCatching ctrl-C, cleaning temporary files!";
	&cleanup;
	die;
}

sub cleanup 
{
	remove_tree($tmpdir);
}

sub advertise
{
	use Term::ANSIColor;
	print color 'bold yellow on_blue';
	disp($scriptname);
	print color 'bold green';
	disp("Read profile plot generator for NGS experiments... Copyright: Panagiotis Moulos (moulos\@fleming.gr)\n");
	print color 'reset';
}
