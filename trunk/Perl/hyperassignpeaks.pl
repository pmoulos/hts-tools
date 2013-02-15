#!/usr/bin/perl

# Description later...
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 13 - 01 - 2009 (dd - mm - yyyy)
# Last Update : 30 - 01 - 2013 (dd - mm - yyyy)
# Version     : 1.3
#
# Version 1.1 supports multiple outputs (simple stats, pretty shaped or .gff files)
# Version 1.11 has a bug fix in exporting gff-gene files and pretty-gene files
# Version 1.2 can also work with chi-square test instead of only hypergeometric
# Version 1.3 can run without a statistical test in which case, no background file is
# needed. Fixed a distance bug in gff-gene export, where a gene having more than one
# peak had the same distance for both peaks. Removed waitbar, total time-memory loss
#

use strict;
use File::Basename;
use Getopt::Long;
 			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "hyperassignpeaks.pl";
our @peakfile;		   # The peak file(s) (BED format)
our $sigfile;	       # Significant regions to be assigned to peaks (BED format)
our $backfile;		   # Background regions to be assigned to peaks (BED format)
our @span;		  	   # Upstream and downstream  (default +/-100k from TSS)
our @sbcols;		   # Columns containing unique sig/back region ID and strand (default (4,5))
our @pcols;		       # Columns containing unique peak region ID and mode (default (4,5))
our @out;			   # Output filetype
our $test = "hypgeom"; # Default hypergeometric test
our $pval = 0.05; 	   # Hypergeometric p-value cutoff (default 0.05)
our $header = 0;	   # ALL files contain header line? 
our $silent = 0;	   # Display verbose messages
our $help = 0;		   # Help?

# Check inputs
&checkInputs;

# Check for existence of Math::Cephes in case of chi-square test
if ($test eq "chi2")
{
	eval
	{
		require Math::Cephes;
	};
	if ($@)
	{
		my $butch = "Module Math::Cephes is required to perform selected chi square test. If you\n". 
					"are in Windows and you have ActiveState Perl installed, use the Package Manager\n".
					"to get the module. If you are under Linux, log in as a super user (or use\n".
					"sudo under Ubuntu) and type \"perl -MCPAN -e shell\" (you will possibly have\n".
					"to answer some questions). After this type \"install Math::Cephes\" to\n".
					"install the module. If you don't know how to install the package, contact\n".
					"your system administrator or use another over-representation test.";
		die "\n$butch\n";
	}
	else
	{
		eval "use Math::Cephes qw(:dists)";	
	}
}

# Check for existence of Tie::IxHash::Easy in case of matrix generation
if (@out ~~ /matrix/)
{
	eval
	{
		require Tie::IxHash::Easy;
	};
	if ($@)
	{
		my $murd = "Module Tie::IxHash::Easy is required to continue with the execution. If you\n". 
					"are in Windows and you have ActiveState Perl installed, use the Package Manager\n".
					"to get the module. If you are under Linux, log in as a super user (or use\n".
					"sudo under Ubuntu) and type \"perl -MCPAN -e shell\" (you will possibly have\n".
					"to answer some questions). After this type \"install Math::Cephes\" to\n".
					"install the module. If you don't know how to install the package, contact\n".
					"your system administrator or use another over-representation test.";
		die "\n$murd\n";
	}
	else
	{
		eval "use Tie::IxHash::Easy";	
	}
}

# Bavard
disp("HyperAssignPeaks - A Perl script to assign peaks to genes based on over-representation");
disp("Copyright 2009 - Panagiotis Moulos - All rights reserved\n");
my $date = &now;
disp("$date - Started...\n");

# More bavard... I like information...
disp("Genomic span from TSS : ($span[0],$span[1])");
disp("Over-representation test : $test");
disp("p-value threshold : $pval") if ($test ne "none");
disp("Chosen output(s) : ",join("\, ",@out),"\n");

# General indices
my ($i,$j,$k);
# Some info about gff output
my $gffreq = my $pdataout = 0;
foreach my $o (@out)
{
	$gffreq = 1 if ($o =~ /gff/);
	$pdataout = 1 if ($o =~ /peakdata/);
}

# Some intialization
my (%sigID,%sigStart,%sigEnd);
my (%backID,%backStart,%backEnd);
my (@all,$chr,$start,$end,$id,$strand);
my $lensig = my $lenback = 0;
my @lenpeak;

# Variable for matrix generation
my %hasPeak;
tie %hasPeak, "Tie::IxHash::Easy";

# Suck in significant region file
open (SIG,$sigfile) or die "\nThe file $sigfile does not exist!\n";;
disp("Reading region file $sigfile...");
my $headerline = <SIG> if ($header);
while (my $line = <SIG>)
{
	next if ($line =~/^chrM/);
	next if ($line =~/rand/);
	$line =~ s/\r|\n$//g;
	@all = split(/\t/,$line);
	$chr = $all[0];
	$start = $all[1];
	$end = $all[2];
	$id = $all[$sbcols[0]];
	$strand = $all[$sbcols[1]];
	if ($strand == 1 || $strand eq "+" || $strand eq "F")
	{
		push(@{$sigStart{$chr}},$start);
		push(@{$sigEnd{$chr}},$end);
	}
	elsif ($strand == -1 || $strand eq "-" || $strand eq "R")
	{
		push(@{$sigStart{$chr}},$end);
		push(@{$sigEnd{$chr}},$start);
	}
	else # Some self-defense
	{
		disp("Improper strand format... Skipping line $. from $sigfile");
		next;
	}
	push(@{$sigID{$chr}},$id);
	$lensig++;

	# Initiate the hash to keep record of peaks for a peak matrix file generation
	if (@out ~~ /matrix/)
	{
		for ($i=0; $i<@peakfile; $i++)
		{
			$hasPeak{$id}{basename($peakfile[$i])} = ();
		}
	}
}
close(SIG);

# Suck in background region file if test is to be performed
if ($test ne "none")
{
	open (BACK,$backfile) or die "\nThe file $backfile does not exist!\n";
	disp("Reading background file $backfile...");
	$headerline = <BACK> if ($header);
	while (my $line = <BACK>)
	{
		next if ($line =~/^chrM/);
		next if ($line =~/rand/);
		$line =~ s/\r|\n$//g;
		@all = split(/\t/,$line);
		$chr = $all[0];
		$start = $all[1];
		$end = $all[2];
		$id = $all[$sbcols[0]];
		$strand = $all[$sbcols[1]];
		if ($strand == 1 || $strand eq "+" || $strand eq "F")
		{
			push(@{$backStart{$chr}},$start);
			push(@{$backEnd{$chr}},$end);
		}
		elsif ($strand == -1 || $strand eq "-" || $strand eq "R")
		{
			push(@{$backStart{$chr}},$end);
			push(@{$backEnd{$chr}},$start);
		}
		else # Some self-defense
		{
			disp("Improper strand format... Skipping line $. from $backfile");
			next;
		}
		push(@{$backID{$chr}},$id);
		$lenback++;
	}
	close(BACK);
}

# Some self-defense... Check uniqueness of gene IDs...
die "\nThe region IDs in regions file $sigfile are not unique! Exiting...\n" if (!&checkUnique(%sigID));
die "\nThe region IDs in background file $backfile are not unique! Exiting...\n" if (!&checkUnique(%backID) && $test ne "none");

# Read and process peak files
for ($i=0; $i<@peakfile; $i++)
{
    disp("\nReading file $peakfile[$i]...");
    
    my (%peakID,%peakMode,%countSig,%countBack,%allPeakData,%countGenesSig,%countGenesBack);
    my ($pstart,$pend,%peakStarts,%peakEnds) if ($gffreq);
    my (@pall,$pchr,$pmode,$pid);
    my (@starts,@ends,@modes);
    my (@peakchrs,@peakids,@geneids);
    my (@sigallpeakids,@backallpeakids,@sigassgenes,@backassgenes);
    my ($currchr,$currdist,$currpeak,$currgene,$currout,@currgenes);
    my (%genePeaksSig,%peaksGenesSig,%distPeakBased,%distGeneBased,%peakIndex,%geneIndex);#,%genePeaksBack);
    my ($p,$cp,@elems);
    my %finalPeaks;
    
	open(PEAK,$peakfile[$i]) or die "\nThe file $peakfile[$i] does not exist!\n";
	$headerline = <PEAK> if ($header);
	while (my $line = <PEAK>)
	{
		$line =~ s/\r|\n$//g;
		@pall = split(/\t/,$line);
		$pchr = $pall[0];
		$pid = $pall[$pcols[0]];
		$pmode = $pall[$pcols[1]];
		push(@{$peakID{$pchr}},$pid);
		push(@{$peakMode{$pchr}},$pmode);
		if ($gffreq) # Starts and ends required in this case for GFF files
		{
			$pstart = $pall[1];
			$pend = $pall[2];
			push(@{$peakStarts{$pchr}},$pstart);
			push(@{$peakEnds{$pchr}},$pend);
		}
		$allPeakData{$pid} = join("\t",@pall) if ($pdataout);
		$lenpeak[$i]++;
	}
	close(PEAK);
	
	disp("Processing file $peakfile[$i]...");
	die "\nThe region IDs in file $peakfile[$i] are not unique! Exiting...\n" if (!&checkUnique(%peakID));
	
	@peakchrs = keys(%peakMode);
	
	# Do stuff with siginificant file
	if ($test ne "none")
	{
		disp("Associating peaks with genes in significant and background region files :");
		disp("Significant region file : $sigfile");
		disp("Background region file : $backfile");
	}
	else
	{
		disp("Associating peaks with genes in region file : $sigfile");
	}
	foreach $currchr (@peakchrs)
	{	
		disp("Peaks at $currchr...");
		
		@modes = @{$peakMode{$currchr}};
		@peakids = @{$peakID{$currchr}};
		
		if ($sigID{$currchr}) # Could not have significant genes at chromosome
		{
			@geneids = @{$sigID{$currchr}};
			@starts = @{$sigStart{$currchr}};
			@ends = @{$sigEnd{$currchr}};
		
			for ($j=0; $j<@starts; $j++)
			{
				for ($k=0; $k<@modes; $k++)
				{
					$currdist = &dist($starts[$j],$ends[$j],$modes[$k]);
					if ($currdist > $span[0] && $currdist < $span[1])
					{
						push(@{$genePeaksSig{$currchr}{$geneids[$j]}},$peakids[$k]);
						push(@{$peaksGenesSig{$currchr}{$peakids[$k]}},$geneids[$j]);
						push(@{$distGeneBased{$currchr}{$geneids[$j]}},$currdist);
						push(@{$distPeakBased{$currchr}{$peakids[$k]}},$currdist);
						push(@{$peakIndex{$currchr}{$peakids[$k]}},$k);
						push(@{$geneIndex{$currchr}{$geneids[$j]}},$j);
						push(@sigallpeakids,$peakids[$k]);
						push(@sigassgenes,$geneids[$j]);
						#$hasPeak{$geneids[$j]}{basename($peakfile[$i])}++;
						push(@{$hasPeak{$geneids[$j]}{basename($peakfile[$i])}},$peakids[$k]."_".$currdist);
					}
				}
			}
		}

		if ($test ne "none") # If no test performed, no need to do anything with background file
		{
			if ($backID{$currchr}) # Could not have genes at chromosome (unlikely for background...)
			{
				@geneids = @{$backID{$currchr}};
				@starts = @{$backStart{$currchr}};
				@ends = @{$backEnd{$currchr}};
			
				for ($j=0; $j<@starts; $j++)
				{
					for ($k=0; $k<@modes; $k++)
					{
						$currdist = &dist($starts[$j],$ends[$j],$modes[$k]);
						if ($currdist > $span[0] && $currdist < $span[1])
						{
							#push(@{$genePeaksBack{$currchr}{$geneids[$j]}},$peakids[$k]);
							push(@backallpeakids,$peakids[$k]);
							push(@backassgenes,$geneids[$j]);
						}
					}
				}
			}
		}
	}
	
	# Get peak counts in significant and background file
	%countSig = &unique(@sigallpeakids);
	%countBack = &unique(@backallpeakids) if ($test ne "none");
	%countGenesSig = &unique(@sigassgenes);
	%countGenesBack = &unique(@backassgenes) if ($test ne "none");
	
	# Run hypergeometric test
	if ($test eq "hypgeom")
	{
		disp("Running hypergeometric test for each assigned peak...");	
		@elems = keys(%countSig);
		for ($j=0; $j<@elems; $j++)
		{
			$p = abs(1 - &hypergeomCDF($lensig,$lenback - $lensig,$countBack{$elems[$j]},$countSig{$elems[$j]}));
			($p*@elems > 1) ? ($cp = 1) : ($cp = $p*@elems); # Bonferroni type correction
			$finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}/$countBack{$elems[$j]}" if ($p < $pval);
		}
	}
	elsif ($test eq "chi2")
	{
		disp("Running chi-square test for each assigned peak...");	
		@elems = keys(%countSig);
		for ($j=0; $j<@elems; $j++)
		{
			$p = &chisquarecont($countSig{$elems[$j]},$lensig - $countSig{$elems[$j]},$countBack{$elems[$j]},$lenback - $countBack{$elems[$j]});
			($p*@elems > 1) ? ($cp = 1) : ($cp = $p*@elems); # Bonferroni type correction
			$finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}/$countBack{$elems[$j]}" if ($p < $pval);
		}
	}
	elsif ($test eq "none")
	{
		disp("No statistical testing performed...");	
		@elems = keys(%countSig);
		for ($j=0; $j<@elems; $j++)
		{
			$p = "NA";
			$cp = "NA";
			$finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}";
		}
	}
	
	my $sp = keys(%finalPeaks);
	my $ap = @elems;
	my $bp = keys(%countBack) if ($test ne "none");
	my $sag = keys(%countGenesSig);
	my $sbg = keys(%countGenesBack) if ($test ne "none");
	if ($test ne "none")
	{
		disp("\nAssigned peaks in significant list : $ap out of $lenpeak[$i] peaks in $sag out of $lensig genes");
		disp("Assigned peaks in background : $bp out of $lenpeak[$i] peaks in $sbg out of $lenback genes");
		disp("Over-represented at p-value<$pval : $sp out of $ap\n");
	}
	else
	{
		disp("\nAssigned peaks in gene list : $ap out of $lenpeak[$i] peaks in $sag out of $lensig genes");
	}
	
	# Free some memory...
	($sp,$ap,$bp,$sag,$sbg,@sigassgenes,@backassgenes,%countGenesSig,%countGenesBack) = 
	(undef,undef,undef,undef,undef,undef,undef,undef,undef);
	
	# Construct output
	foreach my $opt (@out)
	{
		if ($opt eq "stats")
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			my $co;
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			foreach $co (sort(keys(%finalPeaks)))
			{
				print OUTPUT "$co\t$finalPeaks{$co}\n";
			}
		}
		if ($opt =~ /gff-peak/)
		{ 
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "Chromosome\tProgram\tFeature\tStart\tEnd\tp-value\tStrand\tFrame\t".
						 "GeneID\tPeakID\tCorrected p-value\tEnrichment\tDistance\n" if ($opt eq "gff-peak-db");
			foreach $currchr (@peakchrs)
			{
				foreach $currpeak (keys(%finalPeaks))
				{	
					if ($peaksGenesSig{$currchr}{$currpeak})
					{
						@currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};
						my ($q,$cq,$r) = split(/\t/,$finalPeaks{$currpeak});
						my $pos = 0;
						foreach $currgene (@currgenes)
						{	
							#my $ct = "+";
							#my $cd = ${$distPeakBased{$currchr}{$currpeak}}[$pos];
							#$ct = "-" if ($cd > 0);
							$currout = "$currchr\tHyperAssignPeaks\tPeak\t".
									   "${$peakStarts{$currchr}}[${$peakIndex{$currchr}{$currpeak}}[$pos]]\t".
									   "${$peakEnds{$currchr}}[${$peakIndex{$currchr}{$currpeak}}[$pos]]\t".
									   "$q\t+\t\.\t$currgene\t$currpeak\t$cq\t$r\t".
									   "${$distPeakBased{$currchr}{$currpeak}}[$pos]\n";
							print OUTPUT "$currout";
							$pos++;
						}
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt =~ /gff-gene/)
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "Chromosome\tProgram\tFeature\tStart\tEnd\tp-value\tStrand\tFrame\t".
						 "PeakID\tGeneID\tCorrected p-value\tEnrichment\tDistance\n" if ($opt eq "gff-gene-db");
			foreach $currchr (@peakchrs)
			{
				foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
				{	
					my $pos = 0;
					foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
					{
						#my $pos = 0; # It was a tricky one bug...
						if ($finalPeaks{$currpeak})
						{
							my ($q,$cq,$r) = split(/\t/,$finalPeaks{$currpeak});
							my $s = ${$sigStart{$currchr}}[${$geneIndex{$currchr}{$currgene}}[$pos]];
							my $e = ${$sigEnd{$currchr}}[${$geneIndex{$currchr}{$currgene}}[$pos]];
							my $st = "+";
							if ($s > $e)
							{
								($s,$e) = ($e,$s); # Swap them
								$st = "-";
							}
							$currout = "$currchr\tHyperAssignPeaks\tGene\t$s\t$e\t$q\t$st\t\.\t$currpeak\t".
									   "$currgene\t$cq\t$r\t${$distGeneBased{$currchr}{$currgene}}[$pos]\n";
							print OUTPUT "$currout";
						}
						$pos++;
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "peak")
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			foreach $currchr (@peakchrs)
			{
				foreach $currpeak (keys(%finalPeaks))
				{	
					if ($peaksGenesSig{$currchr}{$currpeak})
					{
						@currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};			
						print OUTPUT "$currpeak\t",join("\ ",@currgenes),"\n";
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "gene")
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			foreach $currchr (@peakchrs)
			{
				foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
				{	
					my @outpeaks;
					foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
					{
						push(@outpeaks,$currpeak) if ($finalPeaks{$currpeak});
					}
					print OUTPUT "$currgene\t",join("\ ",@outpeaks),"\n"  if (@outpeaks);
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "pretty-peak")
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "PeakID/GeneID\tp-value/Distance\tBonferroni p-value/Strand\tEnrichment\n\n";
			foreach $currchr (@peakchrs)
			{
				foreach $currpeak (keys(%finalPeaks))
				{	
					if ($peaksGenesSig{$currchr}{$currpeak})
					{
						print OUTPUT "$currpeak\t$finalPeaks{$currpeak}\n";
						@currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};
						my $pos = 0;
						foreach $currgene (@currgenes)
						{	
							my $ct = "+";
							my $cd = ${$distPeakBased{$currchr}{$currpeak}}[$pos];
							$ct = "-" if ($cd > 0);
							print OUTPUT "$currgene\t$cd\t$ct\t\n";
							$pos++;
						}
						print OUTPUT "\n";
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "pretty-gene")
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "GeneID/PeakID\tStrand/p-value\tBonferroni p-value\tEnrichment\tDistance\n\n";
			foreach $currchr (@peakchrs)
			{
				foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
				{	
					my $s = ${$sigStart{$currchr}}[${$geneIndex{$currchr}{$currgene}}[0]];
					my $e = ${$sigEnd{$currchr}}[${$geneIndex{$currchr}{$currgene}}[0]];
					my $st = "+";
					if ($s > $e)
					{
						($s,$e) = ($e,$s); # Swap them
						$st = "-";
					}
					my (@outpeaks,@dis);
					my $pos = 0;
					foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
					{
						if ($finalPeaks{$currpeak})
						{
							push(@outpeaks,$currpeak);
							push(@dis,${$distGeneBased{$currchr}{$currgene}}[$pos]);		
						}
						$pos++;
					}
					if (@outpeaks)
					{
						print OUTPUT "$currgene\t$st\t\t\t\n";
						for (my $cpo=0; $cpo<@outpeaks; $cpo++)
						{
							print OUTPUT "$outpeaks[$cpo]\t$finalPeaks{$outpeaks[$cpo]}\t$dis[$cpo]\n";
						}
						print OUTPUT "\n";
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "peakdata")
		{
			my $outfile = &createOutputFile($peakfile[$i],$opt);
			disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "$headerline" if ($header);
			my @fkeys = keys(%finalPeaks);
			foreach $currpeak (@fkeys)
			{
				print OUTPUT "$allPeakData{$currpeak}\n";
			} 
		}	
		&printGeneOrPeak($peakfile[$i],"all-peak",%peaksGenesSig) if ($opt eq "all-peak");
		&printGeneOrPeak($peakfile[$i],"all-gene",%genePeaksSig) if ($opt eq "all-gene");
	}
}

&printMatrix(\%hasPeak,"matrix-number") if (@out ~~ /matrix-number/);
&printMatrix(\%hasPeak,"matrix-presence") if (@out ~~ /matrix-presence/);
&printMatrix(\%hasPeak,"matrix-peaks") if (@out ~~ /matrix-peaks/);

$date = &now;
disp("$date - Finished!\n\n");


sub dist
{
	my $s = shift @_;
	my $e = shift @_;
	my $m = shift @_;
	my $d;
	if ($s < $e) # + strand
	{
		#(($m > $s) && ($m < $e)) ? ($d = 0) : ($d = $m - $s);
		$d = $m - $s;
	}
	else # - strand, reversed while reading file
	{
		#(($m < $s) && ($m > $e)) ? ($d = 0) : ($d = $s - $m);
		$d = $s - $m;
	}
	return($d);
}

sub unique
{
	my @list = @_;
	my (%seen,$item);
	foreach $item (@list) 
	{
		$seen{$item}++;
	}
	return(%seen);
}

# Hypergeometric probabilty density function (pdf)
# There are m "bad" and n "good" balls in an urn.
# Pick N of them. The probability of i or more successful selections:
# (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
sub hypergeomPDF
{   
	my ($n,$m,$N,$i) = @_;
	my $loghyp1 = logfact($m) + logfact($n) + logfact($N) + logfact($m+$n-$N);
	my $loghyp2 = logfact($i) + logfact($n-$i) + logfact($m+$i-$N) + logfact($N-$i) + logfact($m+$n);
	return(exp($loghyp1 - $loghyp2));
}

# Hypergeometric cumulative distribution function (cdf)
sub hypergeomCDF
{  
	my ($n,$m,$N,$i) = @_; 
	my @lessthan = (0..$i);
	my $cum = 0; #;-)
	foreach my $j (@lessthan)
	{
		$cum += hypergeomPDF($n,$m,$N,$j);
	}
	return($cum);
}

sub logfact 
{
	return gammaln(shift(@_) + 1.0);
}

sub gammaln 
{
	my $xx = shift;
	my $eps = 1; 
	my @cof = (76.18009172947146, -86.50532032941677,
			   24.01409824083091, -1.231739572450155,
			   0.12086509738661e-2, -0.5395239384953e-5);
	my $y = my $x = $xx;
	my $tmp = $x + 5.5;
	$tmp -= ($x + .5) * log($tmp);
	my $ser = 1.000000000190015;
	for my $j (0..5) 
	{
	 $ser += $cof[$j]/++$y;
	}
	($x = $eps) if ($x == 0);
	return(-$tmp + log(2.5066282746310005*$ser/$x));
}

sub checkUnique
{
	my %h = @_;
	my @vals2check;
	foreach my $k (keys(%h))
	{
		push (@vals2check,$h{$k});
	}
	my $inlen = @vals2check;
	my %ch = &unique(@vals2check);
	my $outlen = keys(%ch);
	($inlen == $outlen) ? (return(1)) : (return(0));
}

sub printGeneOrPeak
{
	my $infile = shift @_;
	my $otype = shift @_;
	my %inhash = @_;
	my ($outchr,$ind,$outhash);
	my (@k,@v);
	my $outfilename = &createOutputFile($infile,$otype);
	disp("Writing output in $outfilename...");
	my @chrs = keys(%inhash);
	open(OUTPUT,">$outfilename");
	foreach $outchr (sort(@chrs))
	{
		$outhash = $inhash{$outchr};
		@k = keys(%$outhash);
		@v = values(%$outhash);
		for ($ind=0; $ind<@k; $ind++)
		{
			print OUTPUT "$k[$ind]\t";
			print OUTPUT join("\ ",@{$v[$ind]});
			print OUTPUT "\n";
		}
	}
	close(OUTPUT);
}

sub printMatrix
{
	my $inhash = $_[0];
	my $type = $_[1];
	my ($row,$column,$colhash);
	#my (@tmp,@v);
	my $outfilename = &createOutputFile($peakfile[0],$type);
	disp("Writing output in $outfilename...");
	my @rows = keys(%$inhash);
	open(OUTPUT,">$outfilename");
	my $headhash = $inhash->{$rows[0]};
	my @headers = keys(%$headhash);
	print OUTPUT "GeneID\t",join("\t",@headers),"\n";
	if ($type =~ m/number/i)
	{
		foreach $row (@rows)
		{
			$colhash = $inhash->{$row};
			my @tmp = keys(%$colhash);
			my @v;
			foreach $column (@tmp)
			{
				(defined($colhash->{$column})) ?
				(push(@v,scalar(@{$colhash->{$column}}))) :
				(push(@v,0));
			}
			print OUTPUT "$row\t",join("\t",@v),"\n";
		}
	}
	elsif ($type =~ m/presence/i)
	{
		foreach $row (@rows)
		{
			$colhash = $inhash->{$row};
			my @tmp = keys(%$colhash);
			my @v;
			foreach $column (@tmp)
			{
				(defined($colhash->{$column})) ? (push(@v,"+")) : (push(@v,"-"));
			}
			print OUTPUT "$row\t",join("\t",@v),"\n";
		}
	}
	elsif ($type =~ m/peaks/i)
	{
		foreach $row (@rows)
		{
			$colhash = $inhash->{$row};
			my @tmp = keys(%$colhash);
			my @v;
			foreach $column (@tmp)
			{
				(defined($colhash->{$column})) ?
				(push(@v,join("; ",@{$colhash->{$column}}))) :
				(push(@v,"NP"));
			}
			print OUTPUT "$row\t",join("\t",@v),"\n";
		}
	}
	close(OUTPUT);
}

# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@peakfile,
    		   "region|r=s" => \$sigfile,
    		   "background|b=s" => \$backfile,
    		   "span|n=i{,}" => \@span,
    		   "idstrand|c=i{,}" => \@sbcols,
    		   "idmode|m=i{,}" => \@pcols,
    		   "test|t=s" => \$test,
    		   "pvalue|p=f" => \$pval,
    		   "outformat|o=s{,}" => \@out,
    		   "header|d" => \$header,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&programUsage;
    	exit;
    }
    $stop .= "--- Please specify peak file(s) ---\n" if (!@peakfile);
    $stop .= "--- Please specify significant region file ---\n" if (!$sigfile);
    $stop .= "--- Please specify background region file ---\n" if (!$backfile && $test ne "none");
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    # Check statistical test
	if ($test ne "hypgeom" && $test ne "chi2" && $test ne "none" && $test ne "auto")
	{
		disp("--test should be \"hypgeom\", \"chi2\", \"none\" or \"auto\"... Using default (hypgeom)");
		$test = "hypgeom";
	}
    # Check if span given
    if (!@span)
    {
    	disp("Search range from TSS not given... Using defaults (-100kbp,100kbp)");
    	@span = (-1e+5,1e+5);
	}
	# Check if id and strand columns given for sig/back files
    if (!@sbcols)
    {
    	disp("Unique ID and strand columns for region and background files not given... Using defaults (4,5)");
    	@sbcols = (3,4);
	}
	else # Proper perl indexing
	{
		$sbcols[0] = $sbcols[0] - 1;
		$sbcols[1] = $sbcols[1] - 1;
	}
	# Check if id and mode columns given for peak files
    if (!@pcols)
    {
    	disp("Unique ID and mode columns for peak files not given... Using defaults (4,5)");
    	@pcols = (3,4);
	}
	else # Proper perl indexing
	{
		$pcols[0] = $pcols[0] - 1;
		$pcols[1] = $pcols[1] - 1;
	}
	# Check proper output format
	if (@out)
	{
		foreach my $c (@out)
		{
			if ($c ne "stats" && $c ne "gff-peak" && $c ne "gff-gene" && $c ne "peak" &&  
				$c ne "gene" && $c ne "all-peak" && $c ne "all-gene" && $c ne "pretty-peak" && 
				$c ne "pretty-gene" && $c ne "gff-peak-db" && $c ne "gff-gene-db" && $c ne "peakdata"
				&& $c ne "matrix-number" && $c ne "matrix-presence" && $c ne "matrix-peaks")
			{
				my $msg = "WARNING! --outformat options should be \"gff-peak\", \"gff-gene\", \"peak\",\n".
						  "\"gene\", \"all-peak\", \"all-gene\", \"pretty-peak\", \"pretty-gene\", \"gff-peak-db\",\n".
						  "\"gff-gene-db\", \"peakdata\", \"stats\", \"matrix-number\", \"matrix-presence\" or ".
						  "\"matrix-peaks\" \nUsing default (\"gff-peak\")...";
				disp($msg);
				@out = ("gff-peak");
			}
		}
	}
	else
	{
		disp("Output format not given... Using default (gff-peak)");
		@out = ("gff-peak");
	}
}

sub chisquarecont
{
	# Expected values
	my ($a,$b,$c,$d) = @_;
	my $tot = $a + $b + $c + $d;
	my $A = (($a + $b)*($a + $c))/$tot;
	my $B = (($a + $b)*($b + $d))/$tot;
	my $C = (($c + $d)*($a + $c))/$tot;
	my $D = (($c + $d)*($b + $d))/$tot;
	
	# Chi square statistic
	my $x2 = (($a - $A)**2)/$A + (($b - $B)**2)/$B + (($c - $C)**2)/$C +
	         (($d - $D)**2)/$D;

	return(1-chdtr(1,$x2));
}

sub createOutputFile
{
	my $in = shift @_;
	my $type = shift @_;
	my $ext;
	my ($base,$dir) = fileparse($in,'\.[^.]*');
	if ($type =~/gff/)
	{
		($type =~/db/) ? ($ext = ".txt") : ($ext = ".gff");
	}
	else { $ext = ".txt" }
	if ($type =~ /matrix/)
	{
		return($dir."gene-peak-$type.txt");
	}
	else
	{
		return($dir.$base."_ASSIGNED_".$type."$ext");
	}
}

sub now
{
	my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
	$year += 1900;
	$month++;
	$month = "0".$month if (length($month)==1);
	$day = "0".$day if (length($day)==1);
	$hour = "0".$hour if (length($hour)==1);
	$min = "0".$min if (length($min)==1);
	$sec = "0".$sec if (length($sec)==1);
	return($day."/".$month."/".$year." ".$hour.":".$min.":".$sec);
}

sub disp
{
	print "\n@_" if (!$silent);
}

sub programUsage 
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";
	
$scriptname
A perl program to assign ChIP-Seq peaks to a set of regulated genes based on
the gene-peak distances between the genes of the regulated set and the gene-
peak distances between the genes of a background set (e.g. the whole genome).
In order to find peaks associatevely enriched in the set of the regulated
genes as compared to the background set, the program uses the hypergeometric
test on the population of peaks at a p-value threshold. There can be multiple
output format, all of them containing the gene-peak associations in different
formats. The hypergeometric method is NOT verified and there are other methods
out there that may perform better. Use at your own risk. The tools works very
nicely to calculate peak-gene distances with a very nice output.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input peakfile(s) --region regfile --background backfile [OPTIONS]

--- Required ---
  --input|i  file(s)	Peak BED file(s) containing a column with a
  			UNIQUE peak ID and a column with the peak mode (the
  			point with the highest tag pile-up) or a location that
  			the user thinks as the point of the peak from which 
  			the distance to the genes will be calculated.
  --region|r  file	A BED file with the set of regulated genes,
  			containing a column with a UNIQUE gene (or whatever
  			region) ID and a column with the gene strand.
--- Optional ---
  --background|b		A BED file with the set of background 
  			genes, containing a column with a UNIQUE gene (or 
  			whatever region) ID and a column with the gene strand.
  			Required when running a statistical test.
  --span|d		Use this option to set the genomic span (distance
  			upstream and downstream from TSS) into which the program
  			will look for any peaks. It should be two values (e.g.
  			--span -50000 50000) and defaults to (-100000,100000).
  --idstrand|t		The columns in BOTH the gene files where their
			unique IDs and strands are. You should provide two values
			(e.g. --idstrand 4 5) where the first denotes the unique
			ID column and the second the strand column. It defaults
			to (4,5).
  --idmode|m		The columns in the peak files where their unique
			IDs and modes are. You should provide two values (e.g
			--idmode 4 5) where the first denotes the unique ID
			column and the second the mode column. It defaults to (4,5).
  --test|t		What over-representation statistical test to perform.
			Can be one of hypgeom for hypergeometric test, chi2 for
			chi-square test, auto for automatic selection and none for
			no testing. Defaults to hypgeom.
  --pvalue|p		The hypergeometric test p-value threshold. It
			defaults to 0.05.
  --outformat|o		Use this option to determine which output format
			filetype(s) you wish to retrieve.	Possible choices are:
   				"stats" for retrieving the significantly associated 
				peaks with their p-values, Bonferroni corrected p-values
				and enrichment ratios.
   				"gff-peak" for retrieving a peak-based gff file which
   				contains additional columns with peak ids, distances
   				and enrichment ratios. The score column is the p-value.
				"gff-gene" for similar to "gff-peak" but gene-based.
				"gff-peak-db" for same as "gff-peak" but with a header,
				suitable for incorporating to a database.
				"gff-gene-db" for same as "gff-gene" but with a header,
				suitable for incorporating to a database.
				"peak" for a simple file which contains the significantly
				associated peak IDs in the first column and a list of 
				associated genes in the second column.
				"gene" for similar to "peak" but gene based.
				"all-peak" for a simple file which contains ALL (based on
				distance) associated peak IDs in the first column and a 
				list of associated genes in the second column.
				"all-gene" for similar to "all-peak" but gene based.
				"pretty-peak" for retrieving a more human-readable format
				quite self-explicating (please see output).
				"pretty-gene" similar to "pretty-peak" but gene-based
				(please see output).
				"peakdata" for retrieving only the assigned peaks from
				the original peak file.
				"matrix-number" to retrieve a spreadsheet-like file where 
				rows correspond to genes (or the --region file) and columns
				correspond to peak files. The cell (i,j) contains the number
				of peaks in peak file j assigned to gene i.
				"matrix-presence" to retrieve a spreadsheet-like file where 
				rows correspond to genes (or the --region file) and columns
				correspond to peak files. The cell (i,j) contains "+" if peak
				in peak file j assigned to gene i, "-" otherwise.
				"matrix-peaks" to retrieve a spreadsheet-like file where 
				rows correspond to genes (or the --region file) and columns
				correspond to peak files. The cell (i,j) contains the peaks
				in peak file j assigned to gene i, "NP" otherwise.
			Example: --outformat stats gff-peak pretty-gene matrix
  --header|e		Use this option if you have a header line in
			your input files.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.
	
The main output of the program is up to nine files with information on gene-peak
association.

END
	print $usagetext;
	exit;
}


## Old chunks
#sub printMatrix
#{
#	my %inhash = @_;
#	my ($row,$colhash);
#	my @v;
#	my $outfilename = &createOutputFile($peakfile[0],"matrix");
#	disp("Writing output in $outfilename...");
#	my @rows = keys(%inhash);
#	print Dumper(\@rows);
#	open(OUTPUT,">$outfilename");
#	my $headhash = $inhash{$rows[0]};
#	my @headers = keys(%$headhash);
#	print OUTPUT "GeneID\t",join("\t",@headers),"\n";
#	foreach $row (@rows)
#	{
#		$colhash = $inhash{$row};
#		@v = values(%$colhash);
#		print OUTPUT "$row\t",join("\t",@v),"\n";
#	}
#	#close(OUTPUT);
#}
#
############################### PROTOTYPE SECTION ############################## 
#sub intersect (\@\@)
#{
#	my ($arr1,$arr2) = @_;
#	my %h1 = map{$_ => 1} @{$arr1};
#	my %h2 = map{$_ => 1} @{$arr2};
#	my @inter = grep($h1{$_},@{$arr2});
#	my $len1 = @{$arr1};
#	my $len2 = @{$arr2};
#	my $len3 = @inter;
#	print "\n$len1\t$len2\t$len3\n";
#	
#	return(@inter);
#}
############################## END PROTOTYPE SECTION ############################
