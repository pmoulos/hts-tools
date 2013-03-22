#!/usr/bin/perl

# Description later...
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 15 - 01 - 2009 (dd - mm - yyyy)
# Last Update : 12 - 06 - 2009 (dd - mm - yyyy)
# Version     : 2.2
#
# Version 1.1 introduces overlaping based on both-way extension from the peak mode for a
# number of pbs and the --both switch
# Version 2.0 sorts files with headers (Unix/Linux only) and has extra options for returning
# distance distributions (longest, shortest and from peak modes) together with pairs of
# intersected locations and also allows gaps between non-intersected regions so as to find
# distances of let's say closey binding factors which do not overlap
# Version 2.1 fixes certain bugs of overlap accuracy and displays summary output
# Version 2.2 does not require the peak mode column to be given in the --extend option. It
# has its own options instead (--mode). Also extension can be done in different lengths
# upstream and downstream

#################### SPECIAL SECTION ####################
sub outsort 
{
  my @one = split(/\t/,$a);
  my @two = split(/\t/,$b);
  return 1 if ($one[0] > $two[0]);
  return 0 if ($one[0] == $two[0]);
  return -1 if ($one[0] < $two[0]);
}
################## END SPECIAL SECTION ##################

use strict;
use Getopt::Long;
use File::Basename;
use POSIX qw(floor
 			 ceil);
 			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "intersectbed.pl";
our $fileA;		   	# Input bedfile A 
our $fileB;			# Input bedfile B
our $sort = 0;		# Sort input files? (Necessary for binary search)
our @percent = ();	# Overlap percentage (not the default)
our $any = 1;		# Any overlap (default)
our @extend;	 	# Extend peak modes and use this instead of whole width
our $mode;			# The column where peak mode is in case of extend or distance calculation
our $autoxtend = 0; # Extend peaks based on half of median peak length, peak length column
our $both = 0;		# When overlaping with percent, overlap based on both tags or the first?
our $exact = 0;		# Percentage intersections exactly as the UCSC Table Browser?
our $npass = 3;		# Number of passes for dynamic binary search
our $agap = 0;      # Acceptable gap to say that two peaks from two conditions are close
our @out;		 	# Output files
our $multi;	# Called from outside to multiply intersect?
our $header = 0;	# Files contain header?
our $waitbar = 0;  	# Do not use waitbar
our $silent = 0;  	# Display verbose messages
our $help = 0;	   	# Help?

# Check inputs
&checkInputs;

# Start work... bed files have standard format (1st columns, chr, start, end) so sorting
# is easy even if containing additional columns...
my $windows = 0; # Windows or Linux?
my $ofileA = $fileA;
my $ofileB = $fileB; # Keep original filenames (verbose purposes)

if ($sort)
{
	# Simulate try and catch
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		if ($header)
		{
			disp("Sorting file $fileA...");
            my $widA = "$$";
            my $sortcmdA = "awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $fileA > /tmp/tempA.in$widA";
            `$sortcmdA`;
            #`awk 'NR==1; NR > 1 {print $0 | " sort -k1,1 -k2g,2"}' $fileA > /tmp/tempA.in$$ `;
            $fileA = "/tmp/tempA.in$widA";
            disp("Sorting file $fileB...");
            my $widB = "$$";
            my $sortcmdB = "awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $fileB > /tmp/tempB.in$widB";
            `$sortcmdB`;
            #`awk 'NR==1; NR > 1 {print $0 | " sort -k1,1 -k2g,2"}' $fileB > /tmp/tempB.in$$ `;
            $fileB = "/tmp/tempB.in$widB";

		}
		else
		{	
			disp("Sorting file $fileA...");
			`sort -k1,1 -k2g,2 $fileA > /tmp/tempA.in$$ `;
			$fileA = "/tmp/tempA.in$$";
			disp("Sorting file $fileB...");
			`sort -k1,1 -k2g,2 $fileB > /tmp/tempB.in$$ `;
			$fileB = "/tmp/tempB.in$$";
		}
	}
	else # We are in Windows... package required
	{
		eval
		{
			require File::Sort;
		};
		if ($@)
		{
			my $killer = "Module File::Sort is required to continue with file sorting. If you\n".
						 "have ActiveState Perl installed, use the Package Manager to get the\n".
						 "module. If you don't know how to install the package, sort the files\n".
						 "using another tool like Excel or contact your system administrator.";
			die "\n$killer\n\n";
		}
		else # Else sort the bed files according to chromosome and start position
		{
			if ($header) # Cannot find a solution for direct sorting in Windows... sorry :-(
			{
				my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
						   "messing up data. Please sort files outside $scriptname first (e.g. using\n".
						   "Excel or something similar.";
				die "\n$dmsg\n\n";
			}
			eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
			disp("Sorting region file $fileA...");
			sort_file(
			{
				I => $fileA,
				o => "tempA.tmp",
				k => ['1,1','2n,2'],
				t => "\t"
			});
			$fileA = "tempA.tmp";
			disp("Sorting region file $fileB...");
			sort_file(
			{
				I => $fileB,
				o => "tempB.tmp",
				k => ['1,1','2n,2'],
				t => "\t"
			});
			$fileB = "tempB.tmp";
			$windows = 1;
		}
	}
}

# Bavard a little... I like information...
disp("Type of overlap: any overlap, percentage overlap ignored...") if ($any && @percent);
disp("Type of overlap: any overlap") if ($any && !@percent);
disp("Type of overlap: $percent[0]% overlap") if ($percent[0] && !$percent[1] && !$any);
disp("Extending region modes upstream $extend[0] and downstream $extend[1] bps") if (@extend && !$autoxtend);
disp("Region modes extension on each side will be auto-calculated...") if ($autoxtend);
if ($percent[1])
{
	disp("Multiple overlap percentages... Running in batch mode to determine overlapping distributions...");
	disp("No overlapping output files will be produced..."); 
}

# Save some memory...
my $ovA = my $ovB = my $oA = my $oB = my $op = my $np = 0;
foreach my $opt (@out)
{
	$ovA = 1 if ($opt eq "overlapA");
	$ovB = 1 if ($opt eq "overlapB");
	$oA = 1 if ($opt eq "onlyA");
	$oB = 1 if ($opt eq "onlyB");
	$op = 1 if ($opt eq "overpairs");
	$np = 1 if ($opt eq "nonpairs");
}

if (($op || $np) && !$mode)
{
	disp("overpairs and nonpairs output formats can only be given with --mode option... Ignoring.");
	$op = 0;
	$np = 0;
}
disp("Retrieving distances between non-overlapping regions if distance <= $agap bps") if (($op || $np) && $agap);
disp(" ");

# Get number of lines for fileA (mostly for waitbar, I know, inefficient but I like waitbars ;-)
my $linesA = &countLines($fileA) if ($waitbar);

# Suck in fileB
my ($chr,@rest);
my (%hashB,%onlyB);
open(INB,$fileB) or die "\nThe file $fileB does not exist!\n";
disp("Reading file $ofileB...");
my $headerB = <INB> if ($header);
while (my $line = <INB>)
{
	next if ($line =~/^chrM/);
	next if ($line =~/rand/);
	$line =~ s/\r|\n$//g;
	($chr,@rest) = split(/\t/,$line);
	push(@{$hashB{$chr}},join("\t",@rest));
	#push(@{$onlyB{$chr}},join("\t",@rest));
	$onlyB{$chr}{join("\t",@rest)}++ if ($oB);
}
close(INB);

# Doesn't work!!! Creates rather reference... Seen also in some forums...
# my %onlyB = %hashB;

##########################################################################################

# Do we run in batch mode in order to determine distributions?
if (!$any && $percent[1])
{	
	my $i;
	my @distributions;
	my ($cchr,@crest);
	my ($bsr,$ci,$bsf,$n,@currvals);
	my $linB = &countLines($fileB);
	$linB-- if ($header);

	if (@extend || $autoxtend)
	{	
		if ($autoxtend) # We have to suck in the whole peak file once...
		{
			open (INA,$fileA) or die "\nThe file $fileA does not exist!\n";;
			disp("Reading file $ofileA and processing overlaps...");
			my $headerA = <INA> if ($header);
			my (@lines,@medmodes);
			&waitbarInit if ($waitbar);
			while (my $line = <INA>)
			{
				&waitbarUpdate($.,$linesA) if ($waitbar);
				next if ($line =~/^chrM/);
				next if ($line =~/rand/);
				$line =~ s/\r|\n$//g;
				push(@lines,$line);
				@crest = split(/\t/,$line);
				push(@medmodes,$crest[2] - $crest[1]); # BED format
			}
			close(INA);
			my $med = &median((@medmodes,&getLengths(%hashB)));
			@extend = (int($med/2),int($med/2));
			disp("Median region length is $med bps. Extending each region mode $extend[1] bps on each side...\n");
			for ($i=0;$i<@percent;$i++)
			{	
				disp("Overlap percentage: $percent[$i]");
				my (%overlapA,%overlapB,%onlyA,%conlyB);
				my @overpairs;
				# Make a hard copy of onlyB hash
				foreach my $ock (keys(%onlyB))
				{
					foreach my $ick (keys(%{$onlyB{$ock}})) 
					{
						$conlyB{$ock}{$ick} = $onlyB{$ock}{$ick};
					}
				}
				my $countopic = 0;
				&waitbarInit if ($waitbar);
				foreach my $l (@lines)
				{
					$countopic++;
					&waitbarUpdate($countopic,$linesA-1) if ($waitbar);
					($cchr,@crest) = split(/\t/,$l);
					@currvals = @{$hashB{$cchr}} if ($hashB{$cchr}); 
					$n = 0;
					$bsf = 0;
					while ($n < $npass) 
					{
						($ci,$bsr) = &binSearchPercentCenter($crest[$mode],$mode,@extend,$percent[$i]/100,@currvals);
						if ($bsr) # Found in overlap, put into overlap hash of both files
						{
							$overlapA{$cchr}{join("\t",@crest)}++;
							$overlapB{$cchr}{$bsr}++ if ($ovB);
							delete $conlyB{$cchr}{$bsr} if ($oB);
							if ($op)
							{
								my $cl = join("\t",@crest);
								my @ds = &distsCenter($cl,$bsr,$mode,@extend);
								push(@overpairs,$ds[1]);
							}
							$bsf++;
						}
						if ($ci)
						{
							splice(@currvals,$ci,1); # Remove it from areas else it will be found again
							$n++;
						} else { last; }
					}

					if (!$bsf) 
					{
						$onlyA{$cchr}{join("\t",@crest)}++ if ($oA);
					}
				}
				
				push(@{$distributions[$i]},&countHoH(%overlapA)) if ($ovA);
				push(@{$distributions[$i]},&countHoH(%overlapB)) if ($ovB);
				push(@{$distributions[$i]},&countHoH(%onlyA)) if ($oA);
				push(@{$distributions[$i]},&countHoH(%conlyB)) if ($oB);
				push(@{$distributions[$i]},&mean(@overpairs)) if ($op);
				push(@{$distributions[$i]},&median(@overpairs)) if ($op);
			}
		}
		else
		{
			for ($i=0;$i<@percent;$i++)
			{
				
				disp("Overlap percentage: $percent[$i]");
				my (%overlapA,%overlapB,%onlyA,%conlyB);
				my @overpairs;
				# Make a hard copy of onlyB hash
				foreach my $ock (keys(%onlyB))
				{
					foreach my $ick (keys(%{$onlyB{$ock}})) 
					{
						$conlyB{$ock}{$ick} = $onlyB{$ock}{$ick};
					}
				}		
				open (INA,$fileA) or die "\nThe file $fileA does not exist!\n";;
				disp("Reading file $ofileA and processing overlaps...");
				&waitbarInit if ($waitbar);
				my $headerA = <INA> if ($header);
				while (my $line = <INA>)
				{
					&waitbarUpdate($.,$linesA) if ($waitbar);	
					next if ($line =~/^chrM/);
					next if ($line =~/rand/);
					$line =~ s/\r|\n$//g;
					($cchr,@crest) = split(/\t/,$line);
					@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
					$n = 0;
					$bsf = 0;
					while ($n < $npass) 
					{
						($ci,$bsr) = &binSearchPercentCenter($crest[$mode],$mode,@extend,$percent[$i]/100,@currvals);
						if ($bsr) # Found in overlap, put into overlap hash of both files
						{
							$overlapA{$cchr}{join("\t",@crest)}++;
							$overlapB{$cchr}{$bsr}++ if ($ovB);
							delete $conlyB{$cchr}{$bsr} if ($oB);
							if ($op)
							{
								my $cl = join("\t",@crest);
								my @ds = &distsCenter($cl,$bsr,$mode,@extend);
								push(@overpairs,$ds[1]);
							}
							$bsf++;
						}
						if ($ci)
						{
							splice(@currvals,$ci,1);
							$n++;
						} else { last; }
					}
					
					if (!$bsf) 
					{
						$onlyA{$cchr}{join("\t",@crest)}++ if ($oA);
					}
				}
				close(INA);
				
				push(@{$distributions[$i]},&countHoH(%overlapA)) if ($ovA);
				push(@{$distributions[$i]},&countHoH(%overlapB)) if ($ovB);
				push(@{$distributions[$i]},&countHoH(%onlyA)) if ($oA);
				push(@{$distributions[$i]},&countHoH(%conlyB)) if ($oB);
				push(@{$distributions[$i]},&mean(@overpairs)) if ($op);
				push(@{$distributions[$i]},&median(@overpairs)) if ($op);
			}
		}
	}
	else
	{
		for ($i=0;$i<@percent;$i++)
		{
			disp("Overlap percentage: $percent[$i]");
			my (%overlapA,%overlapB,%onlyA,%conlyB);
			my @overpairs;
			# Make a hard copy of onlyB hash
			foreach my $ock (keys(%onlyB))
			{
				foreach my $ick (keys(%{$onlyB{$ock}})) 
				{
					$conlyB{$ock}{$ick} = $onlyB{$ock}{$ick};
				}
			}
			
			my $cntovA = my $cntovB = my $cntonA = 0;
			my $cntonB = $linB;
			
			open (INA,$fileA) or die "\nThe file $fileA does not exist!\n";;
			disp("Reading file $ofileA and processing overlaps...");
			&waitbarInit if ($waitbar);
			my $headerA = <INA> if ($header);
			while (my $line = <INA>)
			{
				&waitbarUpdate($.,$linesA) if ($waitbar);	
				next if ($line =~/^chrM/);
				next if ($line =~/rand/);
				$line =~ s/\r|\n$//g;
				($cchr,@crest) = split(/\t/,$line);
				@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
				$n = 0;
				$bsf = 0;
				while ($n < $npass) 
				{
					if ($exact)
					{
						($both) ? (($ci,$bsr) = &binSearchPercentExactBoth($crest[0],$crest[1],$percent[$i]/100,@currvals)) :
						(($ci,$bsr) = &binSearchPercentExact($crest[0],$crest[1],$percent[$i]/100,@currvals));
					}
					else
					{
						($both) ? (($ci,$bsr) = &binSearchPercentBoth($crest[0],$crest[1],$percent[$i]/100,@currvals)) :
						(($ci,$bsr) = &binSearchPercent($crest[0],$crest[1],$percent[$i]/100,@currvals));
					}
					
					if ($bsr) # Found in overlap, put into overlap hash of both files
					{
						$overlapA{$cchr}{join("\t",@crest)}++;
						$overlapB{$cchr}{$bsr}++ if ($ovB);
						delete $conlyB{$cchr}{$bsr} if ($oB);
						if ($op)
						{
							my $cl = join("\t",@crest);
							my @ds = &distsEvery($cl,$bsr,$mode);
							push(@overpairs,$ds[1]);
						}
						$bsf++;
					}
					if ($ci)
					{
						splice(@currvals,$ci,1);
						$n++;
					} else { last; }
				}
				if (!$bsf) 
				{
					$onlyA{$cchr}{join("\t",@crest)}++ if ($oA);
				}
			}
			close(INA);
			
			push(@{$distributions[$i]},&countHoH(%overlapA)) if ($ovA);
			push(@{$distributions[$i]},&countHoH(%overlapB)) if ($ovB);
			push(@{$distributions[$i]},&countHoH(%onlyA)) if ($oA);
			push(@{$distributions[$i]},&countHoH(%conlyB)) if ($oB);
			push(@{$distributions[$i]},&mean(@overpairs)) if ($op);
			push(@{$distributions[$i]},&median(@overpairs)) if ($op);
		}
	}
	
	my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
	$year += 1900;
	$month++;
	my $fdn = "distrib_".$day.$month.$year.$hour.$min.$sec;
	disp("Writing output in $fdn.txt");
	open(OUTDISTRIB,">$fdn.txt");
	my ($cdo,@disthead);
	push(@disthead,"Overlap $fileA") if ($ovA);
	push(@disthead,"Overlap $fileB") if ($ovB);
	push(@disthead,"Only $fileA") if ($oA);
	push(@disthead,"Only $fileB") if ($oB);
	push(@disthead,"Mean distance\tMedian distance") if ($op);
	print OUTDISTRIB join("\t",@disthead),"\n";
	foreach $cdo (@distributions)
	{
		print OUTDISTRIB join("\t",@{$cdo}),"\n";
	}
	close(OUTDISTRIB);
	
	# Remove garbage
	(!$windows) ? (`rm /tmp/temp*.in$$ `) : (`del /f /q temp*.tmp `) if ($sort);
	disp("Finished!\n\n");
	exit;
}

##########################################################################################

# Do job with lines of fileA
my ($cchr,@crest);
my ($bsr,$ci,$bsf,$n,@currvals);
my (%overlapA,%overlapB,%onlyA);
my (@overpairs,@nonpairs);
	
open (INA,$fileA) or die "\nThe file $fileA does not exist!\n";;
disp("Reading file $ofileA and processing overlaps...");
my $headerA = <INA> if ($header);
&waitbarInit if ($waitbar);
if (@extend || $autoxtend)
{
	if ($autoxtend) # We have to suck in the whole peak file once...
	{
		my (@lines,@medmodes);
		while (my $line = <INA>)
		{
			&waitbarUpdate($.,$linesA) if ($waitbar);
			next if ($line =~/^chrM/);
			next if ($line =~/rand/);
			$line =~ s/\r|\n$//g;
			push(@lines,$line);
			@crest = split(/\t/,$line);
			push(@medmodes,$crest[2] - $crest[1]); # BED format
		}
		my $med = &median((@medmodes,&getLengths(%hashB)));
		@extend = (int($med/2),int($med/2));
		disp("Median region length is $med bps. Extending each region mode $extend[1] bps on each side...\n");
		my $countopic = 0;
		&waitbarInit if ($waitbar);
		foreach my $l (@lines)
		{
			$countopic++;
			&waitbarUpdate($countopic,$linesA-1) if ($waitbar);
			($cchr,@crest) = split(/\t/,$l);
			
			# Perform binary search according to user choices...
			# A specific chromosome might not exist in one of the two files...
			@currvals = @{$hashB{$cchr}} if ($hashB{$cchr}); 
			$n = 0;
			$bsf = 0;
			while ($n < $npass) 
			{
				($any) ? (($ci,$bsr) = &binSearchAnyCenter($crest[$mode],$mode,@extend,@currvals)) :
				(($ci,$bsr) = &binSearchPercentCenter($crest[$mode],$mode,@extend,$percent[0]/100,@currvals));
				if ($bsr) # Found in overlap, put into overlap hash of both files
				{
					$overlapA{$cchr}{join("\t",@crest)}++;
					#push(@{$overlapA{$cchr}},join("\t",@crest));
					$overlapB{$cchr}{$bsr}++ if ($ovB);
					# Remove from fileB hash so as what remains will be the only B.
					#splice(@{$onlyB{$cchr}},$ci,1);
					delete $onlyB{$cchr}{$bsr} if ($oB);
					if ($op)
					{
						my $cl = join("\t",@crest);
						my @ds = &distsCenter($cl,$bsr,$mode,@extend);
						push(@overpairs,"$cchr\t$cl\t$cchr\t$bsr\t$ds[0]\t$ds[1]\t$ds[2]");
					}
					$bsf++;
				}
				if ($ci)
				{
					splice(@currvals,$ci,1); # Remove it from areas else it will be found again
					$n++;
				} else { last; }
			}
			
			# Not found in any of the binary search passes
			if (!$bsf) 
			{
				$onlyA{$cchr}{join("\t",@crest)}++ if ($oA);
				#push(@{$onlyA{$cchr}},join("\t",@crest)) if ($oA);
				if ($np)
				{
					my $cl = join("\t",@crest);
					my @ds = &distsCenter($cl,$bsr,$mode,@extend);
					if ($agap)
					{
						push(@nonpairs,"$cchr\t$cl\t$cchr\t$bsr\t$ds[0]\t$ds[1]\t$ds[2]") if (abs($ds[1]) <= $agap);
					}
				}
			}
		}
	}
	else
	{
		while (my $line = <INA>)
		{
			&waitbarUpdate($.,$linesA) if ($waitbar);	
			next if ($line =~/^chrM/);
			next if ($line =~/rand/);
			$line =~ s/\r|\n$//g;
			($cchr,@crest) = split(/\t/,$line);

			# Perform binary search according to user choices...
			# A specific chromosome might not exist in one of the two files...
			@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
			$n = 0;
			$bsf = 0;
			while ($n < $npass) 
			{
				($any) ? (($ci,$bsr) = &binSearchAnyCenter($crest[$mode],$mode,@extend,@currvals)) :
				(($ci,$bsr) = &binSearchPercentCenter($crest[$mode],$mode,@extend,$percent[0]/100,@currvals));
				if ($bsr) # Found in overlap, put into overlap hash of both files
				{
					$overlapA{$cchr}{join("\t",@crest)}++;
					#push(@{$overlapA{$cchr}},join("\t",@crest));
					$overlapB{$cchr}{$bsr}++ if ($ovB);
					# Remove from fileB hash so as what remains will be the only B.
					delete $onlyB{$cchr}{$bsr} if ($oB);
					if ($op)
					{
						my $cl = join("\t",@crest);
						my @ds = &distsCenter($cl,$bsr,$mode,@extend);
						push(@overpairs,"$cchr\t$cl\t$cchr\t$bsr\t$ds[0]\t$ds[1]\t$ds[2]");
					}
					$bsf++;
				}
				if ($ci)
				{
					splice(@currvals,$ci,1);
					$n++;
				} else { last; }
			}
			
			# Not found in any of the binary search passes
			if (!$bsf) 
			{
				$onlyA{$cchr}{join("\t",@crest)}++ if ($oA);
				#push(@{$onlyA{$cchr}},join("\t",@crest)) if ($oA);
				if ($np)
				{
					my $cl = join("\t",@crest);
					my @ds = &distsCenter($cl,$bsr,$mode,@extend);
					if ($agap)
					{
						push(@nonpairs,"$cchr\t$cl\t$cchr\t$bsr\t$ds[0]\t$ds[1]\t$ds[2]") if (abs($ds[1]) <= $agap);
					}
				}
			}
		}
	}
}
else
{
	while (my $line = <INA>)
	{
		&waitbarUpdate($.,$linesA) if ($waitbar);	
		next if ($line =~/^chrM/);
		next if ($line =~/rand/);
		$line =~ s/\r|\n$//g;
		($cchr,@crest) = split(/\t/,$line);

		# Perform binary search according to user choices...
		# A specific chromosome might not exist in one of the two files...
		@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
		$n = 0;
		$bsf = 0;
		while ($n < $npass) 
		{
			if ($any)
			{
				($ci,$bsr) = &binSearchAny($crest[0],$crest[1],@currvals);
			}
			else
			{
				if ($exact)
				{
					($both) ? (($ci,$bsr) = &binSearchPercentExactBoth($crest[0],$crest[1],$percent[0]/100,@currvals)) :
					(($ci,$bsr) = &binSearchPercentExact($crest[0],$crest[1],$percent[0]/100,@currvals));
				}
				else
				{
					($both) ? (($ci,$bsr) = &binSearchPercentBoth($crest[0],$crest[1],$percent[0]/100,@currvals)) :
					(($ci,$bsr) = &binSearchPercent($crest[0],$crest[1],$percent[0]/100,@currvals));
				}
			}
			if ($bsr) # Found in overlap, put into overlap hash of both files
			{
				$overlapA{$cchr}{join("\t",@crest)}++;
				$overlapB{$cchr}{$bsr}++ if ($ovB);
				# Remove from fileB hash so as what remains will be the only B.
				delete $onlyB{$cchr}{$bsr} if ($oB);
				if ($op)
				{
					my $cl = join("\t",@crest);
					my @ds = &distsEvery($cl,$bsr,$mode);
					push(@overpairs,"$cchr\t$cl\t$cchr\t$bsr\t$ds[0]\t$ds[1]\t$ds[2]");
				}
				$bsf++;
			}
			if ($ci)
			{
				splice(@currvals,$ci,1);
				$n++;
			} else { last; }
		}
		
		# Not found in any of the binary search passes
		if (!$bsf) 
		{
			$onlyA{$cchr}{join("\t",@crest)}++ if ($oA);
			#push(@{$onlyA{$cchr}},join("\t",@crest)) if ($oA);
			if ($np)
			{
				my $cl = join("\t",@crest);
				my @ds = &distsEvery($cl,$bsr,$mode);
				if ($agap)
				{
					push(@nonpairs,"$cchr\t$cl\t$cchr\t$bsr\t$ds[0]\t$ds[1]\t$ds[2]") if (abs($ds[1]) <= $agap);
				}
			}
		}
	}
}
close(INA);

# Display stats
disp(" ") if (!$waitbar);
if (!$silent)
{
	my ($lA,$lB);
	($waitbar) ? ($lA = $linesA) : ($lA = &countLines($fileA));
	$lB = &countLines($fileB);
	if ($header)
	{
		$lA--;	
		$lB--;
	}
	if ($ovA)
	{
		my $covA = &countHoH(%overlapA);
		disp("$covA out of $lA regions from $ofileA overlap with regions from $ofileB");
	}
	if ($ovB)
	{
		my $covB = &countHoH(%overlapB);
		disp("$covB out of $lB regions from $ofileB overlap with regions from $ofileA");
	}
	if ($oA)
	{
		my $coA = &countHoH(%onlyA);
		disp("$coA out of $lA regions exist only in $ofileA");
	}
	if ($oB)
	{
		my $coB = &countHoH(%onlyB);
		disp("$coB out of $lB regions exist only in $ofileB");
	}
	disp(" ");
} 

# Good... print outputs...
disp("Writing output...");
my ($bA,$bB);
$bA = fileparse($ofileA,'\..*?') if ($ovA || $oA);
$bB = fileparse($ofileB,'\..*?') if ($ovB || $oB);
# A and B overlap elements with elements of file A
&printComplexOutput($ofileA,$ofileB,$headerA,"OVERLAP_FROM_$bA",%overlapA) if ($ovA); 
# A and B overlap elements with elements of file B
&printComplexOutput($ofileA,$ofileB,$headerB,"OVERLAP_FROM_$bB",%overlapB) if ($ovB);
# A only elements
&printComplexOutput($ofileA,$ofileB,$headerA,"ONLY_$bA",%onlyA) if ($oA);
# B only elements
&printComplexOutput($ofileA,$ofileB,$headerB,"ONLY_$bB",%onlyB) if ($oB);
# Distances of overlaping
&printArray($ofileA,$ofileB,"OVERDIST",@overpairs) if ($op);
# Distances of non-overlaping
&printArray($ofileA,$ofileB,"NONDIST",@nonpairs) if ($np);

# Remove garbage
(!$windows) ? (`rm /tmp/temp*.in$$ `) : (`del /f /q temp*.tmp `) if ($sort);
disp("Finished!\n\n") if (!$multi);


sub binSearchAny
{
	my $start = shift @_;
	my $end = shift @_;
	my @areas = @_;
	my ($ind,$currstart,$currend);
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end) ||
			($currstart < $start && $currend < $end && $start < $currend) ||
			($currstart > $start && $currend > $end && $end > $currstart))
		{
			return($ind,$areas[$ind]);
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub binSearchPercent
{
	my $start = shift @_;
	my $end = shift @_;
	my $p = shift @_;
	my @areas = @_;
	my ($ind,$currstart,$currend);
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end))
		{
			return($ind,$areas[$ind]);
		}
		elsif ($currstart < $start && $currend < $end && $start < $currend)
		{
			(($currend - $start) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart > $start && $currend > $end && $end > $currstart)
		{
			(($end - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) :
			(return($ind,0));
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub binSearchAnyCenter
{
	my $mode = shift @_;
	my $mpos = shift @_;
	my $dxval = shift @_;
	my $uxval = shift @_;
	my @areas = @_;
	my ($ind,$start,$end,$currstart,$currend,@arr);
	$start = $mode - $dxval;
	$end = $mode + $uxval;
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		@arr = split(/\t/,$areas[$ind]);
		$currstart = $arr[$mpos] - $dxval;
		$currend = $arr[$mpos] + $uxval;
		if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end) ||
			($currstart < $start && $currend < $end && $start < $currend) ||
			($currstart > $start && $currend > $end && $end > $currstart))
		{
			return($ind,$areas[$ind]);
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub binSearchPercentCenter
{
	my $mode = shift @_;
	my $mpos = shift @_;
	my $dxval = shift @_;
	my $uxval = shift @_;
	my $p = shift @_;
	my @areas = @_;
	my ($ind,$start,$end,$currstart,$currend,@arr);
	$start = $mode - $dxval;
	$end = $mode + $uxval;
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		@arr = split(/\t/,$areas[$ind]);
		$currstart = $arr[$mpos] - $dxval;
		$currend = $arr[$mpos] + $uxval;
		if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end))
		{
			return($ind,$areas[$ind]);
		}
		elsif ($currstart < $start && $currend < $end && $start < $currend)
		{
			(($currend - $start) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart > $start && $currend > $end && $end > $currstart)
		{
			(($end - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) :
			(return($ind,0));
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub binSearchPercentBoth
{
	my $start = shift @_;
	my $end = shift @_;
	my $p = shift @_;
	my @areas = @_;
	my ($ind,$currstart,$currend,$diff);
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end))
		{
			return($ind,$areas[$ind]);
		}
		elsif ($currstart < $start && $currend < $end && $start < $currend)
		{
			$diff = $currend - $start;
			($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			(return($ind,$areas[$ind])) : (return($ind,0));
		}
		elsif ($currstart > $start && $currend > $end && $end > $currstart)
		{
			$diff = $end - $currstart;
			($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			(return($ind,$areas[$ind])) : (return($ind,0));
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub binSearchPercentExact
{
	my $start = shift @_;
	my $end = shift @_;
	my $p = shift @_;
	my @areas = @_;
	#print "\n$start\t$end\t$p";
	my ($ind,$currstart,$currend);
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		if ($currstart >= $start && $currend <= $end) # TagB <= TagA
		{
			(($currend - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart <= $start && $currend >= $end) # TagB >= TagA
		{
			(($end - $start) >= $p*($currend - $currstart)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart < $start && $currend < $end && $start < $currend)
		{
			(($currend - $start) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart > $start && $currend > $end && $end > $currstart)
		{
			(($end - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) :
			(return($ind,0));
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub binSearchPercentExactBoth
{
	my $start = shift @_;
	my $end = shift @_;
	my $p = shift @_;
	my @areas = @_;
	my ($ind,$currstart,$currend,$diff);
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		if ($currstart >= $start && $currend <= $end) # TagB <= TagA
		{
			(($currend - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart <= $start && $currend >= $end) # TagB >= TagA
		{
			(($end - $start) >= $p*($currend - $currstart)) ? (return($ind,$areas[$ind])) : 
			(return($ind,0));
		}
		elsif ($currstart < $start && $currend < $end && $start < $currend)
		{
			$diff = $currend - $start;
			($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			(return($ind,$areas[$ind])) : (return($ind,0));
		}
		elsif ($currstart > $start && $currend > $end && $end > $currstart)
		{
			$diff = $end - $currstart;
			($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			(return($ind,$areas[$ind])) : (return($ind,0));
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub distsEvery
{
	my $A = shift @_;
	my $B = shift @_;
	my $ei = shift @_;
	my ($sa,$ea,$ma,$sb,$eb,$mb,$ds,$da,$dc);
	my (@aA,@aB);
	
	@aA = split(/\t/,$A);
	@aB = split(/\t/,$B);
	($sa,$ea,$ma) = ($aA[0],$aA[1],$aA[$ei]);
	($sb,$eb,$mb) = ($aB[0],$aB[1],$aB[$ei]);
	
	if ($sa < $sb && $ea < $eb && $ea < $sb)
	{
		$ds = $ea - $sb;
		$da = $ma - $mb if ($ei);
		$dc = $sa - $eb;
	}
	elsif ($sa > $sb && $ea > $eb && $sa > $eb)
	{
		$ds = $sa - $eb;
		$da = $ma - $mb if ($ei);
		$dc = $ea - $sb;
	}
	elsif ($sa < $sb && $ea < $eb && $sb < $ea)
	{
		$ds = $sb - $ea;
		$da = $ma - $mb if ($ei);
		$dc = $sa - $eb;
	}
	elsif ($sa > $sb && $ea > $eb && $eb > $sa)
	{
		$ds = $sa - $eb;
		$da = $ma - $mb if ($ei);
		$dc = $ea - $sb;
	}
	elsif (($sa <= $sb && $ea >= $eb) || ($sa >= $sb && $ea <= $eb))
	{
		$ds = 0;
		$da = $ma - $mb if ($ei);
		$dc = 0;
	}
	
	$ds = 0 if (!$ds);
	$da = 0 if (!$da);
	$dc = 0 if (!$dc);
	
	return (($ds,$da,$dc));
}

sub distsCenter
{
	my $A = shift @_;
	my $B = shift @_;
	my $ei = shift @_;
	my $exu = shift @_;
	my $exd = shift @_;
	my ($sa,$ea,$ma,$sb,$eb,$mb,$ds,$da,$dc);
	my (@aA,@aB);
	
	@aA = split(/\t/,$A);
	@aB = split(/\t/,$B);
	($sa,$ea,$ma) = ($aA[$ei] - $exu,$aA[$ei] + $exd,$aA[$ei]);
	($sb,$eb,$mb) = ($aB[$ei] - $exu,$aB[$ei] + $exd,$aB[$ei]);
	
	if ($sa < $sb && $ea < $eb && $ea < $sb)
	{
		$ds = $ea - $sb;
		$da = $ma - $mb;
		$dc = $sa - $eb;
	}
	elsif ($sa > $sb && $ea > $eb && $sa > $eb)
	{
		$ds = $sa - $eb;
		$da = $ma - $mb;
		$dc = $ea - $sb;
	}
	elsif ($sa < $sb && $ea < $eb && $sb < $ea)
	{
		$ds = $sb - $ea;
		$da = $ma - $mb;
		$dc = $sa - $eb;
	}
	elsif ($sa > $sb && $ea > $eb && $eb > $sa)
	{
		$ds = $sa - $eb;
		$da = $ma - $mb;
		$dc = $ea - $sb;
	}
	
	$ds = 0 if (!$ds);
	$da = 0 if (!$da);
	$dc = 0 if (!$dc);
	
	return (($ds,$da,$dc));
}

sub mean 
{
    my ($curr,$sum,$n);
	foreach $curr (@_)
    {
    	$sum+= $curr;
    	$n++;
	}
    return ($sum/$n);
}
			
sub median 
{
    my @pole = sort { $a <=> $b } @_;
    my $ret;

    if((@pole % 2) == 1)
    {
    	$ret = $pole[((@pole+1)/2) - 1];
    } else 
    {
        $ret = ($pole[(@pole/2) - 1] + $pole[@pole/2])/2;
    }

    return $ret;
}

# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("inputA|a=s" => \$fileA,
    		   "inputB|b=s" => \$fileB,
    		   "sort|r" => \$sort,
    		   "percent|p=f{,}" => \@percent,
    		   "any|y" => \$any,
    		   "extend|e=i{,}" => \@extend,
    		   "mode|m=i" => \$mode,
    		   "autoextend|u" => \$autoxtend,
    		   "both|t" => \$both,
    		   "exact|c" => \$exact,
    		   "pass|n=i" => \$npass,
    		   "gap|g=i" => \$agap,
    		   "output|o=s{,}" => \@out,
    		   "multi|m" => \$multi,
    		   "header|d" => \$header,
    		   "waitbar|w" => \$waitbar,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&programUsage;
    	exit;
    }
    $stop .= "--- Please specify input files ---\n" if (!($fileA || $fileB));
    # Check if percent in proper range
    if ($percent[0] =~ /\d\:\d+/) 
    {
    	my ($s,$e) = split(":",$percent[0]);
    	@percent = ($s..$e);
	}
	foreach my $cpp (@percent)
	{
		if ($cpp < 0 || $cpp > 100)
		{
			$stop .= "--- Overlap percentage should be a value between 0 and 100 ---\n";
			last;
		}
	}
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($npass < 0)
	{
		disp("--pass option should be >0... Using default (3)");
		$npass = 3;
	}
    # Check gap
    disp("--gap option parameter should be >=0... Using default (10000)") if ($agap < 0);
    # Check what is given on extend
    if (@extend)
    {
    	if (!$extend[1])
    	{
    		disp("--extend parameter has one argument... Assuming same extension in both sides.");
    		$extend[1] = $extend[0];
    	}
    }
    if (@extend && $autoxtend)
    {
    	disp("--extend and --autoextend options cannot be given together... Ignoring --autoextend.");
    	$autoxtend = 0;
	}
    $mode -= 2 if ($mode);
    # Check proper output format
    if (@out)
    {
		foreach my $c (@out)
		{
			if ($c ne "overlapA" && $c ne "overlapB" && $c ne "onlyA" && $c ne "onlyB" && $c ne "overpairs" &&
			    $c ne "nonpairs")
			{
				my $msg = "WARNING! --output options should be one or more of \"overlapA\", \"overlapB\",".
				          " \"onlyA\", \"onlyB\", \"overpairs\" or \"nonpairs\",\n".
						  "Using default (\"overlapA\")...";
				disp($msg);
				@out = ("overlapA");
			}
		}
	}
	else
	{
		if (!$percent[1])
		{
			disp("Output file type not given... Using default (overlapA)") ;
			@out = ("overlapA");
		}
	}
}

# Waitbar initiation
# 1st argument: the length of the waitbar
sub waitbarInit
{
	my $initlen = $_[0];
	$initlen = 50 if ($initlen eq "");
	my $printlen = ' 'x$initlen;
	print "\nProgress\n";
	print "|$printlen|\n";
	print("|");
}

# Update waitbar (floor and ceil from POSIX module)
# 1st argument: the current index, 2nd argument: the total length
# 3rd argument: the length of the waitbar  (optional)
sub waitbarUpdate
{
	my $curr = $_[0];
	my $tot = $_[1];
	my $waitbarlen = $_[2];
	$waitbarlen = 50 if ($waitbarlen eq "");
	my $step;
	if ($tot > $waitbarlen)
	{
		$step = ceil($tot/$waitbarlen);
		print "#" if ($curr%$step == 0);
	}
	else
	{
		$step = floor($waitbarlen/$tot);
		print "#" x $step;
	}
	if ($curr == $tot)
	{
		my $rem = $waitbarlen - floor($tot/$step);
		($rem != 0) ? (print "#" x $rem."|\n") : print "|\n";
	}
}

#sub printSimpleOutput
#{
	#my $infileA = shift @_;
	#my $infileB = shift @_;
	#my $he = shift @_;
	#my $otype = shift @_;
	#my %inhash = @_;
	#my ($outchr,$ind,@tdata);
	#my $outfilename = &createOutputFile($infileA,$infileB,$otype);
	#my @chrs = keys(%inhash);
	#open(OUTPUT,">$outfilename");
	#print OUTPUT "$he" if ($he);
	#foreach $outchr (sort(@chrs))
	#{		
		#@tdata = @{$inhash{$outchr}};
		#for ($ind=0; $ind<@tdata; $ind++)
		#{
			#print OUTPUT "$outchr\t$tdata[$ind]\n";
		#}
	#}
	#close(OUTPUT);
#}

sub printComplexOutput
{
	my $infileA = shift @_;
	my $infileB = shift @_;
	my $he = shift @_;
	my $otype = shift @_;
	my %inhash = @_;
	my ($outchr,$ind,$outhash,@k);
	my $outfilename = &createOutputFile($infileA,$infileB,$otype);
	my @chrs = keys(%inhash);
	open(OUTPUT,">$outfilename");
	print OUTPUT "$he" if ($he);
	foreach $outchr (sort(@chrs))
	{
		$outhash = $inhash{$outchr};
		@k = sort outsort keys(%$outhash);
		for ($ind=0; $ind<@k; $ind++)
		{
			print OUTPUT "$outchr\t$k[$ind]\n";
		}
	}
	close(OUTPUT);
}

sub printArray
{
	my $infileA = shift @_;
	my $infileB = shift @_;
	my $otype = shift @_;
	my @inarr = @_;
	my $outfilename = &createOutputFile($infileA,$infileB,$otype);
	open(OUTPUT,">$outfilename");
	print OUTPUT join("\n",@inarr);
	close(OUTPUT);
}

sub countHoH
{
	my %inhash = @_;
	my ($c,$h,$i,$count);
	foreach $c (keys(%inhash))
	{
		$h = $inhash{$c};
		$count+=keys(%$h); 
	}
	$count = 0 if (!$count);
	return($count);
}

sub getLengths
{
	my %inhash = @_;
	my @lens;
	my ($c,$l,@t);
	foreach $c (keys(%inhash))
	{
		foreach $l (@{$inhash{$c}})
		{
			@t = split(/\t/,$l);
			#push(@lens,$t[$autoxtend-1]);
			push(@lens,$t[1] - $t[0]);
		}
	}
	return(@lens);
}

sub createOutputFile
{
	my $inA = shift @_;
	my $inB = shift @_;
	my $type = shift @_;
	my ($baseA,$dirA,$extA) = fileparse($inA,'\..*?');
	my $baseB = fileparse($inB,'\..*?');
	($multi) ? (return($dirA.$baseA.$baseB)) :
	(return($dirA.$baseA."_".$baseB."_".$type.$extA));
}

sub countLines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
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
A perl program to calculate intersection between bed files (can contain as
many additional columns as desired, as long as the first 3 are chromosome,
start, end). It requires as input only the bed files and returns the overlap
regions between the first and the second file. However, if fileA and fileB
are intersected, the final file containing name OVERLAP...fileA output will 
contain regions from fileA that overlap with regions of fileB and the final 
file containing name OVERLAP...B will contain those regions from fileB.
The ONLY files contain regions only in fileA or only in fileB, when performing 
the intersection of fileA and fileB. The actions performed are similar to 
those of the UCSC Table Browser, only faster, in the sense that a simple command 
is required for multiple outputs with many options. Intersection can be done at 
specific overlap percentages or any overlap (like Table Browser). Extension from
a point (e.g. the peak mode) towards both directions is also possible as long as 
a column with this point is given in both files and the appropriate option is 
used. The user has the option to retrieve only certain of the four available 
output file types.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --inputA fileA --inputB fileB [OPTIONS]

--- Required ---
  --inputA|a  file  First input file
  --inputB|b  file  Second input file
--- Optional ---
  --sort|r		Use this option to sort the input files first. NECESSARY
  			if they are not sorted beforehand as the script uses a
  			binary search algorithm which requires sorted input.
  --percent|p		The overlap percentage between the two files. It  
			should be a value between 0 and 100. Alternatively, the
			overlap percentage can be a series of values. In this case
			the program will run as a batch procedure, calculating the
			number of overlapping regions for the determined output types
			for each given percentage, creating and writing as output the
			distribution of the number of overlapping regions. In this
			case, no other output files are produced. Instead of a series
			of numbers, it can be an element in the format a:b where a is
			the starting percentage and b the ending. A list between a and
			b sith interval 1 will be autogeneratied. This options does 
			not work with --any option. Examples are --percent 50, 
			--percent 10 20 30 40 50, --percent 1:100.
  --any|y		Use this switch for returning tags with any overlap
			between fileA and fileB. When used, --percent is ignored.
  --extend|e		Use this option to supply the program the upstream and
			downstream extension length. It should be two numbers, the first
			the number of bps to extend upstream of peak center and the second
			to extend downstream. Use this option in combination with --mode
			option. For example --extend 100 200.
  --autoextend|u	Use this switch to calculate half of the median peak 
			length which will then be used for extending from peak mode. Use
			this switch in combination with --mode option to provide the column
			with the peak mode.
  --mode|m		Use this option to supply the program the column in both
			fileA and fileB that contains the peak mode or any point from
			which extension to left and right will start. E.g. --mode 4.
			This parameter has to be provided when --extend and/or
			--autoextend are used. Also when overpairs is chosen as an
			output format.
  --both|t		Normally, when overlaping fileA with fileB at a percentage
			level, if the tag of the fileA overlaps with tag(s) from fileB,
			this percentage will be calculated based on tag from fileA.
			That is, if the required percentage is p and the overlap is d,
			then if d < p*length(tagA) then there is no overlap. However, if
			tagB is smaller than tagA it could d > p*length(tagB), so there
			is p percent overlap between tagA and tagB when using the length
			of tagB. This option allows for such comparisons by checking
			both cases and returning a positive overlap if one of the above
			two is true. This swicth is useless when used with --extend
			option and ignored.
  --exact|c		Use this switch to calculate overlaps based on exact 
  			genomic locations rather than when a region is "totally"
  			included inside the other, when performing overlaps based
  			on percentages. For example if overlaping fileA with fileB,
  			if a region in fileA starts at 300 and ends at 800 and the
  			overlaping region in fileB starts at 500 and ends at 700,
  			if using the --exact swicth, a 50% percent overlap will not
  			return this region as positive because (700-500)<0.5*(800-300)
  			even if tagB is totally included in tagA. This switch forces
  			the program to calculate percent overlaps exactly as the 
  			UCSC Table Browser. Can also be used in combination with the
  			--both switch.
  --pass|n		Use this option to supply the program the number of times
			that the dynamic binary search algorithm will search regions
			from the fileB for each region of the fileA. One pass returns
			at maximum one hit because if the algorithm finds a hit from
			fileB, it will exit. The number of passes determines how many
			times the algorithm will search for overlapping regions of
			fileA and fileB according to specified overlapping crteria. 
			Use a larger number of passes for when regions from fileA are
			likely to overlap many regions from fileB (e.g. when fileA has
			a large peak which could correspond to more than one peak in
			fileB). It defaults to 3.
  --gap|g		Use this option to retrieve distances between non-
  			overlaping (according to specified criteria) regions 
  			that is not larger than the number specified with the
  			--gap option. Can be useful with onlyA and nonpairs
  			output options. The default gap is set to 0 so the option
  			is not used. If you wish to use that option you should
  			provide a positive integer, for example --gap 10000.
  --output|o		Use this option to determine which intersection output
			filetypes you wish to retrieve.	Possible choices are:
   				"overlapA" for retrieving regions from fileA 
				that overlap with fileB.
   				"overlapB" for retrieving those regions that 
				overlap was found with fileA regions when using 
				fileA as first input (ATTENTION! NOT REGIONS FROM
				fileB THAT OVERLAP WITH fileA).
				"onlyA" for retrieving regions from fileA that
				DO NOT overlap with regions in fileB.
				"onlyB" for retrieving those regions that overlap
				was NOT found with fileA when using fileA as first
				input (ATTENTION! SAME BEHAVIOUR AS "overlapB"
				choice).
  --header|d		Use this option if you have a header line in
			your input files. It will also be written in your output
			files. Should be the same! Defaults to no header.
  --waitbar|w		Use this option if you wish to display a simple  
			progress bar while selecting and printing the final files.   
			For small files it is probably useless as the program finishes
			very quickly.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.

Usage examples:

perl intersectbed.pl --inputA A.bed --inputB B.bed --sort --any --output overlapA 
onlyA --waitbar

perl intersectbed.pl --inputA A.bed --inputB B.bed --percent 50 --extend 4 100 --both
--output overlapA onlyA onlyB --waitbar
	
The main output of the program is up to four files in BED format containing also
any additional data columns.

END
	print $usagetext;
	exit;
}

# Correct autoextend (only -1 as chromosome is not separated this time)
#$autoxtend-- if ($autoxtend);

# Template for R venn diagrams
#library(VennDiagrams)
#overrideTriple=TRUE
#draw.triple.venn(
	#area1=100,
	#area2=80,
	#area3=60,
	#n12=30,
	#n23=25,
	#n13=20,
	#n123=10,
	#category=c("One","Two","Three"),
	#fill=c("red","green","blue"),
	#lty="blank",
	#cex=2,
	#cat.cex=2,
	#cat.col=c("red","green","blue"),
	#cat.fontfamily=c("arial","arial","arial")
#)

