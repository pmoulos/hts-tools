#!/usr/bin/perl

# Description later...
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 13 - 05 - 2009 (dd - mm - yyyy)
# Last Update : 11 - 12 - 2009 (dd - mm - yyyy)
# Version     : 1.2
#
# Version 1.1 assigns a proper scoring scheme in scanning output bed conversion so that it
# can be used in combination with the useScore attribute of UCSC genome bropwser and the
# conversion creates proper (accroding to browser format) 6-column bed files.
# Version 1.2 supports scanning with MotifScanner (Thijs et al., 2001)

# TODO List
# Incorporate support for MotifLocator
# Put an option to get background from same chromosome
# Automate background generation (somewhen...)

use strict;
use File::Basename;
use Getopt::Long;
use Switch;
use Cwd;
use POSIX qw(floor
 			 ceil);
 			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "scanthemall.pl";
our @seqfile;		      # The peak-region-whatever file(s) (FASTA format)
our $motfile;		      # A file with PWMs 
our $backfile;	          # A possible background regions file to choose from
our @cntfile;		      # Original peak files, containing the centers in case of peak ids in fasta
our @cntcol;		      # Columns with peak IDs and peak centers in cntfile
our @range;		  	      # Range of cutoffs to use with p53Scan
our $scanner = "p53scan"; # Which scanner to use
our $fpr = 0.05; 	   	  # Allowed false positive rate to be used with p53Scan
our $times = 10;	   	  # How many times larger background?
our $length = 400;	      # Length for background sequences
our @spacer;		      # Spacer length for motif
our $besthit = 1;      	  # How many best matches to return (default 1)
our @output;		      # gff? bed? stats? report?
our $unistats = 0;        # In stats files, return unique numbers?
our $justscan = 0;		  # Just scan given only one threshold, no background
our $waitbar = 0;	      # Use waitbar?
our $header = 0;	      # Do peak files for bed generation contain headers?	
our $silent = 0;	      # Display verbose messages
our $help = 0;		      # Help?

# Check inputs
&checkInputs;

# Kavouries... tetoia ora tetoia logia...
my ($ogff,$obed,$ostat);
our $olog;
foreach my $o (@output)
{
	$ogff = 1 if ($o eq "gff");
	$obed = 1 if ($o eq "bed");
	$ostat = 1 if ($o eq "stats");
	$olog = 1 if ($o eq "log");
}

# Start logging if chosen
our @log;
push(@log,"##########     scanthemall.pl logfile     ##########\n") if ($olog);

# Bavard
my $date = &now("display");
disp("$date - Started...\n");

# More bavard... I like information...
disp("Selected scanner : $scanner");
($justscan) ? (disp("Cutoff range : $justscan")) : (disp("Cutoff range : ",$range[0],"to",$range[$#range]));
disp("False Positive Rate : $fpr");
disp("Chosen output(s) : ",join("\, ",@output),"\n");

# If bed output chosen and cntfile given, read them and merge them to database
my %cnthash;
if ($obed && @cntfile)
{
	disp("Merging files containing peak centers...");
	foreach my $c (@cntfile)
	{
		my ($line,@columns);
		open(CNT,"$c");
		$line = <CNT> if ($header);
		while ($line = <CNT>)
		{
			@columns = split(/\t/,$line);
			$cnthash{$columns[$cntcol[0]]} = $columns[$cntcol[1]];
		}
		close(CNT);
	}
}
elsif ($obed && !(@cntfile || @cntcol))
{
	disp("WARNING! BED output requested but no files containing peak centers are given...");
	disp("Proceeding only with co-ordinates as will can be derived from sequence IDs...");
	disp("Assuming sequence ID format : chr:start-end or chr:start:end...\n");
}

my ($i,$j,$k); # General indices
my @stats; # Stats array in case

# Parse the motifs file, split it into different files
my @mfiles = &parseMotifs($motfile,$scanner);

# Check what is going on with spacers if scanner is p53scan
if ($scanner eq "p53scan")
{
	if (@spacer)
	{
		my $nmf = @mfiles;
		die "\nNumber of spacer lengths is different from number of given motifs! Exiting...\n" if ($nmf != @spacer);
	}
	else # Generate spacers vector if not given (0-length spacers)
	{
		for ($i=0; $i<@mfiles; $i++)
		{
			push(@spacer,0);
		}
	}
}

# If running the full process, not just scanning given only one threshold
if (!$justscan)
{
	# Quick and dirty checking about tabular file
	my $chk = &dirtyChkTab($backfile);
	my $backbackup = $backfile;
	if ($chk)
	{
		disp("Background file $backfile does not appear to be a tabular sequence file...");
		disp("Checking if $backfile is in FASTA format...");
		my $chkchk = &dirtyChkFASTA($backfile);
		die "\nBackground file $backfile does not appear to be in FASTA format either! Exiting...\n" if ($chkchk);
		disp("$backfile is in FASTA format. Converting to tabular format...");
		$backbackup = $backfile;
		$backfile = &fasta2tab($backfile);
	}

	# Count number of background sequences
	disp("Counting sequences in background file $backfile...");
	my $nback = &countLines($backfile); 
	disp("Background file $backfile contains $nback sequences.");

	# Construct background index file for quick reference
	my ($bbi,$bbd) = fileparse($backfile,'\..*?');
	$bbi = $bbd.$bbi.".idx";
	if (! -e $bbi)
	{
		disp("Constructing index file for background file $backfile...");
		open(BACK,"< $backfile");
		open(INDEX,"+>$bbi") or die "Can't open $bbi for read/write: $!\n";
		buildIndex(*BACK,*INDEX);
		close(BACK);
		close(INDEX);
		disp("Index file for background constructed in $bbi.");
	}
	else
	{
		disp("Index file $bbi for background file $backfile already exists. Proceeding...\n");
	}

	if ($scanner =~ /MotifScanner/i)
	{
		disp("Constructing Markov background model to be used with MotifScanner...");
		&constructMSBackground($backbackup);
	}

	# Do job for peak-region-whatever files
	for ($i=0; $i<@seqfile; $i++)
	{
		disp("Processing $seqfile[$i]...");
		my $chkat1 = 0;
		my $chk = &dirtyChkFASTA($seqfile[$i]);
		if ($chk)
		{
			disp("File $seqfile[$i] does not appear to be a FASTA file. Proceeding to next...");
			next;
		}
		else { $chkat1 = 1; }
		
		# Open, count lines, call background etc.
		disp("Counting sequences in file $seqfile[$i]...");
		my $nseqs = &countFASTA($seqfile[$i]);
		disp("File $seqfile[$i] has $nseqs sequences.");
		
		# Generate adjusted background for each input file
		my $n = $times*$nseqs;
		disp("Generating $n background sequences to be used for scanning of $seqfile[$i]...");
		my $bfile = &getRandomSeq($backfile,$bbi,$nback,$length,$n);
		disp("$n background sequences written in file $bfile...\n");
		
		# For each motif, calculate threshold and scan according to selected scanner
		switch($scanner)
		{
			case /p53scan/i
			{
				my ($sl,$cutoff);
				for ($j=0; $j<@mfiles; $j++)
				{
					disp("Calculating score cutoff based on sequences in $bfile for motif $mfiles[$j] using p53scan...");
					disp("Output will be written in bgcut.gff.\n");
					# Use range vector
					for  ($k=0; $k<@range; $k++)
					{
						disp("Now testing threshold $range[$k]...\n");
						`python p53scan.py -i $bfile -c $range[$k] -s $spacer[$j] -p $mfiles[$j] -n $besthit > bgcut.gff `;
						if ($besthit > 1) 
						{
							($unistats) ? ($sl = &countUniLines("bgcut.gff")) : ($sl = &countLines("bgcut.gff"));
						}
						else { $sl = &countLines("bgcut.gff"); }
						if ($sl <= $fpr*$n)
						{
							$cutoff = $range[$k];
							last;
						}
						disp("FPR $fpr not reached... $sl matches out of $n sequences remain in background.");
					}
					if ($k == @range)
					{
						disp ("Last number of sequences remaining in background : $sl");
						disp ("No cutoff within given range satisfies FPR criteria... Skipping to next.\n");
						$cutoff = 0;
						$stats[$i][$j] = 0 if ($ostat);
						next;
					} 
					else 
					{ 
						disp("Cutoff for FPR $fpr determined at $cutoff.");
						disp("$sl matches out of $n sequences remain in background.\n");
					}
					disp("Scanning $seqfile[$i] for motif $mfiles[$j] using p53scan... Cutoff: $cutoff.\n");
					my $mbase = fileparse($mfiles[$j],'\..*?');
					my $currout = &createOutputFile($seqfile[$i],"output",$mbase);
					`python p53scan.py -i $seqfile[$i] -c $cutoff -s $spacer[$j] -p $mfiles[$j] -n $besthit > $currout `;
					&convert2bed($currout,%cnthash) if ($obed);
					my $nmatch;
					if ($besthit > 1) 
					{
						($unistats) ? ($nmatch = &countUniLines($currout)) : ($nmatch = &countLines($currout));
					}
					else { $nmatch = &countLines($currout); }
					$stats[$i][$j] = $nmatch if ($ostat);
					disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
					if (!$ogff)
					{
						($^O !~ /MSWin/) ? (`rm $currout `) : (`del /f /q $currout `);
					}
				}
			}
			
			case /MotifScanner/i
			{
				my ($sl,$cutoff);
				for ($j=0; $j<@mfiles; $j++)
				{
					disp("Calculating score cutoff based on sequences in $bfile for motif $mfiles[$j] using MotifScanner...");
					disp("Output will be written in bgcut.gff.\n");
					# Use range vector
					for  ($k=0; $k<@range; $k++)
					{
						disp("Now testing p-value threshold $range[$k]...\n");
						($^O !~ /MSWin/) ? (`./MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $range[$k] -s 0 -o bgcut.gff `) :
						(`MotifScanner -f $bfile -b MSmodel.bkg -m $mfiles[$j] -p $range[$k] -s 0 -o bgcut.gff `);
						($unistats) ? ($sl = &countUniLines("bgcut.gff")) : ($sl = &countLines("bgcut.gff"));
						$sl--; # One line header of MotifScanner output
						if ($sl <= $fpr*$n) 
						{
							$cutoff = $range[$k];
							last;
						}
						disp("FPR $fpr not reached... $sl matches out of $n sequences remain in background.");
					}
					if ($k == @range)
					{
						disp ("Last number of sequences remaining in background : $sl");
						disp ("No cutoff within given range satisfies FPR criteria... Skipping to next.\n");
						$cutoff = 0;
						$stats[$i][$j] = 0 if ($ostat);
						next;
					} 
					else 
					{ 
						disp("Cutoff for FPR $fpr determined at $cutoff.");
						disp("$sl matches out of $n sequences remain in background.\n");
					}
					disp("Scanning $seqfile[$i] for motif $mfiles[$j]... Cutoff: $cutoff.\n");
					my $mbase = fileparse($mfiles[$j],'\..*?');
					my $currout = &createOutputFile($seqfile[$i],"output",$mbase);
					($^O !~ /MSWin/) ? (`./MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $cutoff -s 0 -o $currout `) :
					(`MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $cutoff -s 0 -o $currout `);
					&convert2bed($currout,%cnthash) if ($obed);
					my $nmatch;
					($unistats) ? ($nmatch = &countUniLines($currout)) : ($nmatch = &countLines($currout));
					$nmatch--; # One line header of MotifScanner output
					$stats[$i][$j] = $nmatch if ($ostat);
					disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
					if (!$ogff)
					{
						($^O !~ /MSWin/) ? (`rm $currout `) : (`del /f /q $currout `);
					}
				}
			}
		}
		disp(" ");
	}
}
else # Just scan the sequence files using one defined threshold
{
	# Do job for peak-region-whatever files
	for ($i=0; $i<@seqfile; $i++)
	{
		disp("Processing $seqfile[$i]...");
		my $chkat1 = 0;
		my $chk = &dirtyChkFASTA($seqfile[$i]);
		if ($chk)
		{
			disp("File $seqfile[$i] does not appear to be a FASTA file. Proceeding to next...");
			next;
		}
		else { $chkat1 = 1; }

		# For each motif, calculate threshold and scan according to selected scanner
		switch($scanner)
		{
			case /p53scan/i
			{
				for ($j=0; $j<@mfiles; $j++)
				{
					disp("Scanning $seqfile[$i] for motif $mfiles[$j] using p53scan... Cutoff: $justscan.\n");
					my $mbase = fileparse($mfiles[$j],'\..*?');
					my $currout = &createOutputFile($seqfile[$i],"output",$mbase);
					`python p53scan.py -i $seqfile[$i] -c $justscan -s $spacer[$j] -p $mfiles[$j] -n $besthit > $currout `;
					&convert2bed($currout,%cnthash) if ($obed);
					my $nmatch;
					if ($besthit > 1) 
					{
						($unistats) ? ($nmatch = &countUniLines($currout)) : ($nmatch = &countLines($currout));
					}
					else { $nmatch = &countLines($currout); }
					$stats[$i][$j] = $nmatch if ($ostat);
					disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
					if (!$ogff)
					{
						($^O !~ /MSWin/) ? (`rm $currout `) : (`del /f /q $currout `);
					}
				}
			}
			
			case /MotifScanner/i
			{
				for ($j=0; $j<@mfiles; $j++)
				{
					disp("Scanning $seqfile[$i] for motif $mfiles[$j]... Cutoff: $justscan.\n");
					my $mbase = fileparse($mfiles[$j],'\..*?');
					my $currout = &createOutputFile($seqfile[$i],"output",$mbase);
					($^O !~ /MSWin/) ? (`./MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o $currout `) :
					(`MotifScanner -f $seqfile[$i] -b MSmodel.bkg -m $mfiles[$j] -p $justscan -s 0 -o $currout `);
					&convert2bed($currout,%cnthash) if ($obed);
					my $nmatch;
					($unistats) ? ($nmatch = &countUniLines($currout)) : ($nmatch = &countLines($currout));
					$nmatch--; # One line header of MotifScanner output
					$stats[$i][$j] = $nmatch if ($ostat);
					disp("$nmatch matches found in $seqfile[$i]. Output written in $currout.\n");
					if (!$ogff)
					{
						($^O !~ /MSWin/) ? (`rm $currout `) : (`del /f /q $currout `);
					}
				}
			}
		}
		disp(" ");
	}
}

# Print stats if requested
if ($ostat)
{
	my $outstat = &createOutputFile(" ","stats");
	open(STATS,">$outstat");
	print STATS "\t",join("\t",@mfiles),"\n";
	for ($i=0; $i<@seqfile; $i++)
	{
		print STATS "$seqfile[$i]";
		for ($j=0; $j<@mfiles; $j++)
		{
			print STATS "\t$stats[$i][$j]";
		}
		print STATS "\n";
	}
	close(STATS);
}

$date = &now("display");
disp("$date - Finished!\n\n");

if ($olog)
{
	my $logfile = &createOutputFile(" ","log");
	open(LOG,">$logfile");
	print LOG join("\n",@log),"\n";
	close(LOG);
}


# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@seqfile,
    		   "motif|m=s" => \$motfile,
    		   "background|b=s" => \$backfile,
    		   "scanner|a=s" => \$scanner,
    		   "center|c=s{,}" => \@cntfile,
    		   "colext|x=i{,}" => \@cntcol,
    		   "range|r=s{,}" => \@range,
    		   "fpr|p=f" => \$fpr,
    		   "times|n=i" => \$times,
    		   "length|l=i" => \$length,
			   "spacer|r=i{,}" => \@spacer,
			   "besthit|e=i" => \$besthit,
    		   "output|o=s{,}" => \@output,
    		   "uniquestats|u" => \$unistats,
    		   "justscan|j=f" => \$justscan,
    		   "waitbar|w" => \$waitbar,
    		   "header|d" => \$header,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&programUsage;
    	exit;
    }
    $stop .= "--- Please specify input sequence file(s) ---\n" if (!@seqfile);
    $stop .= "--- Please specify motifs file ---\n" if (!$motfile);
    $stop .= "--- Please specify background sequences file ---\n" if (!$backfile && !$justscan);
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    # Check the scanner
    if ($scanner !~ /p53scan/i && $scanner !~ /MotifScanner/i)
    {
    	disp("The scanner should be one of \"p53scan\" or \"MotifScanner\". Using default (p53scan)...");
    	$scanner = "p53scan";
	}
	# Check range - increase per real number, very simple expression, use with caution
	if (@range)
	{
		if (@range == 1 && $range[0] =~ /\d+\:\d+(\.\d*)?\:\d+/)
		{
			my ($s,$inc,$e) = split(":",$range[0]);
			if (($s > 1 || $inc > 1 || $e > 1) && $scanner =~ /MotifScanner/i)
			{
				disp("Range for MotifScanner should be <1. Using default (0.2:0.01:0.7)...");
				@range = &rangeVector(0.2,0.7,0.01);
			}
			else { @range = &rangeVector($s,$e,$inc); }
		}
		elsif (@range == 1 && $range[0] =~ /\d+\:\d+/) # Increase per 1 if p53scan or per 0.1 if MotifScanner
		{
			my ($s,$e) = split(":",$range[0]);
			switch($scanner)
			{
				case /p53scan/i
				{
					@range = ($s..$e);
				}
				case /MotifScanner/i
				{
					if ($s > 1 || $e > 1)
					{
						disp("Range for MotifScanner should be <1. Using default (0.2:0.1:0.7)...");
						@range = &rangeVector(0.5,0.01,0.01);
					}
					else { @range = &rangeVector($s,$e,0.1); }
				}
			}
		}
	}
	else
	{
		if (!$justscan)
		{
			disp("Cutoff range not given. Using default (8..13)...");
			@range = (8..13);
		} else { $range[0] = $justscan; }
	}
	disp("Spacer(s) not given. Default spacing (0) will be used for all motifs...") if (!@spacer && $scanner =~ /p53scan/i);
	# Check fpr
	$stop .= "--- False Positive Rate should be a value between 0 and 1 ---\n" if ($fpr<0 || $fpr>1);
	if (@output)
    {
		foreach my $c (@output)
		{
			if ($c ne "gff" && $c ne "bed" && $c ne "stats" && $c ne "log")
			{
				my $msg = "WARNING! --output options should be one or more of \"gff\", \"bed\", \"stats\" or \"log\"\n".
						  "Using default (\"gff\")...";
				disp($msg);
				@output = ("gff");
			}
		}
	}
	if (@cntcol)
    {
    	my $l = @cntcol;
    	if ($l != 3)
    	{
    		disp("ID, center columns and downstream extension not given properly... Using defaults (1,2,200)...");
    		@cntcol = (0,1,200);
		}
		else
		{
    		$cntcol[0]--;
    		$cntcol[1]--;
		}
    }
    if ($besthit < 0)
	{
		disp("--besthit option should be >0... Using default (1)");
		$besthit = 1;
	}
}

sub getRandomSeq
{
	my ($tabfasta,$itsindex,$itslen,$length,$num) = @_;
	my ($line,$id,$seq,$currlen,$start,$end);
	my $count = my $safeswitch = 0;
	my $BIG = 1e+6;
	srand;
	my $outfile = &createOutputFile(" ","sequence");
	&waitbarInit(50) if ($waitbar);
	open(RANDSEQ,">$outfile");
	open(FASTA,"<$tabfasta");
	open(INDEX,"$itsindex");
	while ($count < $num && $safeswitch < $BIG)
	{
		# Safeswitch in case too many sequences have small lengths, process shouldn't
		# take forever to complete...
		$safeswitch++;
		$line = getIndexedLine(*FASTA,*INDEX,int(rand($itslen)));
		($id,$seq) = split(/\t/,$line);
		$currlen = length($seq);
		next if ($currlen < $length);
		if ($currlen == $length)
		{
			($start,$end) = (1,$length);
			$id .= "_".$start."-".$end;
			&writeSeq(*RANDSEQ,$id,$seq);
			$count++;
		}
		if ($currlen > $length)
		{
			# Restrict the random index generation so that we don't go beyond the sequence end
			$start = int(rand($currlen - $length));
			$end = $start + $length;
			$id .= "_".$start."-".$end;
			&writeSeq(*RANDSEQ,$id,substr($seq,$start,$length));
			$count++;
		}
		&waitbarUpdate($count,$num,50) if ($waitbar);
	}
	close(RANDSEQ);
	close(FASTA);
	close(INDEX);
	if ($safeswitch >= $BIG)
	{
		disp("Sequence fetching discontinued... $count sequences fetched in total in $outfile...");
		disp("Probably the FASTA file you supplied to get random sequences from contains too many short sequences.");
		disp("Try again with larger sequences or smaller length.");
	}
	return($outfile);
}

sub writeSeq
{
	my $file = shift @_;
	my $id = shift @_;
	my $seq = shift @_;
	$id = ">".$id if ($id !~ /^>/);
	print $file "$id\n";
	while ($seq) # Write sequences of length 100
	{
		my $wseq = substr($seq,0,100,"");
		print $file "$wseq\n";
	}
}

# Usage: buildIndex(*DATAHANDLE,*INDEXHANDLE)
sub buildIndex 
{
    my $datafile = shift @_;
    my $indexfile = shift @_;
    my $offset = 0;
    while (<$datafile>) 
    {
        print $indexfile pack("N",$offset);
        $offset = tell($datafile);
    }
}

# Usage: getIndexedLine(*DATAHANDLE,*INDEXHANDLE, $LINENUMBER)
# Returns line or undef if LINENUMBER was out of range
sub getIndexedLine 
{
    my $datafile = shift @_;
    my $indexfile = shift @_;
    my $linenumber = shift @_;
    my $size;    # Size of an index entry
    my $ioffset; # Offset into the index of the entry
    my $entry;   # Index entry
    my $doffset; # Offset into the data file

    $size = length(pack("N",0));
    $ioffset = $size * ($linenumber - 1);
    seek($indexfile,$ioffset,0) or return;
    read($indexfile,$entry,$size);
    $doffset = unpack("N",$entry);
    seek($datafile,$doffset,0);
    return scalar(<$datafile>);
}

sub constructMSBackground
{
	my $infile = shift @_;
	if (-e "MSmodel.bkg")
	{
		disp("Background model MSmodel.bkg already exists. Proceeding...");
		return 0;
	}
	my ($base,$dir,$ext) = fileparse($infile,'\..*?');
	my $chkfa = &dirtyChkFASTA($infile);
	if ($chkfa && ! -e $dir.$base.".fa")
	{
		disp("Background file $infile does not appear to be in FASTA format...");
		disp("Checking if $infile is in tabular format...");
		my $chktab = &dirtyChkTab($infile);
		die "\nBackground file $infile does not appear to be in tabular format either! Exiting...\n" if ($chktab);
		disp("$infile is in tabular format. Converting to FASTA format...");
		my $infile = &tab2fasta($infile);
	}
	else { $infile = $dir.$base.".fa"; }
	disp("Creating Markov model background from $infile... Output will be written in MSmodel.bkg\n");
	if ($^O !~ /MSWin/)
	{
		`./CreateBackgroundModel -f $infile -b MSmodel.bkg > temp.std `;
	 	`rm temp.std `;
	}
	else
	{
		`CreateBackgroundModel -f \"$infile\" -b MSmodel.bkg > temp.std `;
		`del /f /q temp.std `;
	}
	if (! -e "MSmodel.bkg")
	{
		disp("Markov background model not created!");
		disp("Probably program CreateBackgroundModel cannot be located on the system! Exiting...");
		die;
	}
}

sub parseMotifs
{
	my $infile = shift @_;
	my $scn = shift @_;
	my ($base,$dir) = fileparse($infile,'\..*?');
	my ($c,$line,$om,$f,@fhs,@outnames);
	open(MOTIF,"$infile");
	
	switch($scn)
	{
		case /p53scan/i
		{
			$line = <MOTIF>;
			if ($line =~ /^#/) # Case of known motif names
			{
				$c = 0; # Reset counter
				$line = <MOTIF>;
				die "\nMotif file (p53scan) does not appear to have the correct syntax!\n" if ($line !~ /^>/);
				seek(MOTIF,0,0);
				while ($line = <MOTIF>)
				{
					if ($line =~ /^#/) # Signal to open new file
					{
						$line =~ s/\r|\n$//g;
						$line =~ s/^#//g;
						$om = &createOutputFile($infile,$line);
						push(@outnames,$om);
						local *PWM;
						open(PWM,">$om");
						$fhs[$c] = *PWM;
						$c++;
					}
					else 
					{ 
						$f = $fhs[$c-1];
						print $f $line;
					}
				}
				# Close the split motifs
				foreach my $fh (@fhs) { close($fh); }
			}
			elsif ($line =~ /^>/) # Case of unknown motif names
			{
				$c = 0; # Reset file counter
				my $inc = 0; # > counter
				seek(MOTIF,0,0);
				while($line = <MOTIF>)
				{
					if ($line =~ /^>/ && $inc%2 == 0) # Signal to open new file
					{
						$om = &createOutputFile($infile,"Motif_$c"."_");
						push(@outnames,$om);
						local *PWM;
						open(PWM,">$om");
						$fhs[$c] = *PWM;
						$f = $fhs[$c];
						print $f $line;
						$inc++;
						$c++;
					}
					elsif  ($line =~ /^>/ && $inc%2 != 0)
					{
						$f = $fhs[$c-1];
						print $f $line;
						$inc++;
					}
					else
					{ 
						$f = $fhs[$c-1];
						print $f $line;
					}	
				}
				# Close the split motifs
				foreach my $fh (@fhs) { close($fh); }
			}
			else { die "\nMotif file does not appear to have the correct syntax!\n"; }
		}
		case /MotifScanner/i
		{
			# Some checking
			my $fline = <MOTIF>;
			$fline =~ s/\r|\n$//g;
			die "\nMotif file (MotifScanner) does not appear to have the correct syntax (line 1)!\n" if ($fline !~ /^#INCLUSive/i);
			$c = 0; # Reset counter
			my %motifs;
			while ($line = <MOTIF>)
			{
				$line =~ s/\r|\n$//g;
				if ($line =~ /^#ID/) # Signal to open new file
				{ 
					my @sep = split("=",$line);
					$sep[1] =~ s/^\s+|\s+$//g;
					my $name = $sep[1];
					$om = &createOutputFile($infile,$name);
					push(@outnames,$om);
					local *PWM;
					open(PWM,">$om");
					$fhs[$c] = *PWM;
					push(@{$motifs{$c}},$fline);
					$c++;
				}
				next if ($line =~ /^#W/); # Skip, we will create afterwards...
				push(@{$motifs{$c-1}},$line);		
			}
			foreach my $k (sort { $a <=> $b } keys(%motifs))
			{
				my @fileconts = @{$motifs{$k}};
				my $mw = $#fileconts - 1;
				splice(@fileconts,2,0,"#W = $mw");
				$f = $fhs[$k];
				print $f join("\n",@fileconts);
				print $f "\n";
				close($f)
			}
		}
	}
	
	return(@outnames);
}

# Convert the gff output from p53scan to bed format
# The first column of the gff (that is the peak/region ID) MUST contain coordinates
# information in the form chr:start-end (track2fasta) or chr:start:end (format I used).
# WARNING! If the fasta files used for scanning have been generated with a program like
# track2fasta, then the bed co-ordinates for each occurence can be correctly generated.
# If the sequence id's in the fasta files correspond to peak ids rather than exact
# sequence locations, another file with peak ids and peak centers must be provided
# Converts to 6 column bed files
# Converts also the motif score in gff file to the 0-1000 scale of UCSC genome browser
# so that motif strength can be visualized by color
# This is done by linear conversion of the form new_value = a*old_value + b and by solving
# the special case of a 2x2 linear system (since we know several of the values):
# min(score)*a + b = 0
# max(score)*a + b = 1000
sub convert2bed
{
	my $f2c = shift @_;
	my %ch = @_;
	my ($base,$dir) = fileparse($f2c,'\..*?');
	my ($line,@lines,@scores,@content,@locs,@newcoord);
	my $bedout = $dir.$base.".bed";
	# In order to determine the coefficients of linear conversion we have to suck in all
	# gff file, the hard way...
	open(F2C,$f2c);
	while ($line = <F2C>)
	{
		$line =~ s/\r|\n$//g;
		push(@lines,$line);
		@content = split(/\t/,$line);
		push(@scores,$content[5]);
	}
	close(F2C);
	# If the scanner was MotifScanner/MotifLocator, there is one xtra line...
	my $tort = shift @scores if (!$scores[0]);
	# Get min, max score
	my ($min,$max) = &minmax(@scores);
	# Get coefficients
	my ($a,$b) = &naiveSolveTwo($min,$max);
	open(BED,">$bedout");
	if (%ch)
	{
		foreach $line (@lines)
		{
			@content = split(/\t/,$line);
			@locs = split(":",$content[0]);
			if ($#locs == 1) # track2fasta format
			{
				my $joined = pop(@locs);
				my @splitted = split("-",$joined);
				push(@locs,@splitted);
			}
			# Shift coordinates to motif occurence
			my $strand = "+";
			$strand = "-" if ($content[6] == -1 || $content[6] eq "R");
			@newcoord = ($locs[0],
						 $ch{$content[0]} - $cntcol[2] + $content[3],
						 $ch{$content[0]} - $cntcol[2] + $content[4],
						 $content[0],$a*$content[5] + $b,$strand);
			print BED join("\t",@newcoord),"\n";
		}
	}
	else
	{
		foreach $line (@lines)
		{
			@content = split(/\t/,$line);
			@locs = split(":",$content[0]);
			if ($#locs == 1) # track2fasta format
			{
				my $joined = pop(@locs);
				my @splitted = split("-",$joined);
				push(@locs,@splitted);
			}
			# Shift coordinates to motif occurence
			my $strand = "+";
			$strand = "-" if ($content[6] == -1 || $content[6] eq "R");
			@newcoord = ($locs[0],$locs[1] + $content[3],$locs[1] + $content[4],
						 $content[0],$a*$content[5] + $b,$strand);
			print BED join("\t",@newcoord),"\n";
		}
	}
	close(BED);
}

sub naiveSolveTwo
{
	my ($min,$max) = @_;
	my $eps = 0.000001;
	my $div;
	(abs($min-$max) < $eps) ? ($div = $eps) : ($div = $max/$min);
	$min = $eps if (!$min);
	my $y = 1000/(1-($div));
	my $x = -$y/$min;
	return($x,$y);
}

sub createOutputFile
{
	use Switch; # TEMPORARY!
	my $in = shift @_;
	my $type = shift @_;
	my $subtype = shift @_;
	my $cdir = getcwd;
	my $date = &now("file");
	my $sep = &getSysSep;
	switch($type)
	{
		case /sequence/i
		{
			return($cdir.$sep."randseq"."$date".".fa");
		}
		case /stats/i
		{
			return($cdir.$sep."stats"."$date".".txt");
		}
		case /log/i
		{
			return($cdir.$sep."log"."$date".".txt");
		}
		case /output/i
		{
			my $base = fileparse($in,'\..*?');
			return($cdir.$sep.$base."_".$subtype.".gff");
		}
		else # Motif file
		{
			my ($base,$dir) = fileparse($in,'\..*?');
			return($dir.$type."motif.pwm");
		}
	}
}

sub fasta2tab
{
	my $infile = shift @_;
	my $currid;
	open(INPUT,$infile) or die "\nThe file $infile does not exist!\n";
	my ($base,$dir) = fileparse($infile,'\..*?');
	my $outfile = $dir.$base.".tab";
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

sub tab2fasta
{
	my $infile = shift @_;
	my @conts;
	open(INPUT,$infile) or die "\nThe file $infile does not exist!\n";
	my ($base,$dir) = fileparse($infile,'\..*?');
	my $outfile = $dir.$base.".fa";
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

sub countLines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
}

sub countFASTA
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totfa=0;
	$totfa += tr/\>/\>/ while sysread(IN,$_,2**16);
	close(IN);
	return $totfa;
}

sub countUniLines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my ($l,@li,%lh);
	while ($l = <IN>)
	{
		@li = split(/\t/,$l);
		$lh{$li[0]}++;
	}
	close(IN);
	my $totlines = keys(%lh);
	return $totlines;
}

sub rangeVector
{
	my ($s,$e,$inc) = @_;
	my @rout;
	my $counter = 0;
	my $temp = $s;
	# Create output vector according to order
	if ($s < $e)
	{
		while ($temp <= $e)
		{
			$rout[$counter] = $temp;
			$counter++;
			$temp += $inc;
		}
	}
	else
	{
		while ($temp >= $e)
		{
			$rout[$counter] = $temp;
			$counter++;
			$temp -= $inc;
		}
	}
	return(@rout);
}

sub dirtyChkFASTA
{
	my $in = shift @_;
	open(FCHK,$in) or die "\nThe file $in does not exist!\n";
	my $chk = 0;
	$chk = 1 if (<FCHK> !~ /^>/);
	close(FCHK);
	return($chk);
}

sub dirtyChkTab
{
	my $in = shift @_;
	open(TCHK,$in) or die "\nThe file $in does not exist!\n";
	my $chk = 0;
	my @test = split(/\t/,<TCHK>);
	my $len = @test;
	$chk = 1 if ($len < 2);
	close(TCHK);
	return($chk);
}

sub waitbarInit
{
	my $initlen = shift @_;
	$initlen = 50 if (!$initlen);
	my $printlen = " "x$initlen;
	print "\nProgress\n";
	print "|$printlen|\n";
	print("|");
}

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
		my $rem;
		($tot > $waitbarlen) ? ($rem = $waitbarlen - floor($tot/$step)) : 
		($rem = $waitbarlen - $tot*$step);
		($rem != 0) ? (print "#" x $rem."|\n") : print "|\n";
	}
}

sub now
{
	my $type = shift @_;
	my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
	$year += 1900;
	$month++;
	$month = "0".$month if (length($month)==1);
	$day = "0".$day if (length($day)==1);
	$hour = "0".$hour if (length($hour)==1);
	$min = "0".$min if (length($min)==1);
	$sec = "0".$sec if (length($sec)==1);
	if ($type eq "file")
	{
		return($day.$month.$year.$hour.$min.$sec);
	}
	elsif ($type eq "display")
	{
		return($day."/".$month."/".$year." ".$hour.":".$min.":".$sec);
	}
}

sub minmax
{
	my @s = sort { $a <=> $b } @_;
	return($s[0],$s[$#s]);
}

sub getSysSep
{
	($^O !~ /MSWin/) ? (return("/")) : (return("\\"))
}

sub disp
{
	my @in = @_;
	print "\n@in" if (!$silent);
	push(@log,"@in") if ($olog);
}

sub programUsage
{
	my $usagetext = << "END";

$scriptname
Help will be soon available.

END
	print $usagetext;
	exit;
}

#options are:
#--input : as many fasta files as you want to scan
#--motif : the input motifs file
#--background : a background sequences file formatted in TABULAR format
#--range : range of cutoffs for p53scan (for example if range=5:0.1:15 then it will perform searches for cutoffs 5 to 15 raised each time by 0.1)
#--fpr : The false positive rate (0.05)
#--times : how many times should be the background larger
#--length : length of sequences to select from background
#--spacer : as many spacers as the motifs included in the motif file (e.g. --spacer 3 0 0 0 )
#--output : types of output for hits, gff for gff file, bed for bed file, stats for file with numbers of hits formatted for all input seq files and motifs, log for a log file, can be all of them
#--waitbar : obvious :)
#--silent : no output messaged
#--help : not yet written

# Old BED conversion
#sub convert2bed
#{
	#my $f2c = shift @_;
	#my %ch = @_;
	#my ($base,$dir) = fileparse($f2c,'\..*?');
	#my ($line,@content,@locs,@newcoord);
	#my $bedout = $dir.$base.".bed";
	#open(F2C,$f2c);
	#open(BED,">$bedout");
	#if (%ch)
	#{
		#while ($line = <F2C>)
		#{
			#$line =~ s/\r|\n$//g;
			#@content = split(/\t/,$line);
			#@locs = split(":",$content[0]);
			#if ($#locs == 1) # track2fasta format
			#{
				#my $joined = pop(@locs);
				#my @splitted = split("-",$joined);
				#push(@locs,@splitted);
			#}
			## Shift coordinates to motif occurence
			#my $strand = "+";
			#$strand = "-" if ($content[6] == -1 || $content[6] eq "R");
			#@newcoord = ($locs[0],
						 #$ch{$content[0]} - $cntcol[2] + $content[3],
						 #$ch{$content[0]} - $cntcol[2] + $content[4],
						 #$content[0],$content[5],$strand);
			#print BED join("\t",@newcoord),"\n";
		#}
	#}
	#else
	#{
		#while ($line = <F2C>)
		#{
			#$line =~ s/\r|\n$//g;
			#@content = split(/\t/,$line);
			#@locs = split(":",$content[0]);
			#if ($#locs == 1) # track2fasta format
			#{
				#my $joined = pop(@locs);
				#my @splitted = split("-",$joined);
				#push(@locs,@splitted);
			#}
			## Shift coordinates to motif occurence
			#my $strand = "+";
			#$strand = "-" if ($content[6] == -1 || $content[6] eq "R");
			#@newcoord = ($locs[0],$locs[1] + $content[3],$locs[1] + $content[4],
						 #$content[0],$content[5],$strand);
			#print BED join("\t",@newcoord),"\n";
		#}
	#}
	#close(F2C);
	#close(BED);
#}
