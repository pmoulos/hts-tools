#!/usr/bin/perl

# Program to join reads in bed file(s) with their genomic regions in a genome and count
# presences of reads (uniquely) in each gene. Simplifies the operation on genomic intervals
# of Galaxy and solw by performing all operations locally and unified. Additionally, it
# addresses the problem of partial overlaping between a read and a genomic regions by 3
# possible ways:
# a. By providing an overlaping coefficient between 0 and 1, the program adds a read into
#	 a genomic region if the part of the read within the boundaries (start or end) of the
#	 genomic region is >= coefficient*read length. For example, if a read has length 200
#	 and the overlaping coefficient is 0.9, then 0.9*200 = 180 bases should be within the
#	 region boundaries to be included in that region. This should provide more accuracy
#	 when having large reads. It can be used in the case of small reads too, but lucks
#	 automation.
# b. By using a linear probabilistic score. In this case, if a read is located half inside
#    a region, it is added up to that region, else a score is calculated for the read based
#    on the formula: score = #bps inside region/(length of read/2). This score lies between
#	 0 and 1 and it is compared to a uniform random number p between 0 and 1. If p<score
#    then this tag is added up to the genomic region else discarded. This scoring scheme
#    should be used with reads of small length and should be avoided in the case of large
#	 reads because of causing possible under-representation in promoter regions, in case
#	 one wishes to summarize reads in gene regions. However it could perform well when
#	 studying promoter or binding regions.
# c. By using an exponential probabilistic score. This case, works as case (b) but the
#    scoring function is score = exp(-c^2/#bps inside region) where c is constant. The
#	 constant c determines determines the steepness of the score. The lower the value is,
#	 the higher is the probability to include in the region tags whose larger part lies
#	 outside the region. This scoring scheme should be used with reads of small length and
#    for reasons similar to those of (b) but would maybe perform better when used in
#	 promoter or binding regions.
# Apart from the bed files, the program requires also a bed file (3 first columns should
# be chromosome, start, end and the 4th a UNIQUE region identification, e.g. Ensembl ID, 
# could also contain further information in extra columns) which contains the genomic 
# regions of interest, in which number of reads should be summarized. The program returns 
# bed file(s) that contain the provided genomic regions with any other information together 
# with read counts per region.
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 26 - 05 - 2008 (dd - mm - yyyy)
# Last Update : 28 - 02 - 2013 (dd - mm - yyyy)
# Version     : 2.1
#
# CHANGELOG
#
# Version 1.1 implements a binary search algorithm to look for the proper regions in the
# regions hash instead of greedily searching through all of them. Additionally, it checks
# when the chromosome changes in the bed file so much time is saved from calling the "keys"
# and "values" functions to get hash elements.
# Version 1.2 allows the user to choose whether to discard or not counts in very small
# genomic regions.
# Version 1.21 introduces a bug fix in the binary search
# Version 2.0 is an almost complete rewrite of the original program, more efficient, more
# modularized and OO ready. It also works in paraller with multiple input files and the
# output is only one matrix-like file. The current implementation is slightly more memory
# consuming than the previous one, but faster, parallelizable and can easily return a table
# of counts and several stats with minimum effort.
# Version 2.1 allows for automatic download of region files in the case of Ensembl genes
# or exons for human, mouse, rat, fruitlfy and zebrafish

# TODO: Check how we can use gtf2tree.pl from the ngsplot package (https://code.google.com/p/ngsplot/)
# TODO: Make a note that when we are multi-coring, a large regions file (e.g. an organism
#		genes) takes longer because the hash is initialized for every file (so if we have
#		a lot of files like the tcga case it takes some time).

#################### SPECIAL SECTION ####################
sub cosort
{
	my @one = split("-",$a);
	my @two = split("-",$b);
	return 1 if ($one[0] > $two[0]);
	return -1 if ($one[0] < $two[0]);
	if ($one[0] == $two[0])
	{
		return 1 if ($one[1] > $two[1]);
		return 0 if ($one[1] == $two[1]);
		return -1 if ($one[1] < $two[1]);
	}
}
################## END SPECIAL SECTION ##################

use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
use File::Spec;
use File::Path qw(make_path remove_tree);
use Parallel::Loops;
use LWP::UserAgent;

use constant MAXCORES => 12;
use constant REMOTE_HOST => "genome-mysql.cse.ucsc.edu";
use constant REMOTE_USER => "genome";

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# On Ctrl-C or die, do cleanup
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "readcounter.pl";
our @infile; # Input file(s)
our $regionfile;	# Genomic regions file
our $sort = 0;	 # Sort files
our $percentage = 0.95; # Percentage of tag inside genomic region
our $lscore = 0; # Linear scoring to determine if tag partially inside genomic region included
our $escore = 0; # Exponential scoring to determine if tag partially inside genomic region included
our $c = 3; # Default constant to be used in the exponential scoring scheme
our $countsmall = 0;	# Count tags in genomic regions that are smaller than the tag length?
our $splitper = 0; # Split genomic regions in sub-regions of xbps and count individually
our $ncore = 1; # Number of cores for parallel calculations
our $stats = 0; # Return some statistics about tag distributions per area
our $output = ""; # The output file, can be auto
our $source; # Source in case of downloading data
our $waitbar = 0;	# Do not use waitbar
our $silent = 0; # Display verbose messages
our $help = 0; # Help?

# Advertise
&advertise;

# Check inputs
&check_inputs;

&try_module("Tie::IxHash::Easy");
&try_module("Term::ProgressBar") if ($waitbar);

# Global variables
our $tmpdir = File::Temp->newdir();

my @originfile = @infile; # Keep original filenames (verbose purposes)

# Read regions file in memory and initialize the hash that will host read counts per gene
# per file. They operate on bed files. Finally, write the output.
#FIXME Until given ... when statement becomes standard since Switch does not work properly anymore!
if ($regionfile =~ m/human-gene|mouse-gene|rat-gene|fly-gene|zebrafish-gene/i)
{
	$regionfile = &fetch_ensembl_genes($regionfile);
	$regionfile = &sort_ensembl_genes($regionfile);
}
elsif ($regionfile =~ m/human-exon|mouse-exon|rat-exon|fly-exon|zebrafish-exon/i)
{
	$regionfile = &fetch_ensembl_exons($regionfile);
	$regionfile = &sort_ensembl_exons($regionfile);
}
elsif ($regionfile =~ m/human-(5|3)utr|mouse-(5|3)utr|rat-(5|3)utr|fly-(5|3)utr|zebrafish-(5|3)utr/i)
{
	($regionfile =~ m/5utr/i) ? ($regionfile = &fetch_ensembl_utr($regionfile,5)) :
		($regionfile = &fetch_ensembl_utr($regionfile,3));
	$regionfile = &sort_ensembl_exons($regionfile);
}
elsif ($regionfile =~ m/human-cds|mouse-cds|rat-cds|fly-cds|zebrafish-cds/i)
{
	$regionfile = &fetch_ensembl_cds($regionfile);
	$regionfile = &sort_ensembl_exons($regionfile);
}
else { ($regionfile,@infile) = &sort_inputs($regionfile,@infile) if ($sort); }

my ($chromosome,$gencounts,$splitcounts,$theHeader) = &read_region_file($regionfile,\@originfile);
($gencounts,$splitcounts) = &count_all_reads($chromosome,$gencounts,$splitcounts,\@infile,\@originfile);
&write_reads($chromosome,$gencounts,$splitcounts,\@originfile,$theHeader);

# Remove garbage
&cleanup;

disp("Finished!\n\n");


#################### MAIN SUBROUTINES ####################

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@infile,
    		   "region|r=s" => \$regionfile,
    		   "sort|a" => \$sort,
    		   "percent|p=f" => \$percentage,
    		   "lscore|l" => \$lscore,
    		   "escore|e" => \$escore,
    		   "constant|c=f" => \$c,
    		   "small|m" => \$countsmall,
               "split|t=i" => \$splitper,
               "stats|z" => \$stats,
               "output|o=s" => \$output,
               "ncore|n=i" => \$ncore,
               "source|u=s" => \$source,
    		   "waitbar|w" => \$waitbar,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&program_usage;
    	exit;
    }
    $stop .= "--- Please specify input file(s) ---\n" if (!@infile);
    $stop .= "--- Please specify region file ---\n" if (!$regionfile);
    $stop .= "--- The supported genomes are organism-type, where organism is human, mouse, rat, fly or zebrafish and type is gene, exon, 5utr, 3utr or cds! ---\n"
		if ( !-f $regionfile && $regionfile !~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i);
    if ($stop)
    {
            print STDERR "\n$stop\n";
            print STDERR "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($lscore && $escore) # If both given, use exponential scoring
    {
    	disp("You chose both linear and exponential scoring. Only exponential will be used...");
    	$lscore = 0 ;
    }
    if ($stats && !$splitper)
    {
    	disp("You can't calculate area statistics without splitting to sub-areas. Option deactivated...");
    	$stats = 0 ;
    }
    if ($ncore)
    {
		eval "&try_module(\"Parallel::ForkManager\")";
		if ($@)
		{
			disp("Module Parallel::Loops not found, proceeding with one core...");
			$ncore = 1;
		}
		else { use Parallel::ForkManager; }
		if ($ncore > MAXCORES)
		{
			my $c = MAXCORES;
			disp("The maximum number of cores allowed is $c...");
			$ncore = MAXCORES;
		}
	}
	else { $ncore = 1; }
	if ($output eq "auto")
	{
		my ($base,$dir,$ext) = fileparse($infile[0],'\.[^.]*');
		$output = File::Spec->catfile($dir,$base."_REGIONCOUNTSTATS".$ext);
	}
	if ($regionfile =~ m/human-(gene|exon|(5|3)utr|cds)|mouse-(gene|exon|(5|3)utr|cds)|rat-(gene|exon|(5|3)utr|cds)|fly-(gene|exon|(5|3)utr|cds)|zebrafish-(gene|exon|(5|3)utr|cds)/i)
	{
		my %sources = ("ucsc" => "UCSC","refseq" => "RefSeq","ensembl" => "Ensembl");
		if ($source)
		{
			$source = lc($source);
			if ($source ~~ keys(%sources))
			{
				disp("Selected template regions source: ",$sources{$source});
			}
			else
			{
				disp("Source for template region files not given or is not well-defined! Using default (ucsc)...");
				$source = "ucsc";
			}
		}
	}
}

sub split_area
{
	my $start = shift @_;
	my $end = shift @_;
	my $splitlen = shift @_;
	
	my $regstart;
	my $regend = $start - 1;
	my @subareas;
	
	while ($regend < $end)
	{
		$regstart = $regend + 1;
		$regend = $regstart + $splitlen - 1;
		push(@subareas,$regstart."\t".$regend);
	}
	if ($regend > $end)
	{
		pop(@subareas);
		$regend = $end;
		push(@subareas,$regstart."\t".$regend);
	}
	
	return @subareas;
}

sub bin_search_loc
{
	my $start = shift @_;
	my $end = shift @_;
	my @areas = @_;
	my ($ind,$currstart,$currend,$center);
	my ($l,$u) = (0,$#areas);
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		$center = $start + &round(($end - $start)/2); # Location of center of the tag
		if ($currstart < $center && $currend > $center)
		{
			return (1,$ind);
		}
		elsif ($currstart == $center) # Randomly assign to one window
		{
			(int(10*rand(1)) < 5) ? (return(1,$ind-1)) : (return(1,$ind));
		}
		elsif ($currend == $center)
		{
			(int(10*rand(1)) < 5) ? (return(1,$ind)) : (return(1,$ind+1));
		}
		else
		{
			$u = $ind - 1 if ($center <= $currstart);
            $l = $ind + 1 if ($center >= $currend);
		}
	}
	return (0,-1);
}

sub read_region_file
{
	my $regionfile = shift @_;
	my @originfile = @{shift @_};
	@originfile = @infile unless(@originfile);
	
	my ($theHeader,$chr,$start,$end,$uid,$score,$strand,$rest);
	my @arest;

	my (%chromosome,%gencounts,%splitcounts);
	tie %chromosome, "Tie::IxHash::Easy";
	tie %gencounts, "Tie::IxHash::Easy";

	my ($regionlines,$rprog,$f);
	$regionlines = &count_lines($regionfile) if ($waitbar);

	my ($regstart,$regend,$regconts,$regcount,$regline);
	my @tempconts;

	disp("Reading genomic regions file...");
	disp("...also splitting genomic regions in areas of $splitper base pairs...") if ($splitper);
	if ($waitbar)
	{
		print "\n";
		$rprog = &waitbar_init($regionlines);
	}
	open(REGIONS,$regionfile) or die "\nThe file $regionfile does not exist!\n";
	$regline = <REGIONS>;
	$theHeader = &decide_header($regline);
	seek(REGIONS,0,0) if (!$theHeader);
	while ($regline = <REGIONS>)
	{
		$rprog->update($.) if ($waitbar);
		$regline =~ s/\r|\n$//g; # Make sure to remove carriage returns
		
		($chr,$start,$end,$uid,@arest) = split(/\t/,$regline);
		$rest = join("\t",@arest);
		($rest eq "") ? ($chromosome{$chr}{$uid} = $start."\t".$end) :
			($chromosome{$chr}{$uid} = $start."\t".$end."\t".$rest);
		
		foreach $f (@originfile)
		{
			($ncore > 1) ? ($gencounts{basename($f)}{$chr}{$uid} = 0) :
				($gencounts{$chr}{$uid}{basename($f)} = 0);

			if ($splitper)
			{
				$regcount = 0;
				$regend = $start - 1;
				while ($regend < $end)
				{
					$regstart = $regend + 1;
					$regend = $regstart + $splitper - 1;
					push(@tempconts,$regstart." ".$regend."\t");
					($ncore > 1) ? ($splitcounts{basename($f)}{$chr}{$uid}{$regcount} = 0) :
						($splitcounts{$chr}{$uid}{basename($f)}{$regcount} = 0);
					$regcount++;
				}
			}
		}
	}
	close(REGIONS);

	return(\%chromosome,\%gencounts,\%splitcounts,$theHeader);
}

sub count_all_reads
{
	my $chromosome = shift @_;
	my $gencounts = shift @_;
	my $splitcounts = shift @_;
	my @infile = @{shift @_};
	my @originfile = @{shift @_};
	@originfile = @infile unless(@originfile);
	
	my @base;
	foreach my $of (@originfile)
	{
		#my $bb = fileparse($of,'\.[^.]*');
		my $bb = basename($of);
		push(@base,$bb);
	}

	if ($ncore == 1)
	{
		for (my $i=0; $i<@infile; $i++)
		{
			($gencounts,$splitcounts) = &count_reads($chromosome,$gencounts,$splitcounts,$infile[$i],$originfile[$i]);
		}
	}
	else
	{
		my (%vhash1,%vhash2);
		foreach (my $i=0; $i<@infile; $i++)
		{
			$vhash1{$base[$i]} = $infile[$i];
			$vhash2{$base[$i]} = $originfile[$i];
		}
		#my $pl = Parallel::Loops->new($ncore);
		#$pl->share($chromosome);
		#$pl->share($gencounts);
		#$pl->share($splitcounts);
		#$pl->share(\%vhash1);
		#$pl->share(\%vhash2);
		#$pl->foreach(\@base, sub {
			#($gencounts->{$_},$splitcounts->{$_}) = &count_reads_multi($chromosome,$gencounts->{$_},$splitcounts->{$_},$vhash1{$_},$vhash2{$_});
		#});
		#my %pidh;
		my $pm = new Parallel::ForkManager($ncore);
		$pm -> run_on_finish(sub {
			my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data_structure_reference) = @_;
			if (defined($data_structure_reference)) {
				my @hr = @{$data_structure_reference};
				while(my($chr,$genehash) = each(%{$hr[0]}))
				{		
					while(my($genename,$count) = each(%{$genehash}))
					{
						$gencounts->{$ident}->{$chr}->{$genename} = $count;
					}
				}
				if ($splitper)
				{
					while(my($chr,$genehash) = each(%{$hr[1]}))
					{		
						while(my($genename,$idhash) = each(%{$genehash}))
						{
							while(my($id,$scount) = each(%{$idhash}))
							{
								$splitcounts->{$ident}->{$chr}->{$genename}->{$id} = $scount;
							}
						}
					}
				}
			} else {
				print STDERR qq|No message received from child process $pid!\n|;
			}
		});
		foreach my $ff (@base)
		{
			my $pid = $pm->start($ff) and next;
			my @arr = &count_reads_multi($chromosome,$gencounts->{$ff},$splitcounts->{$ff},$vhash1{$ff},$vhash2{$ff});
			$pm->finish(0,\@arr);
		}
		$pm->wait_all_children;
	}
	
	return($gencounts,$splitcounts);
}

sub count_reads
{
	my $chromosome = shift @_;
	my $gencounts = shift @_;
	my $splitcounts = shift @_;
	my $infile = shift @_;
	my $originfile = shift @_;
	$originfile = $infile unless($originfile);

	my (@k,@v);
	my ($gene,$coords,$diff,$score,$p,$filename,$bprog,$numberlines);
	my ($tmp1,$tmp2);
	my ($bedchr,$bedstart,$bedend,$bedrest);
	my ($currhash,$currchr,$currstart,$currend,$currrest);
	my (@regarr,@bsr);

	$numberlines= &count_lines($infile) if ($waitbar);
	
	disp("Reading and processing file $originfile...");

	if ($waitbar)
	{
		print "\n";
		$bprog = &waitbar_init($numberlines,"Counting...","yes");
	}
	
	$filename = basename($originfile);
	open(BEDFILE,$infile);
	my $nextchr = "";

	while (<BEDFILE>)
	{
		$bprog->update($.) if ($waitbar);
		chomp $_;
		($bedchr,$bedstart,$bedend,$bedrest) = split(/\t/,$_);

		if ($bedchr ne $nextchr)
		{
			$currhash = $chromosome->{$bedchr};
			@k = keys(%$currhash);
			@v = values(%$currhash);
			$nextchr = $bedchr;
		}

		# Start binary search
		my $i;
		my ($l,$u) = (0,$#k);
		while ($l <= $u)
		{
			$i = int(($l + $u)/2);
			my $gene = $k[$i];
			my $coords = $v[$i];
			($currstart,$currend,$currrest) = split(/\t/,$coords);

			# Conditions that satisfy our search
			# Security, too small region, tag does not fit?
			# Mark that there is something there if user wishes and proceed
			if ($currstart >= $bedstart && $currend <= $bedend)
			{
				$gencounts->{$bedchr}->{$gene}->{$filename}++ if ($countsmall);
				if ($splitper && $countsmall) # If splitting of regions
				{
					@regarr = &split_area($currstart,$currend,$splitper);
					@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
					$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
				}
				last;
			}
			# Ideal case... whole tag inside gene... no problem
			elsif ($currstart <= $bedstart && $currend >= $bedend)
			{
				$gencounts->{$bedchr}->{$gene}->{$filename}++;
				if ($splitper) # If splitting of regions
				{
					@regarr = &split_area($currstart,$currend,$splitper);
					@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
					$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
				}
				last;
			}
			# Part of tag outside after region
			elsif ($currstart < $bedstart && $currend < $bedend && $bedstart < $currend)
			{
				$diff = $currend - $bedstart;
				if ($lscore || $escore) # Linear or exponential scoring
				{
					if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside no need for further check
					{
						$gencounts->{$bedchr}->{$gene}->{$filename}++;
						if ($splitper) # If splitting of regions
						{
							@regarr = &split_area($currstart,$currend,$splitper);
							@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
							$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
						}
					}
					else # Use probabilistic scoring scheme
					{
						$score = ($currend - $bedstart)/(($bedend - $bedstart)/2) if ($lscore); # Linear
						$score = exp(-$c**2/($currend - $bedstart)) if ($escore); # Exponential
						$p = rand();
						if ($p < $score)
						{
							$gencounts->{$bedchr}->{$gene}->{$filename}++; 
							if ($splitper) # If splitting of regions
							{
								@regarr = &split_area($currstart,$currend,$splitper);
								@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
								$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
							}
						}
					}
				}
				else # Simple overlap
				{
					if ($diff >= $percentage*($bedend - $bedstart))
					{
						$gencounts->{$bedchr}->{$gene}->{$filename}++;
						if ($splitper) # If splitting of regions
						{
							@regarr = &split_area($currstart,$currend,$splitper);
							@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
							$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
						}
					}
				}
				last;
			}
			# Part of tag outside before region
			elsif ($currstart > $bedstart && $currend > $bedend && $bedend > $currstart)
			{
				$diff = $bedend - $currstart;
				if ($lscore || $escore) # Linear or exponential scoring
				{
					if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside
					{
						$gencounts->{$bedchr}->{$gene}->{$filename}++;
						if ($splitper) # If splitting of regions
						{
							@regarr = &split_area($currstart,$currend,$splitper);
							@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
							$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
						}
					}
					else # Use probabilistic scoring scheme
					{
						$score = ($bedend - $currstart)/(($bedend - $bedstart)/2) if ($lscore); # Linear
						$score = exp(-$c**2/($bedend - $currstart)) if ($escore); # Exponential
						$p = rand();
						if ($p < $score)
						{
							$gencounts->{$bedchr}->{$gene}->{$filename}++;
							if ($splitper) # If splitting of regions
							{
								@regarr = &split_area($currstart,$currend,$splitper);
								@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
								$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
							}
						}
					}
				}
				else # Simple overlap
				{
					if ($diff >= $percentage*($bedend - $bedstart))
					{
						$gencounts->{$bedchr}->{$gene}->{$filename}++;
						if ($splitper) # If splitting of regions
						{
							@regarr = &split_area($currstart,$currend,$splitper);
							@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
							$splitcounts->{$bedchr}->{$gene}->{$filename}->{$bsr[1]}++ if ($bsr[0]);
						}
					}
				}
				last;
			}
			else # If none of them met, reduce the searching subset
			{
				$u = $i - 1 if ($bedend <= $currstart);
				$l = $i + 1 if ($bedstart >= $currend);
			}

		}
	}
	close(BEDFILE);

	return($gencounts,$splitcounts);
}

sub count_reads_multi
{
	my $chromosome = shift @_;
	my $gencounts = shift @_;
	my $splitcounts = shift @_;
	my $infile = shift @_;
	my $originfile = shift @_;
	$originfile = $infile unless($originfile);

	my (@k,@v);
	my ($gene,$coords,$diff,$score,$p,$filename,$bprog,$numberlines);
	my ($bedchr,$bedstart,$bedend,$bedrest);
	my ($currhash,$currchr,$currstart,$currend,$currrest);
	my (@regarr,@bsr);

	$numberlines= &count_lines($infile) if ($waitbar);
	
	disp("Reading and processing file $originfile...");

	if ($waitbar)
	{
		print "\n";
		$bprog = &waitbar_init($numberlines,"Counting...","yes");
	}
	
	$filename = basename($originfile);
	open(BEDFILE,$infile);
	my $nextchr = "";

	while (<BEDFILE>)
	{
		$bprog->update($.) if ($waitbar);
		chomp $_;
		($bedchr,$bedstart,$bedend,$bedrest) = split(/\t/,$_);

		if ($bedchr ne $nextchr)
		{
			$currhash = $chromosome->{$bedchr};
			@k = keys(%$currhash);
			@v = values(%$currhash);
			$nextchr = $bedchr;
		}

		#use Data::Dumper;
		#print Dumper(\@v);

		# Start binary search
		#my %found;
		#my $found;
		my $n = 0;
		#while ($n < 3)
		#{
			my $i;
			#$found = 0;
			my ($l,$u) = (0,$#k);
			while ($l <= $u)
			{
				$i = int(($l + $u)/2);
				my $gene = $k[$i];
				my $coords = $v[$i];
				($currstart,$currend,$currrest) = split(/\t/,$coords);

				# Conditions that satisfy our search
				# Security, too small region, tag does not fit?
				# Mark that there is something there if user wishes and proceed
				if ($currstart >= $bedstart && $currend <= $bedend)
				{
					#next if $found{$i};
					$gencounts->{$bedchr}->{$gene}++ if ($countsmall);
					if ($splitper && $countsmall) # If splitting of regions
					{
						@regarr = &split_area($currstart,$currend,$splitper);
						@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
						$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
					}
					#$found{$i}++; # Find a way not to splice the array of coordinates, check if the found value is the same...
					#$found = $i;
					last;
				}
				# Ideal case... whole tag inside gene... no problem
				elsif ($currstart <= $bedstart && $currend >= $bedend)
				{
					#next if $found{$i};
					$gencounts->{$bedchr}->{$gene}++;
					#print "\t$filename $gencounts->{$bedchr}->{$gene}";
					if ($splitper) # If splitting of regions
					{
						@regarr = &split_area($currstart,$currend,$splitper);
						@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
						$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
					}
					#$found{$i}++;
					#$found = $i;
					last;
				}
				# Part of tag outside after region
				elsif ($currstart < $bedstart && $currend < $bedend && $bedstart < $currend)
				{
					#next if $found{$i};
					$diff = $currend - $bedstart;
					if ($lscore || $escore) # Linear or exponential scoring
					{
						if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside no need for further check
						{
							$gencounts->{$bedchr}->{$gene}++;
							if ($splitper) # If splitting of regions
							{
								@regarr = &split_area($currstart,$currend,$splitper);
								@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
								$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
							}
						}
						else # Use probabilistic scoring scheme
						{
							$score = ($currend - $bedstart)/(($bedend - $bedstart)/2) if ($lscore); # Linear
							$score = exp(-$c**2/($currend - $bedstart)) if ($escore); # Exponential
							$p = rand();
							if ($p < $score)
							{
								$gencounts->{$bedchr}->{$gene}++; 
								if ($splitper) # If splitting of regions
								{
									@regarr = &split_area($currstart,$currend,$splitper);
									@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
									$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
								}
							}
						}
					}
					else # Simple overlap
					{
						if ($diff >= $percentage*($bedend - $bedstart))
						{
							$gencounts->{$bedchr}->{$gene}++;
							if ($splitper) # If splitting of regions
							{
								@regarr = &split_area($currstart,$currend,$splitper);
								@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
								$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
							}
						}
					}
					#$found{$i}++;
					#$found = $i;
					last;
				}
				# Part of tag outside before region
				elsif ($currstart > $bedstart && $currend > $bedend && $bedend > $currstart)
				{
					#next if $found{$i};
					$diff = $bedend - $currstart;
					if ($lscore || $escore) # Linear or exponential scoring
					{
						if ($diff > ($bedend - $bedstart)/2) # If half of the tag inside
						{
							$gencounts->{$bedchr}->{$gene}++;
							if ($splitper) # If splitting of regions
							{
								@regarr = &split_area($currstart,$currend,$splitper);
								@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
								$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
							}
						}
						else # Use probabilistic scoring scheme
						{
							$score = ($bedend - $currstart)/(($bedend - $bedstart)/2) if ($lscore); # Linear
							$score = exp(-$c**2/($bedend - $currstart)) if ($escore); # Exponential
							$p = rand();
							if ($p < $score)
							{
								$gencounts->{$bedchr}->{$gene}++;
								if ($splitper) # If splitting of regions
								{
									@regarr = &split_area($currstart,$currend,$splitper);
									@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
									$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
								}
							}
						}
					}
					else # Simple overlap
					{
						if ($diff >= $percentage*($bedend - $bedstart))
						{
							$gencounts->{$bedchr}->{$gene}++;
							if ($splitper) # If splitting of regions
							{
								@regarr = &split_area($currstart,$currend,$splitper);
								@bsr = &bin_search_loc($bedstart,$bedend,@regarr);
								$splitcounts->{$filename}->{$bedchr}->{$gene}->{$bsr[1]}++ if ($bsr[0]);
							}
						}
					}
					#$found{$i}++;
					#$found = $i;
					last;
				}
				else # If none of them met, reduce the searching subset
				{
					$u = $i - 1 if ($bedend <= $currstart);
					$l = $i + 1 if ($bedstart >= $currend);
				}
			}
			#if ($found)
			#{
				##print STDERR "\nfound is $found!";
				#splice(@k,$found,1); # Remove it from areas else it will be found again
				#splice(@v,$found,1);
				#$n++;
				##print STDERR "\nI found it in $n!";
			#} else { last; }
		}
	#}
	close(BEDFILE);

	#print Dumper($gencounts);

	#return($gencounts->{$filename},$splitcounts->{$filename});
	return($gencounts,$splitcounts);
}

sub write_reads
{
	my $chromosome = shift @_;
	my $gencounts = shift @_;
	my $splitcounts = shift @_;
	my @originfile = @{shift @_};
	my $theHeader = shift @_;
	@originfile = @infile unless(@originfile);
	
	my (@outrest,@distr,@filenames,@outchrs,@outgeneids);
	my ($outgene,$outcount,$files,$filename,$b,$retrhash,$outhead,$outsubhead);
	my ($outchr,$outgeneid,$outcounts,$outstart,$outend,$outrest,$finalout);
	my ($mean,$median,$stdev,$mad,$jdistr);
	
	my @base;
	foreach my $of (@originfile)
	{
		my $bb = fileparse($of,'\.[^.]*');
		push(@base,$bb);
	}

	my $out;
	if ($output)
	{
		open(OUTPUT,">$output");
		$out = *OUTPUT;
	}
	else { $out = *STDOUT; }

	if ($output)
	{
		($stats) ? (disp("Writing reads and calculating statistics per genomic regions in $output...")) :
			disp("Writing reads per genomic regions for file in $output...");
	}
	else
	{
		($stats) ? (disp("Writing reads and calculating statistics per genomic regions to standard output...\n")) :
			disp("Writing reads per genomic regions for file to standard output...\n");
	}

	my ($regionlines,$fprog);
	if ($waitbar)
	{
		print "\n";
		$regionlines = &count_lines($regionfile);
		$fprog = &waitbar_init($regionlines,"Writing...");
	}
	$outcount = 1;

	if ($theHeader) # Add a header based on what we did plus additional data given
	{
		$outhead = $theHeader;
		if (scalar @base == 1)
		{
			$outhead .= "\t".$base[0];
			$outhead .= "\tSub-area Counts/$splitper bp" if ($splitper); # Return also sub-area distributions
			$outhead .= "\tMean Counts/$splitper bp\tMedian Counts/$splitper bp\tStDev Counts\tMAD Counts" if ($stats);
			print $out "$outhead\n";
		}
		else
		{
			foreach my $b (@base)
			{
				$outhead .= "\t".$b;
				$outhead .= "\t" if ($splitper);
				$outhead .= "\t"x4 if ($stats);
			}
			# We must also print a subheader in case of reporting stats
			$outsubhead = "Counts";
			$outsubhead .= "\tSub-area Counts/$splitper bp" if ($splitper);
			$outsubhead .= "\tMean Counts/$splitper bp\tMedian Counts/$splitper bp\tStDev Counts\tMAD Counts" if ($stats);
			$outsubhead = "$outsubhead"x(scalar @base);
			$outsubhead = ("\t"x(scalar split(/\t/,$theHeader))).$outsubhead;
			print $out "$outhead\n";
			print $out "$outsubhead\n" if ($splitper);
		}
	}

	if ($ncore == 1)
	{
		while(($outchr,$outgene) = each(%$gencounts))
		{
			while (($outgeneid,$files) = each(%$outgene))
			{
				$outcount++;
				$fprog->update($outcount) if ($waitbar);
				
				($outstart,$outend,@outrest) = split(/\t/,$chromosome->{$outchr}->{$outgeneid});
				if (scalar @outrest != 0)
				{
					$outrest = join("\t",@outrest);
					$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid."\t".$outrest;
				}
				else
				{
					$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid;
				}

				while (($filename,$outcounts) = each(%$files))
				{
					if ($splitper) # Retrieve calculated distributions per region
					{
						$retrhash = $splitcounts->{$outchr}->{$outgeneid}->{$filename};
						@distr = values(%$retrhash);
						$jdistr = join(" ",@distr);
						if ($stats) # Calculate stats
						{
							$mean = &mean(@distr);
							$median = &median(@distr);
							$stdev = &stdev(@distr);
							$mad = &mad(@distr);
						}
					}
					$finalout .= "\t".$outcounts;
					$finalout .= "\t".$jdistr if ($splitper);
					$finalout .= "\t".$mean."\t".$median."\t".$stdev."\t".$mad if ($stats);
				}
				print $out "$finalout\n";
			}
		}
	}
	else
	{
		@outchrs = keys(%$chromosome);
		@filenames = keys(%$gencounts);
		foreach my $outchr (@outchrs)
		{
			@outgeneids = keys(%{$chromosome->{$outchr}});
			foreach $outgeneid (@outgeneids)
			{
				$outcount++;
				$fprog->update($outcount) if ($waitbar);

				($outstart,$outend,@outrest) = split(/\t/,$chromosome->{$outchr}->{$outgeneid});
				if (scalar @outrest != 0)
				{
					$outrest = join("\t",@outrest);
					$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid."\t".$outrest;
				}
				else
				{
					$finalout = $outchr."\t".$outstart."\t".$outend."\t".$outgeneid;
				}

				foreach $filename (@filenames)
				{
					$outcounts = $gencounts->{$filename}->{$outchr}->{$outgeneid};
					if ($splitper) # Retrieve calculated distributions per region
					{
						$retrhash = $splitcounts->{$filename}->{$outchr}->{$outgeneid};
						@distr = values(%$retrhash);
						$jdistr = join(" ",@distr);
						if ($stats) # Calculate stats
						{
							$mean = &mean(@distr);
							$median = &median(@distr);
							$stdev = &stdev(@distr);
							$mad = &mad(@distr);
						}
					}
					$finalout .= "\t".$outcounts;
					$finalout .= "\t".$jdistr if ($splitper);
					$finalout .= "\t".$mean."\t".$median."\t".$stdev."\t".$mad if ($stats);
				}
				print $out "$finalout\n";
			}
		}
	}
	close(OUTPUT) if ($output);
}


#################### HELP SUBROUTINES ####################

sub round
{
	my $number = shift;
	return int($number + .5*($number <=> 0));
}

sub mean 
{
  my $result;
  foreach (@_) { $result += $_ ;}
  return $result/@_;
}

sub stdev 
{
  my $mean = mean(@_);
  my @elemsquared;
  foreach (@_)
  {
    push (@elemsquared,($_**2));
  }
  return sqrt(mean(@elemsquared) - ($mean**2));
}	

sub median 
{
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

sub mad
{
	my @absdiff;
	my $med = median(@_);
	foreach (@_)
	{
		push(@absdiff,($_ - $med));
	}
	return abs(median(@absdiff));
}

sub sort_inputs
{
	my $regionfile = shift @_;
	my @infile = @_;
	my $tmpfile;
	
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		for (my $i=0; $i<@infile; $i++)
		{
			disp("Sorting bed file $infile[$i]...");
			$tmpfile = File::Spec->catfile($tmpdir,"temp"."$i".".in$$");
			`sort -k1,1 -k2g,2 $infile[$i] > $tmpfile `;
			$infile[$i] = $tmpfile;
		}
		disp("Sorting region file $regionfile...");
		$tmpfile = File::Spec->catfile($tmpdir,"tempreg.in$$");
		`sort -k1,1 -k2g,2 $regionfile > $tmpfile `;
		$regionfile = $tmpfile;
	}
	else # We are in Windows... package required
	{
		&tryModule("File::Sort","sort_file");
		eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
		for (my $i=0; $i<@infile; $i++)
		{
			disp("Sorting file $infile[$i]...");
			$tmpfile = File::Spec->catfile($tmpdir,"temp"."$i".".tmp");
			sort_file(
			{
				I => $infile[$i],
				o => $tmpfile,
				k => ['1,1','2n,2'],
				t => "\t"
			});
			$infile[$i] = $tmpfile;
		}
		disp("Sorting region file $regionfile...");
		$tmpfile = File::Spec->catfile($tmpdir,"tempreg.tmp");
		sort_file(
		{
			I => $regionfile,
			o => $tmpfile,
			k => ['1,1','2n,2'],
			t => "\t"
		});
		$regionfile = $tmpfile;
	}

	return($regionfile,@infile);
}

sub sort_one
{
	my $file = shift @_;
	my $tmpfile;
	
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		disp("Sorting file $file...");
		$tmpfile = File::Spec->catfile($tmpdir,"temptss.in$$");
		`sort -k1,1 -k2g,2 $file > $tmpfile `;
		$file = $tmpfile;
	}
	else # We are in Windows... package required
	{
		&tryModule("File::Sort","sort_file");
		eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
		disp("Sorting file $file...");
		$tmpfile = File::Spec->catfile($tmpdir,"temptss.tmp");
		sort_file(
		{
			I => $file,
			o => $tmpfile,
			k => ['1,1','2n,2'],
			t => "\t"
		});
		$file = $tmpfile;
	}

	return($file);
}

sub sort_ensembl_genes
{
	my $infile = shift @_;
	my $tmpfile;
	
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		disp("Sorting bed file $infile...");
		$tmpfile = File::Spec->catfile($tmpdir,"temp".".in$$");
		`awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $infile > $tmpfile `;
		$infile = $tmpfile;
	}
	else # We are in Windows... package required not able to sort file with header...
	{
		my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
				    "messing up data. Please sort files outside $scriptname first (e.g. using\n".
				    "Excel or something similar.";
		die "\n$dmsg\n\n";
		#&tryModule("File::Sort","sort_file");
		#eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
		#disp("Sorting file $infile...");
		#$tmpfile = File::Spec->catfile($tmpdir,"temp".".tmp");
		#sort_file(
		#{
			#I => $infile,
			#o => $tmpfile,
			#k => ['1,1','2n,2'],
			#t => "\t"
		#});
		#$infile = $tmpfile;
	}

	return($infile);
}

sub sort_ensembl_exons
{
	my $infile = shift @_;
	my $tmpfile;
	
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		disp("Sorting bed file $infile...");
		$tmpfile = File::Spec->catfile($tmpdir,"temp".".in$$");
		`awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2 -k3g,3 -u\"}' $infile > $tmpfile `;
		$infile = $tmpfile;
	}
	else # We are in Windows... package required not able to sort file with header...
	{
		my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
				    "messing up data. Please sort files outside $scriptname first (e.g. using\n".
				    "Excel or something similar.";
		die "\n$dmsg\n\n";
		#&tryModule("File::Sort","sort_file");
		#eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
		#disp("Sorting file $infile...");
		#$tmpfile = File::Spec->catfile($tmpdir,"temp".".tmp");
		#sort_file(
		#{
			#I => $infile,
			#o => $tmpfile,
			#k => ['1,1','2n,2'],
			#t => "\t"
		#});
		#$infile = $tmpfile;
	}

	return($infile);
}

sub try_module
{
	my $module = shift @_;
	my @fun = @_;
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

sub decide_header
{
	my $line = $_[0];
	$line =~ s/\r|\n$//g;
	my @cols = split(/\t/,$line);
	if ($cols[0] =~ m/^chr/ && $cols[1] =~ m/\d+/ && $cols[2] =~ m/\d+/)
	{
		return(0); # Does not contain a header, is proper bed line
	}
	else
	{
		return($line);
	}
}

sub fetch_ensembl_genes
{
	my $sp = &format_species($_[0]);
	my $xml = &get_xml_genes_query($sp);
	my $path="http://www.biomart.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
	my $ua = LWP::UserAgent->new;
	my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
	my $regs = File::Spec->catfile($tmpdir,"regions.txt");
	my $response;

	# We need to put data in a temporary file because it's scrambled by asynchronicity
	disp("Querying Biomart...");
	$ua->request($request,sub {
		my ($data,$response) = @_;
		if ($response->is_success) {
			print $tmpfh "$data";
		}
		else {
			warn ("Problems with the web server: ".$response->status_line);
		}
	},1000);

	seek($tmpfh,0,SEEK_SET);
	my %strand = ( 1 => "+", -1 => "-" );
	open(REGS,">$regs");
	print REGS "chromosome\tstart\tend\tensembl_id\tgc_content\tstrand\tgene_name\n";
	while (my $line = <$tmpfh>)
	{
		next if ($line !~ m/^[0-9XY]/);
		$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\n";
	}
	close(REGS);
	close($tmpfh);

	return($regs);
}

sub fetch_ensembl_exons
{
	my $sp = &format_species($_[0]);
	my $xml = &get_xml_exons_query($sp);
	my $path="http://www.biomart.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
	my $ua = LWP::UserAgent->new;
	my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
	my $regs = File::Spec->catfile($tmpdir,"regions.txt");
	my $response;

	# We need to put data in a temporary file because it's scrambled by asynchronicity
	disp("Querying Biomart...");
	$ua->request($request,sub {   
		my ($data,$response) = @_;
		if ($response->is_success) {
			print $tmpfh "$data";
		}
		else {
			warn ("Problems with the web server: ".$response->status_line);
		}
	},1000);

	seek($tmpfh,0,SEEK_SET);
	my %strand = ( 1 => "+", -1 => "-" );
	open(REGS,">$regs");
	print REGS "chromosome\tstart\tend\tensembl_exon_id\tensembl_gene_id\tstrand\tgene_name\n";
	while (my $line = <$tmpfh>)
	{
		next if ($line !~ m/^[0-9XY]/);
		$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\n";
	}
	close(REGS);
	close($tmpfh);

	return($regs);
}

sub fetch_ensembl_utr
{
	my $sp = &format_species($_[0]);
	my $utr = $_[1];
	my $xml = &get_xml_utr_query($sp,$utr);
	my $path="http://www.biomart.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
	my $ua = LWP::UserAgent->new;
	my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
	my $regs = File::Spec->catfile($tmpdir,"regions.txt");
	my $response;

	# We need to put data in a temporary file because it's scrambled by asynchronicity
	disp("Querying Biomart...");
	$ua->request($request,sub {
		my ($data,$response) = @_;
		if ($response->is_success) {
			print $tmpfh "$data";
		}
		else {
			warn ("Problems with the web server: ".$response->status_line);
		}
	},1000);

	seek($tmpfh,0,SEEK_SET);
	my %strand = ( 1 => "+", -1 => "-" );
	open(REGS,">$regs");
	print REGS "chromosome\tstart\tend\tensembl_id\ttranscript_count\tstrand\tname\n";
	while (my $line = <$tmpfh>)
	{
		next if ($line !~ m/^[0-9XY]/);
		$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		next if (!$cols[1]);
		print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\n";
	}
	close(REGS);
	close($tmpfh);

	return($regs);
}

sub fetch_ensembl_cds
{
	my $sp = &format_species($_[0]);
	my $xml = &get_xml_cds_query($sp);
	my $path="http://www.biomart.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
	my $ua = LWP::UserAgent->new;
	my $tmpfh = File::Temp->new(DIR => $tmpdir,SUFFIX => ".ens");
	my $regs = File::Spec->catfile($tmpdir,"regions.txt");
	my $response;

	# We need to put data in a temporary file because it's scrambled by asynchronicity
	disp("Querying Biomart...");
	$ua->request($request,sub {   
		my ($data,$response) = @_;
		if ($response->is_success) {
			print $tmpfh "$data";
		}
		else {
			warn ("Problems with the web server: ".$response->status_line);
		}
	},1000);

	seek($tmpfh,0,SEEK_SET);
	my %strand = ( 1 => "+", -1 => "-" );
	open(REGS,">$regs");
	print REGS "chromosome\tstart\tend\tensembl_exon_id\tensembl_gene_id\tstrand\tgene_name\n";
	while (my $line = <$tmpfh>)
	{
		next if ($line !~ m/^[0-9XY]/);
		$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		next if (!$cols[1]);
		print REGS "chr$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$strand{$cols[5]}\t$cols[6]\n";
	}
	close(REGS);
	close($tmpfh);

	return($regs);
}

sub get_xml_genes_query
{
	my $species = $_[0];
	my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
			  "<!DOCTYPE Query>".
			  "<Query virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
			  "<Dataset name = \"$species\" interface = \"default\" >\n".
			  "<Attribute name = \"chromosome_name\" />".
			  "<Attribute name = \"start_position\" />".
			  "<Attribute name = \"end_position\" />".
			  "<Attribute name = \"ensembl_gene_id\" />".
			  "<Attribute name = \"percentage_gc_content\" />".
			  "<Attribute name = \"strand\" />".
			  "<Attribute name = \"external_gene_id\" />".
			  "</Dataset>\n".
			  "</Query>\n";
	return($xml);
}

sub get_xml_exons_query
{
	my $species = $_[0];
	my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
		  "<!DOCTYPE Query>".
		  "<Query virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
		  "<Dataset name = \"$species\" interface = \"default\" >\n".
		  "<Attribute name = \"chromosome_name\" />".
		  "<Attribute name = \"exon_chrom_start\" />".
		  "<Attribute name = \"exon_chrom_end\" />".
		  "<Attribute name = \"ensembl_exon_id\" />".
		  "<Attribute name = \"ensembl_gene_id\" />".
		  "<Attribute name = \"strand\" />".
		  "<Attribute name = \"external_gene_id\" />".
		  "</Dataset>\n".
		  "</Query>\n";
	return($xml);
}

sub get_xml_utr_query
{
	my ($species,$utr) = @_;
	my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
		"<!DOCTYPE Query>\n".
		"<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
		"<Dataset name = \"$species\" interface = \"default\" >\n".
		"<Attribute name = \"chromosome_name\" />\n".
		"<Attribute name = \"".$utr."_utr_start\" />\n".
		"<Attribute name = \"".$utr."_utr_end\" />\n".
		"<Attribute name = \"ensembl_gene_id\" />\n".
		"<Attribute name = \"transcript_count\" />\n".
		"<Attribute name = \"strand\" />\n".
		"<Attribute name = \"external_gene_id\" />\n".
		"</Dataset>\n".
		"</Query>\n";
	return($xml);
}

sub get_xml_cds_query
{
	my $species = $_[0];
	my $xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".
		"<!DOCTYPE Query>\n".
		"<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n".
		"<Dataset name = \"$species\" interface = \"default\" >\n".
		"<Attribute name = \"chromosome_name\" />\n".
		"<Attribute name = \"genomic_coding_start\" />\n".
		"<Attribute name = \"genomic_coding_end\" />\n".
		"<Attribute name = \"ensembl_exon_id\" />\n".
		"<Attribute name = \"ensembl_gene_id\" />\n".
		"<Attribute name = \"strand\" />\n".
		"<Attribute name = \"external_gene_id\" />\n".
		"</Dataset>\n".
		"</Query>\n";
	return($xml);
}

sub format_species
{
	# Impossible otherwise!
	use v5.14;
	my ($org,$source) = @_;
	given($source)
	{
		when(/ucsc|refseq/)
		{
			given($org)
			{
				when(/human/) { return("hg19"); }
				when(/mouse/) { return("mm9"); }
				when(/rat/) { return("rn5"); }
				when(/fly/) { return("dm3"); }
				when(/zebrafish/) { return("danRer7"); }
			}
		}
		when(/ensembl/)
		{
			given($org)
			{
				when(/human/) { return("hsapiens_gene_ensembl"); }
				when(/mouse/) { return("mmusculus_gene_ensembl"); }
				when(/rat/) { return("rnorvegicus_gene_ensembl"); }
				when(/fly/) { return("dmelanogaster_gene_ensembl"); }
				when(/zebrafish/) { return("drerio_gene_ensembl"); }
			}
		}
	}
	#if ($_[0] =~ m/human/) { return("hsapiens_gene_ensembl"); }
	#elsif ($_[0] =~ m/mouse/) { return("mmusculus_gene_ensembl"); }
	#elsif ($_[0] =~ m/rat/) { return("rnorvegicus_gene_ensembl"); }
	#elsif ($_[0] =~ m/fly/) { return("dmelanogaster_gene_ensembl"); }
	#elsif ($_[0] =~ m/zebrafish/) { return("drerio_gene_ensembl"); }
}

sub waitbar_init
{
	my $title;
	($_[1]) ? ($title = $_[1]) : ($title = "Counting...");
	my $progress = Term::ProgressBar->new({name => $title,count => $_[0],fh => \*STDOUT,ETA => 'linear'});
	$progress->minor(0) if (!$_[2]);
	return($progress);
}

sub waitbar_update { }

sub count_lines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
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

sub advertise
{
	#use Term::ANSIColor;
	#print color 'bold yellow on_blue';
	disp($scriptname);
	#print color 'bold green';
	disp("Short sequence read counting in genomic regions... Copyright: Panagiotis Moulos (moulos\@fleming.gr)\n");
	#print color 'reset';
}

# Print the program usage
sub program_usage
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";

A perl program to join reads in bed file(s) with their genomic regions
in a genome and count presences of reads (uniquely) in each gene. It
simplifies the operation on genomic intervals of Galaxy and solw by
performing all operations locally and unified. Additionally, it
addresses the problem of partial overlaping between a read and a
genomic regions by 3 possible ways: an overlaping coefficient, a
linear probabilitistic score and an exponential probabilistic score.
Please see program description in the header of the script for further
details and file format descriptions.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file1 [file2, file3, ..., filen] --region file [OPTIONS]

--- Required ---
  --input|i  file(s)  Input file(s). The input bed files can be 3 or 6
  			column bedfiles.
  --region|r file Genomic regions file. The first 3 columns should
  			be of the same structure as in bed files. The 4th column
  			should contain a UNIQUE identifier for each region (e.g.
  			Ensembl IDs). The rest columns can contain additional data
  			about the regions. Instead of a local file, it can be one
  			of the following, which will automatically download and use
  			genomic annotations from the latest version of Ensembl:
  			"human-gene" for Ensembl homo sapiens gene co-ordinates
  			"human-exon" for Ensembl homo sapiens exon co-ordinates
  			"mouse-gene" for Ensembl mus musculus gene co-ordinates
  			"mouse-exon" for Ensembl mus musculus exon co-ordinates
  			"rat-gene" for Ensembl rattus norvegicus gene co-ordinates
  			"rat-exon" for Ensembl rattus norvegicus exon co-ordinates
  			"fly-gene" for Ensembl drosophila melanogaster gene co-ordinates
  			"fly-exon" for Ensembl drosophila melanogaster exon co-ordinates
  			"zebrafish-gene" for Ensembl danio rerio gene co-ordinates
  			"zebrafish-exon" for Ensembl danio rerio exon co-ordinates
--- Optional ---
  --source|u		Use this option to set the online data source in
			the case of selecting one of the prefefined region templates
			with --region. Can be one of "ucsc", "refseq" or "ensembl".
			Default to "ensembl".
  --sort|a		Use this option if you wish to sort the input files
  			first. It is obligatory if the files are not already
  			sorted. If you do not use it and the files are not sorted,
  			results will be incorrect if the program will not crash.
  			If you wish to sort the files manually, be sure to sort
  			firstly by chromosome (1st column) and then by region
  			start. In Linux systems, the command should look like
  			this: 'sort -k1,1 -k2g,2 inputfile > outputfile'. In
  			Windows systems, it would be better to use a spreadsheet
  			program like Excel to perform sorting by columns. Both
  			the input bed files and the region file should be sorted.
  			If the structure of the region file is as described (1st
  			column: chromosome, 2nd column: start 3rd column: end,
  			4th column: unique ID etc.) the sorting command is
  			exactly the same as in common bed files (3 or 6 columns).
  --percent|p		Use this option to provide the overlaping
  			coefficient according to which tags that partially fall
  			outside provided genomic regions will be assigned to
  			genomic regions. Defaults to 0.95 and is the default
  			algorithm for assigning tags that fall partially inside
  			provided genomic regions.
  --lscore|l		Use this option to use the linear probabilistic
  			scoring scheme to address the problem of partial tag
  			overlap with the provided genomic regions. Please see
  			header of the script file for description. Defaults to 0
  			(not use).
  --escore|e		Use this option to use the exponential probabilistic
  			scoring scheme to address the problem of partial tag
  			overlap with the provided genomic regions. Please see
  			header of the script file for description. Defaults to 0
  			(not use).
  --constant|c		Use this option to provide the constant for the
  			exponential scoring scheme (see description of --escore
  			option). Defaults to 3.
  --small|m		Use this option if you wish to take into consideration
  			genomic regions in your genomic regions file which are
  			smaller than the tag length (rare but may occur, especially
  			with customized genomic regions). Defaults to 0 (not use).
  --split|t		Use this option if you wish to further split your genomic
  			regions in smaller areas. In this case, tags will be counted
  			per area and the distribution of counts per area will be
  			returned as a column in the output file. After --split,
  			you can provide the length of sub-areas. The default is
  			1000 and splits the regions per 1kb, unless the regions are
  			smaller. In this case, the area will consist of only one
  			sub-area of length equal to that of the area. Note also
  			that when using this option, tags that are found inside
  			sub-areas are not assigned to those sub-areas based on
  			scoring schemes (options --percent, --lscore and --escore)
  			but tags are assigned based on the location of their center.
  --stats|z		Use this option to also return basic statistics of
			counts in the windows used returned by using --split.
  --ncore|n		If the machine has multicore processor(s) and the
			package Parallel::ForkManager is installed, you can use
			parallel processing. Default is 1 and can go up to 12.
  --output|o		A file to write the output to. If "auto", then it
			generates an automatic filename in the folder where the
			input files are. If not provided, output is written to STDOUT.
  --waitbar|w		Use this option if you wish to display a progress
			bar using Term::ProgressBar module. Caution, as for large
			files it explodes the processing time.
  --silent|s		Use this option if you want to turn informative
  			messages off.
  --help|h		Display this help text.

The program returns bed file(s) that contain the provided genomic regions
with any other information together with read counts per region.

Package dependecies:
	Tie::IxHash::Easy (mandatory)
	Term::ProgressBar (optional)
	File::Sort (optional)

END
	print $usagetext;
	exit;
}

#sub get_file_parts()
#{
	#my @in = @_;
	#my (@base,@dir,@ext);
	#for (my $i=0; $i<@in; $i++)
	#{
		#($base[$i],$dir[$i],$ext[$i]) = fileparse($in[$i],'\.[^.]*');
	#}
	#return(\@base,\@dir,\@ext);
#}

#given($regionfile)
#{
	#when(/human-gene|mouse-gene|rat-gene|fly-gene|zebrafish-gene/i)
	#{
		#$regionfile = &fetch_ensembl_genes($regionfile);
		#$regionfile = &sort_ensembl_genes($regionfile);
	#}
	#when (/human-exon|mouse-exon|rat-exon|fly-exon|zebrafish-exon/i)
	#{
		#$regionfile = &fetch_ensembl_exons($regionfile);
		#$regionfile = &sort_ensembl_exons($regionfile);
	#}
	#when(/human-(5|3)utr|mouse-(5|3)utr|rat-(5|3)utr|fly-(5|3)utr|zebrafish-(5|3)utr/i)
	#{
		#($regionfile =~ m/5utr/i) ? ($regionfile = &fetch_ensembl_utr($regionfile,5)) :
			#($regionfile = &fetch_ensembl_utr($regionfile,3));
		#$regionfile = &sort_ensembl_exons($regionfile);
	#}
	#when(/human-cds|mouse-cds|rat-cds|fly-cds|zebrafish-cds/i)
	#{
		#$regionfile = &fetch_ensembl_cds($regionfile);
		#$regionfile = &sort_ensembl_exons($regionfile);
	#}
	#default
	#{
		#($regionfile,@infile) = &sort_inputs($regionfile,@infile) if ($sort);
	#}
#}
