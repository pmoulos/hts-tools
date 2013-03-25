=head1 NAME

HTS::Tools::Intersect - Intersect small genomic regions files, usually results from peak callers.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

A pure perl module to calculate intersection between bed files (can contain as many additional columns 
as desired, as long as the first 3 are chromosome, start, end). It requires as input only the bed files 
and returns the overlap regions between the first and the second file. However, if fileA and fileB are 
intersected, the final file containing name OVERLAP...fileA output will  contain regions from fileA 
that overlap with regions of fileB and the final file containing name OVERLAP...B will contain those 
regions from fileB. The _ONLY files contain regions only in fileA or only in fileB, when performing  the 
intersection of fileA and fileB. The actions performed are similar to those of the UCSC Table Browser, 
and BEDTools intersecting functions, but quite slower because of the slower internal structures and search, 
algorithm. However, this module can be safely used when the user wants to maintain several precalculated 
statistics, e.g. the average number of reads per specified window size under a peak region, something
which is not currently possible using BEDTools. Intersection can be done at specific overlap percentages 
or any overlap (like Table Browser, default, 1bp). In addition, one can specify an array of percentages
so that various overlaps and overlap percentages can be calculated in a batch mode to answer questions
like 'is there an overlap saturation between my two peak sets?' Extension from a point (e.g. the peak mode) 
towards both directions is also possible as long as a column with this point (the peak 'summit') is given 
in both files and the appropriate option is used. The user has the option to retrieve only certain of 
the four available output file types.

    use HTS::Tools::Intersect;
	my %params1 = (
		'inputA' => 'normal_nfkb_peaks.txt',
		'inputB' => 'disease_nfkb_peaks.txt',
		'any' => 1
	)
    my $intersecter = HTS::Tools::Count->new(\%params1);
    $intersecter->run;
    
    my %params2 = (
		'inputA' => 'normal_nfkb_peaks.txt',
		'inputB' => 'disease_nfkb_peaks.txt',
		'percent' => 0.5,
		'extend' => [4,100],
		'both' => 1,
		'output' => ['overlapA','onlyA','onlyB'],
		'waitbar' => 1
	)
    $intersecter->set_params(\%params2);
    $intersecter->run;
    
The acceptable parameters are as follows:

=over 4

=item I<inputA> B<(required)>

First input BED-like file

=item I<inputB> B<(required)>

Second input BED-like file

=item I<sort> B<(optional)>

Use this option to sort the input files first. NECESSARY if they are not sorted beforehand as the script 
uses a binary search algorithm which requires sorted input.

=item I<percent> B<(optional)>

The overlap percentage between the two files. It should be a value between 0 and 100. Alternatively, the
overlap percentage can be a series of values. In this case the program will run as a batch procedure, 
calculating the number of overlapping regions for the determined output types for each given percentage, 
creating and writing as output the distribution of the number of overlapping regions. In this case, no 
other output files are produced. Instead of a series of numbers, it can be an element in the format a:b 
where a is the starting percentage and b the ending. A list between a and b sith interval 1 will be 
autogeneratied. This options does not work with --any option. Examples are I<percent> => 50, I<percent> => 
10 20 30 40 50, I<percent> => 1:100.

=item I<any> B<(optional)> 

Use this switch for returning regions with any (at least 1bp) overlap between fileA and fileB. When used, 
I<percent> is ignored. This is the default (1).

=item I<extend> B<(optional)>

Use this option to supply the program the upstream and downstream extension length. It should be two 
numbers, the first the number of bps to extend upstream of peak center and the second to extend downstream. 
Use this option in combination with I<mode> option. For example I<extend> => [100,200].

=item I<autoextend> B<(optional)>

Use this switch to calculate half of the median peak length which will then be used for extending from 
peak mode. Use this switch in combination with I<mode> option to provide the column with the peak mode
in the input files

=item I<mode> B<(optional)>

Use this option to supply the program the column in both fileA and fileB that contains the peak mode 
or any point from which extension to left and right will start. E.g. I<mode> => 4. This parameter has 
to be provided when I<extend> => 1 and/or I<autoextend> => 1 are used. Also when overpairs is chosen 
as an output format.

=item I<both> B<(optional)>

Normally, when overlaping fileA with fileB at a percentage level, if the tag of the fileA overlaps with 
tag(s) from fileB, this percentage will be calculated based on tag from fileA. That is, if the required 
percentage is p and the overlap is d, then if d < p*length(regionA) then there is no overlap. However, if
region B is smaller than regionA it could d > p*length(regionB), so there is p percent overlap between 
region A and region B when using the length of region B. This option allows for such comparisons by checking
both cases and returning a positive overlap if one of the above two is true. This swicth is useless when 
used with I<extend> option and ignored.

=item I<exact> B<(optional)>

Use this switch to calculate overlaps based on exact genomic locations rather than when a region is 
"totally" included inside the other, when performing overlaps based on percentages. For example if 
overlaping fileA with fileB, if a region in fileA starts at 300 and ends at 800 and the overlaping region
in fileB starts at 500 and ends at 700, if using the I<exact> swicth, a 50% percent overlap will not 
return this region as positive because (700-500)<0.5*(800-300) even if region B is totally included in 
region A. This switch forces the program to calculate percent overlaps exactly as the UCSC Table Browser.
Can also be used in combination with the I<both> switch.

=item I<pass> B<(optional)>

Use this option to supply the program the number of times that the dynamic binary search algorithm will 
search regions from the fileB for each region of the fileA. One pass returns at maximum one hit because 
if the algorithm finds a hit from fileB, it will exit. The number of passes determines how many times 
the algorithm will search for overlapping regions of fileA and fileB according to specified overlapping 
crteria. Use a larger number of passes for when regions from fileA are likely to overlap many regions 
from fileB (e.g. when fileA has a large peak which could correspond to more than one peak in fileB). 
It defaults to 3.

=item I<gap> B<(optional)>

Use this option to retrieve distances between non-overlaping (according to specified criteria) regions 
that is not larger than the number specified with the I<gap> option. Can be useful with onlyA and nonpairs
output options. The default gap is set to 0 so the option is not used. If you wish to use that option 
you should provide a positive integer, for example I<gap> => 10000.

=item I<output> B<(optional)>

Use this option to determine which intersection output filetypes you wish to retrieve.	Possible choices 
are: "overlapA" for retrieving regions from fileA that overlap with fileB. "overlapB" for retrieving 
those regions that overlap was found with fileA regions when using  fileA as first input (ATTENTION! 
NOT REGIONS FROM fileB THAT OVERLAP WITH fileA). "onlyA" for retrieving regions from fileA that DO NOT 
overlap with regions in fileB. "onlyB" for retrieving those regions that overlap was NOT found with fileA 
when using fileA as first input (ATTENTION! SAME BEHAVIOUR AS "overlapB" choice).

=item I<header> B<(optional)>

Use this option if you have a header line in your input files. It will also be written in your output
files. Should be the same! Defaults to no header.

=item I<waitbar> B<(optional)>

Use this option if you wish to display a simple progress bar while running the procedures. For small 
files it is probably useless as the program finishes very quickly.

=item <silent>

Use this option if you want to turn informative messages off.

=back  			

=head1	OUTPUT

The main output of the module is up to four files in BED format containing also any additional data 
columns.

=head1 SUBROUTINES/METHODS

=cut

package HTS::Tools::Intersect;

use v5.10;
use strict;
use warnings FATAL => 'all';

our $MODNAME = "HTS::Tools::Intersect";
our $VERSION = '0.01';
our $AUTHOR = "Panagiotis Moulos";
our $EMAIL = "moulos\@fleming.gr";
our $DESC = "Advanced intersection of genomic region sets.";

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

use Carp;
use File::Basename;
use File::Temp;
use File::Spec;

use lib '/media/HD4/Fleming/dev/HTS-Tools/lib';
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

The HTS::Tools::Intersect object constructor. It accepts a set of parameters that are required to run
the counter and get the output.

	my $intersecter = HTS::Tools::Intersect->new({'inputA' => 'myfileA.bed','inputB' => 'myfileB.bed'});

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
	$helper->advertise($MODNAME,$VERSION,$AUTHOR,$EMAIL,$DESC);

	# Validate the input parameters
	my $checker = HTS::Tools::Paramcheck->new();
	$checker->set("tool","intersect");
	$checker->set("params",$params);
	$params = $checker->validate;

	# After validating, bless and initialize
	bless($self,$class);
	$self->init($params);
	return($self);
}

=head2 init

HTS::Tools::Intersect object initialization method. NEVER use this directly, use new instead.

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

The HTS::Tools::Intersect run subroutine. It runs the interssecter with the given parameters in the 
constructor.

	$intersecter->run;
	
=cut

sub run
{
	my $self = shift @_;
	my $ofileA = $fileA;
	my $ofileB = $fileB; # Keep original filenames (verbose purposes)

	if ($self->get($sort))
	{
		($fileA,$fileB) = $self->sort_input($self->get("inputA"),$self->get("inputB"));
	}
	
	# Copy some memory-less variables to avoid rewriting the whole thing...
	my $any = $self->get("any");
	my $mode = $self->get("mode");
	my $agap = $self->get("agap");
	my @percent = @{$self->get("percent")};
	my @extend = @{$self->get("extend")};
	my @out = @{$self->get("output")};
	
	# Bavard a little... I like information...
	$helper->$helper->disp("Type of overlap: any overlap, percentage overlap ignored...") if ($any && @percent);
	$helper->$helper->disp("Type of overlap: any overlap") if ($any && !@percent);
	$helper->$helper->disp("Type of overlap: $percent[0]% overlap") if ($percent[0] && !$percent[1] && !$any);
	$helper->$helper->disp("Extending region modes upstream $extend[0] and downstream $extend[1] bps") if (@extend && !$autoxtend);
	$helper->$helper->disp("Region modes extension on each side will be auto-calculated...") if ($autoxtend);
	if ($percent[1])
	{
		$helper->$helper->disp("Multiple overlap percentages... Running in batch mode to determine overlapping distributions...");
		$helper->$helper->disp("No overlapping output files will be produced..."); 
	}

	# Save some more memory...
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
		$helper->$helper->disp("overpairs and nonpairs output formats can only be given with a peak summit column! Ignoring...");
		$op = 0;
		$np = 0;
	}
	$helper->$helper->disp("Retrieving distances between non-overlapping regions if distance <= $agap bps") if (($op || $np) && $agap);
	$helper->$helper->disp(" ");

	# Get number of lines for fileA (mostly for waitbar, I know, inefficient but I like waitbars ;-)
	my $linesA = $helper->count_lines($fileA) if ($self->get("waitbar"));

	# Suck in fileB
	my $chr;
	my @rest;
	my (%hashB,%onlyB);
	open(INB,$fileB) or croak "\nThe file $fileB does not exist!\n";
	$helper->$helper->disp("Reading file $ofileB...");
	my $headerB = <INB> if ($self->get("header"));
	while (my $line = <INB>)
	{
		next if ($line =~/^chrM/);
		next if ($line =~/rand|hap|chrU/);
		$line =~ s/\r|\n$//g;
		($chr,@rest) = split(/\t/,$line);
		push(@{$hashB{$chr}},join("\t",@rest));
		$onlyB{$chr}{join("\t",@rest)}++ if ($oB);
	}
	close(INB);
	
	##########################################################################################

	# Do we run in batch mode in order to determine distributions?
	if (!$any && $percent[1])
	{	
		my $i;
		my @distributions;
		my ($cchr,@crest);
		my ($bsr,$ci,$bsf,$n,@currvals);
		my $linB = $helper->count_lines($fileB);
		$linB-- if ($header);

		if (@extend || $autoxtend)
		{	
			if ($autoxtend) # We have to suck in the whole peak file once...
			{
				open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";;
				$helper->$helper->disp("Reading file $ofileA and processing overlaps...");
				my $headerA = <INA> if ($self->get("header"));
				my (@lines,@medmodes);
				$helper->waitbar_init if ($self->get("waitbar"));
				while (my $line = <INA>)
				{
					$helper->waitbar_update($.,$linesA) if ($self->get("waitbar"));
					next if ($line =~/^chrM/);
					next if ($line =~/rand|hap|chrU/);
					$line =~ s/\r|\n$//g;
					push(@lines,$line);
					@crest = split(/\t/,$line);
					push(@medmodes,$crest[2] - $crest[1]); # BED format
				}
				close(INA);
				my $med = $helper->median((@medmodes,$self->get_lengths(\%hashB)));
				@extend = (int($med/2),int($med/2));
				$helper->$helper->disp("Median region length is $med bps. Extending each region mode $extend[1] bps on each side...\n");
				for ($i=0; $i<@percent; $i++)
				{	
					$helper->disp("Overlap percentage: $percent[$i]");
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
					$helper->waitbar_init if ($self->get("waitbar"));
					foreach my $l (@lines)
					{
						$countopic++;
						$helper->waitbar_update($countopic,$linesA-1) if ($self->get("waitbar"));
						($cchr,@crest) = split(/\t/,$l);
						@currvals = @{$hashB{$cchr}} if ($hashB{$cchr}); 
						$n = 0;
						$bsf = 0;
						while ($n < $npass) 
						{
							($ci,$bsr) = $self->bin_search_percent_center($crest[$mode],$mode,@extend,$percent[$i]/100,@currvals);
							if ($bsr) # Found in overlap, put into overlap hash of both files
							{
								$overlapA{$cchr}{join("\t",@crest)}++;
								$overlapB{$cchr}{$bsr}++ if ($ovB);
								delete $conlyB{$cchr}{$bsr} if ($oB);
								if ($op)
								{
									my $cl = join("\t",@crest);
									my @ds = $self->dists_center($cl,$bsr,$mode,@extend);
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
					
					push(@{$distributions[$i]},$helper->count_hoh(\%overlapA)) if ($ovA);
					push(@{$distributions[$i]},$helper->count_hoh(\%overlapB)) if ($ovB);
					push(@{$distributions[$i]},$helper->count_hoh(\%onlyA)) if ($oA);
					push(@{$distributions[$i]},$helper->count_hoh(\%conlyB)) if ($oB);
					push(@{$distributions[$i]},$helper->mean(@overpairs)) if ($op);
					push(@{$distributions[$i]},$helper->median(@overpairs)) if ($op);
				}
			}
			else
			{
				for ($i=0; $i<@percent; $i++)
				{
					
					$helper->$helper->disp("Overlap percentage: $percent[$i]");
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
					open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";;
					$helper->$helper->disp("Reading file $ofileA and processing overlaps...");
					$helper->waitbar_init if ($self->get("waitbar"));
					my $headerA = <INA> if ($self->get("header"));
					while (my $line = <INA>)
					{
						$helper->waitbar_update($.,$linesA) if ($self->get("waitbar"));	
						next if ($line =~/^chrM/);
						next if ($line =~/rand|hap|chrU/);
						$line =~ s/\r|\n$//g;
						($cchr,@crest) = split(/\t/,$line);
						@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
						$n = 0;
						$bsf = 0;
						while ($n < $npass) 
						{
							($ci,$bsr) = $self->bin_search_percent_center($crest[$mode],$mode,@extend,$percent[$i]/100,@currvals);
							if ($bsr) # Found in overlap, put into overlap hash of both files
							{
								$overlapA{$cchr}{join("\t",@crest)}++;
								$overlapB{$cchr}{$bsr}++ if ($ovB);
								delete $conlyB{$cchr}{$bsr} if ($oB);
								if ($op)
								{
									my $cl = join("\t",@crest);
									my @ds = $self->dists_center($cl,$bsr,$mode,@extend);
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
					
					push(@{$distributions[$i]},$helper->count_hoh(\%overlapA)) if ($ovA);
					push(@{$distributions[$i]},$helper->count_hoh(\%overlapB)) if ($ovB);
					push(@{$distributions[$i]},$helper->count_hoh(\%onlyA)) if ($oA);
					push(@{$distributions[$i]},$helper->count_hoh(\%conlyB)) if ($oB);
					push(@{$distributions[$i]},$helper->mean(@overpairs)) if ($op);
					push(@{$distributions[$i]},$helper->median(@overpairs)) if ($op);
				}
			}
		}
		else
		{
			for ($i=0; $i<@percent; $i++)
			{
				$helper->$helper->disp("Overlap percentage: $percent[$i]");
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
				
				open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";;
				$helper->$helper->disp("Reading file $ofileA and processing overlaps...");
				$helper->waitbar_init if ($self->get("waitbar"));
				my $headerA = <INA> if ($self->get("header"));
				while (my $line = <INA>)
				{
					$helper->waitbar_update($.,$linesA) if ($self->get("waitbar"));	
					next if ($line =~/^chrM/);
					next if ($line =~/rand|hap|chrU/);
					$line =~ s/\r|\n$//g;
					($cchr,@crest) = split(/\t/,$line);
					@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
					$n = 0;
					$bsf = 0;
					while ($n < $npass) 
					{
						if ($exact)
						{
							($both) ? (($ci,$bsr) = $self->bin_search_percent_exact_both($crest[0],$crest[1],$percent[$i]/100,@currvals)) :
							(($ci,$bsr) = $self->bin_search_percent_exact($crest[0],$crest[1],$percent[$i]/100,@currvals));
						}
						else
						{
							($both) ? (($ci,$bsr) = $self->bin_search_percent_both($crest[0],$crest[1],$percent[$i]/100,@currvals)) :
							(($ci,$bsr) = $self->bin_search_percent($crest[0],$crest[1],$percent[$i]/100,@currvals));
						}
						
						if ($bsr) # Found in overlap, put into overlap hash of both files
						{
							$overlapA{$cchr}{join("\t",@crest)}++;
							$overlapB{$cchr}{$bsr}++ if ($ovB);
							delete $conlyB{$cchr}{$bsr} if ($oB);
							if ($op)
							{
								my $cl = join("\t",@crest);
								my @ds = $self->dists_every($cl,$bsr,$mode);
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
				
				push(@{$distributions[$i]},$helper->count_hoh(\%overlapA)) if ($ovA);
				push(@{$distributions[$i]},$helper->count_hoh(\%overlapB)) if ($ovB);
				push(@{$distributions[$i]},$helper->count_hoh(\%onlyA)) if ($oA);
				push(@{$distributions[$i]},$helper->count_hoh(\%conlyB)) if ($oB);
				push(@{$distributions[$i]},$helper->mean(@overpairs)) if ($op);
				push(@{$distributions[$i]},$helper->median(@overpairs)) if ($op);
			}
		}
		
		my $time = $helper->now("machine");
		my $fdn = "distrib_".$time;
		$helper->$helper->disp("Writing output in $fdn.txt");
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
		$helper->cleanup;
		$helper->$helper->disp("Finished!\n\n");
		exit;
	}

	##########################################################################################
	
	# Do job with lines of fileA
	my ($cchr,@crest);
	my ($bsr,$ci,$bsf,$n,@currvals);
	my (%overlapA,%overlapB,%onlyA);
	my (@overpairs,@nonpairs);
		
	open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";;
	$helper->$helper->disp("Reading file $ofileA and processing overlaps...");
	my $headerA = <INA> if ($self->get("header"));
	$helper->waitbar_init if ($self->get("waitbar"));
	if (@extend || $autoxtend)
	{
		if ($autoxtend) # We have to suck in the whole peak file once...
		{
			my (@lines,@medmodes);
			while (my $line = <INA>)
			{
				$helper->waitbar_update($.,$linesA) if ($self->get("waitbar"));
				next if ($line =~/^chrM/);
				next if ($line =~/rand|hap|chrU/);
				$line =~ s/\r|\n$//g;
				push(@lines,$line);
				@crest = split(/\t/,$line);
				push(@medmodes,$crest[2] - $crest[1]); # BED format
			}
			my $med = $helper->median((@medmodes,$self->get_lengths(\%hashB)));
			@extend = (int($med/2),int($med/2));
			$helper->disp("Median region length is $med bps. Extending each region mode $extend[1] bps on each side...\n");
			my $countopic = 0;
			$helper->waitbar_init if ($self->get("waitbar"));
			foreach my $l (@lines)
			{
				$countopic++;
				$helper->waitbar_update($countopic,$linesA-1) if ($self->get("waitbar"));
				($cchr,@crest) = split(/\t/,$l);
				
				# Perform binary search according to user choices...
				# A specific chromosome might not exist in one of the two files...
				@currvals = @{$hashB{$cchr}} if ($hashB{$cchr}); 
				$n = 0;
				$bsf = 0;
				while ($n < $npass) 
				{
					($any) ? (($ci,$bsr) = $self->bin_search_any_center($crest[$mode],$mode,@extend,@currvals)) :
					(($ci,$bsr) = $self->bin_search_percent_center($crest[$mode],$mode,@extend,$percent[0]/100,@currvals));
					if ($bsr) # Found in overlap, put into overlap hash of both files
					{
						$overlapA{$cchr}{join("\t",@crest)}++;
						$overlapB{$cchr}{$bsr}++ if ($ovB);
						# Remove from fileB hash so as what remains will be the only B.
						delete $onlyB{$cchr}{$bsr} if ($oB);
						if ($op)
						{
							my $cl = join("\t",@crest);
							my @ds = $self->dists_center($cl,$bsr,$mode,@extend);
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
					if ($np)
					{
						my $cl = join("\t",@crest);
						my @ds = $self->dists_center($cl,$bsr,$mode,@extend);
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
				$helper->waitbar_update($.,$linesA) if ($self->get("waitbar"));	
				next if ($line =~/^chrM/);
				next if ($line =~/rand|hap|chrU/);
				$line =~ s/\r|\n$//g;
				($cchr,@crest) = split(/\t/,$line);

				# Perform binary search according to user choices...
				# A specific chromosome might not exist in one of the two files...
				@currvals = @{$hashB{$cchr}} if ($hashB{$cchr});
				$n = 0;
				$bsf = 0;
				while ($n < $npass) 
				{
					($any) ? (($ci,$bsr) = $self->bin_search_any_center($crest[$mode],$mode,@extend,@currvals)) :
					(($ci,$bsr) = $self->bin_search_percent_center($crest[$mode],$mode,@extend,$percent[0]/100,@currvals));
					if ($bsr) # Found in overlap, put into overlap hash of both files
					{
						$overlapA{$cchr}{join("\t",@crest)}++;
						$overlapB{$cchr}{$bsr}++ if ($ovB);
						# Remove from fileB hash so as what remains will be the only B.
						delete $onlyB{$cchr}{$bsr} if ($oB);
						if ($op)
						{
							my $cl = join("\t",@crest);
							my @ds = $self->dists_center($cl,$bsr,$mode,@extend);
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
					if ($np)
					{
						my $cl = join("\t",@crest);
						my @ds = $self->dists_center($cl,$bsr,$mode,@extend);
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
			$helper->waitbar_update($.,$linesA) if ($self->get("waitbar"));	
			next if ($line =~/^chrM/);
			next if ($line =~/rand|hap|chrU/);
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
					($ci,$bsr) = $self->bin_search_any($crest[0],$crest[1],@currvals);
				}
				else
				{
					if ($exact)
					{
						($both) ? (($ci,$bsr) = $self->bin_search_percent_exact_both($crest[0],$crest[1],$percent[0]/100,@currvals)) :
						(($ci,$bsr) = $self->bin_search_percent_exact($crest[0],$crest[1],$percent[0]/100,@currvals));
					}
					else
					{
						($both) ? (($ci,$bsr) = $self->bin_search_percent_both($crest[0],$crest[1],$percent[0]/100,@currvals)) :
						(($ci,$bsr) = $self->bin_search_percent($crest[0],$crest[1],$percent[0]/100,@currvals));
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
						my @ds = $self->dists_every($cl,$bsr,$mode);
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
					my @ds = $self->dists_every($cl,$bsr,$mode);
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
	$helper->disp(" ") if (!$waitbar);
	if (!$silent)
	{
		my ($lA,$lB);
		($waitbar) ? ($lA = $linesA) : ($lA = $helper->count_lines($fileA));
		$lB = $helper->count_lines($fileB);
		if ($header)
		{
			$lA--;	
			$lB--;
		}
		if ($ovA)
		{
			my $covA = $helper->count_hoh(\%overlapA);
			$helper->disp("$covA out of $lA regions from $ofileA overlap with regions from $ofileB");
		}
		if ($ovB)
		{
			my $covB = $helper->count_hoh(\%overlapB);
			$helper->disp("$covB out of $lB regions from $ofileB overlap with regions from $ofileA");
		}
		if ($oA)
		{
			my $coA = $helper->count_hoh(\%onlyA);
			$helper->disp("$coA out of $lA regions exist only in $ofileA");
		}
		if ($oB)
		{
			my $coB = $helper->count_hoh(\%onlyB);
			$helper->disp("$coB out of $lB regions exist only in $ofileB");
		}
		$helper->disp(" ");
	} 

	# Good... print outputs...
	$helper->disp("Writing output...");
	my ($bA,$bB);
	$bA = fileparse($ofileA,'\..*?') if ($ovA || $oA);
	$bB = fileparse($ofileB,'\..*?') if ($ovB || $oB);
	# A and B overlap elements with elements of file A
	$self->print_complex_output($ofileA,$ofileB,$headerA,"OVERLAP_FROM_$bA",\%overlapA) if ($ovA); 
	# A and B overlap elements with elements of file B
	$self->print_complex_output($ofileA,$ofileB,$headerB,"OVERLAP_FROM_$bB",\%overlapB) if ($ovB);
	# A only elements
	$self->print_complex_output($ofileA,$ofileB,$headerA,"ONLY_$bA",\%onlyA) if ($oA);
	# B only elements
	$self->print_complex_output($ofileA,$ofileB,$headerB,"ONLY_$bB",\%onlyB) if ($oB);
	# Distances of overlaping
	$self->print_array($ofileA,$ofileB,"OVERDIST",@overpairs) if ($op);
	# Distances of non-overlaping
	$self->print_array($ofileA,$ofileB,"NONDIST",@nonpairs) if ($np);

	# Remove garbage
	$helper->cleanup;
	$helper->disp("Finished!\n\n") if (!$multi);
}

=head2 bin_search_any

Binary search algorithm for any overlap between genomic regions. Internal use.

	$intersecter->bin_search_any($start,$end,@candidate_areas);
	
=cut

sub bin_search_any
{
	my ($self,$start,$end,@areas) = @_;
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

=head2 bin_search_percent

Binary search algorithm for percent overlap between genomic regions. Internal use.

	$intersecter->bin_search_percent($start,$end,$percentage,@candidate_areas);
	
=cut

sub bin_search_percent
{
	my ($self,$start,$end,$p,@areas) = @_;
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

=head2 bin_search_any_center

Binary search algorithm for any overlap between genomic regions using their centers. Internal use.

	$intersecter->bin_search_any_center($mode,$position,$downstream,$upstream,@candidate_areas);
	
=cut

sub bin_search_any_center
{
	my ($self,$mode,$mpos,$dxval,$uxval,@areas) = @_;
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

=head2 bin_search_any_center

Binary search algorithm for percentage overlap between genomic regions using their centers. Internal use.

	$intersecter->bin_search_percent_center($mode,$position,$downstream,$upstream,$percentage,@candidate_areas);
	
=cut

sub bin_search_percent_center
{
	my ($self,$mode,$mpos,$dxval,$uxval,$p,@areas) = @_;
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

=head2 bin_search_percent_both

Binary search algorithm for percentage overlap between genomic regions for the "both" case. Internal use.

	$intersecter->bin_search_any_center($start,$end,$percentage,@candidate_areas);
	
=cut

sub bin_search_percent_both
{
	my ($self,$start,$end,$p,@areas) = @_;
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

=head2 bin_search_percent_exact

Binary search algorithm for percentgae overlap between genomic regions for the "exact" case. Internal use.

	$intersecter->bin_search_percent_exact($start,$end,$percentage,@candidate_areas);
	
=cut

sub bin_search_percent_exact
{
	my ($self,$start,$end,$p,@areas) = @_;
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

=head2 bin_search_percent_both

Binary search algorithm for percentage overlap between genomic regions for the "exact" and "both" case.
Internal use.

	$intersecter->bin_search_percent_both($start,$end,$percentage,@candidate_areas);
	
=cut

sub bin_search_percent_exact_both
{
	my ($self,$start,$end,$p,@areas) = @_;
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

=head2 dists_every

Distance calculation subroutine. Internal use.

	$intersecter->dists_every($A,$B,$ei);
	
=cut

sub dists_every
{
	my ($self,$A,$B,$ei) = @_;
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

=head2 dists_every

Distance calculation subroutine using centers. Internal use.

	$intersecter->dists_center($A,$B,$ei,$up,$down);
	
=cut

sub dists_center
{
	my ($self,$A,$B,$ei,$exu,$exd) = @_;
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

=head2 get_lengths

Get lengths of input genomic regions. Internal use.

	$intersecter->get_lengths(%input_hash);

=cut

sub get_lengths
{
	my ($self,$inhash) = @_;
	my @lens;
	my ($c,$l,@t);
	foreach $c (keys(%$inhash))
	{
		foreach $l (@{$inhash->{$c}})
		{
			@t = split(/\t/,$l);
			push(@lens,$t[1] - $t[0]);
		}
	}
	return(@lens);
}

=head2 print_complex_output

Module specific output printing function. Internal use.

	$intersecter->print_complex_output($A,$B,$output_type,$header,$the_hash);

=cut

sub print_complex_output
{
	my ($self,$infileA,$infileB,$he,$otype,$inhash) = @_;
	my ($outchr,$ind,$outhash,@k);
	my $outfilename = $self->create_output_file($infileA,$infileB,$otype);
	my @chrs = keys(%$inhash);
	open(OUTPUT,">$outfilename");
	print OUTPUT "$he" if ($he);
	foreach $outchr (sort(@chrs))
	{
		$outhash = $inhash->{$outchr};
		@k = sort outsort keys(%$outhash);
		for ($ind=0; $ind<@k; $ind++)
		{
			print OUTPUT "$outchr\t$k[$ind]\n";
		}
	}
	close(OUTPUT);
}

=head2 print_array

Module specific output printing function. Internal use.

	$intersecter->print_array($A,$B,$output_type,@array);

=cut

sub print_array
{
	my ($self,$infileA,$infileB,$otype,@inarr) = @_;
	my $outfilename = $self->create_output_file($infileA,$infileB,$otype);
	open(OUTPUT,">$outfilename");
	print OUTPUT join("\n",@inarr);
	close(OUTPUT);
}

=head2 create_output_file

Create the name of the output file according to output type. Internal use.

	$intersecter->create_output_file($A,$B,$output_type);

=cut

sub create_output_file
{
	my ($self,$inA,$inB,$type) = @_;
	my ($baseA,$dirA,$extA) = fileparse($inA,'\..*?');
	my $baseB = fileparse($inB,'\..*?');
	($multi) ? (return(File::Spec->catfile($dirA,$baseA.$baseB))) :
	(return(File::Spec->catfile($dirA,$baseA."_".$baseB."_".$type.$extA)));
}

=head2 sort_inputs

Input file sorting function. Internal use.

=cut

sub sort_inputs
{
	my ($self,$fileA,$fileB) = @_;
	my $tmpdir = $self->get("tmpdir");
	my $tmpfileA = File::Spec->catfile($tmpdir,"tempA.in$$");
	my $tmpfileB = File::Spec->catfile($tmpdir,"tempB.in$$");
	
	# Simulate try and catch
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		if ($header)
		{
			$helper->disp("Sorting file $fileA...");
            my $sortcmdA = "awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $fileA > $tmpfileA";
            `$sortcmdA`;
            $fileA = $tmpfileA;
            $helper->disp("Sorting file $fileB...");
            my $sortcmdB = "awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $fileB > $tmpfileB";
            `$sortcmdB`;
            $fileB = $tmpfileB;
		}
		else
		{	
			$helper->disp("Sorting file $fileA...");
			`sort -k1,1 -k2g,2 $fileA > $tmpfileA `;
			$fileA = $tmpfileA;
			$helper->disp("Sorting file $fileB...");
			`sort -k1,1 -k2g,2 $fileB > $tmpfileB `;
			$fileB = $tmpfileB;
		}
	}
	else # We are in Windows... package required
	{
		my $status = $helper->try_module("File::Sort");
		
		# Will die if module does not exist
		if ($header) # Cannot find a solution for direct sorting in Windows... sorry :-(
		{
			my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
					   "messing up data. Please sort files outside $scriptname first (e.g. using\n".
					   "Excel or something similar.";
			croak "\n$dmsg\n\n";
		}
		eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
		$helper->disp("Sorting region file $fileA...");
		sort_file(
		{
			I => $fileA,
			o => File::Spec->catfile($tmpdir,"tempA.tmp"),
			k => ['1,1','2n,2'],
			t => "\t"
		});
		$fileA = File::Spec->catfile($tmpdir,"tempA.tmp");
		$helper->disp("Sorting region file $fileB...");
		sort_file(
		{
			I => $fileB,
			o => File::Spect->catfile($tmpdir,"tempB.tmp"),
			k => ['1,1','2n,2'],
			t => "\t"
		});
		$fileB = File::Spec->catfile($tmpdir,"tempB.tmp");
	}
	
	return($fileA,$fileB)
}

=head2 change_params

Massively change the parameters of an HTS::Tools::Intersect object.

	$intersecter->change_params({'input' => 'another_file','region' => 'mouse-exon'})
	$intersecter->run;
	
=cut

sub change_params
{
	my ($self,$params) = @_;
	
	# Validate the new parameters 
	my $checker = HTS::Tools::Paramcheck->new();
	$checker->set("tool","intersect");
	$checker->set("params",$params);
	$params = $checker->validate;
	
	# If validator does not complain, change the parameters
	while (my ($name,$value) = each(%$params))
	{
		$self->set($name,$value);
	}
	return($self);
}

=head2 get

HTS::Tools::Intersect object getter

	my $param_value = $count->get('param_name')
=cut

sub get
{
	my ($self,$name) = @_;
	return($self->{$name});
}

=head2 set

HTS::Tools::Intersect object setter

	$intersecter->set('param_name','param_value');
	
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

    perldoc HTS::Tools::Intersect


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

1; # End of HTS::Tools::Intersect
