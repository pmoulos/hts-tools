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
and BEDTools intersecting functions, at similar speed and more convenience! This module can be safely
used when the user wants to maintain several precalculated statistics, e.g. the average number of reads
per specified window size under a peak region, something which is not currently possible using BEDTools.
Intersection can be done at specific overlap percentages or any overlap (like Table Browser, default, 1bp).
In addition, one can specify an array of percentages so that various overlaps and overlap percentages can
be calculated in a batch mode to answer questions like 'is there an overlap saturation between my two peak
sets?' Extension from a point (e.g. the peak mode) towards both directions is also possible as long as a
column with this point (the peak 'summit') is given in both files and the appropriate option is used.
The user has the option to retrieve only certain of the four available output file types. The user has the
option to retrieve only certain of the four available output file types. Stats are also displayed and the
distances between overlapping and limited non-overlapping regions can be returned.

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
    $intersecter->change_params(\%params2);
    $intersecter->run;
    
The acceptable parameters are as follows:

=over 4

=item I<inputA> B<(required)>

First input BED-like file

=item I<inputB> B<(required)>

Second input BED-like file

=item I<sort> B<(optional)>

Use this option to sort the input files first. This is not necessary but the module runs slighlty faster
if the files are sorted beforehand.

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

=item I<reportonce> B<(optional)>

Use this option to report only once regions from one file that may overlap multiple times from regions in
the other file (this happens for example with BEDTools). Such situations may arise when for example there
is a broad peak in fileA which is split in two peaks in fileB.

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
when using fileA as first input (ATTENTION! SAME BEHAVIOUR AS "overlapB" choice). Two more possible choices
are "overpairs" and "nonpairs". These will return a file with concatenated regions from fileA and fileB
(similar to the BEDPE format of BEDTools) with additional columns that represent:
a. in the case of "overpairs", statistics about the distances of regions found to be overlapping between
fileA and fileB (the statistics depend on the I<extend>, I<autoextend>, I<percent> and I<mode> options.
For example, if I<mode> and I<extend> are given, the distances are from the centers of the regions while
in any other case, all the distances from the start, end and center of regions are reported.
b. in the case of "nonpairs", the distances from non-overlapping regions are reported, if the non-overlapping
regions fall in a range of I<gap> upstream and downstream (where I<gap> must be specified). The maximum
number of non-overlapping regions is defined by the I<maxud> parameter.

=item I<maxud>

Use this option to define the maximum number of non-overlapping regions that will be reported using the
"nonpairs" option of I<output> parameter, within an upstream/downstream region defined by the I<gap> input
parameter.

=item I<keeporder> B<(optional)>

Use this parameter if you want to force the lines of the output files to be in the same order (e.g. sorted 
per chromosome or gene name) as the input files. This is accomplished through the module Tie::IxHash::Easy
which must be present in your machine. If the module is not present, the I<keeporder> option is deactivated.
Keep in mind that maintaining the order requires slighlty more memory during runtime.

=item I<dryrun> B<(optional)>

Use this option if you wish to do a "dry-run", that is just display statistics about chosen overlaps and
not write any output files.

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
use IntervalTree;

#use lib 'D:/Software/hts-tools/HTS-Tools/lib';
use lib '/media/HD4/Fleming/hts-tools/HTS-Tools/lib';
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
	$helper->set_logger($params->{"log"}) if (defined($params->{"log"}));
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
	
	# Copy some memory-less variables to avoid rewriting the whole thing...
	my $fileA = $self->get("inputA");
	my $fileB = $self->get("inputB");
	my $any = $self->get("any");
	my $mode = $self->get("mode");
	my $agap = $self->get("gap");
	my $maxud = $self->get("maxud");
	my $autoxtend = $self->get("autoextend");
	my $exact = $self->get("exact");
	my $both = $self->get("both");
	my $multi = $self->get("multi");
	my @percent = (defined($self->get("percent"))) ? (@{$self->get("percent")}) : (());
	my @extend = (defined($self->get("extend"))) ? (@{$self->get("extend")}) : (());
	my @out = (defined($self->get("output"))) ? (@{$self->get("output")}) : (());
	my $waitbar = $self->get("waitbar");
	my $ofileA = $fileA;
	my $ofileB = $fileB; # Keep original filenames (verbose purposes)

	if ($self->get("sort"))
	{
		($fileA,$fileB) = $self->sort_input($fileA,$fileB);
	}
	
	# Bavard a little... I like information...
	$helper->disp("Type of overlap: any overlap, percentage overlap ignored...") if ($any && @percent);
	$helper->disp("Type of overlap: any overlap") if ($any && !@percent);
	$helper->disp("Type of overlap: $percent[0]% overlap") if ($percent[0] && !$percent[1] && !$any);
	$helper->disp("Extending region modes upstream $extend[0] and downstream $extend[1] bps") if (@extend && !$autoxtend);
	$helper->disp("Region modes extension on each side will be auto-calculated...") if ($autoxtend);
	if ($percent[1])
	{
		$helper->disp("Multiple overlap percentages... Running in batch mode to determine overlapping distributions...");
		$helper->disp("No overlapping output files will be produced..."); 
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

	$helper->disp("Retrieving distances between non-overlapping regions if distance <= $agap bps") if (($op || $np) && $agap);
	$helper->disp(" ");

	# Get number of lines for fileA (mostly for waitbar, I know, inefficient but I like waitbars ;-)
	my $linesA = $helper->count_lines($fileA) if ($waitbar);

	# Suck in fileB
	my ($chromtreeB,$headerB) = $self->read_input($fileB,$ofileB);

	# General varianles to be used
	my (%overlapA,%overlapB,%onlyA,%onlyB,%saveB);
	my (@overpairs,@nonpairs,@ds,@rest);
	my ($chr,$start,$end,$tree,$result,$cit,$line,$headerA,$ups,$downs,$i,$j);
	my %strands = $self->strand_hash;

	if ($self->get("keeporder"))
	{
		tie %overlapA, "Tie::IxHash::Easy";
		tie %overlapB, "Tie::IxHash::Easy";
		tie %onlyA, "Tie::IxHash::Easy";
		tie %onlyB, "Tie::IxHash::Easy";
		tie %saveB, "Tie::IxHash::Easy";
	}
	
	##########################################################################################

	# Do we run in batch mode in order to determine distributions?
	if (!$any && $percent[1])
	{
		my @distributions;

		if (@extend || $autoxtend)
		{	
			if ($autoxtend) # We have to suck in the whole peak file once...
			{
				my @medmodes;
				my ($chromtreeA,$headerA) = $self->read_input($fileA,$ofileA);
				while(($chr,$tree) = each(%$chromtreeA))
				{
					$tree->traverse(sub {
						push(@medmodes,$_[0]->{"interval"}->{"end"} - $_[0]->{"interval"}->{"start"});
					});
				}
				my $med = $helper->median((@medmodes,$self->get_lengths($chromtreeB)));
				@extend = (int($med/2),int($med/2));
				$helper->disp("Median region length is $med bps. Extending each region mode $extend[1] bps on each side...\n");
				for ($i=0; $i<@percent; $i++)
				{	
					$helper->disp("Overlap percentage: $percent[$i]");
					while(($chr,$tree) = each($chromtreeA))
					{
						$helper->disp("Traversing chromosome $chr...");
						$tree->traverse(sub {
							(!$chromtreeB->{$chr}) ? ($result->[0] = 0) :
							($result = $self->search_percent_center($_[0]->{"interval"},$mode,@extend,$percent[$i]/100,$chromtreeB->{$chr}));
							
							if ($result->[0])
							{
								# We don't put an if for overlapA because we always start from this point
								# TODO: Add a report-once option in case a region from one input overlaps multiple regions in the other input
								$overlapA{$chr} = IntervalTree->new() if (!$overlapA{$chr});
								for ($j=0; $j<scalar @$result; $j++)
								{
									$overlapA{$chr}->insert_interval(
										IntervalTree::Interval->new(
											$_[0]->{"interval"}->{"start"},
											$_[0]->{"interval"}->{"end"},{
												"id" => (!defined($_[0]->{"interval"}->{"value"}->{"id"})) ?
													$chr.":".$_[0]->{"interval"}->{"start"}."-".$_[0]->{"interval"}->{"end"} :
													($_[0]->{"interval"}->{"value"}->{"id"}),
												"rest" => $_[0]->{"interval"}->{"value"}->{"rest"}
											},$chr,
											(!defined($_[0]->{"interval"}->{"value"}->{"strand"})) ? (undef) :
												($_[0]->{"interval"}->{"value"}->{"strand"})));

									if ($op)
									{
										@ds = $self->dists_center($_[0]->{"interval"},$result->[$j],$mode,@extend);
										#push(@overpairs,$self->node2text($_[0]->{"interval"})."\t".$result->[$j]."\t".$ds[1]);
										push(@overpairs,$ds[1]);
									}
								}

								if ($ovB)
								{
									$overlapB{$chr} = IntervalTree->new() if (!$overlapB{$chr});
									for ($j=0; $j<scalar @$result; $j++)
									{
										$overlapB{$chr}->insert_interval(
											IntervalTree::Interval->new(
												$result->[$j]->{"start"},$result->[$j]->{"end"},{
													"id" => (!defined($result->[$j]->{"value"}->{"id"})) ?
														$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"} :
														($result->[$j]->{"value"}->{"id"}),
													"rest" => $result->[$j]->{"value"}->{"rest"}
												},$chr,
												(!defined($result->[$j]->{"value"}->{"strand"})) ? (undef) :
													($result->[$j]->{"value"}->{"strand"})));
									}
								}

								if ($oB)
								{
									for ($j=0; $j<scalar @$result; $j++)
									{
										$saveB{$chr}{$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"}}++;
									}
								}
							}
							else # We don't find an overlap, we report onlyA and onlyB if wished
							{
								if ($oA)
								{
									$onlyA{$chr} = IntervalTree->new() if (!$onlyA{$chr});									
									$onlyA{$chr}->insert_interval(
										IntervalTree::Interval->new(
											$_[0]->{"interval"}->{"start"},$_[0]->{"interval"}->{"end"},{
												"id" => (!defined($_[0]->{"interval"}->{"value"}->{"id"})) ?
													$chr.":".$_[0]->{"interval"}->{"start"}."-".$_[0]->{"interval"}->{"end"} :
													($_[0]->{"interval"}->{"value"}->{"id"}),
												"rest" => $_[0]->{"interval"}->{"value"}->{"rest"}
											},$chr,
											(!defined($_[0]->{"interval"}->{"value"}->{"strand"})) ? (undef) :
												($_[0]->{"interval"}->{"value"}->{"strand"})));
								}
							}
						});
					}

					if ($oB)
					{
						$helper->disp("Retrieving onlyB regions...");
						%onlyB = $self->make_onlyB_tree(\%saveB);
					}
					
					push(@{$distributions[$i]},$self->chrom_itree_size(\%overlapA)) if ($ovA);
					push(@{$distributions[$i]},$self->chrom_itree_size(\%overlapB)) if ($ovB);
					push(@{$distributions[$i]},$self->chrom_itree_size(\%onlyA)) if ($oA);
					push(@{$distributions[$i]},$self->chrom_itree_size(\%onlyB)) if ($oB);
					push(@{$distributions[$i]},$helper->mean(@overpairs)) if ($op);
					push(@{$distributions[$i]},$helper->median(@overpairs)) if ($op);

					%overlapA = ();
					%overlapB = ();
					%onlyA = ();
					%onlyB = ();
					%saveB = ();
					@overpairs = ();
					@nonpairs = ();
					@ds = ();
				}
			}
			else
			{
				for ($i=0; $i<@percent; $i++)
				{
					$helper->disp("Overlap percentage: $percent[$i]");
					open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";
					$helper->disp("Reading file $ofileA and processing overlaps...");
					$helper->waitbar_init if ($waitbar);
					$line = <INA>;
					$headerA = $helper->decide_header($line);
					seek(INA,0,0) if (!$headerA);
					while ($line = <INA>)
					{
						$helper->waitbar_update($.,$linesA) if ($waitbar);	
						next if ($line =~/^chrM/);
						next if ($line =~/rand|hap|chrU/);
						$line =~ s/\r|\n$//g;
						($chr,$start,$end,@rest) = split(/\t/,$line);
						next if (!$chromtreeB->{$chr});

						$cit = IntervalTree::Interval->new(
							$start,$end,{
							"id" => $chr.":".$start."-".$end,
							"rest" => join("\t",@rest)
						},$chr,
						(grep {$_ eq $rest[2]} keys(%strands)) ? ($strands{$rest[2]}) : (undef));
							
						$result = $self->search_percent_center($cit,$mode,@extend,$percent[$i]/100,$chromtreeB->{$chr});

						if ($result->[0])
						{
							$overlapA{$chr} = IntervalTree->new() if (!$overlapA{$chr});
							for ($j=0; $j<scalar @$result; $j++)
							{
								$overlapA{$chr}->insert_interval(
									IntervalTree::Interval->new(
										$cit->{"start"},$cit->{"end"},{
											"id" => (!defined($cit->{"value"}->{"id"})) ?
												$chr.":".$cit->{"start"}."-".$cit->{"end"} :
												($cit->{"value"}->{"id"}),
											"rest" => $cit->{"value"}->{"rest"}
										},$chr,
										(!defined($cit->{"value"}->{"strand"})) ? (undef) :
											($cit->{"value"}->{"strand"})));

								if ($op)
								{
									@ds = $self->dists_center($cit,$result->[$j],$mode,@extend);
									push(@overpairs,$ds[1]);
								}
							}

							if ($ovB)
							{
								$overlapB{$chr} = IntervalTree->new() if (!$overlapB{$chr});
								for ($j=0; $j<scalar @$result; $j++)
								{
									$overlapB{$chr}->insert_interval(
										IntervalTree::Interval->new(
											$result->[$j]->{"start"},$result->[$j]->{"end"},{
												"id" => (!defined($result->[$j]->{"value"}->{"id"})) ?
													$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"} :
													($result->[$j]->{"value"}->{"id"}),
												"rest" => $result->[$j]->{"value"}->{"rest"}
											},$chr,
											(!defined($result->[$j]->{"value"}->{"strand"})) ? (undef) :
												($result->[$j]->{"value"}->{"strand"})));
								}
							}

							if ($oB)
							{
								for ($j=0; $j<scalar @$result; $j++)
								{
									$saveB{$chr}{$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"}}++;
								}
							}
						}
						else # We don't find an overlap, we report onlyA and onlyB if wished
						{
							if ($oA)
							{
								$onlyA{$chr} = IntervalTree->new() if (!$onlyA{$chr});									
								$onlyA{$chr}->insert_interval(
									IntervalTree::Interval->new(
										$cit->{"start"},$cit->{"end"},{
											"id" => (!defined($cit->{"value"}->{"id"})) ?
												$chr.":".$cit->{"start"}."-".$cit->{"end"} :
												($cit->{"value"}->{"id"}),
											"rest" => $cit->{"value"}->{"rest"}
										},$chr,
										(!defined($cit->{"value"}->{"strand"})) ? (undef) :
											($cit->{"value"}->{"strand"})));
							}
						}
					}
					close(INA);

					if ($oB)
					{
						$helper->disp("Retrieving onlyB regions...");
						%onlyB = $self->make_onlyB_tree(\%saveB);
					}
					
					push(@{$distributions[$i]},$self->chrom_itree_size(\%overlapA)) if ($ovA);
					push(@{$distributions[$i]},$self->chrom_itree_size(\%overlapB)) if ($ovB);
					push(@{$distributions[$i]},$self->chrom_itree_size(\%onlyA)) if ($oA);
					push(@{$distributions[$i]},$self->chrom_itree_size(\%onlyB)) if ($oB);
					push(@{$distributions[$i]},$helper->mean(@overpairs)) if ($op);
					push(@{$distributions[$i]},$helper->median(@overpairs)) if ($op);
					
					%overlapA = ();
					%overlapB = ();
					%onlyA = ();
					%onlyB = ();
					%saveB = ();
					@overpairs = ();
					@nonpairs = ();
					@ds = ();
				}
			}
		}
		else
		{
			for ($i=0; $i<@percent; $i++)
			{
				$helper->disp("Overlap percentage: $percent[$i]");
				open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";
				$helper->disp("Reading file $ofileA and processing overlaps...");
				$helper->waitbar_init if ($waitbar);
				$line = <INA>;
				$headerA = $helper->decide_header($line);
				seek(INA,0,0) if (!$headerA);
				while ($line = <INA>)
				{
					$helper->waitbar_update($.,$linesA) if ($waitbar);	
					next if ($line =~/^chrM/);
					next if ($line =~/rand|hap|chrU/);
					$line =~ s/\r|\n$//g;
					($chr,$start,$end,@rest) = split(/\t/,$line);
					next if (!$chromtreeB->{$chr}); # A chromosome might not exist in a peak file

					$cit = IntervalTree::Interval->new(
						$start,$end,{
						"id" => $chr.":".$start."-".$end,
						"rest" => join("\t",@rest)
					},$chr,
					(grep {$_ eq $rest[2]} keys(%strands)) ? ($strands{$rest[2]}) : (undef));
					
					if ($exact)
					{
						($both) ? ($result = $self->search_percent_exact_both($cit,$percent[$i]/100,$chromtreeB->{$chr})) :
						($result = $self->search_percent_exact($cit,$percent[$i]/100,$chromtreeB->{$chr}));
					}
					else
					{
						($both) ? ($result = $self->search_percent_both($cit,$percent[$i]/100,$chromtreeB->{$chr})) :
						($result = $self->search_percent($cit,$percent[$i]/100,$chromtreeB->{$chr}));
					}
					
					if ($result->[0])
					{
						$overlapA{$chr} = IntervalTree->new() if (!$overlapA{$chr});
						for ($j=0; $j<scalar @$result; $j++)
						{
							$overlapA{$chr}->insert_interval(
								IntervalTree::Interval->new(
									$cit->{"start"},$cit->{"end"},{
										"id" => (!defined($cit->{"value"}->{"id"})) ?
											$chr.":".$cit->{"start"}."-".$cit->{"end"} :
											($cit->{"value"}->{"id"}),
										"rest" => $cit->{"value"}->{"rest"}
									},$chr,
									(!defined($cit->{"value"}->{"strand"})) ? (undef) :
										($cit->{"value"}->{"strand"})));

							if ($op)
							{
								@ds = $self->dists_every($cit,$result->[$j],$mode);
								push(@overpairs,$ds[1]);
							}
						}

						if ($ovB)
						{
							$overlapB{$chr} = IntervalTree->new() if (!$overlapB{$chr});
							for ($j=0; $j<scalar @$result; $j++)
							{
								$overlapB{$chr}->insert_interval(
									IntervalTree::Interval->new(
										$result->[$j]->{"start"},$result->[$j]->{"end"},{
											"id" => (!defined($result->[$j]->{"value"}->{"id"})) ?
												$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"} :
												($result->[$j]->{"value"}->{"id"}),
											"rest" => $result->[$j]->{"value"}->{"rest"}
										},$chr,
										(!defined($result->[$j]->{"value"}->{"strand"})) ? (undef) :
											($result->[$j]->{"value"}->{"strand"})));
							}
						}

						if ($oB)
						{
							for ($j=0; $j<scalar @$result; $j++)
							{
								$saveB{$chr}{$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"}}++;
							}
						}
					}
					else # We don't find an overlap, we report onlyA and onlyB if wished
					{
						if ($oA)
						{
							$onlyA{$chr} = IntervalTree->new() if (!$onlyA{$chr});									
							$onlyA{$chr}->insert_interval(
								IntervalTree::Interval->new(
									$cit->{"start"},$cit->{"end"},{
										"id" => (!defined($cit->{"value"}->{"id"})) ?
											$chr.":".$cit->{"start"}."-".$cit->{"end"} :
											($cit->{"value"}->{"id"}),
										"rest" => $cit->{"value"}->{"rest"}
									},$chr,
									(!defined($cit->{"value"}->{"strand"})) ? (undef) :
										($cit->{"value"}->{"strand"})));
						}
					}
				}
				close(INA);

				if ($oB)
				{
					$helper->disp("Retrieving onlyB regions...");
					%onlyB = $self->make_onlyB_tree(\%saveB);
				}
				
				push(@{$distributions[$i]},$self->chrom_itree_size(\%overlapA)) if ($ovA);
				push(@{$distributions[$i]},$self->chrom_itree_size(\%overlapB)) if ($ovB);
				push(@{$distributions[$i]},$self->chrom_itree_size(\%onlyA)) if ($oA);
				push(@{$distributions[$i]},$self->chrom_itree_size(\%onlyB)) if ($oB);
				push(@{$distributions[$i]},$helper->mean(@overpairs)) if ($op);
				push(@{$distributions[$i]},$helper->median(@overpairs)) if ($op);

				%overlapA = ();
				%overlapB = ();
				%onlyA = ();
				%onlyB = ();
				%saveB = ();
				@overpairs = ();
				@nonpairs = ();
				@ds = ();
			}
		}

		my $co = 0;
		my $time = $helper->now("machine");
		my $fdn = "distrib_".$time;
		$helper->disp("Writing output in $fdn.txt");
		open(OUTDISTRIB,">$fdn.txt");
		my ($cdo,@disthead);
		push(@disthead,"Percentage");
		push(@disthead,"Overlap ".basename($fileA)) if ($ovA);
		push(@disthead,"Overlap ".basename($fileB)) if ($ovB);
		push(@disthead,"Only ".basename($fileA)) if ($oA);
		push(@disthead,"Only ".basename($fileB)) if ($oB);
		push(@disthead,"Mean distance\tMedian distance") if ($op);
		print OUTDISTRIB join("\t",@disthead),"\n";
		foreach $cdo (@distributions)
		{
			print OUTDISTRIB $percent[$co]."\t";
			print OUTDISTRIB join("\t",@{$cdo}),"\n";
			$co++;
		}
		close(OUTDISTRIB);
		
		# Remove garbage
		$helper->cleanup;
		$helper->disp("Finished!\n\n");
		exit;
	}

	##########################################################################################

	if (@extend || $autoxtend)
	{
		if ($autoxtend) # We have to suck in the whole peak file once...
		{
			my @medmodes;
			my ($chromtreeA,$headerA) = $self->read_input($fileA,$ofileA);
			while(($chr,$tree) = each(%$chromtreeA))
			{
				$tree->traverse(sub {
					push(@medmodes,$_[0]->{"interval"}->{"end"} - $_[0]->{"interval"}->{"start"});
				});
			}
			my $med = $helper->median((@medmodes,$self->get_lengths($chromtreeB)));
			@extend = (int($med/2),int($med/2));
			$helper->disp("Median region length is $med bps. Extending each region mode $extend[1] bps on each side...\n");

			while(($chr,$tree) = each($chromtreeA))
			{
				$helper->disp("Traversing chromosome $chr...");
				$tree->traverse(sub {
					if ($chromtreeB->{$chr})
					{
						($any) ? ($result = $self->search_any_center($_[0]->{"interval"},$mode,@extend,$chromtreeB->{$chr})) :
						($result = $self->search_percent_center($_[0]->{"interval"},$mode,@extend,$percent[0]/100,$chromtreeB->{$chr}));
					} else { $result->[0] = 0; }
					
					if ($result->[0])
					{
						$overlapA{$chr} = IntervalTree->new() if (!$overlapA{$chr});
						for ($j=0; $j<scalar @$result; $j++)
						{
							$overlapA{$chr}->insert_interval(
								IntervalTree::Interval->new(
									$_[0]->{"interval"}->{"start"},$_[0]->{"interval"}->{"end"},{
										"id" => (!defined($_[0]->{"interval"}->{"value"}->{"id"})) ?
											$chr.":".$_[0]->{"interval"}->{"start"}."-".$_[0]->{"interval"}->{"end"} :
											($_[0]->{"interval"}->{"value"}->{"id"}),
										"rest" => $_[0]->{"interval"}->{"value"}->{"rest"}
									},$chr,
									(!defined($_[0]->{"interval"}->{"value"}->{"strand"})) ? (undef) :
										($_[0]->{"interval"}->{"value"}->{"strand"})));

							if ($op)
							{
								@ds = $self->dists_center($_[0]->{"interval"},$result->[$j],$mode,@extend);
								push(@overpairs,$chr."\t".$self->node2text($_[0]->{"interval"})."\t".$chr."\t".
									$self->node2text($result->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
							}
						}

						if ($ovB)
						{
							$overlapB{$chr} = IntervalTree->new() if (!$overlapB{$chr});
							for ($j=0; $j<scalar @$result; $j++)
							{
								$overlapB{$chr}->insert_interval(
									IntervalTree::Interval->new(
										$result->[$j]->{"start"},$result->[$j]->{"end"},{
											"id" => (!defined($result->[$j]->{"value"}->{"id"})) ?
												$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"} :
												($result->[$j]->{"value"}->{"id"}),
											"rest" => $result->[$j]->{"value"}->{"rest"}
										},$chr,
										(!defined($result->[$j]->{"value"}->{"strand"})) ? (undef) :
											($result->[$j]->{"value"}->{"strand"})));
							}
						}

						if ($oB)
						{
							for ($j=0; $j<scalar @$result; $j++)
							{
								$saveB{$chr}{$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"}}++;
							}
						}
					}
					else # We don't find an overlap, we report onlyA and onlyB if wished
					{
						if ($oA)
						{
							$onlyA{$chr} = IntervalTree->new() if (!$onlyA{$chr});									
							$onlyA{$chr}->insert_interval(
								IntervalTree::Interval->new(
									$_[0]->{"interval"}->{"start"},$_[0]->{"interval"}->{"end"},{
										"id" => (!defined($_[0]->{"interval"}->{"value"}->{"id"})) ?
											$chr.":".$_[0]->{"interval"}->{"start"}."-".$_[0]->{"interval"}->{"end"} :
											($_[0]->{"interval"}->{"value"}->{"id"}),
										"rest" => $_[0]->{"interval"}->{"value"}->{"rest"}
									},$chr,
									(!defined($_[0]->{"interval"}->{"value"}->{"strand"})) ? (undef) :
										($_[0]->{"interval"}->{"value"}->{"strand"})));
						}

						if ($np)
						{
							$ups = $chromtreeB->{$chr}->upstream_of_interval($cit,$maxud,$agap);
							$downs = $chromtreeB->{$chr}->downstream_of_interval($cit,$maxud,$agap);
							if ($ups)
							{
								for ($j=0; $j<scalar @$ups; $j++)
								{
									@ds = $self->dists_every($cit,$ups->[$j],$mode);
									push(@nonpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".
										$self->node2text($ups->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
								}
							}
							if ($downs)
							{
								for ($j=0; $j<scalar @$downs; $j++)
								{
									@ds = $self->dists_every($cit,$downs->[$j],$mode);
									push(@nonpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".
										$self->node2text($downs->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
								}
							}
						}
					}
				});
			}
		}
		else
		{
			$helper->disp("Reading file $ofileA and processing overlaps...");
			$helper->waitbar_init if ($waitbar);
			open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";
			$line = <INA>;
			$headerA = $helper->decide_header($line);
			seek(INA,0,0) if (!$headerA);
			while ($line = <INA>)
			{
				$helper->waitbar_update($.,$linesA) if ($waitbar);
				next if ($line =~/^chrM/);
				next if ($line =~/rand|hap|chrU/);
				$line =~ s/\r|\n$//g;
				($chr,$start,$end,@rest) = split(/\t/,$line);
				next if (!$chromtreeB->{$chr});

				$cit = IntervalTree::Interval->new(
					$start,$end,{
					"id" => $chr.":".$start."-".$end,
					"rest" => join("\t",@rest)
				},$chr,
				(grep {$_ eq $rest[2]} keys(%strands)) ? ($strands{$rest[2]}) : (undef));

				($any) ? ($result = $self->search_any_center($cit,$mode,@extend,$chromtreeB->{$chr})) :
				($result = $self->search_percent_center($cit,$mode,@extend,$percent[0]/100,$chromtreeB->{$chr}));

				if ($result->[0])
				{
					$overlapA{$chr} = IntervalTree->new() if (!$overlapA{$chr});
					for ($j=0; $j<scalar @$result; $j++)
					{
						$overlapA{$chr}->insert_interval(
							IntervalTree::Interval->new(
								$cit->{"start"},$cit->{"end"},{
									"id" => (!defined($cit->{"value"}->{"id"})) ?
										$chr.":".$cit->{"start"}."-".$cit->{"end"} :
										($cit->{"value"}->{"id"}),
									"rest" => $cit->{"value"}->{"rest"}
								},$chr,
								(!defined($cit->{"value"}->{"strand"})) ? (undef) :
									($cit->{"value"}->{"strand"})));

						if ($op)
						{
							@ds = $self->dists_center($cit,$result->[$j],$mode,@extend);
							push(@overpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".$self->node2text($result->[$j])."\t".$ds[1]);
						}
					}

					if ($ovB)
					{
						$overlapB{$chr} = IntervalTree->new() if (!$overlapB{$chr});

						for ($j=0; $j<scalar @$result; $j++)
						{
							$overlapB{$chr}->insert_interval(
								IntervalTree::Interval->new(
									$result->[$j]->{"start"},$result->[$j]->{"end"},{
										"id" => (!defined($result->[$j]->{"value"}->{"id"})) ?
											$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"} :
											($result->[$j]->{"value"}->{"id"}),
										"rest" => $result->[$j]->{"value"}->{"rest"}
									},$chr,
									(!defined($result->[$j]->{"value"}->{"strand"})) ? (undef) :
										($result->[$j]->{"value"}->{"strand"})));
						}
					}

					if ($oB)
					{
						for ($j=0; $j<scalar @$result; $j++)
						{
							$saveB{$chr}{$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"}}++;
						}
					}
				}
				else # We don't find an overlap, we report onlyA and onlyB if wished
				{
					if ($oA)
					{
						$onlyA{$chr} = IntervalTree->new() if (!$onlyA{$chr});									
						$onlyA{$chr}->insert_interval(
							IntervalTree::Interval->new(
								$cit->{"start"},$cit->{"end"},{
									"id" => (!defined($cit->{"value"}->{"id"})) ?
										$chr.":".$cit->{"start"}."-".$cit->{"end"} :
										($cit->{"value"}->{"id"}),
									"rest" => $cit->{"value"}->{"rest"}
								},$chr,
								(!defined($cit->{"value"}->{"strand"})) ? (undef) :
									($cit->{"value"}->{"strand"})));
					}

					if ($np)
					{
						$ups = $chromtreeB->{$chr}->upstream_of_interval($cit,$maxud,$agap);
						$downs = $chromtreeB->{$chr}->downstream_of_interval($cit,$maxud,$agap);
						if ($ups)
						{
							for ($j=0; $j<scalar @$ups; $j++)
							{
								@ds = $self->dists_every($cit,$ups->[$j],$mode);
								push(@nonpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".
									$self->node2text($ups->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
							}
						}
						if ($downs)
						{
							for ($j=0; $j<scalar @$downs; $j++)
							{
								@ds = $self->dists_every($cit,$downs->[$j],$mode);
								push(@nonpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".
									$self->node2text($downs->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
							}
						}
					}
				}
			}
			close(INA);
		}
	}
	else
	{
		$helper->disp("Reading file $ofileA and processing overlaps...");
		$helper->waitbar_init if ($waitbar);
		open (INA,$fileA) or croak "\nThe file $fileA does not exist!\n";
		$line = <INA>;
		$headerA = $helper->decide_header($line);
		seek(INA,0,0) if (!$headerA);
		while ($line = <INA>)
		{
			$helper->waitbar_update($.,$linesA) if ($waitbar);
			next if ($line =~ m/chrM|rand|hap|chrU/);
			$line =~ s/\r|\n$//g;
			($chr,$start,$end,@rest) = split(/\t/,$line);
			next if (!$chromtreeB->{$chr});

			$cit = IntervalTree::Interval->new(
				$start,$end,{
				"id" => $chr.":".$start."-".$end,
				"rest" => join("\t",@rest)
			},$chr,
			(grep {$_ eq $rest[2]} keys(%strands)) ? ($strands{$rest[2]}) : (undef));

			if ($any)
			{
				$result = $self->search_any($cit,$chromtreeB->{$chr});
			}
			else
			{
				if ($exact)
				{
					($both) ? ($result = $self->search_percent_exact_both($cit,$percent[0]/100,$chromtreeB->{$chr})) :
					($result = $self->search_percent_exact($cit,$percent[0]/100,$chromtreeB->{$chr}));
				}
				else
				{
					($both) ? ($result = $self->search_percent_both($cit,$percent[0]/100,$chromtreeB->{$chr})) :
					($result = $self->search_percent($cit,$percent[0]/100,$chromtreeB->{$chr}));
				}
			}

			if ($result->[0])
			{
				$overlapA{$chr} = IntervalTree->new() if (!$overlapA{$chr});
				for ($j=0; $j<scalar @$result; $j++)
				{
					$overlapA{$chr}->insert_interval(
						IntervalTree::Interval->new(
							$cit->{"start"},$cit->{"end"},{
								"id" => (!defined($cit->{"value"}->{"id"})) ?
									$chr.":".$cit->{"start"}."-".$cit->{"end"} :
									($cit->{"value"}->{"id"}),
								"rest" => $cit->{"value"}->{"rest"}
							},$chr,
							(!defined($cit->{"value"}->{"strand"})) ? (undef) :
								($cit->{"value"}->{"strand"})));

					if ($op)
					{
						@ds = $self->dists_every($cit,$result->[$j],$mode);
						push(@overpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".$self->node2text($result->[$j])."\t".$ds[1]);
					}
				}

				if ($ovB)
				{
					$overlapB{$chr} = IntervalTree->new() if (!$overlapB{$chr});
					for ($j=0; $j<scalar @$result; $j++)
					{
						$overlapB{$chr}->insert_interval(
							IntervalTree::Interval->new(
								$result->[$j]->{"start"},$result->[$j]->{"end"},{
									"id" => (!defined($result->[$j]->{"value"}->{"id"})) ?
										$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"} :
										($result->[$j]->{"value"}->{"id"}),
									"rest" => $result->[$j]->{"value"}->{"rest"}
								},$chr,
								(!defined($result->[$j]->{"value"}->{"strand"})) ? (undef) :
									($result->[$j]->{"value"}->{"strand"})));
					}
				}

				if ($oB)
				{
					for ($j=0; $j<scalar @$result; $j++)
					{
						$saveB{$chr}{$chr.":".$result->[$j]->{"start"}."-".$result->[$j]->{"end"}}++;
					}
				}
			}
			else # We don't find an overlap, we report onlyA and onlyB if wished
			{
				if ($oA)
				{
					$onlyA{$chr} = IntervalTree->new() if (!$onlyA{$chr});									
					$onlyA{$chr}->insert_interval(
						IntervalTree::Interval->new(
							$cit->{"start"},$cit->{"end"},{
								"id" => (!defined($cit->{"value"}->{"id"})) ?
									$chr.":".$cit->{"start"}."-".$cit->{"end"} :
									($cit->{"value"}->{"id"}),
								"rest" => $cit->{"value"}->{"rest"}
							},$chr,
							(!defined($cit->{"value"}->{"strand"})) ? (undef) :
								($cit->{"value"}->{"strand"})));
				}

				if ($np)
				{
					$ups = $chromtreeB->{$chr}->upstream_of_interval($cit,$maxud,$agap);
					$downs = $chromtreeB->{$chr}->downstream_of_interval($cit,$maxud,$agap);
					if ($ups)
					{
						for ($j=0; $j<scalar @$ups; $j++)
						{
							@ds = $self->dists_every($cit,$ups->[$j],$mode);
							push(@nonpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".
								$self->node2text($ups->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
						}
					}
					if ($downs)
					{
						for ($j=0; $j<scalar @$downs; $j++)
						{
							@ds = $self->dists_every($cit,$downs->[$j],$mode);
							push(@nonpairs,$chr."\t".$self->node2text($cit)."\t".$chr."\t".
								$self->node2text($downs->[$j])."\t".$ds[0]."\t".$ds[1]."\t".$ds[2]);
						}
					}
				}
			}
		}
		close(INA);	
	}

	if ($oB)
	{
		$helper->disp("Retrieving onlyB regions...");
		%onlyB = $self->make_onlyB_tree(\%saveB);
	}

	$helper->disp(" ") if (!$waitbar);

	if (!$self->get("dryrun"))
	{
		$helper->disp("Writing output...");
		my ($bA,$bB);
		$bA = fileparse($ofileA,'\.[^.]*') if ($ovA || $oA);
		$bB = fileparse($ofileB,'\.[^.]*') if ($ovB || $oB);
		# A and B overlap elements with elements of file A
		$self->print_itree_output($ofileA,$ofileB,"OVERLAP_FROM_$bA",\%overlapA,$headerA) if ($ovA); 
		# A and B overlap elements with elements of file B
		$self->print_itree_output($ofileA,$ofileB,"OVERLAP_FROM_$bB",\%overlapB,$headerB) if ($ovB);
		# A only elements
		$self->print_itree_output($ofileA,$ofileB,"ONLY_$bA",\%onlyA,$headerA) if ($oA);
		# B only elements
		$self->print_itree_output($ofileA,$ofileB,"ONLY_$bB",\%onlyB,$headerB) if ($oB);
		# Distances of overlaping
		$self->print_array($ofileA,$ofileB,"OVERDIST",$headerA,$headerB,@overpairs) if ($op);
		# Distances of non-overlaping
		$self->print_array($ofileA,$ofileB,"NONDIST",$headerA,$headerB,@nonpairs) if ($np);
	}

	# Display stats
	if (!$self->get("silent"))
	{
		$helper->disp("\n--- STATS ---");
		my ($lA,$lB);
		($waitbar) ? ($lA = $linesA) : ($lA = $helper->count_lines($fileA));
		$lB = $helper->count_lines($fileB);
		$lA-- if ($headerA);
		$lB-- if ($headerB);
		if ($ovA || $op)
		{
			my $covA = $self->chrom_itree_size(\%overlapA);
			$helper->disp("$covA out of $lA regions from ".basename($ofileA)." overlap with regions from ".basename($ofileB));
		}
		if ($ovB)
		{
			my $covB = $self->chrom_itree_size(\%overlapB);
			$helper->disp("$covB out of $lB regions from ".basename($ofileB)." overlap with regions from ".basename($ofileA));
		}
		if ($oA)
		{
			my $coA = $self->chrom_itree_size(\%onlyA);
			$helper->disp("$coA out of $lA regions exist only in ".basename($ofileA));
		}
		if ($oB)
		{
			my $coB = $self->chrom_itree_size(\%onlyB);
			$helper->disp("$coB out of $lB regions exist only in ".basename($ofileB));
		}
		$helper->disp(" ");
	}

	# Remove garbage
	$helper->cleanup;
	$helper->disp("Finished!\n\n") if (!$multi);
}

=head2 make_onlyB_tree

Make an interval tree out of a hash storing intervals to exclude from an original file. This is just a hack as the current IntervalTree
implementation in pure Perl does not have a remove function. Internal use.

	$intersecter->make_onlyB_tree(\%onlyB);

=cut

sub make_onlyB_tree
{
	my ($self,$bhash) = @_;
	my ($chr,$start,$end,$tree,$he,$line);
	my @rest;
	my %onlyB;
	my %strands = $self->strand_hash;
	
	open(IN,$self->get("inputB"));
	$line = <IN>;
	$he = $helper->decide_header($line);
	seek(IN,0,0) if (!$he);
	while ($line = <IN>)
	{
		next if ($line =~ m/chrM|rand|chrU|hap/i);
		$line =~ s/\r|\n$//g; # Make sure to remove carriage returns
		($chr,$start,$end,@rest) = split(/\t/,$line);
		next if ($bhash->{$chr}->{$chr.":".$start."-".$end}); # Skip overlapping intervals, as stored in $bhash

		$onlyB{$chr} = IntervalTree->new() if (!$onlyB{$chr});
		$onlyB{$chr}->insert_interval(
			IntervalTree::Interval->new(
				$start,$end,{
					"id" => $chr.":".$start."-".$end,
					"rest" => join("\t",@rest)
				},$chr,
				(grep {$_ eq $rest[2]} keys(%strands)) ? ($strands{$rest[2]}) : (undef)));
	}
	close(IN);

	return(%onlyB);
}

=head2 read_input

Read an input file and create an Interval Tree. Internal use.

	$intersecter->read_input($file);

=cut

sub read_input
{
	my ($self,$file,$ofile) = @_;
	$ofile = $file unless ($ofile);
	my ($line,$header,$chr,$start,$end);
	my @rest;
	my %chromosome;
	my %strands = $self->strand_hash;

	tie %chromosome, "Tie::IxHash::Easy" if ($self->get("keeporder"));
	
	open(INPUT,$file) or croak "\nThe file $file does not exist!\n";
	$helper->disp("Reading file $ofile...");
	$line = <INPUT>;
	$header = $helper->decide_header($line);
	seek(INPUT,0,0) if (!$header);
	while ($line = <INPUT>)
	{
		next if ($line =~ m/chrM|rand|chrU|hap/i);
		$line =~ s/\r|\n$//g; # Make sure to remove carriage returns
		($chr,$start,$end,@rest) = split(/\t/,$line);
		
		$chromosome{$chr} = IntervalTree->new() if (!$chromosome{$chr});

		# Fill the IntervalTree
		$chromosome{$chr}->insert_interval(
			IntervalTree::Interval->new(
				$start,$end,{
					"id" => $chr.":".$start."-".$end,
					"rest" => join("\t",@rest)
				},$chr,
				(grep {$_ eq $rest[2]} keys(%strands)) ? ($strands{$rest[2]}) : (undef)));
				
		# If we have to define the "onlyB" structure, make a copy of the above, following the old style	
	}
	close(INPUT);
	
	return(\%chromosome,$header);
}

=head2 print_itree_output

Module specific output printing function. Internal use.

	$intersecter->print_itree_output($A,$B,$output_type,$header,$the_hash);

=cut

sub print_itree_output
{
	my ($self,$infileA,$infileB,$otype,$chromhash,$he) = @_;
	my ($chr,$tree);
	my %seen;
	my $outfilename = $self->create_output_file($infileA,$infileB,$otype);
	open(OUTPUT,">$outfilename");
	print OUTPUT "$he\n" if ($he);
	if ($self->get("reportonce"))
	{
		while(($chr,$tree) = each(%$chromhash))
		{
			$tree->traverse(sub {
				print OUTPUT $chr."\t".$self->node2text($_[0]->{"interval"})."\n"
					if (!$seen{$chr}{$_[0]->{"interval"}->{"start"}.$_[0]->{"interval"}->{"end"}}++);
				$seen{$chr}{$_[0]->{"interval"}->{"start"}.$_[0]->{"interval"}->{"end"}}++;
			});
		}
	}
	else
	{
		while(($chr,$tree) = each(%$chromhash))
		{
			$tree->traverse(sub {
				print OUTPUT $chr."\t".$self->node2text($_[0]->{"interval"})."\n";
			});
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
	my ($self,$infileA,$infileB,$otype,$hA,$hB,@inarr) = @_;
	my $outfilename = $self->create_output_file($infileA,$infileB,$otype);
	open(OUTPUT,">$outfilename");
	print OUTPUT $hA."\t".$hB."\tdistance\n" if ($hA && $hB);
	print OUTPUT join("\n",@inarr);
	close(OUTPUT);
}

=head2 search_any

Binary search algorithm for any overlap between genomic regions. Internal use.

	$intersecter->search_any($start,$end,@candidate_areas);
	
=cut

sub search_any
{
	my ($self,$node,$tree) = @_;
	my ($ind,$start,$end,$currstart,$currend,$overstru,$i);
	my @result;

	$start = $node->{"start"};
	$end = $node->{"end"};
	$overstru = $tree->find($start,$end);
	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{	
			$currstart = $overstru->[$i]->{"start"};
			$currend = $overstru->[$i]->{"end"};
			
			if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end) ||
			($currstart < $start && $currend < $end && $start < $currend) ||
			($currstart > $start && $currend > $end && $end > $currstart))
			{
				push(@result,$overstru->[$i]);
			}
		}
		return(\@result);
	}
	return(\(0)); # Return back the node if not found
}

=head2 search_percent

Binary search algorithm for percent overlap between genomic regions. Internal use.

	$intersecter->search_percent($start,$end,$percentage,@candidate_areas);
	
=cut

sub search_percent
{
	my ($self,$node,$p,$tree) = @_;
	my ($ind,$start,$end,$currstart,$currend,$overstru,$i);
	my @result;

	$start = $node->{"start"};
	$end = $node->{"end"};
	$overstru = $tree->find($start,$end);
	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{
			$currstart = $overstru->[$i]->{"start"};
			$currend = $overstru->[$i]->{"end"};
			
			if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end))
			{
				push(@result,$overstru->[$i]);
			}
			elsif ($currstart < $start && $currend < $end && $start < $currend)
			{
				push(@result,$overstru->[$i]) if (($currend - $start) >= $p*($end - $start));
			}
			elsif ($currstart > $start && $currend > $end && $end > $currstart)
			{
				push(@result,$overstru->[$i]) if (($end - $currstart) >= $p*($end - $start));
			}
		}
		return(\@result);
	}
	return(\(0));
}

=head2 search_any_center

Binary search algorithm for any overlap between genomic regions using their centers. Internal use.

	$intersecter->search_any_center($mode,$position,$downstream,$upstream,@candidate_areas);
	
=cut

sub search_any_center
{
	my ($self,$node,$mpos,$uxval,$dxval,$tree) = @_;
	my ($start,$end,$center,$currstart,$currend,$overstru,$i,@arr,@marr);
	my @result;
	
	@marr = split(/\t/,$node->{"value"}->{"rest"});
	if ($mpos > 0)
	{
		$start = $marr[$mpos] - $uxval;
		$end = $marr[$mpos] + $dxval;
	}
	else
	{
		$center = $node->{"start"} + $helper->round(($node->{"end"} - $node->{"start"})/2);
		$start = $center - $uxval;
		$end = $center + $dxval;
	}
		
	$overstru = $tree->find($start,$end);

	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{
			@arr = split(/\t/,$overstru->[$i]->{"value"}->{"rest"}); 
			if ($mpos > 0)
			{
				$currstart = $arr[$mpos] - $dxval;
				$currend = $arr[$mpos] + $uxval;
			}
			else
			{
				$center = $overstru->[$i]->{"start"} + $helper->round(($overstru->[$i]->{"end"} - $overstru->[$i]->{"start"})/2);
				$currstart = $center - $dxval;
				$currend = $center + $uxval;
			}
			if (($currstart >= $start && $currend <= $end) ||
				($currstart <= $start && $currend >= $end) ||
				($currstart < $start && $currend < $end && $start < $currend) ||
				($currstart > $start && $currend > $end && $end > $currstart))
			{
				push(@result,$overstru->[$i]);
			}
		}
		return(\@result);
	}
	return(\(0));
}

=head2 search_any_center

Binary search algorithm for percentage overlap between genomic regions using their centers. Internal use.

	$intersecter->search_percent_center($mode,$position,$downstream,$upstream,$percentage,@candidate_areas);
	
=cut

sub search_percent_center
{
	my ($self,$node,$mpos,$uxval,$dxval,$p,$tree) = @_;
	my ($start,$end,$center,$currstart,$currend,$overstru,$i,@arr,@marr);
	my @result;

	@marr = split(/\t/,$node->{"value"}->{"rest"});
	if ($mpos > 0)
	{
		$start = $marr[$mpos] - $uxval;
		$end = $marr[$mpos] + $dxval;
	}
	else
	{
		$center = $node->{"start"} + $helper->round(($node->{"end"} - $node->{"start"})/2);
		$start = $center - $uxval;
		$end = $center + $dxval;
	}
	
	$overstru = $tree->find($start,$end);
	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{
			@arr = split(/\t/,$overstru->[$i]->{"value"}->{"rest"});
			if ($mpos > 0)
			{
				$currstart = $arr[$mpos] - $dxval;
				$currend = $arr[$mpos] + $uxval;
			}
			else
			{
				$center = $overstru->[$i]->{"start"} + $helper->round(($overstru->[$i]->{"end"} - $overstru->[$i]->{"start"})/2);
				$currstart = $center - $dxval;
				$currend = $center + $uxval;
			}
			if (($currstart >= $start && $currend <= $end) ||
				($currstart <= $start && $currend >= $end))
			{
				push(@result,$overstru->[$i]);
			}
			elsif ($currstart < $start && $currend < $end && $start < $currend)
			{
				push(@result,$overstru->[$i]) if (($currend - $start) >= $p*($end - $start));
			}
			elsif ($currstart > $start && $currend > $end && $end > $currstart)
			{
				push(@result,$overstru->[$i]) if (($end - $currstart) >= $p*($end - $start));
			}
		}
		return(\@result);
	}
	return(\(0));
}

=head2 search_percent_both

Binary search algorithm for percentage overlap between genomic regions for the "both" case. Internal use.

	$intersecter->search_any_center($start,$end,$percentage,@candidate_areas);
	
=cut

sub search_percent_both
{
	my ($self,$node,$p,$tree) = @_;
	my ($start,$end,$currstart,$currend,$diff,$overstru,$i);
	my @result;

	$start = $node->{"start"};
	$end = $node->{"end"};
	
	$overstru = $tree->find($start,$end);
	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{
			$currstart = $overstru->[$i]->{"start"};
			$currend = $overstru->[$i]->{"end"};
			
			if (($currstart >= $start && $currend <= $end) ||
				($currstart <= $start && $currend >= $end))
			{
				push(@result,$overstru->[$i]);
			}
			elsif ($currstart < $start && $currend < $end && $start < $currend)
			{
				$diff = $currend - $start;
				push(@result,$overstru->[$i]) if ($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart));
			}
			elsif ($currstart > $start && $currend > $end && $end > $currstart)
			{
				$diff = $end - $currstart;
				push(@result,$overstru->[$i]) if ($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart));
			}
		}
		return(\@result);
	}
	return(\(0));
}

=head2 search_percent_exact

Binary search algorithm for percentgae overlap between genomic regions for the "exact" case. Internal use.

	$intersecter->search_percent_exact($start,$end,$percentage,@candidate_areas);
	
=cut

sub search_percent_exact
{
	my ($self,$node,$p,$tree) = @_;
	my ($start,$end,$currstart,$currend,$overstru,$i);
	my @result;

	$start = $node->{"start"};
	$end = $node->{"end"};
	$overstru = $tree->find($start,$end);
	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{
			$currstart = $overstru->[$i]->{"start"};
			$currend = $overstru->[$i]->{"end"};

			if ($currstart >= $start && $currend <= $end) # TagB <= TagA
			{
				push(@result,$overstru->[$i]) if (($currend - $currstart) >= $p*($end - $start));
			}
			elsif ($currstart <= $start && $currend >= $end) # TagB >= TagA
			{
				push(@result,$overstru->[$i]) if(($end - $start) >= $p*($currend - $currstart));
			}
			elsif ($currstart < $start && $currend < $end && $start < $currend)
			{
				push(@result,$overstru->[$i]) if (($currend - $start) >= $p*($end - $start));
			}
			elsif ($currstart > $start && $currend > $end && $end > $currstart)
			{
				push(@result,$overstru->[$i]) if (($end - $currstart) >= $p*($end - $start));
			}
		}
		return(\@result);
	}
	return(\(0));
}

=head2 search_percent_both

Binary search algorithm for percentage overlap between genomic regions for the "exact" and "both" case.
Internal use.

	$intersecter->search_percent_both($start,$end,$percentage,@candidate_areas);
	
=cut

sub search_percent_exact_both
{
	my ($self,$node,$p,$tree) = @_;
	my ($start,$end,$currstart,$currend,$diff,$overstru,$i);
	my @result;

	$start = $node->{"start"};
	$end = $node->{"end"};
	$overstru = $tree->find($start,$end);
	if ($overstru)
	{
		for ($i=0; $i< scalar @$overstru; $i++)
		{
			$currstart = $overstru->[$i]->{"start"};
			$currend = $overstru->[$i]->{"end"};

			if ($currstart >= $start && $currend <= $end) # TagB <= TagA
			{
				push(@result,$overstru->[$i]) if (($currend - $currstart) >= $p*($end - $start));
			}
			elsif ($currstart <= $start && $currend >= $end) # TagB >= TagA
			{
				push(@result,$overstru->[$i]) if (($end - $start) >= $p*($currend - $currstart));
			}
			elsif ($currstart < $start && $currend < $end && $start < $currend)
			{
				$diff = $currend - $start;
				push(@result,$overstru->[$i])  if ($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart));
			}
			elsif ($currstart > $start && $currend > $end && $end > $currstart)
			{
				$diff = $end - $currstart;
				push(@result,$overstru->[$i]) if ($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart));
			}
		}
		return(\@result);
	}
	return(\(0));
}

=head2 dists_every

Distance calculation subroutine using Interval Tree nodes. Internal use.

	$intersecter->dists_every($A,$B,$ei);
	
=cut

sub dists_every
{
	my ($self,$A,$B,$ei) = @_;
	my ($sa,$ea,$ma,$sb,$eb,$mb,$ds,$da,$dc,$ca,$cb);
	my (@aA,@aB);

	
	@aA = split(/\t/,$A->{"value"}->{"rest"});
	@aB = split(/\t/,$B->{"value"}->{"rest"});
	if ($ei > 0)
	{
		$ca = $aA[$ei];
		$cb = $aB[$ei];
	}
	else
	{
		$ca = $A->{"start"} + $helper->round(($A->{"end"} - $A->{"start"})/2);
		$cb = $B->{"start"} + $helper->round(($B->{"end"} - $B->{"start"})/2);
	}

	($sa,$ea,$ma) = ($A->{"start"},$A->{"end"},$ca);
	($sb,$eb,$mb) = ($B->{"start"},$B->{"end"},$cb);
	
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
	elsif (($sa <= $sb && $ea >= $eb) || ($sa >= $sb && $ea <= $eb))
	{
		$ds = 0;
		$da = $ma - $mb;
		$dc = 0;
	}
	
	$ds = 0 if (!$ds);
	$da = 0 if (!$da);
	$dc = 0 if (!$dc);
	
	return (($ds,$da,$dc));
}

=head2 dists_every

Distance calculation subroutine using centers and Interval Tree nodes. Internal use.

	$intersecter->dists_center($A,$B,$ei,$up,$down);
	
=cut

sub dists_center
{	
	my ($self,$A,$B,$ei,$exu,$exd) = @_;
	my ($sa,$ea,$ma,$sb,$eb,$mb,$ds,$da,$dc,$ca,$cb);
	my (@aA,@aB);
	
	@aA = split(/\t/,$A->{"value"}->{"rest"});
	@aB = split(/\t/,$B->{"value"}->{"rest"});
	if ($ei > 0)
	{
		$ca = $aA[$ei];
		$cb = $aB[$ei];
	}
	else
	{
		$ca = $A->{"start"} + $helper->round(($A->{"end"} - $A->{"start"})/2);
		$cb = $B->{"start"} + $helper->round(($B->{"end"} - $B->{"start"})/2);
	}
	
	($sa,$ea,$ma) = ($ca - $exu,$ca + $exd,$ca);
	($sb,$eb,$mb) = ($cb - $exu,$cb + $exd,$cb);
	
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

=head2 node2text

Collapse a node of a genomic interval tree to text. Internal use.

	$intersecter->node2text($itree_node);

=cut

sub node2text
{
	my ($self,$node,$delim) = @_;
	$delim = "\t" unless ($delim);
	return($node->{"start"}.$delim.$node->{"end"}.$delim.$node->{"value"}->{"rest"});
}

=head2 chrom_itree_size

Get the size of a chromosome hash of Interval Trees. Internal use.

	$intersecter->chrom_itree_size(\%chrom_tree_hash);

=cut

sub chrom_itree_size
{
	my ($self,$ctree) = @_;
	my ($c,$t);
	my %seen;
	my $size = 0;
	if ($self->get("reportonce"))
	{
		while(($c,$t) = each($ctree))
		{
			$t->traverse(sub {
				$size++ if (!$seen{$c}{$_[0]->{"interval"}->{"start"}.$_[0]->{"interval"}->{"end"}});
				$seen{$c}{$_[0]->{"interval"}->{"start"}.$_[0]->{"interval"}->{"end"}}++;
			});
		}
	}
	else
	{
		while(($c,$t) = each($ctree))
		{
			$t->traverse(sub { $size++; });
		}
	}
	return($size);
}

=head2 get_lengths

Get lengths of input genomic regions as hash of Interval Trees. Internal use.

	$intersecter->get_lengths(\%chrom_tree_hash);

=cut

sub get_lengths
{
	my ($self,$ctree) = @_;
	my @lens;
	my ($c,$t);
	while(($c,$t) = each($ctree))
	{
		$t->traverse(sub {
			push(@lens,$_[0]->{"interval"}->{"end"} - $_[0]->{"interval"}->{"start"});
		});
	}
	return(@lens);
}

=head2 strand_hash

Initiate a hash with strand representations. Internal use.

	$intersecter->strand_hash;

=cut

sub strand_hash
{
	my $self = shift @_;
	return(("+" => 1,"-" => -1,"1" => 1,"-1" => -1,"F" => 1,"R" => -1));
}

=head2 create_output_file

Create the name of the output file according to output type. Internal use.

	$intersecter->create_output_file($A,$B,$output_type);

=cut

sub create_output_file
{
	my ($self,$inA,$inB,$type) = @_;
	my ($baseA,$dirA,$extA) = fileparse($inA,'\.[^.]*');
	my $baseB = fileparse($inB,'\.[^.]*');

	($self->get("multi")) ? (return(File::Spec->catfile($dirA,$baseA.$baseB))) :
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
		if ($self->get("header"))
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
		if ($self->get("header")) # Cannot find a solution for direct sorting in Windows... sorry :-(
		{
			my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
					   "messing up data. Please sort files outside $MODNAME first (e.g. using\n".
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

	my $param_value = $count->get('param_name');

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

############################### LEGACY SUBROUTINES, USING BINARY SEARCH ############################

# Each binary search subroutine should be changed to accept as input the reference to the interval tree 
# corresponding to the chromosome of fileB, replacing @currvals with $chromosome->{$currchromtree} and 
# the rest of the inputs (starts, ends, modes etc.) as they are. The code performing binary search in 
# each function will be replaced by the $tree->find($start,$end) function and the results will be handled 
# as in ::Count with an additional for loop to handle multiple overlaps. The $npass will be removed of 
# course! 
# The whole idea must also change: the overlapA, overlapB, onlyA, onlyB hashes must become overlap trees.
# We begin by constructing an interval tree for fileB. If we want to run in batch mode (the first part of
# the script) we also construct an interval tree for fileA. During the run, the onlyA and onlyB interval
# trees are constructed dynamically and in the end we traverse and print the results.
# A temporary solution until we have a "remove" function in IntervalTree is that the search function returns
# also the node apart from 0 if there is no hit in the tree.

#sub bin_search_any
#{
	#my ($self,$start,$end,@areas) = @_;
	#my ($ind,$currstart,$currend);
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#($currstart,$currend) = split(/\t/,$areas[$ind]);
		#if (($currstart >= $start && $currend <= $end) ||
		    #($currstart <= $start && $currend >= $end) ||
			#($currstart < $start && $currend < $end && $start < $currend) ||
			#($currstart > $start && $currend > $end && $end > $currstart))
		#{
			#return($ind,$areas[$ind]);
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}

#sub bin_search_percent
#{
	#my ($self,$start,$end,$p,@areas) = @_;
	#my ($ind,$currstart,$currend);
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#($currstart,$currend) = split(/\t/,$areas[$ind]);
		#if (($currstart >= $start && $currend <= $end) ||
		    #($currstart <= $start && $currend >= $end))
		#{
			#return($ind,$areas[$ind]);
		#}
		#elsif ($currstart < $start && $currend < $end && $start < $currend)
		#{
			#(($currend - $start) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart > $start && $currend > $end && $end > $currstart)
		#{
			#(($end - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) :
			#(return($ind,0));
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}

#sub bin_search_any_center
#{
	#my ($self,$mode,$mpos,$dxval,$uxval,@areas) = @_;
	#my ($ind,$start,$end,$currstart,$currend,@arr);
	#$start = $mode - $dxval;
	#$end = $mode + $uxval;
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#@arr = split(/\t/,$areas[$ind]);
		#$currstart = $arr[$mpos] - $dxval;
		#$currend = $arr[$mpos] + $uxval;
		#if (($currstart >= $start && $currend <= $end) ||
		    #($currstart <= $start && $currend >= $end) ||
			#($currstart < $start && $currend < $end && $start < $currend) ||
			#($currstart > $start && $currend > $end && $end > $currstart))
		#{
			#return($ind,$areas[$ind]);
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}

#sub bin_search_percent_center
#{
	#my ($self,$mode,$mpos,$dxval,$uxval,$p,@areas) = @_;
	#my ($ind,$start,$end,$currstart,$currend,@arr);
	#$start = $mode - $dxval;
	#$end = $mode + $uxval;
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#@arr = split(/\t/,$areas[$ind]);
		#$currstart = $arr[$mpos] - $dxval;
		#$currend = $arr[$mpos] + $uxval;
		#if (($currstart >= $start && $currend <= $end) ||
		    #($currstart <= $start && $currend >= $end))
		#{
			#return($ind,$areas[$ind]);
		#}
		#elsif ($currstart < $start && $currend < $end && $start < $currend)
		#{
			#(($currend - $start) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart > $start && $currend > $end && $end > $currstart)
		#{
			#(($end - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) :
			#(return($ind,0));
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}

#sub bin_search_percent_both
#{
	#my ($self,$start,$end,$p,@areas) = @_;
	#my ($ind,$currstart,$currend,$diff);
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#($currstart,$currend) = split(/\t/,$areas[$ind]);
		#if (($currstart >= $start && $currend <= $end) ||
		    #($currstart <= $start && $currend >= $end))
		#{
			#return($ind,$areas[$ind]);
		#}
		#elsif ($currstart < $start && $currend < $end && $start < $currend)
		#{
			#$diff = $currend - $start;
			#($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			#(return($ind,$areas[$ind])) : (return($ind,0));
		#}
		#elsif ($currstart > $start && $currend > $end && $end > $currstart)
		#{
			#$diff = $end - $currstart;
			#($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			#(return($ind,$areas[$ind])) : (return($ind,0));
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}

#sub bin_search_percent_exact
#{
	#my ($self,$start,$end,$p,@areas) = @_;
	#my ($ind,$currstart,$currend);
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#($currstart,$currend) = split(/\t/,$areas[$ind]);
		#if ($currstart >= $start && $currend <= $end) # TagB <= TagA
		#{
			#(($currend - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart <= $start && $currend >= $end) # TagB >= TagA
		#{
			#(($end - $start) >= $p*($currend - $currstart)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart < $start && $currend < $end && $start < $currend)
		#{
			#(($currend - $start) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart > $start && $currend > $end && $end > $currstart)
		#{
			#(($end - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) :
			#(return($ind,0));
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}

#sub bin_search_percent_exact_both
#{
	#my ($self,$start,$end,$p,@areas) = @_;
	#my ($ind,$currstart,$currend,$diff);
	#my ($l,$u) = (0,$#areas);
	#$u = 1 if ($#areas == 0); # Kavourmadies...
	#while ($l <= $u)
	#{
		#$ind = int(($l + $u)/2);
		#($currstart,$currend) = split(/\t/,$areas[$ind]);
		#if ($currstart >= $start && $currend <= $end) # TagB <= TagA
		#{
			#(($currend - $currstart) >= $p*($end - $start)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart <= $start && $currend >= $end) # TagB >= TagA
		#{
			#(($end - $start) >= $p*($currend - $currstart)) ? (return($ind,$areas[$ind])) : 
			#(return($ind,0));
		#}
		#elsif ($currstart < $start && $currend < $end && $start < $currend)
		#{
			#$diff = $currend - $start;
			#($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			#(return($ind,$areas[$ind])) : (return($ind,0));
		#}
		#elsif ($currstart > $start && $currend > $end && $end > $currstart)
		#{
			#$diff = $end - $currstart;
			#($diff >= $p*($end - $start) || $diff >= $p*($currend - $currstart)) ? 
			#(return($ind,$areas[$ind])) : (return($ind,0));
		#}
		#else
		#{
			#$u = $ind - 1 if ($end <= $currstart);
            #$l = $ind + 1 if ($start >= $currend);
		#}
	#}
	#return(0);
#}
