#!/usr/bin/perl

# The script takes as input from 2 up to 5 input files and calculates all possible intersections,
# outputing a Venn diagram using the R package VennDiagram. For an explanation of all available
# options, see intersectbed.pl. The --ouptut option is suppressed because the only acceptable
# output for a multi-intersect is the intersection, that is the overlapA option in intersectbed.
# If you wish for more details, use intersectbed.pl directly. If you choose to use BEDTools
# intersectBed for speed, keep in mind that additional columns in the region files will be
# lost and the BEDTools have to be installed.

use strict;
use Getopt::Long;
use File::Copy;
use File::Spec;
use File::Temp;
use File::Basename;
use File::Path qw(make_path remove_tree);

use constant INTERSECT => "perl /media/HD4/Fleming/hts-tools/Perl/intersectbed.pl ";

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# On Ctrl-C or die, do cleanup
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "multisectbed.pl";
our @input; # Input files
our $name; # A project name
our $nodraw = 0; # Draw venn by default;
our $export =0; # Do not export individual areas by default
our $figtype = "pdf";
our $usebedtools;
our $bedtoolspath;
our $sort;
our $percent;
our $any;
our @extend;
our $mode;
our $autoxtend;
our $both;
our $exact;
our $npass;
our $agap;
our $outdir;
our $header;
our $waitbar;
our $silent;
our $help;

# Check inputs
&check_inputs;

# Global temp dir
our $tmpdir = File::Temp->newdir();

# Record progress...
my $date = &now;
disp("\n$date - Started...");
# Run the intersections pipeline...
&run_intersections;
# Run time...
$date = &now;
disp("\n$date - Finished!\n\n");

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@input,
			   "name|v=s" => \$name,
			   "nodraw|d" => \$nodraw,
			   "export|x" => \$export,
			   "figformat|f=s" => \$figtype,
			   "outdir|o=s" => \$outdir,
			   "usebedtools|b" => \$usebedtools,
			   "bedtoolspath|j" => \$bedtoolspath,
    		   "sort|r" => \$sort,
    		   "percent|p=f" => \$percent,
    		   "any|y" => \$any,
    		   "extend|e=i{,}" => \@extend,
    		   "mode|m=i" => \$mode,
    		   "autoextend|u" => \$autoxtend,
    		   "both|t" => \$both,
    		   "exact|c" => \$exact,
    		   "pass|n=i" => \$npass,
    		   "gap|g=i" => \$agap,
    		   "header|d" => \$header,
    		   "waitbar|w" => \$waitbar,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&program_usage;
    	exit;
    }
    $stop .= "--- Please specify input files ---\n" if (!@input);
	$stop .= "--- Number of input files must be between 2 and 5 ---\n" if (scalar @input < 2 || scalar @input > 5);
	if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if (!$outdir)
    {
		disp("Output directory not given! Assuming the current...");
		$outdir = ".";
	}
	if (!$name)
    {
		disp("Name for the Venn diagram figure not given! It will be auto-generated...");
		$name = (scalar @input)."_venn_".&now("machine");
	}
    $mode -= 2 if ($mode);
    if ($figtype ne "png" && $figtype ne "jpg" && $figtype ne "bmp" && $figtype ne "pdf" && $figtype ne "ps")
	{
		disp("--figformat must be one of png or pdf! Using pdf...");
		$figtype = "pdf";
	}
	if (!$bedtoolspath)
	{
		disp("--bedtoolspath not provided! Assuming /opt/NGSTools/BEDTools/bin...");
		$bedtoolspath = "/opt/NGSTools/BEDTools/bin";
	}
	# We do not check intersectbed.pl specific inputs, they are passed there and it checks
}

sub run_intersections
{
	my ($a,$b,$k,$s,$p,$cmd);
	
	# Create initial conditions for running
	my $alias = &copy_targets;
	my $pairs = &construct_run_pairs;
	my $optargs = &construct_optargs;

	# Run the actual intersections
	foreach $k (keys(%$pairs))
	{
		$a = File::Spec->catfile($tmpdir,${$pairs->{$k}}[0]);
		$b = File::Spec->catfile($tmpdir,${$pairs->{$k}}[1]);
		(!$usebedtools) ? ($cmd = INTERSECT."--inputA $a --inputB $b --multi$optargs") :
		($cmd = File::Spec->catfile($bedtoolspath,"intersectBed")." -a $a -b $b -wa > ".File::Spec->catfile($tmpdir,"${$pairs->{$k}}[0]"."${$pairs->{$k}}[1]"));
		disp("\nIntersecting ${$pairs->{$k}}[0] and ${$pairs->{$k}}[1]...");
		disp("The command is:");
		disp($cmd);
		system($cmd);
	}

	# Get the the areas and their length
	my $areas = &get_areas;
	my $counts = &count_areas($areas);
	if (!$nodraw)
	{
		$s = &create_venn($counts);
		$p = &run_R($s);
	}
	if ($export)
	{
		&export_areas($areas,$alias,$s,$p);
	}
}

sub run_R
{
	my $sf = $_[0];
	my @f;
	print "\n";
	my $status = system("Rscript --vanilla $sf");
	if ($status) # Dirty hack to recreate venn without overrideTriple
	{
		disp("Known problem in VennDiagram package detected... Re-running... The output will not be proportional...");
		open(IN,$sf);
		while(<IN>)
		{
			next if ($_ =~ m/override/);
			push(@f,$_);
		}
		close(IN);
		open(OUT,">$sf");
		print OUT @f;
		close(OUT);
		system("Rscript --vanilla $sf");
	}
	return(File::Spec->catfile($outdir,$name.".".$figtype));
}

sub export_areas
{
	use Archive::Tar;
	my ($areas,$alias,@filelist) = @_;

	# Create the file list with all areas
	use Switch;
	switch(scalar @input)
	{
		case 2 {
			push(@filelist,(
				File::Spec->catfile($tmpdir,$areas->{"area1"}),
				File::Spec->catfile($tmpdir,$areas->{"area2"}),
				File::Spec->catfile($tmpdir,$areas->{"cross.area"})
			));
		}
		case 3 {
			push(@filelist,(
				File::Spec->catfile($tmpdir,$areas->{"area1"}),
				File::Spec->catfile($tmpdir,$areas->{"area2"}),
				File::Spec->catfile($tmpdir,$areas->{"area3"}),
				File::Spec->catfile($tmpdir,$areas->{"n12"}),
				File::Spec->catfile($tmpdir,$areas->{"n23"}),
				File::Spec->catfile($tmpdir,$areas->{"n13"}),
				File::Spec->catfile($tmpdir,$areas->{"n123"})
			));
		}
		case 4 {
			push(@filelist,(
				File::Spec->catfile($tmpdir,$areas->{"area1"}),
				File::Spec->catfile($tmpdir,$areas->{"area2"}),
				File::Spec->catfile($tmpdir,$areas->{"area3"}),
				File::Spec->catfile($tmpdir,$areas->{"area4"}),
				File::Spec->catfile($tmpdir,$areas->{"area4"}),
				File::Spec->catfile($tmpdir,$areas->{"n12"}),
				File::Spec->catfile($tmpdir,$areas->{"n13"}),
				File::Spec->catfile($tmpdir,$areas->{"n14"}),
				File::Spec->catfile($tmpdir,$areas->{"n23"}),
				File::Spec->catfile($tmpdir,$areas->{"n24"}),
				File::Spec->catfile($tmpdir,$areas->{"n34"}),
				File::Spec->catfile($tmpdir,$areas->{"n123"}),
				File::Spec->catfile($tmpdir,$areas->{"n124"}),
				File::Spec->catfile($tmpdir,$areas->{"n134"}),
				File::Spec->catfile($tmpdir,$areas->{"n234"}),
				File::Spec->catfile($tmpdir,$areas->{"n1234"})
			));
		}
		case 5 {
			push(@filelist,(
				File::Spec->catfile($tmpdir,$areas->{"area1"}),
				File::Spec->catfile($tmpdir,$areas->{"area2"}),
				File::Spec->catfile($tmpdir,$areas->{"area3"}),
				File::Spec->catfile($tmpdir,$areas->{"area4"}),
				File::Spec->catfile($tmpdir,$areas->{"n12"}),
				File::Spec->catfile($tmpdir,$areas->{"n13"}),
				File::Spec->catfile($tmpdir,$areas->{"n14"}),
				File::Spec->catfile($tmpdir,$areas->{"n15"}),
				File::Spec->catfile($tmpdir,$areas->{"n23"}),
				File::Spec->catfile($tmpdir,$areas->{"n24"}),
				File::Spec->catfile($tmpdir,$areas->{"n25"}),
				File::Spec->catfile($tmpdir,$areas->{"n34"}),
				File::Spec->catfile($tmpdir,$areas->{"n35"}),
				File::Spec->catfile($tmpdir,$areas->{"n45"}),
				File::Spec->catfile($tmpdir,$areas->{"n123"}),
				File::Spec->catfile($tmpdir,$areas->{"n124"}),
				File::Spec->catfile($tmpdir,$areas->{"n125"}),
				File::Spec->catfile($tmpdir,$areas->{"n134"}),
				File::Spec->catfile($tmpdir,$areas->{"n135"}),
				File::Spec->catfile($tmpdir,$areas->{"n145"}),
				File::Spec->catfile($tmpdir,$areas->{"n234"}),
				File::Spec->catfile($tmpdir,$areas->{"n235"}),
				File::Spec->catfile($tmpdir,$areas->{"n245"}),
				File::Spec->catfile($tmpdir,$areas->{"n345"}),
				File::Spec->catfile($tmpdir,$areas->{"n1234"}),
				File::Spec->catfile($tmpdir,$areas->{"n1235"}),
				File::Spec->catfile($tmpdir,$areas->{"n1245"}),
				File::Spec->catfile($tmpdir,$areas->{"n1345"}),
				File::Spec->catfile($tmpdir,$areas->{"n2345"}),
				File::Spec->catfile($tmpdir,$areas->{"n12345"})
			));
		}
	}

	# Create a README map file
	my @short = ("A","B","C","D","E");
	my @letters = sort keys(%$alias);
	my $readme = "This file contains the association between the input files and the letters ".join(", ",@short[0..$#input])."\n\n";
	foreach my $letter (@letters[0..$#input])
	{
		$readme.= basename($letter).": ".basename($alias->{$letter})."\n";
	}
	$readme.= "\n";
	my $r = File::Spec->catfile($tmpdir,"README.txt");
	open(R,">$r");
	print R $readme;
	close(R);
	push(@filelist,$r);

	# Create the tar
	my $localdir = File::Spec->catdir($outdir,$name);
	mkdir($localdir);
	my @finalist;
	foreach my $f (@filelist)
	{
		copy($f,$localdir);
		push(@finalist,File::Spec->catfile($localdir,basename($f)));
	}
	my $archive = File::Spec->catfile($outdir,$name.".tar.gz");
	Archive::Tar->create_archive($archive,COMPRESS_GZIP,@finalist);
	# Remove the 2nd temporary directory and the original pdf file, since in the tar
	remove_tree($localdir);
	unlink(File::Spec->catfile($outdir,$name.".".$figtype));
	disp("All output areas stored in $archive!");
}

sub create_venn
{
	my ($counts,$otflag) = @_;
	my $script="require(VennDiagram)\n".&open_Rgraphics($outdir,$figtype);
	use Switch;
	switch(scalar @input)
	{
		case 2 {
			$script.=
				"v <- draw.pairwise.venn(\n".
				"\tarea1=$counts->{\"area1\"},\n".
				"\tarea2=$counts->{\"area2\"},\n".
				"\tcross.area=$counts->{\"cross.area\"},\n".
				"\tcategory=c(\"A\",\"B\"),\n".
				"\tfill=c(\"red\",\"green\"),\n".
				"\tlty=\"blank\",\n".
				"\tcex=1,\n".
				"\tcat.cex=2,\n".
				"\tcat.col=c(\"red\",\"green\"),\n".
				"\tcat.fontfamily=rep(\"Bookman\",2)\n".
			")\n";
		}
		case 3 {
			$script.=
				(($otflag) ? ("") : ("overrideTriple=TRUE\n")).
				"v <- draw.triple.venn(\n".
				"\tarea1=$counts->{\"area1\"},\n".
				"\tarea2=$counts->{\"area2\"},\n".
				"\tarea3=$counts->{\"area3\"},\n".
				"\tn12=$counts->{\"n12\"},\n".
				"\tn13=$counts->{\"n13\"},\n".
				"\tn23=$counts->{\"n23\"},\n".
				"\tn123=$counts->{\"n123\"},\n".
				"\tcategory=c(\"A\",\"B\",\"C\"),\n".
				"\tfill=c(\"red\",\"green\",\"blue\"),\n".
				"\tlty=\"blank\",\n".
				"\tcex=1,\n".
				"\tcat.cex=2,\n".
				"\tcat.col=c(\"red\",\"green\",\"blue\"),\n".
				"\tcat.fontfamily=rep(\"Bookman\",3)\n".
				")\n";
		}
		case 4 {
			$script.=
				"v <- draw.quad.venn(\n".
				"\tarea1=$counts->{\"area1\"},\n".
				"\tarea2=$counts->{\"area2\"},\n".
				"\tarea3=$counts->{\"area3\"},\n".
				"\tarea4=$counts->{\"area4\"},\n".
				"\tn12=$counts->{\"n12\"},\n".
				"\tn13=$counts->{\"n13\"},\n".
				"\tn14=$counts->{\"n14\"},\n".
				"\tn23=$counts->{\"n23\"},\n".
				"\tn24=$counts->{\"n24\"},\n".
				"\tn34=$counts->{\"n34\"},\n".
				"\tn123=$counts->{\"n123\"},\n".
				"\tn124=$counts->{\"n124\"},\n".
				"\tn134=$counts->{\"n134\"},\n".
				"\tn234=$counts->{\"n234\"},\n".
				"\tn1234=$counts->{\"n1234\"},\n".
				"\tcategory=c(\"A\",\"B\",\"C\",\"D\"),\n".
				"\tfill=c(\"red\",\"green\",\"blue\",\"orange\"),\n".
				"\tlty=\"blank\",\n".
				"\tcex=1,\n".
				"\tcat.cex=2,\n".
				"\tcat.col=c(\"red\",\"green\",\"blue\",\"orange\"),\n".
				"\tcat.fontfamily=rep(\"Bookman\",4)\n".
			")\n";
		}
		case 5 {
			$script.=
				"v <- draw.quintuple.venn(\n".
				"\tarea1=$counts->{\"area1\"},\n".
				"\tarea2=$counts->{\"area2\"},\n".
				"\tarea3=$counts->{\"area3\"},\n".
				"\tarea4=$counts->{\"area4\"},\n".
				"\tarea5=$counts->{\"area5\"},\n".
				"\tn12=$counts->{\"n12\"},\n".
				"\tn13=$counts->{\"n13\"},\n".
				"\tn14=$counts->{\"n14\"},\n".
				"\tn15=$counts->{\"n15\"},\n".
				"\tn23=$counts->{\"n23\"},\n".
				"\tn24=$counts->{\"n24\"},\n".
				"\tn25=$counts->{\"n25\"},\n".
				"\tn34=$counts->{\"n34\"},\n".
				"\tn35=$counts->{\"n35\"},\n".
				"\tn45=$counts->{\"n45\"},\n".
				"\tn123=$counts->{\"n123\"},\n".
				"\tn124=$counts->{\"n124\"},\n".
				"\tn125=$counts->{\"n125\"},\n".
				"\tn134=$counts->{\"n134\"},\n".
				"\tn135=$counts->{\"n135\"},\n".
				"\tn145=$counts->{\"n145\"},\n".
				"\tn234=$counts->{\"n234\"},\n".
				"\tn235=$counts->{\"n235\"},\n".
				"\tn245=$counts->{\"n245\"},\n".
				"\tn345=$counts->{\"n345\"},\n".
				"\tn1234=$counts->{\"n1234\"},\n".
				"\tn1235=$counts->{\"n1235\"},\n".
				"\tn1245=$counts->{\"n1245\"},\n".
				"\tn1345=$counts->{\"n1345\"},\n".
				"\tn2345=$counts->{\"n2345\"},\n".
				"\tn12345=$counts->{\"n12345\"},\n".
				"\tcategory=c(\"A\",\"B\",\"C\",\"D\",\"E\"),\n".
				"\tfill=c(\"red\",\"green\",\"blue\",\"orange\",\"grey50\"),\n".
				"\tlty=\"blank\",\n".
				"\tcex=1,\n".
				"\tcat.cex=2,\n".
				"\tcat.col=c(\"red\",\"green\",\"blue\",\"orange\",\"grey50\"),\n".
				"\tcat.fontfamily=rep(\"Bookman\",5)\n".
			"\t)\n";
		}
	}
	$script.= &close_Rgraphics;
	#return($script); # This for when I switch to Statistics::R
	my $o = File::Spec->catfile($tmpdir,$name.".R");
	open(O,">$o");
	print O $script;
	close(O);
	return($o);
}

sub construct_run_pairs
{
	use Tie::IxHash::Easy;
	my %uple;
	tie %uple, "Tie::IxHash::Easy";
	use Switch;
	switch(scalar @input)
	{
		case 2 {
			%uple = (
				1 => ["A","B"]
			);
		}
		case 3 {
			%uple = (
				1 => ["A","B"],
				2 => ["A","C"],
				3 => ["B","C"],
				4 => ["AB","C"]
			);
		}
		case 4 {
			%uple = (
				1 => ["A","B"],
				2 => ["A","C"],
				3 => ["A","D"],
				4 => ["B","C"],
				5 => ["B","D"],
				6 => ["C","D"],
				7 => ["AB","C"],
				8 => ["AB","D"],
				9 => ["AC","D"],
				10 => ["BC","D"],
				11 => ["ABC","D"]
			);
		}
		case 5 {
			%uple = (
				1 => ["A","B"],
				2 => ["A","C"],
				3 => ["A","D"],
				4 => ["A","E"],
				5 => ["B","C"],
				6 => ["B","D"],
				7 => ["B","E"],
				8 => ["C","D"],
				9 => ["C","E"],
				10 => ["D","E"],
				11 => ["AB","C"],
				12 => ["AB","D"],
				13 => ["AB","E"],
				14 => ["AC","D"],
				15 => ["AC","E"],
				16 => ["AD","E"],
				17 => ["BC","D"],
				18 => ["BC","E"],
				19 => ["BD","E"],
				20 => ["CD","E"],
				21 => ["ABC","D"],
				22 => ["ABC","E"],
				23 => ["ABD","E"],
				24 => ["ACD","E"],
				25 => ["BCD","E"],
				26 => ["ABCD","E"]
			);
		}
	}
	return(\%uple);
}

sub get_areas
{
	my %areas;
	use Switch;
	switch(scalar @input)
	{
		case 2 {
			%areas = (
				"area1" => "A",
				"area2" => "B",
				"cross.area" => "AB"
			);
		}
		case 3 {
			%areas = (
				"area1" => "A",
				"area2" => "B",
				"area3" => "C",
				"n12" => "AB",
				"n13" => "AC",
				"n23" => "BC",
				"n123" => "ABC"
			);
		}
		case 4 {
			%areas = (
				 "area1" => "A",
				 "area2" => "B",
				 "area3" => "C",
				 "area4" => "D",
				 "n12" => "AB",
				 "n13" => "AC",
				 "n14" => "AD",
				 "n23" => "BC",
				 "n24" => "BD",
				 "n34" => "CD",
				 "n123" => "ABC",
				 "n124" => "ABD",
				 "n134" => "ACD",
				 "n234" => "BCD",
				 "n1234" => "ABCD"
			);
		}
		case 5 {
			%areas = (
				 "area1" => "A",
				 "area2" => "B",
				 "area3" => "C",
				 "area4" => "D",
				 "area5" => "E",
				 "n12" => "AB",
				 "n13" => "AC",
				 "n14" => "AD",
				 "n15" => "AE",
				 "n23" => "BC",
				 "n24" => "BD",
				 "n25" => "BE",
				 "n34" => "CD",
				 "n35" => "CE",
				 "n45" => "DE",
				 "n123" => "ABC",
				 "n124" => "ABD",
				 "n125" => "ABE",
				 "n134" => "ACD",
				 "n135" => "ACE",
				 "n145" => "ADE",
				 "n234" => "BCD",
				 "n235" => "BCE",
				 "n245" => "BDE",
				 "n345" => "CDE",
				 "n1234" => "ABCD",
				 "n1235" => "ABCE",
				 "n1245" => "ABDE",
				 "n1345" => "ACDE",
				 "n2345" => "BCDE",
				 "n12345" => "ABCDE"
			);
		}
	}
	return(\%areas);
}

sub count_areas
{
	my $areas = $_[0];
	my %counts;
	use Switch;
	switch(scalar @input)
	{
		case 2 {
			%counts = (
				"area1" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area1"})),
				"area2" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area2"})),
				"cross.area" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"cross.area"})),
			);
		}
		case 3 {
			%counts = (
				"area1" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area1"})),
				"area2" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area2"})),
				"area3" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area3"})),
				"n12" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n12"})),
				"n23" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n23"})),
				"n13" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n13"})),
				"n123" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n123"}))
			);
		}
		case 4 {
			%counts = (
				"area1" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area1"})),
				"area2" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area2"})),
				"area3" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area3"})),
				"area4" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area4"})),
				"area5" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area4"})),
				"n12" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n12"})),
				"n13" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n13"})),
				"n14" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n14"})),
				"n23" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n23"})),
				"n24" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n24"})),
				"n34" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n34"})),
				"n123" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n123"})),
				"n124" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n124"})),
				"n134" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n134"})),
				"n234" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n234"})),
				"n1234" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n1234"}))
			)
		}
		case 5 {
			%counts = (
				"area1" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area1"})),
				"area2" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area2"})),
				"area3" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area3"})),
				"area4" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area4"})),
				"area5" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"area5"})),
				"n12" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n12"})),
				"n13" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n13"})),
				"n14" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n14"})),
				"n15" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n15"})),
				"n23" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n23"})),
				"n24" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n24"})),
				"n25" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n25"})),
				"n34" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n34"})),
				"n35" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n35"})),
				"n45" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n45"})),
				"n123" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n123"})),
				"n124" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n124"})),
				"n125" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n125"})),
				"n134" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n134"})),
				"n135" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n135"})),
				"n145" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n145"})),
				"n234" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n234"})),
				"n235" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n235"})),
				"n245" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n245"})),
				"n345" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n345"})),
				"n1234" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n1234"})),
				"n1235" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n1235"})),
				"n1245" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n1245"})),
				"n1345" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n1345"})),
				"n2345" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n2345"})),
				"n12345" => &count_lines(File::Spec->catfile($tmpdir,$areas->{"n12345"}))
			);
		}
	}
	if ($header)
	{
		foreach my $k (keys(%counts))
		{
			$counts{$k}--;
		}
	}
	return(\%counts);
}

sub count_lines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	disp("Counting lines of $_[0]...");
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
}

sub construct_optargs
{
	my $str;
	$str.= " --sort" if ($sort);
	$str.= " --percent $percent" if ($percent);
	$str.= " --any" if ($any);
	$str.= " --extend $extend[0] $extend[1]" if (@extend);
	$str.= " --mode $mode" if ($mode);
	$str.= " --autoextend" if ($autoxtend);
	$str.= " --both" if ($both);
	$str.= " --exact" if ($exact);
	$str.= " --pass $npass" if ($npass);
	$str.= " --gap $agap" if ($agap);
	$str.= " --header" if ($header);
	$str.= " --waitbar" if ($waitbar);
	$str.= " --silent" if ($silent);
	$str.= " --output overlapA";	
	return($str);
}

sub copy_targets
{
	my %alias;
	my @short = ("A","B","C","D","E");
	for my $i (0..$#input)
	{
		my $f = File::Spec->catfile($tmpdir,$short[$i]);
		$alias{$f} = $input[$i];
		($usebedtools) ? (&format_for_bedtools($input[$i],$f)) : (copy($input[$i],$f));
	}
	return(\%alias);
}

sub format_for_bedtools
{
	my ($input,$target) = @_;
	my $line;
	my @cols;
	open(IN,$input);
	open(OUT,">$target");
	my $line = <IN> if ($header); # Header to oblivion
	while ($line = <IN>)
	{
		$line =~ s/\r|\n$//g;
		@cols = split("\t",$line);
		print OUT "$cols[0]\t$cols[1]\t$cols[2]\n";
	}
	close(IN);
	close(OUT);
}

sub open_Rgraphics
{
	my ($outdir,$format) = @_;
	my $os;
	use Switch;
	switch($format)
	{
		case /png/i
		{
			$os = "png(file=file.path(\"$outdir\",paste(\"$name\",\".png\",sep=\"\")))\n";
		}
		case /jpg/i
		{
			$os = "jpg(file=file.path(\"$outdir\",paste(\"$name\",\".jpg\",sep=\"\")))\n";
		}
		case /bmp/i
		{
			$os = "bmp(file=file.path(\"$outdir\",paste(\"$name\",\".bmp\",sep=\"\")))\n";
		}
		case /pdf/i
		{
			$os = "pdf(file=file.path(\"$outdir\",paste(\"$name\",\".pdf\",sep=\"\")))\n";
		}
		case /ps/i
		{
			$os = "postscript(file=file.path(\"$outdir\",paste(\"$name\",\".ps\",sep=\"\")))\n";
		}
	}
	return($os);
}

sub close_Rgraphics
{
	return("dev.off()\n");
}

sub now
{
	my $format = shift @_;
	$format = "human" if (!$format);
	my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
	$year += 1900;
	$month++;
	$month = "0".$month if (length($month)==1);
	$day = "0".$day if (length($day)==1);
	$hour = "0".$hour if (length($hour)==1);
	$min = "0".$min if (length($min)==1);
	$sec = "0".$sec if (length($sec)==1);
	($format ne "machine") ? (return($day."/".$month."/".$year." ".$hour.":".$min.":".$sec)) :
	(return($year.$month.$day.$hour.$min.$sec));
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

sub disp
{
	print "\n@_" if (!$silent);
}

sub program_usage 
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";
	
$scriptname
A perl program to calculate multiple intersections between bed files and output
them as well as a Venn diagram. Supports up to 5 different region files (because
of the current feasible visualization of Venn diagrams. The program can optionally
output all the areas and subareas together with the Venn diagram, an R script to
reproduce it and a README map file between original inputs and the output coding.
The program accepts all the intersectbed.pl options which are not explained below.
For a list of them, please see the intersectbed.pl help. If you are not interested
in the contents of the sub-areas but want just a Venn diagram, the program can use
the much faster intersectBed from BEDTools to calculate the numbers, if BEDTools
are installed on your system. This is not the default behavior and can be switched
on. See the options below.

Author : Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname --input fileA fileB [fileC fileD fileE] [OPTIONS]

--- Required ---
  --input|a  files  At least 2 region files and up to 5
--- Optional ---
  --name|v		Use this option to provide a name for the output
  			Venn diagram. it will be autogenerated if none provided.
  --nodraw|d		Do not draw a Venn diagram (why?). Just output the
			area sizes. Defaults to printing the diagram.
  --export|x		Create a name.tar.gz file containing all the output of
			the program (intersection areas, venn diagram and README file).
  --figformat|f		Use this option to tell the program the wished format
			for the output diagram. Can be one of png, jpg, bmp, pdf or ps.
			Defaults to "pdf".
  --usebedtools		Use this option to switch on the use of intersectBed
			from BEDTools. Τhe sub-areas however, will be simple BED3 files
			and any additional data columns will be lost. BEDTools must be
			installed on your system. This implementation is a lot(!) faster
			at the cost of some data loss.
  --bedtoolspath		The path to BEDTools, e.g. /usr/share/bedtools/bin.
			If not provided, a default path is assumed, which might not be
			appropriate for your system.
  --outdir|o	Use this option to provide an output directory for all
			output produced by the program. Defaults to "." (current).
  --help|h		Display this help text.

Usage examples:

perl multisectbed.pl --input A.peaks B.peaks C.peaks --name MyCommonPeaks --export

Package dependecies:
	intersectbed.pl from hts-tools (optional if BEDTools are installed)
	BEDTools (optional if intersectbed.pl is present)
	At least ONE of the above MUST be present.

END
	print $usagetext;
	exit;
}
