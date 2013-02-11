#!/usr/bin/perl

# TODO: Incorporate a step to also check the strands where each read is mapped to. Reads
# 		with same co-ordinates but with different strands should be kept.

use strict;
use Getopt::Long;
use Pod::Usage;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "set2petset.pl";
our $input;
our $output;
our $score;
our $tagsize;
our $dovetol = 10;
our $distol = 500;
our $insert = 0;
our $man = 0;
our $silent = 0;
our $help = 0;

&check_inputs;

my %namehash;
my %isPaired;
my $isize;
my $ddiscordant = 0;
my $sdiscordant = 0;

($input eq "-" || $input eq "stdin") ? (open(INPUT,"-")) : (open(INPUT,$input));
my $out;
if ($output)
{
	open(OUTPUT,">$output");
	$out = *OUTPUT;
}
else { $out = *STDOUT; }

if ($insert)
{
	# Make a simple check to see if re-pairing part of set2petset has been run
	my $reprun = 0;
	my $checkcounter = 0;
	while ((my $line = <INPUT>) && ($checkcounter < 50))
	{
		$checkcounter++;
		#$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		$reprun = 1 if ($cols[3] =~ m/_(P|U)/i);
	}
	die "\nYou must run $scriptname on collapsed PETs first (directly from a BAM file with e.g. bamToBed)!\n\n" if (!$reprun);

	seek(INPUT,0,0);
	while (my $line = <INPUT>)
	{
		$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		my ($chr,$start,$end) = ($cols[0],$cols[1],$cols[2]);
		if ($cols[3] =~ m/_P/i) # Paired, we must split and create the insert
		{
			my $name = $cols[3];
			my $lstart = $start;
			my $lend = $start + $tagsize - 1;
			my $lname = $name;
			$lname =~ s/_P/_I1/;
			my $rstart = $end - $tagsize + 1;
			my $rend = $end;
			my $rname = $name;
			$rname =~ s/_P/_I3/;
			print "$chr\t$lstart\t$lend\t$lname\t$cols[4]\t$cols[5]\n";
			print "$chr\t$rstart\t$rend\t$rname\t$cols[4]\t$cols[5]\n";
			# Now... the insert is possibly dovetailed... we must check this and NOT return it...
			my $istart = $start + $tagsize;
			my $iend = $end - $tagsize;
			my $iname = $name;
			$iname =~ s/_P/_I2/;
			print $out "$chr\t$istart\t$iend\t$iname\t$cols[4]\t$cols[5]\n"
				unless ($iend < $istart || abs($iend - $istart) <= $dovetol); # Dovetail or too close!
		}
		else
		{
			print join("\t",@cols),"\n";
		}
	}
}
else # Just re-pair
{
	while (my $line = <INPUT>)
	{
		$line =~ s/\r|\n$//g;
		my @cols = split(/\t/,$line);
		next if ($cols[4] < $score || $cols[0] =~ m/chrM|rand|chrU|hap/);
		my ($chr,$start,$end,$score,$strand) = ($cols[0],$cols[1],$cols[2],$cols[4],$cols[5]);
		my $mate = substr($cols[3],-1,length($cols[3]));
		my $seqname = substr($cols[3],0,-2);
		if (!$namehash{$seqname})
		{
			$namehash{$seqname}{"chr"} = $chr;
			$namehash{$seqname}{"start"} = $start;
			$namehash{$seqname}{"end"} = $end;
			$namehash{$seqname}{"strand"} = $strand;
			$namehash{$seqname}{"score"} = $score;
			$isPaired{$seqname} = 0;
		}
		else
		{
			if ($chr eq $namehash{$seqname}{"chr"})
			{
				($mate == 2) ? ($isize = abs($start - $namehash{$seqname}{"end"})) :
					($isize = abs($end - $namehash{$seqname}{"start"}));
				if ($isize <= $distol) # Concordant pair
				{
					if ($start < $namehash{$seqname}{"start"})
					{
						$namehash{$seqname}{"start"} = $start;
					}
					if ($end > $namehash{$seqname}{"end"})
					{
						$namehash{$seqname}{"end"} = $end;
					}
					$namehash{$seqname}{"strand"} = "+";
					$isPaired{$seqname} = 1;
				}
				else # Discordant pair, same chromosome
				{
					$namehash{"Discordant_SameChr_$sdiscordant"}{"chr"} = $chr;
					$namehash{"Discordant_SameChr_$sdiscordant"}{"start"} = $start;
					$namehash{"Discordant_SameChr_$sdiscordant"}{"end"} = $end;
					$namehash{"Discordant_SameChr_$sdiscordant"}{"strand"} = $strand;
					$namehash{"Discordant_SameChr_$sdiscordant"}{"score"} = $score;
					$isPaired{"Discordant_SameChr_$sdiscordant"} = 0;
					$sdiscordant++;
				}
			}
			else # Discordant pair, different chromosoe
			{
				$namehash{"Discordant_DiffChr_$ddiscordant"}{"chr"} = $chr;
				$namehash{"Discordant_DiffChr_$ddiscordant"}{"start"} = $start;
				$namehash{"Discordant_DiffChr_$ddiscordant"}{"end"} = $end;
				$namehash{"Discordant_DiffChr_$ddiscordant"}{"strand"} = $strand;
				$namehash{"Discordant_DiffChr_$ddiscordant"}{"score"} = $score;
				$isPaired{"Discordant_DiffChr_$ddiscordant"} = 0;
				$ddiscordant++;
			}
		}
	}
	close(INPUT);

	# Print the case of simple re-pairing
	my $pflag;
	foreach my $name (keys(%namehash))
	{
		($isPaired{$name}) ? ($pflag = "_P") : ($pflag = "_U");
		print $out "$namehash{$name}{\"chr\"}\t".
		"$namehash{$name}{\"start\"}\t$namehash{$name}{\"end\"}\t".
		"$name$pflag\t$namehash{$name}{\"score\"}\t$namehash{$name}{\"strand\"}\n";
	}
}

print STDERR "\n";

#my $pflag = "_U";
#if ($insert)
#{
	#foreach my $name (keys(%namehash))
	#{
		#if ($isPaired{$name})
		#{
			#my $chr = $namehash{$name}{"chr"};
			#my $strand = $namehash{$name}{"strand"};
			#my $score = $namehash{$name}{"score"};
			
			#my $lstart = $namehash{$name}{"start"};
			#my $lend = $namehash{$name}{"start"} + $tagsize - 1;
			#my $lname = $name."_I1";
			#my $rstart = $namehash{$name}{"end"} - $tagsize + 1;
			#my $rend = $namehash{$name}{"end"};
			#my $rname = $name."_I3";
			#my $istart = $namehash{$name}{"start"} + $tagsize;
			#my $iend = $namehash{$name}{"end"} - $tagsize;
			#my $iname = $name."_I2";
			
			#print $out "$chr\t$lstart\t$lend\t$lname\t$score\t$strand\n";
			#print $out "$chr\t$rstart\t$rend\t$rname\t$score\t$strand\n";
			## The insert can be possibly dovetailed... we must check this and NOT return it...
			#print $out "$chr\t$istart\t$iend\t$iname\t$score\t$strand\n"
				#unless ($iend < $istart || abs($iend - $istart) <= $dovetol); # Dovetail or too close!
		#}
		#else
		#{
			#print $out "$namehash{$name}{\"chr\"}\t".
			#"$namehash{$name}{\"start\"}\t$namehash{$name}{\"end\"}\t".
			#"$name$pflag\t$namehash{$name}{\"score\"}\t$namehash{$name}{\"strand\"}\n";
		#}
	#}
#}
#else
#{
	#foreach my $name (keys(%namehash))
	#{
		#($isPaired{$name}) ? ($pflag = "_P") : ($pflag = "_U");
		#print $out "$namehash{$name}{\"chr\"}\t".
		#"$namehash{$name}{\"start\"}\t$namehash{$name}{\"end\"}\t".
		#"$name$pflag\t$namehash{$name}{\"score\"}\t$namehash{$name}{\"strand\"}\n";
	#}
#}

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("input|i=s" => \$input,
			   "score|c=i" => \$score,
			   "tagsize|t=i" => \$tagsize,
			   "insert|r" => \$insert,
			   "distol|n" => \$distol,
			   "dovetol|d" => \$dovetol,
			   "man|m" => \$man,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
	Pod::Usage::pod2usage( -verbose => 1, -exitstatus => 0 ) if ($help);
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ($man);
    $stop .= "--- Please provide an input BED file ---\n" if (!$input);
    $stop .= "--- Please provide a tagsize with the --insert option ---\n" if ($insert && !$tagsize);
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
	if (!$score)
	{
		disp("Quality mapping filter score not given! Will use the default (10)...");
		$score = 10;
	}
	if (!$distol)
	{
		disp("Allowed distance between pairs not given! Will use the default (500)...");
		$distol = 500;
	}
	if (!$dovetol)
	{
		disp("Allowing dovetailed pair length (or minimum insert size) not given! Will use the default (10)...");
		$dovetol = 10;
	}
}

sub zeroFill
{
	my ($n,$len) = @_;
	my $diff = $len - length($n);
	return($n) if $diff <= 0;
	return(('0'x$diff).$n);
}

sub disp
{
	print STDERR "\n@_" if (!$silent);
}

#if ($namehash{$seqname}{"strand"} eq "+")
#{
	#$namehash{$seqname}{"end"} = $end;
	#$namehash{$seqname}{"score"} = $score if ($score > $namehash{$seqname}{"score"});
#}
#elsif ($namehash{$seqname}{"strand"} eq "-")
#{
	#$namehash{$seqname}{"start"} = $start;
	#$namehash{$seqname}{"score"} = $score if ($score > $namehash{$seqname}{"score"});
	#$namehash{$seqname}{"strand"} = "+";
#}

__END__

=pod

=head1 NAME

set2petset.pl - Re-pair split paired-end bed files and get insert sequences

=head1 SYNOPSIS

set2petset.pl --input bedfile OR - [--score phred_filter] [--tagsize tagsize] [--distol distance_between_pairs] [--dovetol minimum_insert_size] [--insert] [--silent] [--help]

Examples:

=over 4

=item perl set2petset.pl --input collapsed_pet.bed --output repaired.bed

=item bamToBed -i bamfile.bam | perl set2petset.pl -i - -t 35 -c 8 -n 400 -d 20 -r > repaired_with_insert.bed

=back

=head1 DESCRIPTION

A perl program to re-pair paired-end reads in bed format, applying phred filter scores when present,
categorizing the reads to Paired (P), Unpaired (U) or Discordant (Discordant), and also allowing the
use to retrieve the insert reads between a pair of reads given some constraints (e.g. minimum insert
size or maximum distance between pairs). Note that set2petset serves two functionalities, defined by
the combination of parameters. If --insert is not chosen, the program will return a list of re-paired
(very long reads) with their insert size. If executed with --insert and --tagsize, it will return a
set of single-end reads characterized according to their position and also the insert sequence.

=head1 ARGUMENTS

=over 4

=item input B<(required)>

--input or -i

An input bed file or - for reading from stdin

=item score B<(optional)>

--score or -c

The phred score filter to use. Must be present on the 5th bed column. Defaults to 10.

=item insert B<(optional)>

--insert or -r

A switch to indicate that the insert size should be returned with the output. Defaults to off.

=item tagsize B<(required with --insert)>

--tagsize or -t

The read length in the single-end reads input bed file.

=item distol B<(optional)>

--distol or -n

An integer allowed distance between paired-end reads so that they are considered as a proper pair. Defaults to 500.

=item dovetol B<(optional)>

--dovetol or -d

The minimum insert read length. Useful when reads have been mapped too close (<10bp) or are dovetailed. Defaults to 10.

=item silent B<optional>

--silent or -s

Do not display verbose messages.

=item help B<(optional)>

--help or -h

Display this help text.

=item man B<(optional)>

--man or -m

Display the full manual of the script.

=back

=head1 AUTHOR

Panagiotis Moulos (L<moulos@fleming.gr)>

=cut
