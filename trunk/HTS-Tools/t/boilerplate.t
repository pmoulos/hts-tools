#!perl -T
use v5.10;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 15;

sub not_in_file_ok {
    my ($filename, %regex) = @_;
    open( my $fh, '<', $filename )
        or die "couldn't open $filename for reading: $!";

    my %violated;

    while (my $line = <$fh>) {
        while (my ($desc, $regex) = each %regex) {
            if ($line =~ $regex) {
                push @{$violated{$desc}||=[]}, $.;
            }
        }
    }

    if (%violated) {
        fail("$filename contains boilerplate text");
        diag "$_ appears on lines @{$violated{$_}}" for keys %violated;
    } else {
        pass("$filename contains no boilerplate text");
    }
}

sub module_boilerplate_ok {
    my ($module) = @_;
    not_in_file_ok($module =>
        'the great new $MODULENAME'   => qr/ - The great new /,
        'boilerplate description'     => qr/Quick summary of what the module/,
        'stub function definition'    => qr/function[12]/,
    );
}

TODO: {
  local $TODO = "Need to replace the boilerplate text";

  not_in_file_ok(README =>
    "The README is used..."       => qr/The README is used/,
    "'version information here'"  => qr/to provide version information/,
  );

  not_in_file_ok(Changes =>
    "placeholder date/time"       => qr(Date/time)
  );

  module_boilerplate_ok('lib/HTS/Tools.pm');
  module_boilerplate_ok('lib/HTS/Tools/Count.pm');
  module_boilerplate_ok('lib/HTS/Tools/Assign.pm');
  module_boilerplate_ok('lib/HTS/Tools/Utils.pm');
  module_boilerplate_ok('lib/HTS/Tools/Normalize.pm');
  module_boilerplate_ok('lib/HTS/Tools/Normalize/Bed.pm');
  module_boilerplate_ok('lib/HTS/Tools/Normalize/Bedgraph.pm');
  module_boilerplate_ok('lib/HTS/Tools/Intersect.pm');
  module_boilerplate_ok('lib/HTS/Tools/Profile.pm');
  module_boilerplate_ok('lib/HTS/Tools/QC.pm');
  module_boilerplate_ok('lib/HTS/Tools/Motifscan.pm');
  module_boilerplate_ok('lib/HTS/Tools/Convert.pm');
  module_boilerplate_ok('lib/HTS/Tools/Fetch.pm');


}

