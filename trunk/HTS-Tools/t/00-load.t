#!perl -T
use v5.10;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 11;

BEGIN {
    use_ok( 'HTS::Tools' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Count' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Assign' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Utils' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Normalize::Bedgraph' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Intersect' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Profile' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::QC' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Motifscan' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Convert' ) || print "Bail out!\n";
    use_ok( 'HTS::Tools::Fetch' ) || print "Bail out!\n";
}

diag( "Testing HTS::Tools $HTS::Tools::VERSION, Perl $], $^X" );
