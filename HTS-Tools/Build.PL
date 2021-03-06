use v5.10;
use strict;
use warnings FATAL => 'all';
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'HTS::Tools',
    license             => 'artistic_2',
    dist_author         => q{Panagiotis Moulos <moulos@fleming.gr>},
    dist_version_from   => 'lib/HTS/Tools.pm',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => 0,
    },
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        #'ABC'              => 1.6,
        #'Foo::Bar::Module' => 5.0401,
    },
    add_to_cleanup     => [ 'HTS-Tools-*' ],
    create_makefile_pl => 'traditional',
);

$builder->add_build_element('yml');
$builder->create_build_script();
