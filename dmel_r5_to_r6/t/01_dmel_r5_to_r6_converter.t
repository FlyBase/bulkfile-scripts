#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Capture::Tiny ':all';

use Test::More tests => 20;

my $script = "$FindBin::Bin/../dmel_r5_to_r6_converter.pl";

ok(-e $script,'Script exists');

my %coordinates = (
    '3RHet:1936890'          => '3R:3585370',
    '3RHet:1964215'          => '3R:3523947',
    '3RHet:1964215..1936890' => '3R:3523947..3585370',
    '2R:1..10000'            => '2R:4112496..4122495',
    'ChrX:330..19020'        => 'X:103943..122633',
    '4:43943..50320'         => '4:23317..29694',
    '3R:3,020,320'           => '3R:7194598',
    'chr4:123020'            => '4:102394',
    '3RHet:1894569..1895349' => '3R:3626910..3627690',
);


while (my ($r5,$r6) = each %coordinates) {
    my ($stdout, $stderr, $exit) = capture {
        open(my $fh, '|-',$script) or die "Can't open converter script."; 
        print $fh "$r5\n";
        close($fh);
    };
    my @lines = grep { $_ !~ /^#/ } split(/\n/,$stdout);

    is(scalar @lines, 1, 'One result for: ' . $r5);
    
    my ($script_r5,$script_r6,$notes) = split(/\t/,$lines[0]);
    is($script_r6,$r6,'Validating conversion for: '. $r5);
}

{
    my ($stdout, $stderr, $exit) = capture {
        open(my $fh, '|-',$script) or die "Can't open converter script."; 
        while (my ($r5,$r6) = each %coordinates) {
            print $fh "$r5\n";
        }
        close($fh);
    };
    my @lines = grep { $_ !~ /^#/ } split(/\n/,$stdout);

    is(scalar @lines,scalar keys %coordinates,'checking bulk run output count');
}
