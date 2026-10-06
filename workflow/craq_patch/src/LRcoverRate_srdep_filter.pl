#!/usr/bin/perl
# Replacement for src/LRcoverRate_srdep_filter.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output; based on the rewrite in asmqc issue #2.
#
# The original puts every zero-depth position of the genome into one hash
# (hundreds of millions of keys for a pea genome) and asks it about +-1500 bp
# around a few dozen long-read candidates. This version reads the candidates
# first and keeps only the zero-depth positions inside a candidate window.
# Assumes samtools depth order: positions ascending within a chromosome block.
use strict;

if (@ARGV != 2) { print "perl $0  SR_dep_file LR_cover_file \n"; exit 1; }
my ($sr_dep, $lrcover_rate) = @ARGV;
my $W = 1500;

my (@cand, %iv);
open my $cr, '<', $lrcover_rate or die "$lrcover_rate: $!";
while (my $l = <$cr>) {
    chomp $l;
    my ($chr, $pos) = (split /\s+/, $l)[0, 1];
    push @cand, [ $l, $chr, $pos ];
    push @{ $iv{$chr} }, [ $pos - $W, $pos + $W ];
}
close $cr;

# each chromosome's windows as sorted, disjoint intervals
for my $chr (keys %iv) {
    my @s = sort { $a->[0] <=> $b->[0] } @{ $iv{$chr} };
    my @m = ([ @{ $s[0] } ]);
    for my $x (@s[1 .. $#s]) {
        if ($x->[0] <= $m[-1][1] + 1) { $m[-1][1] = $x->[1] if $x->[1] > $m[-1][1] }
        else { push @m, [@$x] }
    }
    $iv{$chr} = \@m;
}

my (%dep0, $cur, $ivs, $j);
open my $sd, '<', $sr_dep or die "$sr_dep: $!";
while (my $l = <$sd>) {
    my $t = index($l, "\t");
    my $chr = substr($l, 0, $t);
    if (!defined $cur || $chr ne $cur) { $cur = $chr; $ivs = $iv{$chr}; $j = 0 }
    next unless $ivs && $j < @$ivs;
    my $u = index($l, "\t", $t + 1);
    my $pos = substr($l, $t + 1, $u - $t - 1);
    $j++ while $j < @$ivs && $ivs->[$j][1] < $pos;
    next if $j >= @$ivs || $pos < $ivs->[$j][0];
    my $dep = (split /\s+/, $l)[2];
    $dep0{"$chr--$pos"} = 0 if $dep == 0;
}
close $sd;

for my $cd (@cand) {
    my ($line, $chr, $pos) = @$cd;
    my ($leftnum0, $rightnum0) = (0, 0);
    for my $i (1 .. $W) {
        my $left = $pos - $i;
        my $right = $pos + $i;
        if (defined $dep0{"$chr--$left"}) {
            $leftnum0++;
            if ($leftnum0 >= 10 && $leftnum0 * 10 > $i) { print "$line\n"; last }
        }
        if (defined $dep0{"$chr--$right"}) {
            $rightnum0++;
            if ($rightnum0 >= 10 && $rightnum0 * 10 > $i) { print "$line\n"; last }
        }
    }
}
