#!/usr/bin/perl
# Replacement for src/LReffect_size.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output (stdout and LRout/Nonmap.loc).
#
# Splits each depth line with index/substr instead of split /\s+/, about twice
# as fast on samtools depth output (tab-separated: sequence, position, depth).
# An optional fourth argument replaces the fixed LRout/Nonmap.loc (used when the
# genome is processed in segments).
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

my ($lrdep_file, $avgdep, $depratio, $nonmap) = @ARGV;
$nonmap //= "LRout/Nonmap.loc";
my $low = int($avgdep * $depratio) + 1;
open my $in, '<', ($lrdep_file eq '-' ? '/dev/stdin' : $lrdep_file) or die "$lrdep_file: $!";
open my $out, '>', $nonmap or die "$nonmap: $!";
my %depmin;
while (my $l = <$in>) {
    my $t1 = index($l, "\t");
    my $t2 = index($l, "\t", $t1 + 1);
    my $t3 = index($l, "\t", $t2 + 1);
    my $dep = $t3 < 0 ? substr($l, $t2 + 1) : substr($l, $t2 + 1, $t3 - $t2 - 1);
    chomp $dep;
    if ($dep < $low) { chomp $l; print $out "$l\n" }
    $depmin{ substr($l, 0, $t1) }++ if $dep >= 2;
}
for (sort keys %depmin) { print "$_\t$depmin{$_}\n" }
close $out;
