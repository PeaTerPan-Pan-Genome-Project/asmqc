#!/usr/bin/perl
# Replacement for src/search_dep0.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same argument,
# byte-identical output.
#
# The original stores every zero-depth position of the genome in a nested hash
# (hundreds of millions for a pea genome), sorts each chromosome's positions
# and reports runs of consecutive zero-depth positions spanning >= 150 bp
# (start "+", end "-"). This version collects the runs while streaming and
# keeps only their ends. Output order as in the original: sequences sorted as
# strings, runs by position. Runs are merged after sorting, so the result does
# not depend on the depth file being sorted.
use strict;

my ($infile) = @ARGV;
my (%starts, %ends);    # chr -> run starts / ends, in reading order
my ($cc, $s, $e);

sub close_run { if (defined $cc) { push @{ $starts{$cc} }, $s; push @{ $ends{$cc} }, $e } }

open my $in, '<', $infile or die "$infile: $!";
while (my $l = <$in>) {
    my $t1 = index($l, "\t");
    my $t2 = index($l, "\t", $t1 + 1);
    next if $t1 < 0 || $t2 < 0;
    my $t3 = index($l, "\t", $t2 + 1);
    my $dep = $t3 < 0 ? substr($l, $t2 + 1) : substr($l, $t2 + 1, $t3 - $t2 - 1);
    chomp $dep;
    next unless $dep == 0;
    my ($chr, $pos) = (substr($l, 0, $t1), substr($l, $t1 + 1, $t2 - $t1 - 1));
    if (defined $cc && $chr eq $cc && $pos == $e + 1) { $e = $pos; next }
    if (defined $cc && $chr eq $cc && $pos == $e) { next }
    close_run();
    ($cc, $s, $e) = ($chr, $pos, $pos);
}
close_run();
close $in;

for my $chr (sort keys %starts) {
    my ($st, $en) = ($starts{$chr}, $ends{$chr});
    my @o = sort { $st->[$a] <=> $st->[$b] } 0 .. $#$st;
    my ($rs, $re);
    for my $k (@o, -1) {    # -1: flush the last run
        if ($k >= 0 && defined $re && $st->[$k] <= $re + 1) {
            $re = $en->[$k] if $en->[$k] > $re;
            next;
        }
        if (defined $re && $re > $rs && $re - $rs >= 150) {
            print $chr . "\t" . $rs . "\t+\t0" . "\n";
            print $chr . "\t" . $re . "\t-\t0" . "\n";
        }
        ($rs, $re) = $k >= 0 ? ($st->[$k], $en->[$k]) : ();
    }
}
