#!/usr/bin/perl
# Replacement for src/synthesize_LRbkdep_and_alldep.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output.
#
# The original keeps every per-base depth of a chromosome in a nested hash and
# reads only +-10 bp around each clip site. This version streams the depth
# file once with a ring buffer and evaluates a site as soon as its window has
# been read; lines outside every window are skipped. Kept on purpose, as in the original:
#  - the first depth line of every chromosome is not stored;
#  - a chromosome with no stored line yields nothing;
#  - a site with no position defined on both sides yields nothing ("next"
#    inside avg());
#  - output per chromosome in depth-file order, sites in input order.
# Assumes samtools depth order: each chromosome in one block, positions
# ascending.
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

if (@ARGV != 2) { print "USE: $0 LR_break.depth LR_sort.depth\n"; exit 1; }
my ($bkfile, $depfile) = @ARGV;
my $span = 10;
my $MASK = 63;

my %cand;    # chr -> [[chr, pos, strand, bkdep], ...] in input order
open my $bk, '<', ($bkfile eq '-' ? '/dev/stdin' : $bkfile) or die "$bkfile: $!";
while (my $l = <$bk>) {
    chomp $l;
    my @arr = split /\t/, $l;
    my ($chr, $pos, $stran, $bkdep) = (split /===/, join("===", @arr))[0, 1, 2, 3];
    push @{ $cand{ $arr[0] } }, [ $chr, $pos, $stran, $bkdep ];
}
close $bk;

my ($chr, $stored, @at, @val, @pend, $pi, @out);

sub depth_at {
    my $k = $_[0] & $MASK;
    return (defined $at[$k] && $at[$k] == $_[0]) ? $val[$k] : undef;
}

sub evaluate {
    my ($idx) = @_;
    my ($c, $pos, $stran, $bkdep) = @{ $cand{$chr}[$idx] };
    my ($sl, $sr, $cnt) = (0, 0, 0);
    for my $i (1 .. $span) {
        my $l = depth_at($pos - $i);
        next unless defined $l;
        my $r = depth_at($pos + $i);
        next unless defined $r;
        $sl += $l; $sr += $r; $cnt++;
    }
    return if $cnt == 0;
    my $leftavg = int($sl / $cnt) + 0.1;
    my $rightavg = int($sr / $cnt) + 0.1;
    if ($stran eq "-") { $out[$idx] = "$c\t$pos\t$stran\t$bkdep\t$rightavg\n" }
    if ($stran eq "+") { $out[$idx] = "$c\t$pos\t$stran\t$bkdep\t$leftavg\n" }
}

sub flush_chr {
    return unless defined $chr;
    evaluate($pend[$pi++][1]) while $pi < @pend;
    if ($stored) { for my $o (@out) { print $o if defined $o } }
}

sub start_chr {
    ($chr) = @_;
    $stored = 0; @at = (); @val = (); @out = ();
    my $c = $cand{$chr} || [];
    @pend = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
            map { [ $c->[$_][1] + $span, $_ ] } 0 .. $#$c;
    $pi = 0;
}

open my $dp, '<', ($depfile eq '-' ? '/dev/stdin' : $depfile) or die "$depfile: $!";
while (my $l = <$dp>) {
    my $t1 = index($l, "\t");
    my $c = substr($l, 0, $t1);
    if (!defined $chr || $c ne $chr) {    # first line of a chromosome: not stored
        flush_chr();
        start_chr($c);
        next;
    }
    next unless $cand{$c};
    my $t2 = index($l, "\t", $t1 + 1);
    my $p = substr($l, $t1 + 1, $t2 - $t1 - 1);
    evaluate($pend[$pi++][1]) while $pi < @pend && $pend[$pi][0] < $p;
    # not within any pending window (pend is ordered by pos + span): not needed
    next if $pi >= @pend || $p < $pend[$pi][0] - 2 * $span;
    my $t3 = index($l, "\t", $t2 + 1);
    my $d = $t3 < 0 ? substr($l, $t2 + 1) : substr($l, $t2 + 1, $t3 - $t2 - 1);
    chomp $d;
    my $k = $p & $MASK;
    $at[$k] = $p; $val[$k] = $d; $stored = 1;
}
flush_chr();
close $dp;
