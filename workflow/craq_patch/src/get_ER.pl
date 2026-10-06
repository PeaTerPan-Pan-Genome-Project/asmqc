#!/usr/bin/perl
# Replacement for src/get_ER.pl of CRAQ 1.10 (https://github.com/JiaoLaboratory/CRAQ,
# MIT License, Copyright (c) 2023 JiaoLaboratory; see ../LICENSE.CRAQ), used by
# asmqc M9. Same arguments, byte-identical output; based on the rewrite in
# asmqc issue #2.
#
# The original keeps every per-base depth of a chromosome in a nested hash (up
# to ~750 M entries, ~130 GB for a pea chromosome) and reads only
# +-window*window_extend_num bp around each candidate. This version streams the
# depth file once with a ring buffer, skips lines outside every window and
# evaluates a candidate as soon as its window has been read; window sums grow incrementally (integer depths, so the
# averages are bit-identical). Kept on purpose, as in the original:
#  - the first depth line of every chromosome is not stored;
#  - a chromosome with no stored line yields nothing;
#  - an empty window skips that window size only ("next" inside avg());
#  - output per chromosome in depth-file order, candidates in input order.
# Assumes samtools depth order: each chromosome in one block, positions
# ascending.
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

if (@ARGV != 5) {
    print "USE: $0 pb/ont_depth_file   NGS_coverrate.file   window_size    window_extend_num   threshold     \n";
    exit 1;
}
my ($depfile, $covfile, $window, $win_n, $skewrate) = @ARGV;
my $span = $window * $win_n;
my $ring = 1;
$ring <<= 1 while $ring < 2 * $span + 2;
my $MASK = $ring - 1;

my %cand;    # chr -> [[line, pos, clipnum], ...] in input order
open my $cv, '<', ($covfile eq '-' ? '/dev/stdin' : $covfile) or die "$covfile: $!";
while (my $l = <$cv>) {
    chomp $l;
    my @arr = split /\t/, $l;    # as the original: trailing empty fields dropped
    push @{ $cand{ $arr[0] } }, [ join("\t", @arr), $arr[1], $arr[3] ];
}
close $cv;

my ($chr, $stored, @at, @val, @pend, $pi, @res);

sub depth_at {
    my $k = $_[0] & $MASK;
    return (defined $at[$k] && $at[$k] == $_[0]) ? $val[$k] : undef;
}

sub evaluate {
    my ($idx) = @_;
    my ($line, $pos, $clipnum) = @{ $cand{$chr}[$idx] };
    my $threshold = $clipnum > 1 ? 1 - $skewrate : 0.3 * (1 - $skewrate);
    my ($sl, $sr, $cnt, $i) = (0, 0, 0, 0);
    for my $n (1 .. $win_n) {
        for ($i = $i + 1; $i <= $window * $n; $i++) {
            my $l = depth_at($pos - $i);
            next unless defined $l;
            my $r = depth_at($pos + $i);
            next unless defined $r;
            $sl += $l; $sr += $r; $cnt++;
        }
        $i = $window * $n;
        next if $cnt == 0;
        my $la = $sl / $cnt;
        my $ra = $sr / $cnt;
        if ($la <= 5 and $ra <= 5) {
            if ($la < $threshold * $threshold * $ra or $ra < $la * $threshold * $threshold) {
                $res[$idx] = $line; last;
            }
        }
        if ($la > 5 or $ra > 5) {
            if ($la < $threshold * $ra or $ra < $la * $threshold) { $res[$idx] = $line; last }
        }
    }
}

sub flush_chr {
    return unless defined $chr;
    evaluate($pend[$pi++][1]) while $pi < @pend;
    if ($stored) { for my $l (@res) { print "$l\n" if defined $l } }
}

sub start_chr {
    ($chr) = @_;
    $stored = 0; @at = (); @val = (); @res = ();
    my $c = $cand{$chr} || [];
    # a candidate is ready once a position beyond pos + span has been read
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
