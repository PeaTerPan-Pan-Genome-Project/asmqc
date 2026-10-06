#!/usr/bin/perl
# Replacement for src/get_nonmap_region.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output; extends the rewrite in asmqc issue #3.
#
# The original stores every position of file 1 in a nested hash, adds p..p+3
# for every position of file 2 that is also in file 1 to a second nested hash,
# and sorts all keys of each sequence to print the maximal runs of consecutive
# positions as "seq<TAB>start<TAB>end<TAB>0". CRAQ calls it with SRout/Nonmap.loc
# (every zero-depth short-read position: ~500 M on a pea genome) as file 1,
# once with the same file as file 2 and once with LRout/Nonmap.loc; the hashes
# take ~100 GB.
#
# This version keeps file 1 as merged runs of consecutive positions (not read
# at all when both arguments are the same file: every position is then a
# member) and builds the output directly as merged intervals [p, p+3]. Output
# order as in the original: sequences sorted as strings, runs by start. Input
# order does not matter. Assumes integer positions written as samtools does
# (the original compares them as hash-key strings).
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

my ($f0, $f1) = @ARGV;
my $DIS = 3;
my @s0 = stat($f0);
my @s1 = stat($f1);
my $same = @s0 && @s1 && $s0[0] == $s1[0] && $s0[1] == $s1[1];

# runs of consecutive positions per sequence: [starts], [ends], merged and sorted
sub runs_of {
    my ($list) = @_;    # [[s, e], ...] in reading order
    my @v = sort { $a->[0] <=> $b->[0] } @$list;
    my (@st, @en);
    for my $x (@v) {
        if (@st && $x->[0] <= $en[-1] + 1) { $en[-1] = $x->[1] if $x->[1] > $en[-1] }
        else { push @st, $x->[0]; push @en, $x->[1] }
    }
    return [ \@st, \@en ];
}

sub fields {    # sequence and position, as (split /\s+/)[0, 1] on samtools lines
    my $t = index($_[0], "\t");
    my $u = index($_[0], "\t", $t + 1);
    return (substr($_[0], 0, $t), $u < 0 ? substr($_[0], $t + 1) : substr($_[0], $t + 1, $u - $t - 1));
}

my %member;
if (!$same) {
    my %raw;
    open my $in0, '<', ($f0 eq '-' ? '/dev/stdin' : $f0) or die "$f0: $!";
    while (my $l = <$in0>) {
        chomp $l;
        my ($c, $p) = fields($l);
        my $r = $raw{$c} ||= [];
        if (@$r && $p >= $r->[-1][0] && $p <= $r->[-1][1] + 1) {
            $r->[-1][1] = $p if $p > $r->[-1][1];
        } else {
            push @$r, [ $p, $p ];
        }
    }
    close $in0;
    $member{$_} = runs_of($raw{$_}) for keys %raw;
}

# Is p inside one of the sorted, disjoint intervals (starts $st, ends $en)? $k
# caches the interval reached last: for ascending positions the pointer only
# moves forward; a position before it is looked up by binary search.
sub hit {
    my ($st, $en, $k, $p) = @_;
    my $i = $$k;
    if ($p < $st->[$i]) {
        return 0 if $i == 0 || $p > $en->[$i - 1];    # between two intervals
        my ($lo, $hi) = (0, $i - 1);                  # last start <= p
        while ($lo < $hi) {
            my $mid = ($lo + $hi + 1) >> 1;
            if ($st->[$mid] <= $p) { $lo = $mid } else { $hi = $mid - 1 }
        }
        $i = $lo;
    } else {
        $i++ while $i < $#$st && $en->[$i] < $p;
    }
    $$k = $i;
    return $st->[$i] <= $p && $p <= $en->[$i];
}

my %cur;
sub is_member {
    my ($c, $p) = @_;
    my $m = $member{$c} or return 0;
    return hit(@$m, \($cur{$c} //= 0), $p);
}

my %iv;    # sequence -> [[start, end], ...]
open my $in1, '<', ($f1 eq '-' ? '/dev/stdin' : $f1) or die "$f1: $!";
while (my $l = <$in1>) {
    chomp $l;
    my ($c, $p) = fields($l);
    next if !$same && !is_member($c, $p);
    my ($s, $e) = ($p + 0, $p + $DIS);
    my $list = $iv{$c} ||= [];
    if (@$list && $s >= $list->[-1][0] && $s <= $list->[-1][1] + 1) {
        $list->[-1][1] = $e if $e > $list->[-1][1];
    } else {
        push @$list, [ $s, $e ];
    }
}
close $in1;

for my $c (sort keys %iv) {
    my ($st, $en) = @{ runs_of($iv{$c}) };
    print "$c\t$st->[$_]\t$en->[$_]\t0\n" for 0 .. $#$st;
}
