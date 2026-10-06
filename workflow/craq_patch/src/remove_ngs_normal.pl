#!/usr/bin/perl
# Replacement for src/remove_ngs_normal.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output.
#
# The original stores every position of the zero-depth list (SRout/Nonmap.loc,
# ~500 M on a pea genome) in a nested hash and then looks only +-10 bp around
# each CRE candidate. This version reads the candidates first and keeps only
# the listed positions inside a candidate window. Assumes integer positions
# written as samtools does.
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

my ($covrate, $ngsdep0, $cre_file) = @ARGV;
my $W = 10;

my (@cre, %win);
open my $in3, '<', ($cre_file eq '-' ? '/dev/stdin' : $cre_file) or die "$cre_file: $!";
while (my $l = <$in3>) {
    chomp $l;
    my ($chr, $pos) = (split /\s+/, $l)[0, 1, 2];
    push @cre, [ $l, $chr, $pos ];
    push @{ $win{$chr} }, [ $pos - $W, $pos + $W ];
}
close $in3;
for my $c (keys %win) {    # sorted, disjoint windows per sequence
    my @v = sort { $a->[0] <=> $b->[0] } @{ $win{$c} };
    my (@st, @en);
    for my $x (@v) {
        if (@st && $x->[0] <= $en[-1] + 1) { $en[-1] = $x->[1] if $x->[1] > $en[-1] }
        else { push @st, $x->[0]; push @en, $x->[1] }
    }
    $win{$c} = [ \@st, \@en ];
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

my %hash;
open my $in1, '<', ($covrate eq '-' ? '/dev/stdin' : $covrate) or die "$covrate: $!";
while (<$in1>) {
    chomp;
    my ($chr, $pos, $cov, $depth) = (split /\s+/)[0, 1, 3, 4];
    next if ($depth == 0);
    my $rate = $cov / $depth;
    $hash{$chr}{$pos} = 1 if $rate >= 0.1;
}
close $in1;

my %cur;
open my $in2, '<', ($ngsdep0 eq '-' ? '/dev/stdin' : $ngsdep0) or die "$ngsdep0: $!";
while (my $l = <$in2>) {
    my $t = index($l, "\t");
    my $c = substr($l, 0, $t);
    next unless $win{$c};
    my $u = index($l, "\t", $t + 1);
    my $p = $u < 0 ? substr($l, $t + 1) : substr($l, $t + 1, $u - $t - 1);
    chomp $p;
    $hash{$c}{$p} = 1 if hit(@{ $win{$c} }, \($cur{$c} //= 0), $p);
}
close $in2;

for my $x (@cre) {
    my ($line, $chr, $pos) = @$x;
    for my $posi (($pos - $W) .. ($pos + $W)) {
        if (defined $hash{$chr}{$posi}) { print "$line\n"; last }
    }
}
