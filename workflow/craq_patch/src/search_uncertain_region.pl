#!/usr/bin/perl
# Replacement for src/search_uncertain_region.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments and
# the same output lines.
#
# Prints the zero-depth regions of >= 500 bp (file 1) with no putative error
# (file 2) within 50 bp. The original compares every region with every error
# on the same sequence; this version sorts each sequence's error positions
# and searches them. The original prints the regions in Perl hash order (it
# differs between runs); this version prints them in the order of file 1. The
# only reader, intergrate_uncertain.pl, adds up their lengths.
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

my $dis = 50;
my (@regions, %seen);
open my $in0, '<', ($ARGV[0] eq '-' ? '/dev/stdin' : $ARGV[0]) or die "$ARGV[0]: $!";
while (<$in0>) {
    chomp;
    my ($chr, $s, $e) = (split /\s+/)[0, 1, 2];
    next unless ($e - $s) >= 500;
    my $key = "$chr\t$s\t$e";
    push @regions, $key unless $seen{$key}++;
}
close $in0;

my %err;
open my $in1, '<', ($ARGV[1] eq '-' ? '/dev/stdin' : $ARGV[1]) or die "$ARGV[1]: $!";
while (<$in1>) {
    chomp;
    my ($chr, $error_loc) = (split /\s+/)[0, 1];
    $err{$chr}{$error_loc} = 1;
}
close $in1;
my %sorted = map { $_ => [ sort { $a <=> $b } keys %{ $err{$_} } ] } keys %err;

for my $key (@regions) {
    my ($chr, $s, $e) = (split /\s+/, $key)[0, 1, 2];
    my $v = $sorted{$chr};
    my $common = 0;
    if ($v && @$v) {    # first error position >= s - dis; is it <= e + dis?
        my ($lo, $hi) = (0, scalar @$v);
        while ($lo < $hi) {
            my $mid = ($lo + $hi) >> 1;
            if ($v->[$mid] >= ($s - $dis)) { $hi = $mid } else { $lo = $mid + 1 }
        }
        $common = $lo < @$v && $v->[$lo] <= ($e + $dis);
    }
    print "$chr\t$s\t$e\n" unless $common;
}
