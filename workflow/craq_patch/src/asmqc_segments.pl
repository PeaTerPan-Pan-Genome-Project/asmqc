#!/usr/bin/perl
# asmqc helper for the patched CRAQ 1.10 drivers (not part of CRAQ).
#
# Usage: perl asmqc_segments.pl BAM N OUTPREFIX
# Splits the BAM's sequences (header order) into at most N contiguous
# segments of similar total length and writes OUTPREFIX01.bed, OUTPREFIX02.bed,
# ... (whole sequences). Per-segment outputs concatenated in segment order
# follow header order, as the outputs of a single genome-wide pass do. A new
# segment starts when the next sequence would push the current one beyond
# max(total / N, longest sequence). Prints the BED file names.
use strict;
use warnings;

die "usage: $0 BAM N OUTPREFIX\n" unless @ARGV == 3;
my ($bam, $n, $prefix) = @ARGV;
$n = 1 if $n < 1;
my @sq;
open(my $h, '-|', 'samtools', 'view', '-H', $bam) or die "samtools view -H: $!";
while (<$h>) {
    next unless /^\@SQ\t/;
    my ($sn) = /\tSN:([^\t\n]+)/;
    my ($ln) = /\tLN:(\d+)/;
    push @sq, [ $sn, $ln ];
}
close $h or die "samtools view -H $bam failed\n";
die "$bam: no \@SQ lines\n" unless @sq;

my ($total, $longest) = (0, 0);
for (@sq) { $total += $_->[1]; $longest = $_->[1] if $_->[1] > $longest }
my $target = $total / $n > $longest ? $total / $n : $longest;
my (@seg, $sum);
for my $s (@sq) {
    if (!@seg || $sum + $s->[1] > $target) { push @seg, []; $sum = 0 }
    push @{ $seg[-1] }, $s;
    $sum += $s->[1];
}
my $w = length scalar @seg;
$w = 2 if $w < 2;
for my $i (0 .. $#seg) {
    my $f = sprintf("%s%0*d.bed", $prefix, $w, $i + 1);
    open(my $o, '>', $f) or die "$f: $!";
    print $o "$_->[0]\t0\t$_->[1]\n" for @{ $seg[$i] };
    close $o;
    print "$f\n";
}
