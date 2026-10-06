#!/usr/bin/perl
# asmqc helper for the patched CRAQ 1.10 drivers (not part of CRAQ).
#
# Usage: perl asmqc_region_depth.pl CANDIDATES FLANK BAM SEQS [MINMAPQ]
# Prints the lines of the per-base depth table "samtools depth -a BAM" (with
# MINMAPQ: of "samtools view -q MINMAPQ BAM | samtools depth -a -") that lie
# within FLANK bp of a candidate (columns 1-2 of CANDIDATES), plus position 1
# of each candidate's sequence, in table order. SEQS lists the sequences that
# appear in the full table (written while it was streamed); others are left
# out, as samtools depth -a leaves out sequences without reads. Position 1 is
# included because CRAQ's window scripts never store the first line of a
# sequence: with it, they see the same first line as in the full table.
# Queries go to the indexed BAM; no genome-wide table is written.
use strict;
use warnings;
use File::Temp qw(tempfile);

die "usage: $0 CANDIDATES FLANK BAM SEQS [MINMAPQ]\n" unless @ARGV == 4 || @ARGV == 5;
my ($cand, $flank, $bam, $seqs, $mapq) = @ARGV;

my %present;
open my $sf, '<', $seqs or die "$seqs: $!";
while (<$sf>) { chomp; $present{$_} = 1 if length }
close $sf;

my (%iv);
open my $cf, '<', $cand or die "$cand: $!";
while (my $l = <$cf>) {
    chomp $l;
    my ($chr, $pos) = (split /\t/, $l)[0, 1];
    next unless defined $pos && $present{$chr};
    my $s = $pos - 1 - $flank;
    push @{ $iv{$chr} }, [ $s < 0 ? 0 : $s, $pos + $flank ];
}
close $cf;
exit 0 unless %iv;

my ($bh, $bed) = tempfile("asmqc_region_XXXXXX", SUFFIX => ".bed", TMPDIR => 1, UNLINK => 1);
for my $chr (sort keys %iv) {
    print $bh "$chr\t0\t1\n";
    print $bh "$chr\t$_->[0]\t$_->[1]\n" for @{ $iv{$chr} };
}
close $bh;

my @cmd = ("samtools", "depth", "-a", "-b", $bed);
push @cmd, "-Q", $mapq if defined $mapq;
open(my $dp, '-|', @cmd, $bam) or die "samtools depth: $!";
while (my $l = <$dp>) {
    my $t = index($l, "\t");
    print $l if $present{ substr($l, 0, $t) };
}
close $dp or die "samtools depth failed: " . ($! || $? >> 8) . "\n";
