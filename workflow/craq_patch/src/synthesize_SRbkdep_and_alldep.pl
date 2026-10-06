#!/usr/bin/perl
# Replacement for src/synthesize_SRbkdep_and_alldep.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output.
#
# The original runs a regex and builds two string keys for every line of the
# per-base depth table. This version looks up sequence and position with
# index/substr first and applies the original regex and logic only to lines
# at a clip site.
use strict;

if (@ARGV != 2) {
    print "\nUSE: this script is to combian breakpoint depth and total_depth \n \nperl $0  SR_break.depth SR_sort.depth - > my.cover_rate\n\n";
    exit 1;
}
my ($bk_depth_file, $depth_file) = @ARGV;
my (%bk, %want);
open my $in0, '<', $bk_depth_file or die "$bk_depth_file: $!";
while (my $line = <$in0>) {
    chomp $line;
    if ($line =~ /(\S+)\t(\d+)\t(\S+)\t(\d+)/) {
        $bk{"$1\t$2\t$3"} = $4;
        $want{$1}{$2} = 1 if $3 eq "+" or $3 eq "-";    # only these keys are looked up
    }
}
close $in0;

open my $in1, '<', $depth_file or die "$depth_file: $!";
while (my $line2 = <$in1>) {
    my $t1 = index($line2, "\t");
    next if $t1 < 0;
    my $w = $want{ substr($line2, 0, $t1) } or next;
    my $t2 = index($line2, "\t", $t1 + 1);
    next unless $t2 > 0 && $w->{ substr($line2, $t1 + 1, $t2 - $t1 - 1) };
    chomp $line2;
    if ($line2 =~ /(\S+)\t(\d+)\t(\d+)/) {
        my $key1 = "$1\t$2\t\+";
        my $key2 = "$1\t$2\t\-";
        my $depth = $3;
        if (defined $bk{$key1}) { print "$key1\t$bk{$key1}\t$depth\n"; }
        if (defined $bk{$key2}) { print "$key2\t$bk{$key2}\t$depth\n"; }
    }
}
close $in1;
