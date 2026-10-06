#!/usr/bin/perl
# Replacement for src/synthesize_clipDIcov_and_alldep.pl of CRAQ 1.10
# (https://github.com/JiaoLaboratory/CRAQ, MIT License, Copyright (c) 2023
# JiaoLaboratory; see ../LICENSE.CRAQ), used by asmqc M9. Same arguments,
# byte-identical output.
#
# The original runs a regex and builds two string keys for every line of the
# per-base depth table (4.3 G lines for a pea genome). This version looks up
# sequence and position with index/substr first and applies the original
# regex and logic only to lines at an indel site.
use strict;
# "-" reads standard input, as CRAQ's two-argument open() does.

if (@ARGV != 2) {
    print "\nUSE: this script is to combian breakpoint depth and total_depth \n \nperl $0  break.cov depth - > my.cover_rate\n\n";
    exit 1;
}
my ($bk_depth_file, $depth_file) = @ARGV;
my (%bk, %bk_tmp, %want);
open my $in0, '<', ($bk_depth_file eq '-' ? '/dev/stdin' : $bk_depth_file) or die "$bk_depth_file: $!";
while (my $line = <$in0>) {
    chomp $line;
    my ($chr, $DIpos, $stran, $DInum, $DIlen, $class, $DIsame_num) =
        (split /\s+/, $line)[0, 1, 2, 3, 4, 5, 6];
    next if ($DInum < 3 or $DIsame_num < 0.6);
    $bk{"$chr\t$DIpos\t$stran"} = $DInum;
    $bk_tmp{"$chr\t$DIpos\t$stran"} = "$DIlen\t$class\t$DIsame_num";
    $want{$chr}{$DIpos} = 1 if $stran eq "+";    # only "+" keys are ever looked up
}
close $in0;

open my $in1, '<', ($depth_file eq '-' ? '/dev/stdin' : $depth_file) or die "$depth_file: $!";
while (my $line2 = <$in1>) {
    my $t1 = index($line2, "\t");
    next if $t1 < 0;
    my $w = $want{ substr($line2, 0, $t1) } or next;
    my $t2 = index($line2, "\t", $t1 + 1);
    next unless $t2 > 0 && $w->{ substr($line2, $t1 + 1, $t2 - $t1 - 1) };
    chomp $line2;
    if ($line2 =~ /(\S+)\t(\d+)\t(\d+)/) {
        my $key1 = "$1\t$2\t\+";
        my $depth = $3;
        next if ($depth == 0);
        if (defined $bk{$key1}) {
            next if ($bk{$key1} / $depth < 0.2);
            print "$key1\t$bk{$key1}\t$depth\t$bk_tmp{$key1}\n";
        }
    }
}
close $in1;
