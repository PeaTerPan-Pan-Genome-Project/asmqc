#!/usr/bin/perl
# asmqc helper for the patched CRAQ 1.10 drivers (not part of CRAQ).
#
# Usage: perl asmqc_merge.pl MODE PART... > OUT
# Joins the per-segment outputs of a CRAQ script (parts in segment order) into
# what the script prints for the whole genome:
#   dici    caculate_clipDI_cov.pl: all deletion lines (column 5 "D") of all
#           parts, then all insertion lines, as the script prints them;
#   effsize LReffect_size.pl / SReffect_size.pl: "sequence\tcount" lines,
#           sorted as Perl strings (counts of a sequence split over parts are
#           added, which does not happen with whole-sequence segments);
#   dep0    search_dep0.pl: blocks per sequence, sequences sorted as Perl
#           strings, lines within a sequence in part order.
# Outputs whose order follows the depth table are joined with cat instead.
use strict;
use warnings;

my ($mode, @parts) = @ARGV;
die "usage: $0 dici|effsize|dep0 PART...\n" unless $mode && $mode =~ /^(dici|effsize|dep0)$/;

if ($mode eq 'dici') {
    my (@d, @i);
    for my $p (@parts) {
        open(my $f, '<', $p) or die "$p: $!";
        while (my $l = <$f>) {
            my $k = (split /\t/, $l)[4];
            if (defined $k && $k eq 'I') { push @i, $l } else { push @d, $l }
        }
        close $f;
    }
    print @d, @i;
} elsif ($mode eq 'effsize') {
    my %c;
    for my $p (@parts) {
        open(my $f, '<', $p) or die "$p: $!";
        while (<$f>) { chomp; my ($s, $n) = split /\t/; $c{$s} += $n }
        close $f;
    }
    print "$_\t$c{$_}\n" for sort keys %c;
} else {
    my %blk;
    for my $p (@parts) {
        open(my $f, '<', $p) or die "$p: $!";
        while (my $l = <$f>) {
            my $s = substr($l, 0, index($l, "\t"));
            push @{ $blk{$s} }, $l;
        }
        close $f;
    }
    print @{ $blk{$_} } for sort keys %blk;
}
