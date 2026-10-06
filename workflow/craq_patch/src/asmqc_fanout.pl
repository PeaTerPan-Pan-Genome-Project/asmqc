#!/usr/bin/perl
# asmqc helper for the patched CRAQ 1.10 drivers (not part of CRAQ).
#
# Usage: producer | perl asmqc_fanout.pl 'cmd1' 'cmd2' ...
# Copies stdin, in large blocks, to every command (each run by sh -c) and
# replaces a per-base depth table that several CRAQ scripts would otherwise
# read from disk one after another. Exits non-zero, naming the command, if any
# command fails or stops reading early; a pipeline with a dead consumer cannot
# hang, unlike tee with FIFOs.
use strict;
use warnings;

die "usage: $0 CMD...\n" unless @ARGV;
$SIG{PIPE} = 'IGNORE';
my @out;
for my $cmd (@ARGV) {
    open(my $fh, '|-', $cmd) or die "asmqc_fanout: cannot start '$cmd': $!\n";
    binmode $fh;
    push @out, [ $fh, $cmd ];
}
binmode STDIN;
my $buf;
while (1) {
    my $n = sysread(STDIN, $buf, 1 << 20);
    die "asmqc_fanout: read error: $!\n" unless defined $n;
    last if $n == 0;
    for my $o (@out) {
        my $off = 0;
        while ($off < $n) {
            my $w = syswrite($o->[0], $buf, $n - $off, $off);
            die "asmqc_fanout: '$o->[1]' stopped reading: $!\n" unless defined $w;
            $off += $w;
        }
    }
}
my $fail = 0;
for my $o (@out) {
    next if close($o->[0]);
    my $rc = $! ? "error $!" : "exit " . ($? >> 8) . ($? & 127 ? ", signal " . ($? & 127) : "");
    print STDERR "asmqc_fanout: '$o->[1]' failed ($rc)\n";
    $fail = 1;
}
exit $fail;
