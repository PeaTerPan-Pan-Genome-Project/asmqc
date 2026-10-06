#!/usr/bin/perl
# asmqc helper for the patched CRAQ 1.10 drivers (not part of CRAQ).
#
# Usage: perl asmqc_par.pl N 'cmd1' 'cmd2' ...
# Runs the commands (bash, pipefail) at most N at a time. Exits non-zero,
# naming each failed command, if any command fails; the others still finish.
use strict;
use warnings;

die "usage: $0 N CMD...\n" unless @ARGV >= 1;
my ($n, @cmds) = @ARGV;
$n = 1 if $n < 1;
my (%running, $fail);
for my $cmd (@cmds) {
    reap() while keys %running >= $n;
    my $pid = fork() // die "fork: $!";
    if ($pid == 0) { exec('bash', '-o', 'pipefail', '-c', $cmd) or die "exec: $!" }
    $running{$pid} = $cmd;
}
reap() while %running;
exit($fail ? 1 : 0);

sub reap {
    my $pid = wait();
    return if $pid < 0;
    my $cmd = delete $running{$pid};
    if ($?) {
        print STDERR "asmqc_par: failed (exit " . ($? >> 8) . "): $cmd\n";
        $fail = 1;
    }
}
