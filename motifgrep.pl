#!/usr/bin/env perl

# MotifGrep - grep bio sequence with motif pattern
# Copyright (c) 2021 Kobu.Com. All rights reserved.
# Visit Kobu.Com at https://kobu.com/bioseq/
# Licensed by GNU Public License v3
# https://www.gnu.org/licenses/gpl-3.0.en.html

# MotifGrep searches input sequence for a motif pattern.
# A motif pattern can contain:
# - choice of characters within [..],
# - exclusion with {..}. 
# MotifGrep supports the PROSITE extensions:
# - hyphen (-) delimiter
# - range with e(n) and e(n,m)
# Motifgrep prints a match on each line with its position (one-based).
#
# Example:
# $ cat test.fasta
# >test data
# ATATATCAGAG
# AGCAGAGACC
# $ motifgrep A[TG]C test.fasta
#  5:ATC
# 12:AGC
#
# 2021/01/27 written
# 2021/02/01 readfasta()
# 2021/02/09 amino2nucleo()

use strict;
use warnings;

# support library
use lib "$ENV{BIOSEQ}";
use bioutil qw(readfasta load find);

# constants
use constant { true => 1, false => 0 };

# debug
use Data::Dumper qw(Dumper);

my $debug = 0;

{
    # TODO: build amino code to nucleo triplet table before hand

    my $amino_table = load("amino.txt");
    if ($debug > 1) { print Dumper($amino_table); }

    my $nuc2ami_table = load("nuc2ami.txt");
    if ($debug > 1) { print Dumper($nuc2ami_table); }

    # $nucleos = amino2nucleo($aminos, [$useU])
    # convert known amino codes to corresponding nucleotide patterns
    # eg. N x Y -> UA[UC] xxx AA[UC]
    # it is an error to pass a non-amino character except 'x' which is converted to 'xxx'
    # TODO: support [..] or {..} ? N[YA] will be something like AA[UC](UA[UC]|GCx)
    sub amino2nucleo {
        my ($amino, $useU) = @_;
        my @nucleo =();
        foreach my $code (split(//, $amino)) {
            my $result = undef;
            if ($code eq 'x') {
                $result = 'xxx'; # wildcard triple .. TODO: put this in table?
            }
            else {
                my $name = find($amino_table, $code, 0, 1);
                if ($name) {
                    my $triple = find($nuc2ami_table, $name, 1, 0);
                    if ($triple) { $result = $triple; }
                }
                die "amino2nucleo: not an amino code: $code" unless $result;
            }
            push(@nucleo, $result);
        }
        my $motif = join('', @nucleo);
        if (!$useU) { $motif =~ tr/U/T/; }
        return $motif;
    }
}

sub test_amino2nucleo {
    foreach my $a ("N", "NY", "NXY", "N[YA]") {
        print "$a -> " . amino2nucleo($a) . "\n";
    }
}

# $regexp = motif2regexp($motif)
# convert a motif pattern to a regexp pattern
# x xx -> . ..
# [X] [XY] -> [X] [XY] (no change)
# {X} {XY} -> [^X] [^XY]
# X(n) X(n,m) -> X{m} X{n,m}
# '-' -> just remove
sub motif2regexp {
    my $s = shift;
    $s =~ s/x/\./g; # x -> . (any character)
    $s =~ s/{([^}]+)}/[^$1]/g; # {..} -> [^..] (either of)
    $s =~ s/\(([^)]+)\)/{$1}/g; # (..) -> {..} (except)
    $s =~ s/-//g; # remove hyphen delimiter
    return $s;
}

sub test_motif2regexp {
    my $m = 'xa[bc]d{ef}g-h(3)-i(3,4)x';
    my $r = motif2regexp($m);
    print "$m -> $r\n";
}

# grepline($regexp, $line)
sub grepline {
    my ($regexp, $line) = @_;
    my $len = length($line);
    my $width = length("$len");
    my $n = 1;
    while (true) {
        if ($line =~ /$regexp/p) {
            my $prematch = ${^PREMATCH};
            my $match = ${^MATCH};
            my $postmatch = ${^POSTMATCH};
            $n += length($prematch);
            my $pos = sprintf("%${width}s", $n);
            print "$pos:$match\n";
            $n += length($match);
            $line = $postmatch;
        }
        else {
            # ignore rest in $line
            last;
        }
    }
}

sub test_grepline {
    my $line = "abcdefghijkl0123456789aaadddgggjjj";
    print "$line\n";
    foreach my $r ("abc", "a[ab][^x]", "a{3}d{2,4}ggg") {
        print "[$r]\n";
        grepline($r, $line);
    }
}

sub usage {
    print STDERR "Usage: motifgrep [-n|-nu] motif [path]\n";
    print STDERR "  '-n' option converts amino acid pattern to nucleotide pattern.\n";
    print STDERR "  '-nu' forces use of U instead of T in the converted pattern\n";
    exit 1;
}

sub main {
    my $motif = '';
    my $path = '';
    my $nucleo = false;
    my $useU = false;
    foreach my $arg (@ARGV) {
        if ($arg =~ /^-/) {
            if ($arg eq '-d') { $debug = 1; }
            elsif ($arg eq '-t') { $debug = 2; }
            elsif ($arg =~ /^-n([ut]?)$/) {
                 $nucleo = true;
                 if ($1 eq 'u') { $useU = true; }
            }
            else {
                print STDERR "No such option: $arg\n";
                usage();
            }
        }
        elsif (!$motif) { $motif = $arg; }
        elsif (!$path) { $path = $arg; }
        else { usage; }
    }
    if (!$motif) { usage(); }
    if ($nucleo) {
        my $amino = $motif;
        $motif = amino2nucleo($amino, $useU);
        if ($debug) { print STDERR "amino '$amino' -> "; }
    }
    my $regexp = motif2regexp($motif);
    if ($debug) { print STDERR "motif '$motif' -> regexp '$regexp'\n"; }
    my $line = readfasta($path); # read a single line from a file or stdin
    grepline($regexp, $line);
}

# test
#test_amino2nucleo();
#test_motif2regexp();
#test_grepline();

main();
