#!/usr/bin/env perl

# NucleAmino - convert nucleotide sequence to amino acid sequence
# Copyright (c) 2021 Kobu.Com. All rights reserved.
# Visit Kobu.Com at https://kobu.com/bioseq/
# Licensed by GNU Public License v3
# https://www.gnu.org/licenses/gpl-3.0.en.html

# NucleAmino converts an input nucleotide sequence to an amino acid sequence.
# A triplet of nucleotide codes such as "AAU" and "UAU" corresponds to 
# a type of amino acid that comprise a protein, such as Asparagine and
# Tyrosine, respectively, whose codes are N and Y.
# NucleAmino uses two mapping tables, nuc2ami.txt and amino.txt, to convert
# a triplet to a short amino acid name and then to a single-letter code. 
# Note: Not all part of a nucleotid sequence represents amino acids.
#       One nucleotide triplet corresponds to one amino acid.
#       More than one triplet corresponds to the same amino acid.
# 
# Example:
#  $ cat test.fasta
#   >test data
#   ATATATCAGAG
#   AGCAGAGAGC
#  $ nucleamino test.fasta
#   IYQRAES 
#  $ nucleamino -s2 test.fasta
#   YIREQR
#
# Hint: The '-s' option allows you specify the start position or limit a conversion
# to a portion of the sequence. See usage() at the end of this file.
#
# 2021/01/27 written
# 2021/02/01 readfasta
# 2021/02/02 span option

use strict;
use warnings;

# support library
use lib "$ENV{BIOSEQ}";
use bioutil qw(readfasta segment string load find);

# constants
use constant { true => 1, false => 0 };

# debug
use Data::Dumper qw(Dumper);

my $debug = 0;

{
    # TODO: build nucleo triplet amino code table before hand

    my $nuc2ami_table = load("nuc2ami.txt");
    if ($debug > 1) { print Dumper($nuc2ami_table); }

    my $amino_table = load("amino.txt");
    if ($debug > 1) { print Dumper($amino_table); }

    # $amino = nucle2amino($nucles)
    # return an amino acid code for a nucleotide tuple
    sub nucle2amino {
        my $triple = shift;
        my $code = undef;
        my $name = find($nuc2ami_table, $triple, 0, 1);
        if ($name) {
            $code = find($amino_table, $name, 1, 0);
        }
        if ($debug) {
            print "nucle2amino: $triple -> " .
                ($name ? $name : "???") . " -> " . 
                ($code ? $code : "?") . "\n";
        };
        return $code;
    }
}

sub test_nucle2amino {
    $debug = 1;
    foreach my $nuc ("UUA", "UUG", "UGA", "UCA", "XXX") {
        my $amino = nucle2amino($nuc);
        print "" . ($amino ? $amino : '?') . "\n";
    }
}

# $output = convert($line)
# convert nucleotide sequence to amino acid sequence
sub convert {
    my $line = shift;
    my $len = length($line);
    my $len2 = int($len / 3) * 3; # truncate
    my $outstr = '';
    for (my $i = 0; $i < $len2; $i += 3) {
        my $triple = substr($line, $i, 3);
        my $single = nucle2amino($triple);
        if ($single) {
            $outstr .= $single;
        }
        else {
            if ($debug) { $outstr .= ( '(' . $triple . ')' ); }
            else { $outstr .= 'x'; }
        }
    }
    if ($len != $len2) { print STDERR "convert: line truncated from $len to $len2\n"; }
    return $outstr;
}

sub test_convert {
    $debug = 1;
    my $outstr = convert('AUGCUAGGCUCG');
    print "$outstr\n";
}

# TODO: move to bioutil if others use this

# $segment = cutData($line, $span_spec)
# return a string segment matching the command line span (-s) option
# 'from' and 'to' are one-based:
# if $line is 'hello' you will get:
#  from-to    1-5 -> "hello", 2-4 -> "ell"
#  from+len   1+5 -> "hello", 2+3 -> "ell"
#  from       1   -> "hello", 2   -> "ello"
#  -to        -3  -> "hel"
#  +len       +3  -> "hel"
sub cutData {
    my ($line, $span) = @_;
    my $pos = 0;
    my $len = length($line);
    if ($span =~ /(\d*)([-+]?)(\d*)/) {
        my $n = $1;
        my $op = $2;
        my $m = $3;
        if ($n) { $pos = int($n) - 1; }
        if ($m) {
            if ($op eq '+') { $len = $m; }
            else { $len = int($m) - $pos; } # '-' or missing
        }
    }
    return segment($line, $pos, $len);
}

sub test_cutData {
    foreach my $ss ('1-5', '2-4', '1+5', '2+3', '1', '2', '-3', '+3', '-') {
        my $seg = cutData('hello', $ss);
        print "$ss -> " . string($seg) . "\n";
    }
}

sub usage {
    print STDERR "Usage: nucleamino [-s<span>] [fastafile]\n";
    print STDERR " <span> = N-M where N'th to M'th chars are checked, " .
                 "          N+M where M chars from N'th position are checked.\n";
    exit 1;
}

sub main {
    my $path = '';
    my $span = '';
    foreach my $arg (@ARGV) {
        if ($arg =~ /^-/) {
            if ($arg eq '-d') { $debug = 1; }
            elsif ($arg eq '-t') { $debug = 2; }
            elsif ($arg =~ /^-s(.+)$/) { $span = $1; }
            else {
                print STDERR "No such option: $arg\n";
                usage();
            }
        }
        elsif (!$path) { $path = $arg; }
        else { usage; }
    }
    my $line = readfasta($path); # read a single line from a file or stdin
    if ($span) { $line = string(cutData($line, $span)); }
    $line =~ tr/T/U/; # U -> T; nuc2ami.txt contains U instead of T
    my $outstr = convert($line);
    print "$outstr\n";
}

# test
#test_nucle2amino();
#test_convert();
#test_cutData();

main();
