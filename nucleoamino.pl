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
# 2021/02/14 build internal array from loaded tables
# 2021/02/20 -fasta option

use strict;
use warnings;

# support library
use lib "$ENV{BIOSEQ}";
use bioutil qw(readfasta segment string load find);

# constants
use constant { true => 1, false => 0, AMINOS_PER_LINE => 70 };

# debug
use Data::Dumper qw(Dumper);

my $debug = 0;

{
    my @n3a1_table = (); # nucleotide triplet pattern (n3) to amino code (a1) mapping array

    # $amino_code = nucleo2amino($nucleo_triplet)
    # return an amino acid code for a nucleotide triplet
    sub nucleo2amino {
        my $triple = shift;

        # load tables and build an internal array on the first call
        if (@n3a1_table == 0) {
            my $nuc2ami_table = load("nuc2ami.txt");
            if ($debug > 1) { print STDERR Dumper($nuc2ami_table); }

            my $amino_table = load("amino.txt");
            if ($debug > 1) { print STDERR Dumper($amino_table); }

            foreach my $na (@{$nuc2ami_table}) {
                my $nucleo_pattern = $na->[0];
                my $amino_name = $na->[1];
                my $amino_code = find($amino_table, $amino_name, 1, 0);
                if ($debug) { print STDERR "$nucleo_pattern -> $amino_code\n"; }
                push(@n3a1_table, { n3 => $nucleo_pattern, a1 => $amino_code });
            }
        }

        my $code = undef;
        my $pattern = undef; # debug
        foreach my $na (@n3a1_table) {
            if ($triple =~ m/$na->{n3}/) {
                $pattern = $na->{n3};
                $code = $na->{a1};
                last;
            }
        }
        if ($debug) { print STDERR "nucleo2amino: $triple -> " . ($pattern ? $pattern : '???') . " -> " . ($code ? $code : '?') . "\n"; }
        return $code;
    }
}

sub test_nucleo2amino {
    $debug = 1;
    foreach my $nuc ("UUA", "UUG", "UGA", "UCA", "XXX") {
        my $amino = nucleo2amino($nuc);
        print "" . ($amino ? $amino : '?') . "\n";
    }
}

# $output = convert($line, [$fasta])
# convert nucleotide sequence to amino acid sequence
# if $fasta is true newlines are added; header has been output in main()
sub convert {
    my ($line, $fasta) = @_;
    my $len = length($line);
    my $len2 = int($len / 3) * 3; # truncate
    for (my $i = 0; $i < $len2; $i += 3) {
        my $triple = substr($line, $i, 3);
        my $single = nucleo2amino($triple);
        if ($single) {
            print $single;
        }
        else {
            if ($debug) { print '(' . $triple . ')'; }
            else { print 'x'; }
        }
        if ($fasta && ($i % (AMINOS_PER_LINE * 3) == (AMINOS_PER_LINE * 3) - 3)) {
            print "\n";
        } 
    }
    if ($len != $len2) { print STDERR "convert: line truncated from $len to $len2\n"; }
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
    print STDERR "Usage: nucleamino [-f] [-s<span>] [fastafile]\n";
    print STDERR " -f: output in fasta format\n";
    print STDERR " -s<span> = N-M where N'th to M'th chars are checked,\n" .
                 "            N+M where M chars from N'th position are checked.\n";
    exit 1;
}

sub main {
    my $path = '';
    my $span = '';
    my $fasta = false;
    foreach my $arg (@ARGV) {
        if ($arg =~ /^-/) {
            if ($arg eq '-d') { $debug = 1; }
            elsif ($arg eq '-t') { $debug = 2; }
            elsif ($arg eq '-f') { $fasta = true; }
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
    if ($fasta) { print "> " . ($path ? $path : '(stdin)') . " converted by nucleoamino\n"; }
    convert($line, $fasta);
    print "\n";
}

# test
#test_nucleo2amino();
#test_convert();
#test_cutData();

main();
