#!/usr/bin/env perl

# BioUtil - common functions used by bio-sequence utilities
# Copyright (c) 2021 Kobu.Com. All rights reserved.
# Visit Kobu.Com at https://kobu.com/bioseq/
# Licensed by GNU Public License v3
# https://www.gnu.org/licenses/gpl-3.0.en.html
#
# 2021/02/01 readfasta()

package bioutil;

use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);
@EXPORT_OK = qw(readfasta segment string load find);

# debug
use Data::Dumper qw(Dumper);

# constants
use constant { true => 1, false => 0 };

# \@array = load($filename)
# load a tab-separated file into an array (ref) of arrays (ref)
sub load {
    my $filename = shift;
    my @arr = ();
    open (IN, '<', $ENV{BIOSEQ} . '/' . $filename);
    while(<IN>) {
        if (/^#/ || /^\s*$/) { next; } # skip comment and blank line
        my @items = split(/\s+/);
        push @arr, \@items;
    }
    close(IN);
    return \@arr;
}

# $match = find(\@table, $key, $n, $m)
# search @table and return $m'th value in the line where $key matches $n'th value
sub find {
    my ($table_ref, $key, $n, $m) = @_;
    foreach my $item_ref (@{$table_ref}) {
        my @items = @{$item_ref};
        my $vn = $items[$n];
        # if the value is a pattern do pattern match otherwise exact equality
        if ($vn =~ m/[[{.]/ && $key =~ m/$vn/ || $vn eq $key) { return $items[$m]; }
    }
    return undef;
}

# 'segment' data structure, { str, pos, len }, represents a portion of a string 

# $segment_ref = segment($str, [$pos, [$len]]])
# build a segment hash reference from the raw string, start position and the length.
# defaults: $pos == zero, $len == entire string
sub segment {
    my ($str, $pos, $len) = @_;
    return { str => $str, pos => (defined($pos) ? $pos : 0), len => (defined($len) ? $len : length($str)) };
}

sub test_segment {
    my $str = "hello";
    print Dumper(segment($str));
    print Dumper(segment($str, 1));
    print Dumper(segment($str, 1, 3));
}

# $substr_string = string($segment_ref)
# return substr'd string from a segment reference
sub string {
    my $seg = shift;
    if ($seg->{pos} == 0 && $seg->{len} == length($seg->{str})) { return $seg->{str}; } # original string
    return substr($seg->{str}, $seg->{pos}, $seg->{len});
}

sub test_string {
    my $str = "hello";
    print string(segment($str)) . "\n";
    print string(segment($str, 1)) . "\n";
    print string(segment($str, 1, 3)) . "\n";
}

# $data = readfasta($path)
# read a fasta file and return pure data without header and newlines as a single line
# if $path is missing, standard input is used (if $path is undef, "" or 0)
sub readfasta {
    my $path = shift;
    my $data = '';
    my $fh;

    if ($path) {
        open($fh, '<', $path) or die "readfasta: open: '$path' $!";
    }
    while ( my $line = $path ? <$fh> : <STDIN>) {
        if ($line =~ /^>/) { next; }
        chomp $line;
        $data .= $line;
    }
    if ($path) {
        close $fh or die "readfasta: close: '$path' $!";
    }

    return $data;
}

sub test_readfasta {
    system('echo ">header\nabcdefg\nhijklmn\n" > test.fasta');
    print bioutil::readfasta('test.fasta');
    print "\n\n";
    # no header
    system('echo "abcdefg\nhijklmn\n" > test.fasta');
    print bioutil::readfasta('test.fasta');
    print "\n\n";
    # no trailing newline
    system('echo "abcdefg\nhijklmn" > test.fasta');
    print bioutil::readfasta('test.fasta');
    print "\n\n";
    # stdin
    print bioutil::readfasta();
    print "\n\n";
}

#test_segment();
#test_string();
#test_readfasta();

true; # need to end with a true value
