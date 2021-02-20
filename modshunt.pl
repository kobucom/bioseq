#!/usr/bin/env perl

# ModsHunt - find differences (and identities) in two bio-sequences
# Copyright (c) 2021 Kobu.Com. All rights reserved.
# Visit Kobu.Com at https://kobu.com/bioseq/
# Licensed by GNU Public License v3
# https://www.gnu.org/licenses/gpl-3.0.en.html

# ModsHunt searches two input sequences and shows differences in them
# Core functions used in ModsHunt, lcss and split3, are moved to huntutil.pm.
# This file contains a variety of 'show' functions for comparison result output. 
# The default format is 'single-view' below.
#
# Example:
# $ cat orginal.fasta
#  >original sequence
#  ATATATCAGAGddd
#  AGCArrrGAGAGC
# $ cat variant.fasta
#  >modified sequence
#  ATATATCAGAG
#  AGCARRRGAGAGCiii
# $ modshunt original.fasta variant.fasta
#  ATATATCAGAG{12-ddd}AGCA{19^rrr/RRR}GAGAGC{28+iii}
#
# ModsHunt supports four output formats. 
# See usage() below for what format options are available.
#
# Note: The name 'modshunt' came from the name of a text document diff tool used in our house very long time.
# 'mods' or 'modifications' to a text document corresponds to mutations in a bio-sequence.
#
# 2021/01/28 started
# 2021/01/29 lcss()
# 2021/01/31 split3()
# 2021/02/01 single view
# 2021/02/04 ranges in tab-separated, json
# 2021/02/05 double view
# 2021/02/06 lcss() and split3() moved to huntutil.pm

use strict;
use warnings;

use List::Util qw(min);

# support library
use lib "$ENV{BIOSEQ}";
use bioutil qw(readfasta segment string shorter);
use huntutil qw(split3 rangeName SYNC DIFF INS DEL);
# use huntutil_sl qw(split3 SYNC DIFF INS DEL);
use buffer;

# constants
use constant { true => 1, false => 0, DBL_BUF_LEN => 70 };

# debug
use Data::Dumper qw(Dumper);

my $debug = 0;

{
    # double-view output function uses 'buffer' class implementing space-padded
    # line buffer. See buffer.pm.

    # An output 'line' consists of four 'lanes'
    #  0   4   8   12  16           <- upper, num
    #  AAAABBBBCCCCddddEEEE         <- upper, seq
    #  AAAArrrrCCCC    EEEEiiii     <- lower, seq
    #                  12  16       <- lower, num
    my $upper = { seq => new buffer(DBL_BUF_LEN), num => new buffer(DBL_BUF_LEN) };
    my $lower = { seq => new buffer(DBL_BUF_LEN), num => new buffer(DBL_BUF_LEN) };

    # output all lane buffers then reset them for the next line
    sub flushAll {
        my $un = $upper->{num}->{str}; $un =~ s/ +$//;
        my $us = $upper->{seq}->{str}; $us =~ s/ +$//;
        my $ls = $lower->{seq}->{str}; $ls =~ s/ +$//;
        my $ln = $lower->{num}->{str}; $ln =~ s/ +$//;
        print "$un\n";
        print "$us\n";
        print "$ls\n";
        print "$ln\n\n";
        $upper->{num}->clear();
        $upper->{seq}->clear();
        $lower->{seq}->clear();
        $lower->{num}->clear();
    }

    # setLabel($upper_or_lower, $pos, $capacity, $label)
    # set the label in the number lane
    sub setLabel {
        my ($buf, $pos, $capacity, $label) = @_;
        my $len = length($label);
        if ($len < $capacity) {
            substr($buf->{num}->{str}, $pos, $len, $label);
        }
    }

    # showSegment($longer, $shorter, $seg1, $seg2, [$fill])
    # $longer/$shorter: either of upper or lower langes
    # $seg1: longer sequence, $seg2: shorter sequence 
    # $longer should be upper if $seg1 is src, lower if $seg1 is dst
    # if $fill is defined, $shorter region is filled with that character
    # instead of the sequence data.
    sub showSegment {
        my ($longer, $shorter, $seg1, $seg2, $fill) = @_;
        $shorter->{seq}->{pos} = $longer->{seq}->{pos}; # sync buffer positions
        # current positions in segments
        my $pos = $seg1->{pos};
        my $pos2 = $seg2->{pos};
        # remaining lengths to go
        my $restLen = $seg1->{len};
        my $restLen2 = $seg2->{len}; # same if sync
        my $first = true; # controls number label output
        do {
            # calculate set length for this line
            my $spaceLeft = $longer->{seq}->capacity();
            if ($spaceLeft == 0) {
                flushAll(); # buffers automatically sync'ed
                $spaceLeft = $longer->{seq}->capacity();
            }
            my $len = min($spaceLeft, $restLen);
            my $bufpos = $longer->{seq}->{pos}; # sync shorter to longer
            # show position labels
            if ($first) {
                setLabel($longer, $bufpos, $len, "" . ($pos + 1));   # longer
                setLabel($shorter, $bufpos, $len, "" . ($pos2 + 1)); # shorter
                $first = false;
            }
            # longer side
            substr($longer->{seq}->{str}, $bufpos, $len, substr($seg1->{str}, $pos, $len)); # set portion
            $pos += $len; # adjust position and ...
            $restLen -= $len; # length
            # shorter side
            my $len2 = min($len, $restLen2);
            if (!defined($fill)) { # no fill - output sequence data
                substr($shorter->{seq}->{str}, $bufpos, $len2, substr($seg2->{str}, $pos2, $len2));
            }
            elsif ($fill eq ' ') {
                # already space-padded
            }
            else { # pad with specified fill character
                substr($shorter->{seq}->{str}, $bufpos, $len2, $fill x $len2);
            }
            $restLen2 -= $len2;
            $pos2 += $len2;
            # advance buffer pointer
            $longer->{seq}->{pos} += $len;
        } while ($restLen > 0);
        $shorter->{seq}->{pos} = $longer->{seq}->{pos}; # sync buffer positions
    }

    # showRange($range)
    sub showRange {
        my $r = shift;
        if ($r->{type} eq SYNC) { # src
            showSegment($upper, $lower, $r->{src}, $r->{dst}, '.');
        }
        elsif ($r->{type} eq DIFF) { # both, longer side as main
            if ($r->{dst}->{len} > $r->{src}->{len}) {
                showSegment($lower, $upper, $r->{dst}, $r->{src});
            }
            else {
                showSegment($upper, $lower, $r->{src}, $r->{dst});
            }
        }
        elsif ($r->{type} eq INS) { # dst
            showSegment($lower, $upper, $r->{dst}, $r->{dst}, ' ');
        }
        elsif ($r->{type} eq DEL) { # src
            showSegment($upper, $lower, $r->{src}, $r->{dst}, ' ');
        }
        else {
            die "showRange: bad type $r->{type}\n"; # debug
        }
    }

    # showInDouble(\@ranges)
    # show source and target sequences in parallel (double view format)
    sub showInDouble {
        my $ranges = shift;
        foreach my $r (@{$ranges}) {
            showRange($r);
        }
        flushAll();
    }
}

# showInJson(\@ranges)
# output range list in json format
# [
#  {
#   type: [=^+-],
#   src: {
#    str: "..."
#    pos: n,
#    len: m
#   },
#   dst: {
#    str: "..."
#    pos: n,
#    len: m
#   }
#  },
#  {
#   ...
#  },
#   ...
# ]
sub showInJson {
    my $ranges = shift;
    print "[\n"; # begin array
    foreach my $r (@{$ranges}) {
        print " {\n"; # begin object
        print "  type: \"$r->{type}\",\n"; # type
        print "  src: {\n" . # begin src
              "   str: \"" . string($r->{src}) . "\",\n" .
              "   pos: " . ($r->{src}->{pos} + 1) . ",\n" .
              "   len: $r->{src}->{len}\n" .
              "  }\n"; # end src
        print "  dst: {\n" . # begin dst
              "   str: \"" . string($r->{dst}) . "\",\n" .
              "   pos: " . ($r->{dst}->{pos} + 1) . ",\n" .
              "   len: $r->{dst}->{len}\n" .
              "  }\n"; # end dst
        print " }\n"; # end object
    }
    print "]\n"; # end array
}

# showRanges(\@ranges, $x_option)
# output range list as tab-separated text
# type spos+slen dpos+dlen "src" "dst"
# partial strings are displayed if $x_option specified (debug purpose)
sub showRanges {
    my ($ranges, $x_option) = @_;
    foreach my $r (@{$ranges}) {
        my $type = $r->{type};
        my $src = string($r->{src});
        my $dst = string($r->{dst});
        if ($x_option) {
            $type = rangeName($type);
            $src = shorter($src);
            $dst = shorter($dst);
        }
        print "$type\t" . 
            ($r->{src}->{pos} + 1) . "+$r->{src}->{len}\t" .
            ($r->{dst}->{pos} + 1) . "+$r->{dst}->{len}\t\"" .
            $src . "\"\t\"" . $dst . 
            "\"\n";
    }
}

# delimeter characters
use constant { BGN_CHAR => '{', MID_CHAR => '/', END_CHAR => '}' };

# show(\@ranges)
# show source sequence intermixed with differences (single view format)
sub show {
    my $ranges = shift;
    foreach my $r (@{$ranges}) {
        if ($r->{type} eq SYNC) {
            print string($r->{src}); # src as is
        }
        else {
            print BGN_CHAR . ($r->{src}->{pos} + 1) . $r->{type};
            if ($r->{type} eq DEL) {
                print string($r->{src}); # [N-src]
            }
            elsif ($r->{type} eq INS) {
                print string($r->{dst}); # [N+dst]
            }
            elsif ($r->{type} eq DIFF) {
                print string($r->{src}) . MID_CHAR . string($r->{dst}); # [N^src~dst]
            }
            print END_CHAR;
        }
    }
    print "\n";
}

# $data = getData($arg)
# content of the file is returned as a single-line data.
# while debugging, literal data can be specified on the command-line
sub getData {
    my $arg = shift;
    if (-e $arg || !$debug) { return readfasta($arg); }
    return $arg; # inline data allowed only if $debug and file missing
}

sub usage {
    print STDERR "Usage: modshunt [-f<fmt>] reference_fasta_file target_fasta_file\n";
    print STDERR "  -fs: single view (default)  -ft: ranges in tab-separated text\n" .
                 "  -fd: double view            -fj: ranges in json\n";
    exit 1;
}

sub main {
    my $arg1 = '';
    my $arg2 = '';
    my $fmt = 's';
    foreach my $arg (@ARGV) {
        if ($arg =~ /^-/) {
            if ($arg eq '-d') { $debug = 1; }
            elsif ($arg eq '-t') { $debug = 2; }
            elsif ($arg =~ /^-f([sdtjx])$/) { $fmt = $1; }
            else {
                print STDERR "No such option: $arg\n";
                usage();
            }
        }
        elsif (!$arg1) { $arg1 = $arg; }
        elsif (!$arg2) { $arg2 = $arg; }
        else { usage(); }
    }
    if (!$arg1 || !$arg2) { usage(); }
    my $src = getData($arg1);
    my $dst = getData($arg2);
    my @ranges = split3($src, $dst);
    if ($fmt eq 't') { showRanges(\@ranges); } # ranges in tab-separated
    elsif ($fmt eq 'x') { showRanges(\@ranges, true); } # ditto but shorter strings (debug)
    elsif ($fmt eq 'j') { showInJson(\@ranges); } # ranges in json
    elsif ($fmt eq 'd') { showInDouble(\@ranges); } # double view
    else { show(\@ranges); } # single view (default)
}

main();
