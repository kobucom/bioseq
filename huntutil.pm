#!/usr/bin/env perl

# HuntUtil - ModsHunt's core functions made independent
# Copyright (c) 2021 Kobu.Com. All rights reserved.
# Visit Kobu.Com at https://kobu.com/bioseq/
# Licensed by GNU Public License v3
# https://www.gnu.org/licenses/gpl-3.0.en.html

# ModsHunt searches two input sequences and shows differences in them
# It uses a method called 'split into three' algorithm where: 
# - find the longest common substring (lcss) between the source and dest strings
# - split each string into three parts: left-side, identical (or sync) part and right-side
# - recursively find the longest common substring for left- and right- side of two strings
# - merge all found 'ranges' of identity, diference, insertion and deletions into an array
#
# Due to dependence of 'lcss' function, the comparison result will have the greatest identity and
# smallest difference.
#
# Note: This code uses not-very-appropriate words: 'src' for source and 'dst' for destination,
# meaning an original string and its comparison target, respectively.
#
# TODO: A part that can be made concurrent is marked with PARA
# 
# 2021/01/28 started
# 2021/01/29 lcss()
# 2021/01/31 split3()
# 2021/02/06 huntutil.pm
# 2021/02/09 tested
#  cut-xxxx.data  ( 3823) 10 seconds
#  nuc-xxxx.fasta (30500) 11 minutes 
# 2021/02/18 one-char match by index() instead of regexp
# 2021/02/19 global longest match variable for possible concurrency mode
#  cut-xxxx.data  ( 3823) 5 seconds
#  nuc-xxxx.fasta (30500) 7 minutes 
  
package huntutil;

use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);
@EXPORT_OK = qw(lcss split3 rangeName SYNC DIFF INS DEL);

use List::Util qw(min);

# support library
use lib "$ENV{BIOSEQ}";
use bioutil qw(segment string shorter);

# constants
use constant { true => 1, false => 0 };

# debug
use Data::Dumper qw(Dumper);

my $debug = 0; # 1: debug, 2: trace

### lcss - longest common substring

# 'match range' data structure: { len, spos, dpos }
# A match range represents a pair of src and dst regions with the same length.
# It does not include the reference source and target strings.
# A match associates identical regions in the source and target sequence.
#   abcdEFGhijk -> spos = 4, len = 3 --+--> { len => 3, spos => 4, dpos => 2 }
#   xyEFGz      -> dpos = 2, len = 3 --+
# this structure is used for two purposes:
# - result of already matched regions (len is matched length)
# - candidate regions to check for matching (len is the maximum characters to check)

# $debug_string = dbg_mch($match_range, [$src_seg, $dst_seg])
# where optional $src or $dst are string segments
sub dbg_mch {
    my ($mch, $src, $dst) = @_;
    my $str = "[$mch->{len}, $mch->{spos}, $mch->{dpos}]";
    if ($src) { $str .= " '" . shorter(string($src)) . "'"; }
    if ($dst) { $str .= " '" . shorter(string($dst)) . "'"; }
    return $str;
}

{
    # current longest match range
    # scope of this data structure is per invocation of lcss()
    # ie., shared by subordinate functions called from that instance of lcss()
    # when run sequentially this data can be freely accessed by the subordinate functions
    # when run concurrently this data must be update-exclusive among the subordinate functions
    my $longest_match = { len => 0, spos => 0, dpos => 0 };

    # clear the maximum length to zero
    sub clearLongest {
        $longest_match->{len} = 0;
    }

    # CAUTION: this must be atomic operation if run concurrently
    # setLongest($updated)
    # update the longest data
    # when run concurrently, getLongest() and the following setLongest() with a greater value
    # may overlap among tasks; that's why I check the length before updating below
    sub setLongest {
        my $res = shift;
        if ($res->{len} > $longest_match->{len}) { # avoid race condition
            $longest_match->{len} = $res->{len};
            $longest_match->{spos} = $res->{spos};
            $longest_match->{dpos} = $res->{dpos};
        }
    }

    # return the longest match range
    # Note: the returned value is a 'copy' of the original value, otherwise the caller
    # will see the different value if it calls lcss() later.
    sub getLongest {
        return { len => $longest_match->{len}, spos => $longest_match->{spos}, dpos => $longest_match->{dpos} };
    }

    # return the length of the longest match range
    sub getLongestLength {
        return $longest_match->{len};
    }

    # TODO: use 'segment' data structure in lcss() too

    # tryFixedRange($src, $dst, $candidate);
    # finds the longest match at the fixed start positions up to the length specified in $candidate range
    # it does not return value, instead, it updates the longest match value if it finds the longest.
    # it does nothing if the largest sync part is shorter than the curren longest value it sees.
    # assumption: called only if the candidate length is larger than the longest at the time of call
    sub tryFixedRange {
        my ($src, $dst, $cnd) = @_;
        my $res = { len => 0, spos => $cnd->{spos}, dpos => $cnd->{dpos} };
        my $lower = getLongestLength(); # lower bound

        # check src and dst string at specified positions by increasing check range up to the specified length
        for (my $len = $lower > 0 ? $lower : 1; $len <= $cnd->{len}; $len++) {
            if (substr($src, $cnd->{spos}, $len) eq substr($dst, $cnd->{dpos}, $len)) {
                $res->{len} = $len;
                # Note: why don't update the longest here
                # I chose to do it at the end only once to avoid too many commit ops
            }
            else {
                # no more sync
                last;
            }
        }

        if ($res->{len} > 0) {
            # set the longest down here
            setLongest($res);

            my $str = substr($src, $cnd->{spos}, $res->{len});
            if ($debug) { print STDERR "tryFixedRange: '$str' ($res->{len}) at $cnd->{spos} against $cnd->{dpos}\n"; }
        }
        else {
            if ($debug > 1) { print STDERR "tryFixedRange: none\n"; }
        }
    }

    # tryCharAgnstAll($src, $dst, $spos);
    # this is called for each character in src that is checked against each matching char in dst to
    # find a matching string starting at the positions.
    # this function does not return a value but the subordinate function, tryFixedRange() updates the longest match.
    sub tryCharAgnstAll {
        my ($src, $dst, $spos) = @_;

        my $slen = length($src);
        my $dlen = length($dst);
        my $c = substr($src, $spos, 1);

        if ($debug > 1) { print STDERR "tryCharAgnstAll: check '$c' at $spos\n"; }

        my $cnd = { spos => $spos };

        # for each occurence of one-char match in dst and check the sync part from that position
        for (my $dpos = 0; $dpos < $dlen; $dpos++) {
            $dpos = index($dst, $c, $dpos); # where one-char match found
            if ($dpos == -1) { last; } # no more one-char match
            my $len = min($slen - $spos, $dlen - $dpos);
            my $lower = getLongestLength(); # lower bound
            if ($debug > 1) {
                print STDERR "tryCharAgnstAll: min-len=$len, lower=$lower\n";
            }
            if ($lower > 0 && $len <= $lower) { next; } # skip if shorter
            $cnd->{dpos} = $dpos;
            $cnd->{len} = $len;
            tryFixedRange($src, $dst, $cnd, $lower); # PARA - but not much improvement
        }
    }

    # { len, spos, dpos } = lcss($src, $dst)
    # compares two strings and returns the length, source and dest positions of the longest common substring.
    # $src and $dst are just strings, not segments as used in split3(), and return value is a hash reference.
    sub lcss {
        my ($src, $dst) = @_;
        if ($debug) { print STDERR "lcss2: '" . shorter($src) . "' < '" . shorter($dst) . "'\n"; }

        clearLongest();

        my $slen = length($src);
        my $dlen = length($dst);
        my $res = { len => 0, spos => 0, dpos => 0 };
        
        # for each char in src vs entire dst
        for (my $spos = 0; $spos < $slen ; $spos++) {
            tryCharAgnstAll($src, $dst, $spos); # PARA by divided regions of certain length
        }

        my $longest = getLongest();
        if ($longest->{len} > 0) {
            if ($debug) {
                print STDERR "lcss2: " . dbg_mch($longest, segment($src, $longest->{spos}, $longest->{len})) . "\n";
            }
            return $longest;
        }
        else {
            if ($debug > 1) {
                print STDERR "lcss2: none\n";
            }
            return undef;
        }
    }
}

sub test_lcss {
    $debug = 2;

    my $dst = "The quick brown fox jumps over the lazy dog";

    my @data = (
        # # short diff
        # { src => "a", dst => "b" },
        # { src => "a", dst => "bb" },
        # { src => "aa", dst => "b" },
        # { src => "aa", dst => "bb" },

        # # short sync
        # { src => "a", dst => "a" },
        # { src => "aa", dst => "aa" },
        # { src => "aaa", dst => "aaa" },

        # # first, middle and last
        # { src => "The quick", dst => $dst },
        # { src => "lazy dog", dst => $dst },
        # { src => "jumps over", dst => $dst },

        # # one sync parts
        # { src => "abcDEFxyz", dst => "123DEF789" },
        # { src => "aDEFbc", dst => "123DEF89" },

        # two sync parts
        { src => "abcDEFghiJKLMnop", dst => "12DEF34JKLM89" },
        { src => "xxxXXXxxxXXXXXXxxx", dst => "yyyXXXyyyXXXXXXyyy" }
    );

    foreach my $pair (@data) {
        lcss($pair->{src}, $pair->{dst});
    }
}

### split3 - recursive split-into-three diff

# 'segment' data structure, { str, pos, len }, represents a portion of a string.
# See comment in bioutil.pm for detail.

# debug only
# $debug_string = dbg_seg($segment)
# segment shown as "hello world" or "llo wor"2+7
sub dbg_seg {
    my $seg = shift;
    if ($seg->{pos} == 0 && $seg->{len} == length($seg->{str})) { return "'" . shorter($seg->{str}) . "'"; }
    return "'" . shorter(string($seg)) . "'$seg->{pos}+$seg->{len}";
}

# 'range' data structure, { type, src, dst }, associates a pair of string segments that correspond in any way.
# $src is a segment of a source string.
# $dst is a segment of a target string that matches to the source segment.
# types are:
# - sync range (=)
# - diff range (~: replacement)
# - insert range (+: dst-only; added to dst)
# - delete range (-: src-only; existed in src but deleted)

use constant { SYNC => '=', DIFF => '^', INS => '+', DEL => '-' };

# $type_name = rangeName($type_char)
sub rangeName {
    my $type = shift;
    if ($type eq '=') { return 'sync'; }
    elsif ($type eq '^') { return 'diff'; }
    elsif ($type eq '+') { return 'ins'; }
    elsif ($type eq '-') { return 'del'; }
    else { return '???'; }
}

# debug only
# $debug_string = dbg_rng($range)
# sync: "hello" = "hello"
# diff: "hello" ~ "goodbye"
# ins:  "hello" + 
# del:  "hello" - 
sub dbg_rng {
    my $rng = shift;
    return rangeName($rng->{type}) . ' ' . dbg_seg($rng->{src}) . ' ' . dbg_seg($rng->{dst});
}

{
    # debug
    my $recurCount = 0;

    # $range_array_ref = split3($src_seg, $dst_seg)
    # compare two string segments and return an array of ranges each representing sync, diff, ins or del.
    # Note: the third parameter, $depth, is used for debugging purpose
    sub split3 {
        my ($src, $dst, $depth) = @_;
        if (!defined($depth)) { $depth = 0; } # on initial call
        if (ref($src) eq "") { $src = segment($src); } # if raw string passed, make it segment ref
        if (ref($dst) eq "") { $dst = segment($dst); }

        # debug
        if ($depth > $recurCount) { $recurCount = $depth; }
        # if ($depth > 3) { die "abort recursion at 100\n"; }

        if ($debug) { print STDERR "split3: entry($depth) " . dbg_seg($src) . " < " . dbg_seg($dst) . "\n"; }

        # find the longest common substring
        my $mid = lcss(string($src), string($dst));

        # if totally different, return diff (~)
        if (!$mid) {
            my $diff = { type => DIFF, src => $src, dst => $dst };
            if ($debug) { print STDERR "split3: exit($depth) middle -> " . dbg_rng($diff) . "\n"; }
            return ( $diff );
        }

        # otherwise, sync part (=) exists whether entire or partial
        my $sync = { type => SYNC,
            src => segment($src->{str}, $src->{pos} + $mid->{spos}, $mid->{len}),
            dst => segment($dst->{str}, $dst->{pos} + $mid->{dpos}, $mid->{len}) };
        if ($debug) { print STDERR "split3: middle -> " . dbg_rng($sync) . "\n"; }

        # set the sync part in the middle and try left- and right-side
        my @ranges = ( $sync );

        # handle left-side segments

        my $srcL = segment($src->{str}, $src->{pos}, $mid->{spos});
        my $dstL = segment($dst->{str}, $dst->{pos}, $mid->{dpos});

        if ($srcL->{len} == 0 && $dstL->{len} == 0) { # no left side, go through
            if ($debug) { print STDERR "split3: left -> none\n"; }
        }
        elsif ($srcL->{len} > 0 && $dstL->{len} == 0) { # if only src exists, it's delete (-)
            my $del = { type => DEL, src => $srcL, dst => $dstL };
            if ($debug) { print STDERR "split3: left -> " . dbg_rng($del) . "\n"; }
            unshift(@ranges, $del);
            #@ranges = ( $del, @ranges );
        }
        elsif ($srcL->{len} == 0 && $dstL->{len} > 0) { # if only dst exists, it's insert (+)
            my $ins = { type => INS, src => $srcL, dst => $dstL };
            if ($debug) { print STDERR "split3: left -> " . dbg_rng($ins) . "\n"; }
            unshift(@ranges, $ins);
            #@ranges = ( $ins, @ranges );
        }
        else { # if both src and dst exist, the result of a recursive call to split3()
            if ($debug) { print STDERR "split3: left -> " . dbg_seg($srcL) . " < " . dbg_seg($dstL) . "\n"; }
            unshift(@ranges, split3($srcL, $dstL, $depth + 1)); # PARA - left and right side run parallel
            #@ranges = ( split3($srcL, $dstL, $depth + 1), @ranges );
        }

        # next, handle right-side segments

        my $srcRpos = $src->{pos} + $mid->{spos} + $mid->{len};
        my $dstRpos = $dst->{pos} + $mid->{dpos} + $mid->{len};
        my $srcRlen = $src->{len} - ($mid->{spos} + $mid->{len});
        my $dstRlen = $dst->{len} - ($mid->{dpos} + $mid->{len});
        my $srcR = segment($src->{str}, $srcRpos, $srcRlen);
        my $dstR = segment($dst->{str}, $dstRpos, $dstRlen);

        if ($srcRlen == 0 && $dstRlen == 0) { # no right side, go through
            if ($debug) { print STDERR "split3: right -> none\n"; }
        }
        elsif ($srcRlen > 0 && $dstRlen == 0) { # if only src exists, it's delete (-)
            my $del = { type => DEL, src => $srcR, dst => $dstR };
            if ($debug) { print STDERR "split3: right -> " . dbg_rng($del) . "\n"; }
            push(@ranges, $del);
            #@ranges = ( @ranges, $del );
        }
        elsif ($srcRlen == 0 && $dstRlen > 0) { # if only dst exists, it's insert (+)
            my $ins = { type => INS, src => $srcR, dst => $dstR };
            if ($debug) { print STDERR "split3: right -> " . dbg_rng($ins) . "\n"; }
            push(@ranges, $ins);
            #@ranges = ( @ranges, $ins );
        }
        else { # if both src and dst exist, the result of a recursive call to split3()
            if ($debug) { print STDERR "split3: right -> " . dbg_seg($srcR) . " < ". dbg_seg($dstR) . "\n"; }
            push(@ranges, split3($srcR, $dstR, $depth + 1)); # PARA
            # @ranges = ( @ranges, split3($srcR, $dstR, $depth + 1) ); 
        }

        if ($debug) {
            print STDERR "split3: exit($depth):\n";
            foreach my $r (@ranges) {
                print STDERR dbg_rng($r) . "\n";
            }
            if ($debug > 1) { print STDERR Dumper(\@ranges); }
        }

        # debug
        if ($depth == 0) { print STDERR "recurCount: $recurCount\n"; }

        return @ranges;
    }
}

sub test_split3 {
    $debug = 1;
    my @data = (
        # # total sync
        # { src => segment("aaaa"), dst => segment("aaaa") },

        # # total diff
        # { src => "aaaa", dst => "bbbb" },

        # # sync in the middle
        { src => "aaaaAAAAaaaa", dst => "bbbbAAAAbbbb" },

        # # two syncs in the middle
        # { src => "aAAaaAAAAaaa", dst => "bbAAbbbAAAAbbb" },

        # # insert
        # { src => "AAAA", dst => "aaaaAAAA" },
        # { src => "AAAA", dst => "AAAAaaaa" },
        # { src => "AAAAAAAA", dst => "AAAAaaaaAAAA" },

        # # deletes (reverse of above)
        # { src => "aaaaAAAA", dst => "AAAA" },
        # { src => "AAAAaaaa", dst => "AAAA" },
        # { src => "AAAAaaaaAAAA", dst => "AAAAAAAA" },
        # { src => "bbaAAabbbaAAAAaaabbbb", dst => "aAAaaAAAAaaa" }
    );

    foreach my $pair (@data) {
        split3($pair->{src}, $pair->{dst});
        print "called with '$pair->{src}' < '$pair->{dst}'\n";
    }
}

# test_lcss();
# test_split3();

true; # need to end with a true value
