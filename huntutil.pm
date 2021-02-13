#!/usr/bin/env perl

# HuntUtil - ModsHunt's core functions made independent
# Copyright (c) 2021 Kobu.Com. All rights reserved.
# Visit Kobu.Com at https://kobu.com/bioseq/
# Licensed by GNU Public License v3
# https://www.gnu.org/licenses/gpl-3.0.en.html

# ModsHunt searches two input sequences and shows differences in them
# It uses a method called 'split into three' algorithm where: 
# - find the longest common substring (lcss) between the source and dest strings
# - split each string into three parts: left-side, common part and right-side
# - recursively find the longest common substring for left- and right- side of two strings
# The comparison result will have the greatest identity and smallest difference.
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

package huntutil;

use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);
@EXPORT_OK = qw(lcss split3 SYNC DIFF INS DEL);

use List::Util qw(min);

# support library
use lib "$ENV{BIOSEQ}";
use bioutil qw(segment);

# constants
use constant { true => 1, false => 0 };

# debug
use Data::Dumper qw(Dumper);

my $debug = 0; # 1: debug, 2: trace

### lcss - longest common substring

# { len, spos, dpos } = lcss($src, $dst)
# compares two strings and returns the length, source and dest positions of the longest common substring.
# $src and $dst are just strings, not segments as used in split3(), and return value is a hash reference.
# Note: the $src should not contain any regexp special characters.
sub lcss {
    my ($src, $dst) = @_;
    if ($debug) { print STDERR "lcss: '$src' < '$dst'\n"; }
    my $res = { len => 0, spos => 0, dpos => 0 };
    my $slen = length($src);
    my $dlen = length($dst);

    # for each char in src, find one or more occurrences in dst
    for (my $spos = 0; $spos < $slen; $spos++) { # PARA - for each src char
        my $c = substr($src, $spos, 1); # single char pattern
        my $line = $dst;
        my $dpos = 0;
        while (true) {
            if ($line =~ /$c/p) { # PARA - build list first then do them parallel
                my $prematch = ${^PREMATCH};
                my $match = ${^MATCH};
                my $postmatch = ${^POSTMATCH};
                $dpos += length($prematch);
                my $matchLen = length($match); # == 1
                my $postLen = length($postmatch);
                my $restLen = min($matchLen + $postLen, $slen - $spos);
                if ($debug > 1) { print STDERR "lcss: search '$c' from $spos against $dpos length $restLen\n"; }
                # for each occurrence of single char match, find the longest equal strings starting at that position
                for (my $len = $matchLen; $len <= $restLen; $len++) {
                    if (substr($src, $spos, $len) eq substr($dst, $dpos, $len)) { # increase check length
                        if ($debug > 1) {
                            print STDERR "lcss: matched '" . substr($src, $spos, $len) . "' ($len)\n";
                        }
                        # update the longest match record
                        if ($len > $res->{len}) {
                            if ($debug > 1) { print STDERR "lcss: longest: $len\n"; }
                            $res->{len} = $len;
                            $res->{spos} = $spos;
                            $res->{dpos} = $dpos;
                        }
                    }
                    else {
                        # no need to check the rest
                        last;
                    }
                }
                $dpos += $matchLen;
                $line = $postmatch;
            }
            else {
                # ignore rest in $line
                last;
            }
        }
    }

    if ($debug) {
        print STDERR "lcss: '" . substr($src, $res->{spos}, $res->{len}) . "'$res->{spos}+$res->{len} at $res->{dpos}\n";
    }

    return $res;
}

sub test_lcss {
    $debug = 2;

    # short diff
    # lcss("a", "b");
    # lcss("a", "bb");
    # lcss("aa", "b");
    # lcss("aa", "bb");

    # short sync
    # lcss("a", "a");
    # lcss("aa", "aa");
    # lcss("aaa", "aaa");

    # first, middle and last
    my $dst = "The quick brown fox jumps over the lazy dog";
    # lcss("The quick", $dst);
    # lcss("lazy dog", $dst);
    # lcss("jumps over", $dst);

    # one sync parts
    # lcss("abcDEFxyz", "123DEF789");
    # lcss("aDEFbc", "123DEF89");

    # two sync parts
    # lcss("abcDEFghiJKLMnop", "12DEF34JKLM89");
    lcss("xxxXXXxxxXXXXXXxxx", "yyyXXXyyyXXXXXXyyy");
}

### split3 - recursive split-into-three diff

# modshunt uses 'segment' data structure, { str, pos, len }, to represent
# a portion of a string. See comment in bioutil.pm for detail.

# debug only
# $debug_string = dbg_seg($segment)
# segment shown as "hello world" or "llo wor"2+7
sub dbg_seg {
    my $seg = shift;
    if ($seg->{pos} == 0 && $seg->{len} == length($seg->{str})) { return "'" . $seg->{str} . "'"; }
    return "'" . substr($seg->{str}, $seg->{pos}, $seg->{len}) . "'$seg->{pos}+$seg->{len}" ;
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

# $type_name = rangeType($type_char)
sub rangeType {
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
    return dbg_seg($rng->{src}) . " " . rangeType($rng->{type}) . " " . dbg_seg($rng->{dst});
}

# $range_array_ref = split3($src_seg, $dst_seg)
# compare two string segments and return an array of ranges each representing sync, diff, ins or del.
# Note: the third parameter, $depth, is used for debugging purpose
sub split3 {
    my ($src, $dst, $depth) = @_;
    if (!defined($depth)) { $depth = 0; } # on initial call
    if (ref($src) eq "") { $src = segment($src); } # if raw string passed, make it segment ref
    if (ref($dst) eq "") { $dst = segment($dst); }

    if ($debug) { print STDERR "split3: entry($depth) " . dbg_seg($src) . " < " . dbg_seg($dst) . "\n"; }

    # find the longest common substring
    my $mid = lcss(substr($src->{str}, $src->{pos}, $src->{len}), substr($dst->{str}, $dst->{pos}, $dst->{len}));

    # if totally different, return diff (~)
    if ($mid->{len} == 0) {
        my $diff = { type => DIFF, src => $src, dst => $dst };
        if ($debug) { print STDERR "split3: exit($depth) total diff: " . dbg_rng($diff) . "\n"; }
        return ( $diff );
    }

    # otherwise, sync part (=) exists whether entire or partial
    my $sync = { type => SYNC,
        src => segment($src->{str}, $src->{pos} + $mid->{spos}, $mid->{len}),
        dst => segment($dst->{str}, $dst->{pos} + $mid->{dpos}, $mid->{len}) };
    if ($debug) { print STDERR "split3: mid sync: " . dbg_rng($sync) . "\n"; }

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
        if ($debug) { print STDERR "split3: left -> del: " . dbg_rng($del) . "\n"; }
        unshift(@ranges, $del);
        #@ranges = ( $del, @ranges );
    }
    elsif ($srcL->{len} == 0 && $dstL->{len} > 0) { # if only dst exists, it's insert (+)
        my $ins = { type => INS, src => $srcL, dst => $dstL };
        if ($debug) { print STDERR "split3: left -> ins: " . dbg_rng($ins) . "\n"; }
        unshift(@ranges, $ins);
        #@ranges = ( $ins, @ranges );
    }
    else { # if both src and dst exist, the result of a recursive call to split3()
        if ($debug) { print STDERR "split3: left split3: " . dbg_seg($srcL) . " < " . dbg_seg($dstL) . "\n"; }
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
        if ($debug) { print STDERR "split3: right -> del: " . dbg_rng($del) . "\n"; }
        push(@ranges, $del);
        #@ranges = ( @ranges, $del );
    }
    elsif ($srcRlen == 0 && $dstRlen > 0) { # if only dst exists, it's insert (+)
        my $ins = { type => INS, src => $srcR, dst => $dstR };
        if ($debug) { print STDERR "split3: right -> ins: " . dbg_rng($ins) . "\n"; }
        push(@ranges, $ins);
        #@ranges = ( @ranges, $ins );
    }
    else { # if both src and dst exist, the result of a recursive call to split3()
        if ($debug) { print STDERR "split3: right split3: " . dbg_seg($srcR) . " < ". dbg_seg($dstR) . "\n"; }
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

    return @ranges;
}

sub test_split3 {
    $debug = 1;
    # total sync
    #split3(segment("aaaa"), segment("aaaa"), 0);

    # total diff
    #split3(segment("aaaa"), segment("bbbb"), 0);

    # sync in the middle
    #split3(segment("aaaaAAAAaaaa"), segment("bbbbAAAAbbbb"), 0);

    # two syncs in the middle
    # split3("aAAaaAAAAaaa", "bbAAbbbAAAAbbb");

    #####

    # insert
    # split3("AAAA", "aaaaAAAA");
    # split3("AAAA", "AAAAaaaa");
    # split3("AAAAAAAA", "AAAAaaaaAAAA");

    # delets (reverse of above)
    #split3("aaaaAAAA", "AAAA");
    # split3("AAAAaaaa", "AAAA");
    split3("AAAAaaaaAAAA", "AAAAAAAA");
    # split3("bbaAAabbbaAAAAaaabbbb", "aAAaaAAAAaaa");
}

# test
# test_lcss();
# test_split3();

true; # need to end with a true value
