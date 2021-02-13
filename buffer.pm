#!/usr/bin/env perl

# 'buffer' data structure { str, size, pos } for space-padded line buffer management

package buffer;

use strict;
use warnings;

# support library
use lib "$ENV{BIOSEQ}";

# constants
use constant { true => 1, false => 0 };

# debug
use Data::Dumper qw(Dumper);

# constructor
# $buffer_ref = new buffer($size)
# return a fixed-size space-filled string with position marker initialized to zero
sub new {
    my ($class, $size) = @_;
    if (!$size) { $size = 70; }
    return bless {
        str => ' ' x $size,
        size => $size,
        pos => 0
    }, $class;
}

# clear($buffer_ref)
# clear the buffer to initial space-filled status and reset the position marker to zero
sub clear {
    my $self = shift;
    $self->{str} = ' ' x $self->{size};
    $self->{pos} = 0;
}

# $len = capacity($buffer_ref)
# return the remaining space at the end of the line buffer
sub capacity {
    my $self = shift;
    return $self->{size} - $self->{pos};
}

true; # need to end with a true value
