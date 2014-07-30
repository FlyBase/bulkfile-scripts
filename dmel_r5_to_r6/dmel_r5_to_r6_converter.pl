#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use File::Spec;

=pod

=head1 NAME

dmel_r5_to_r6_converter.pl - A coorindate converter for R5 to R6 for Dmel.

=head1 SYNOPSIS

 dmel_r5_to_r6_converter.pl [options]

 Options:
 --output The file to write R6 coordinates to.
 --input The file to read R5 coordinates from.
 --mapfile The R5 to R6 coordinate mapping source file.
 --help  Print a help message
 --debug Print debug information to STDERR
 --man Show a man page help doc.
    
see L<Options> for full details.

 e.g.
 cat my_r5_coordinates.txt | ./dmel_r5_to_r6_converter.pl > my_r6_coordinates.txt 

 #using bash STDERR redirection
 cat my_r5_coordinates.txt | ./dmel_r5_to_r6_converter.pl > my_r6_coordinates.txt 2> conversion_failures.txt

 ./dmel_r5_to_r6_converter.pl --input my_r5_coordinates.txt --output my_r6_coordinates.txt

 ./dmel_r5_to_r6_converter.pl --input my_r5_coordinates.txt --mapfile /path/to/my/own/mapping/file.tsv > my_r6_coordinates.txt

=head1 DESCRIPTION

This script takes a list of scaffold coordinates for release 5.x (R5)
of the D. melanogaster assembly and converts them to their release 6.x (R6) equivalent.
Input can be supplied via STDIN or a file using the --input option.
If no output file is specified, converted coordinates will be sent to STDOUT
and any failures will be directed to STDERR.

If you specify an output file, all R6 coordinates will be sent to the file
and all failures will be placed alongside it in a file called '<filename>.failed'.

=head2 Options

=over 5

=item --output <file>

The output file used for writing out R6 coordinates.  Defaults to STDOUT.

=item  --input <file>

The input file for R5 coordinates.  Defaults to STDIN.

=item --mapfile <file> 

The R5 to R6 mapping file to use.  Defaults to the mapping file distributed
with this script and placed in the same directory as this executable.
Only use this option if you need to read the mapping file
from another location.

DO NOT USE another mapping file unless you know what you are doing.

=item --help 

Print a help page.

=item --debug

Print debug information to STDERR during conversion.

=item --man

Show help as a man page.

=back

=head2 Coordinate input format

All lines that do not meet these requirements will be ignored.

When indicating a range, if the start position is larger than the stop,
the integers will be swapped and converted.

e.g.

2R:10000..900 becomes 2R:900..10000 before conversion.

Coordinates must be formatted accordingly:

=over 5

=item 1. No more than one coordinate position or coordinate range per line.

=item 2. Scaffold names must use FlyBase approved values.  A 'Chr' prefix will be ignored and removed from output.

=item 3. Location integers must only consist of digits 0-9 and/or commas ','.

=item 4. A range must be indicated by two integers seperated by two periods '..'. 

=item 5. Scaffold names must be separated from coordinate(s) by a colon ':'.

=back

=head3 Examples

  2R:1..10000
  ChrX:330..19020
  4:43943..50320
  3R:3,020,320
  chr4:123020

=head2 Output file format

The output file is a tab delimited file with a header containing counts
of coordinates found, counts of coordinates successfully converted, and 
counts of all failures.  

=over 5

=item Column 1

The original coordinates submitted.

=item Column 2

The converted R6 coordinate position or range.

=item Column 3

A notes field with useful information about the nature of the conversion.

This is used in situations such as inversions, coordinate ranges that span
a region(s) of change, change of scaffold name, and others.

=back

=head2 Error file format

Roughly the same as the output, but without the header.  It may include
skipped lines of varying formats that failed format validation.

See L<Output file format>.

=head1 AUTHOR

=over 5

=item Victor Strelets, FlyBase

=item Josh Goodman, FlyBase

=back


=head1 LICENSE

 Copyright (c) 2014, Indiana University & FlyBase 
 All rights reserved. 

 Redistribution and use in source and binary forms, with or without modification, 
 are permitted provided that the following conditions are met: 

  * Redistributions of source code must retain the above copyright 
    notice, this list of conditions and the following disclaimer. 
  * Redistributions in binary form must reproduce the above copyright 
    notice, this list of conditions and the following disclaimer in the 
    documentation and/or other materials provided with the distribution. 
  * Neither the name of Indiana University, Bloomington nor the names 
    of its contributors may be used to endorse or promote products 
    derived from this software without specific prior written 
    permission. 

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

=cut

my $debug = 0;
my $help  = 0;
my $man   = 0;
my $output;
my $input   = '-';
my $mapfile = File::Spec->catfile($FindBin::Bin,'dmel_r5_to_r6_mapping.tsv');

my $getopt = GetOptions(
    'help|?'    => \$help,
    'debug'     => \$debug,
    'man'       => \$man,
    'output=s'  => \$output,
    'input=s'   => \$input,
    'mapfile=s' => \$mapfile,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my $input_fh;

if ( $input eq '-' ) {
    $input_fh = \*STDIN;
}
else {
    open( $input_fh, '<', $input ) or die "Can't open input file $input: $!\n";
}

my $out_fh = \*STDOUT;
my $failed_fh = \*STDERR;

if ( defined $output ) {
    open( $out_fh, '>', $output ) or die "Can't open output file $output: $!\n";
    open( $failed_fh, '>', $output .'.failed' ) or die "Can't open error file " . $output . '.failed: ' . "$!\n";
}
my @coords = ();
my $skipped_num = 0;

while (<$input_fh>) {
    chomp;
    if ($_ !~ /^\w+:[\d,.]+$/) {
        $skipped_num++;
        print $failed_fh "$_\t#Skipped invalid scaffold location\n";
        next;
    }
    $_ =~ tr/,//;                  #Strip commas
    $_ =~ s/^Chr([^:]+:)/$1/io;    # Strip out 'chr' prefix.
    push( @coords, $_ );
}

close($input_fh);

my $chr_map = setConversionArrays( { version => 5, mapfile => $mapfile } );

my @converted_coords = ();
my @notes            = ();

foreach my $coord (@coords) {
    my ( $converted, $note ) = getConvertedCoord(
        { version => 5, coord => $coord, chr_map => $chr_map } );
    push( @converted_coords, $converted );
    push( @notes,            $note );
}

my $num_found  = scalar @coords;
my $num_failed = scalar grep { /failed/ } @notes;
my $num_conv   = $num_found - $num_failed;

print $out_fh "#Number of properly formatted coordinates found = $num_found\n";
print $out_fh "#Number of skipped lines = $skipped_num\n";
print $out_fh "#Number of coordinates converted = $num_conv\n";
print $out_fh "#Number of failed conversions = $num_failed\n";
print $out_fh "#Original\tConverted\tNotes\n";

for ( my $n = 0 ; $n <= $#coords ; $n++ ) {
    my $coord1 = $coords[$n];
    my ( $chr, $start, $end ) = splitRange($coord1);
    my $coord2 = $converted_coords[$n];
    $coord2 = $chr . ':?..?' unless $coord2 =~ /^\s*[^ :]+:[^?-]*\d/;
    if ($notes[$n] !~ /failed/) {
        print $out_fh $coord1 . "\t" . $coord2 . "\t" . $notes[$n] . "\n";
    }
    else {
        print $failed_fh $coord1 . "\t" . $coord2 . "\t" . $notes[$n] . "\n";
    }
}
close($out_fh);
close($failed_fh);

exit;

#***********************************************************
#
#***********************************************************

sub setConversionArrays {
    my ($args) = @_;

    my @sss = ();
    open( my $map_fh, '<', $args->{mapfile} )
      || die "Cannot open R5 to R6 coorindate mapping data file "
      . $args->{mapfile}
      . ": $1\n";
    while (<$map_fh>) {
        my ( $chr, $start, $end, $newchr, $newstart, $newend, $strand, @rest ) =
          split( /[ \t]+/, $_ );
        ( $start,    $end )    = ( $end,    $start )    if $start > $end;
        ( $newstart, $newend ) = ( $newend, $newstart ) if $newstart > $newend;
        if ( $strand =~ /\-/ ) {
            push( @sss,
                    $chr . "\t" 
                  . $start . "\t" 
                  . $newend . "\t"
                  . ( $end - $start ) . "\t"
                  . $newchr . "\t"
                  . $strand );
        }
        else {
            push( @sss,
                    $chr . "\t" 
                  . $start . "\t"
                  . $newstart . "\t"
                  . ( $end - $start ) . "\t"
                  . $newchr . "\t"
                  . $strand );
        }
    }
    close($map_fh);

    my @sorted_lines = sort mapfile_sort @sss;

    my $mapper  = {};
    my $version = $args->{version};
    foreach (@sorted_lines) {    # must be already sorted GFF-like
        my ( $chr, $start, $newstartornewend, $len, $newchr, $strand ) =
          split( /[ \t]+/, $_ );
        push( @{ $mapper->{ $version . '-' . $chr . '-start' } }, $start )
          ;                      # start of region on old
        push( @{ $mapper->{ $version . '-' . $chr . '-end' } }, $start + $len )
          ;                      # end of region on old

        if ( $strand =~ /\-/ ) {
            push(
                @{ $mapper->{ $version . '-' . $chr . '-diff' } },
                $start + $newstartornewend
            );                   # what to subtract to convert
        }
        else {
            push(
                @{ $mapper->{ $version . '-' . $chr . '-diff' } },
                $start - $newstartornewend
            );                   # what to subtract to convert
        }
        push( @{ $mapper->{ $version . '-' . $chr . '-newchr' } }, $newchr );
          ;                      # chromosome (might be new)
        push( @{ $mapper->{ $version . '-' . $chr . '-inverted' } }, ($strand =~ /\-/) ? 1 : 0 );
          ;                      # inverted area
    }
    return $mapper;
}

#***********************************************************
#
#***********************************************************

sub getConvertedCoord {
    my ($args) = @_;

    my $version = $args->{version};

    my ( $chr, $start, $end, $SingleCoordInput ) = splitRange( $args->{coord} );
    return ( $chr, ' failed to convert' ) unless defined $chr;
    ( $start, $end ) = ( $end, $start ) if $start > $end;

    my $is_inversion = 0;

    my ( $newchr, $newchr1, $newstart, $newend, $shifts, $a, $b, $c ) =
      ( '?', '?', '?', '?', "", "", "", "" );

    ( $newchr, $newstart, $is_inversion ) = convertPosition(
        {
            version      => $version,
            chr          => $chr,
            pos          => $start,
            is_inversion => $is_inversion,
            chr_map      => $args->{chr_map},
        }
    );
    print STDERR "  vars state: $newchr:$newstart..\n" if $debug;

    ( $newchr1, $newend, $is_inversion ) = convertPosition(
        {
            version      => $version,
            chr          => $chr,
            pos          => $end,
            is_inversion => $is_inversion,
            chr_map      => $args->{chr_map},
        }
    );
    print STDERR "  vars state: $newchr1:..$newend\n" if $debug;

    return ( "-", ' failed to convert: maps to more than one scaffold' )
      if $newchr ne $newchr1;

    ( $c, $shifts, $a, $b, $is_inversion ) = convertRange(
        {
            version      => $version,
            chr          => $chr,
            start        => $start,
            end          => $end,
            chr_map      => $args->{chr_map},
            is_inversion => $args->{is_inversion}
        }
    );
    ( $newstart, $newend ) = ( $newend, $newstart ) if (($newstart ne '?' && $newend ne '?') && $newstart > $newend);
    print STDERR "  vars state: $newchr:$newstart..$newend\n" if $debug;

    my $converted =
      $newchr . ':' . $newstart . ( ($SingleCoordInput) ? "" : "..$newend" );
    my $note = " ";    # should not be null

    if ( $newstart eq '?' || $newend eq '?' ) {
        $note = "failed: coordinates fully within region of change";
    }
    elsif ( $shifts == 1 ) {
        $note = "spans 1 region of change";
    }
    elsif ( $shifts > 1 ) {
        $note = "spans $shifts regions of change";
    }
    if ( $newchr ne $chr ) {
        $note .= "; " unless $note eq " ";
        $note .= "different scaffold";
    }
    if ($is_inversion) {
        $note .= "; "
          unless $note eq " ";
        $note .= "inversion";
    }
    print STDERR "  conversion result: $args->{coord} => $converted $note\n"
      if $debug;
    return ( $converted, $note );
}

#
##***********************************************************
##
##***********************************************************
#
sub convertRange

  # inside of stretch only!
{
    my ($args)    = @_;
    my $version   = $args->{version};
    my $chr       = $args->{chr};
    my $start     = $args->{start};
    my $end       = $args->{end};
    my $ChrMapper = $args->{chr_map};

    print STDERR "Decoding range $chr:$start..$end (v$version)\n" if $debug;
    ( $start, $end ) = ( $end, $start ) if $start > $end;

    print STDERR "Setting arrays for $chr (v$version)\n" if $debug;
    my @ChrMapperStart     = @{ $ChrMapper->{ $version . '-' . $chr . '-start' } };
    my @ChrMapperEnd       = @{ $ChrMapper->{ $version . '-' . $chr . '-end' } };
    my @ChrMapperDiff      = @{ $ChrMapper->{ $version . '-' . $chr . '-diff' } };
    my @ChrMapperNewchr    = @{ $ChrMapper->{ $version . '-' . $chr . '-newchr' } };
    my @ChrMapperInversion = @{ $ChrMapper->{ $version . '-' . $chr . '-inverted' } };
    my $nstretches = $#ChrMapperStart;

    print "Scanning array of $nstretches stretches\n" if $debug;
    my $nbad = 0;
    for ( my $n = 0 ; $n <= $nstretches ; $n++ ) {
        print STDERR "  testing stretch "
          . $ChrMapperStart[$n] . ".."
          . $ChrMapperEnd[$n] . "\n"
          if $debug;
        last if ( $ChrMapperStart[$n] > $end );
        next if ( $ChrMapperEnd[$n] < $start );

        #if passed, then this stretch of old release somehow overlaps with requested range
        print STDERR "    overlaps (diff " . $ChrMapperDiff[$n] . ") newchr " . $ChrMapperNewchr[$n] . "\n" if $debug;
        if ($ChrMapperInversion[$n]) {
            $args->{is_inversion} = 1;
            $nbad++;
        }
    }
    return ( $chr, $nbad - 1, '?', '?', $args->{is_inversion} );
}

#
##***********************************************************
##
##***********************************************************
#
sub convertPosition {
    my ($args)    = @_;
    my $version   = $args->{version};
    my $chr       = $args->{chr};
    my $pos       = $args->{pos};
    my $ChrMapper = $args->{chr_map};

    print STDERR "Decoding position $chr:$pos (v$version)\n" if $debug;
    print STDERR "Setting arrays for $chr (v$version)\n"     if $debug;
    my @ChrMapperStart      = @{ $ChrMapper->{ $version . '-' . $chr . '-start' } };
    my @ChrMapperEnd        = @{ $ChrMapper->{ $version . '-' . $chr . '-end' } };
    my @ChrMapperDiff       = @{ $ChrMapper->{ $version . '-' . $chr . '-diff' } };
    my @ChrMapperNewchr     = @{ $ChrMapper->{ $version . '-' . $chr . '-newchr' } };
    my @ChrMapperInversion  = @{ $ChrMapper->{ $version . '-' . $chr . '-inverted' } };

    my $nstretches = $#ChrMapperStart;
    for ( my $n = 0 ; $n <= $nstretches ; $n++ ) {
        print STDERR "  testing "
          . $ChrMapperStart[$n] . '-'
          . $ChrMapperEnd[$n] . "\n"
          if $debug;
        next if ( $ChrMapperEnd[$n] < $pos );
        last if ( $ChrMapperStart[$n] > $pos );
        my $newpos = $pos - $ChrMapperDiff[$n];
        $newpos = -$newpos if $newpos < 0;
        $args->{is_inversion} = 1 if ($ChrMapperInversion[$n]);

        my $newchr = $ChrMapperNewchr[$n];
        print STDERR "$newchr:$pos => $chr:$newpos (fully in "
          . $ChrMapperStart[$n] . '-'
          . $ChrMapperEnd[$n] . ")\n"
          if $debug;
        return ( $newchr, $newpos, $args->{is_inversion} );
    }
    return ( $chr, '?', $args->{is_inversion} );
}

#
##***********************************************************
##
##***********************************************************
#
sub splitRange {
    my $coords           = shift;
    my $SingleCoordInput = 0;

    print STDERR "Splitting $coords\n" if $debug;
    $coords =~ s/,//g;
    if ( $coords =~ /^\s*([^ :]+):([0-9,]+)([\.\-]{2})([0-9,]+)\s*$/ ) {   # range
        my ( $chr, $start, $end ) = ( $1, $2, $4 );
        $chr =~ s/^Chr//i;

        #if( $coords=~/[,]/ ) { $start=~s/[,]//; $end=~s/[,]//; }
        if ( $start > $end ) {
            my $a  = $end;
            $end   = $start;
            $start = $a;
        }
        print STDERR "Split as $chr,$start,$end\n" if $debug;
        return ( $chr, $start, $end, $SingleCoordInput );
    }
    elsif ( $coords =~ /^\s*([^ :]+):([0-9,]+)\s*$/ ) {    # single coord
        $SingleCoordInput = 1;
        my ( $chr, $start, $end ) = ( $1, $2, $2 );
        $chr =~ s/^Chr//i;

        #if( $coords=~/[,]/ ) { $start=~s/[,]//; $end=~s/[,]//; }
        if ( $start > $end ) {
            my $a = $end;
            $end   = $start;
            $start = $a;
        }
        print STDERR "  Split as $chr $start $end\n" if $debug;
        return ( $chr, $start, $end, $SingleCoordInput );
    }
    return;
}

#***********************************************************
#
#***********************************************************

sub mapfile_sort {
    my ( $a1, $a2, @ra ) = split( /\t/, $a );
    my ( $b1, $b2, @rb ) = split( /\t/, $b );
    if ( $a1 ne $b1 ) { return ( $a1 cmp $b1 ); }
    return ( $a2 <=> $b2 );
}

