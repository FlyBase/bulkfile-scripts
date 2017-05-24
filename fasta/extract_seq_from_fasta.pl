#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use v5.10;

use Getopt::Long;
use Pod::Usage;
use FindBin;

=pod

=head1 NAME

extract_seq_from_fasta.pl - Extract select isoforms from the FlyBase FASTA files.

=head1 SYNOPSIS

 extract_seq_from_fasta.pl [options]

 Options:
 --fasta <file>  - The file to read R5 coordinates from
 --longest       - Extract the longest isoforms
 --unique        - Extract the unique isoforms
 --byid <file>   - Extract sequences by ID
 --help          - Print a help message
 --debug         - Print debug information to STDERR
 --man           - Show a man page help doc

 e.g.

 zcat dmel-all-translation.fasta.gz | ./extract_seq_from_fasta.pl --longest > extracted_isoforms.fasta

 ./extract_seq_from_fasta.pl --unique --fasta dmel-all-translation.fasta.gz 

 ./extract_seq_from_fasta.pl --byid myids.txt --fasta dmel-all-translation.fasta 

 see L<Options> for full details.

=head1 DESCRIPTION

  This script extracts various types of sequences from FlyBase FASTA files.

  Currently supported types are:

  * longest
  * unique
  * by ID 

  For longest and unique, the FASTA file used as input must have an FBgn ID in the ID field or the parent field.

  >mygene ID=FBgn1234567;
  or
  >myprotein parent=FBtr1234567,FBgn1234567;

=head2 Options

=over 5

=item  --fasta <file>

The FASTA file to read from.  File can be a plain text or gzip compressed FASTA file.

=item  --longest

  Extract the longest sequence of a gene product (transcript, polypeptide, etc.).

=item  --unique

  Extract the unique sequence of a gene product (transcript, polypeptide, etc.).

=item  --byid <file>

  Extract sequences by their ID contained in <file>.  <file> must have one ID per line.

=item --help 

Print a help page.

=item --debug

Print debug information to STDERR during conversion.

=item --man

Show help as a man page.

=back

=head2 Description of types.

=head3 Longest isoform

Longest transcript or polypeptide isoform of a gene in the given input file.

=head3 Unique sequences

Filters out unique sequences of a gene in the given input file.
Uniqueness is determined by the MD5 checksum of the sequence itself.

=head3 By ID

Extracts FASTA entries by ID from the input file.

=head1 AUTHOR

=over 5

=item Josh Goodman, FlyBase

=back


=head1 LICENSE

 Copyright (c) 2017, Indiana University & FlyBase 
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

my $debug   = 0;
my $help    = 0;
my $man     = 0;
my $fasta   = '-';
my $longest = 0;
my $byid;
my $unique  = 0;

my $getopt = GetOptions(
  'help|?'    => \$help,
  'debug'     => \$debug,
  'man'       => \$man,
  'fasta=s'   => \$fasta,
  'longest'   => \$longest,
  'unique'    => \$unique,
  'byid=s'    => \$byid,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;


my $fasta_fh;
if ( $fasta eq '-' ) {
  $fasta_fh = \*STDIN;
}
else {
  # IF FASTA file ends with .gz extension, pipe it through gunzip first.
  if ($fasta =~ /\.gz$/) {
    open( $fasta_fh, '-|', "gzip -dc $fasta" );
  }
  else {
    open( $fasta_fh, '<', $fasta );
  }
}

pod2usage( -exitstatus => 1, -verbose => 1, -msg => "Must specify at least one type to extract" ) unless ($longest || $unique || $byid);

my $seqs = {};

my $ids;
if ($byid) {
  $ids = read_ids($byid);
}

local $/ = '>';

while (my $buf = <$fasta_fh>) {
  chomp $buf;
  next if ($buf =~ m/^\s*?$/);

  if ($longest) {
    my $record = create_record({ entry=>$buf });
    my $fbgn = $record->{fbgn};

    if (!defined $seqs->{$fbgn}) {
      $seqs->{$fbgn} = $record;
    }
    elsif ($record->{length} >= $seqs->{$fbgn}{length}) {
      $seqs->{$fbgn} = $record;
    }
    else {
      say STDERR "Skipping $fbgn..." if $debug;
    }
  }
  elsif ($unique) {
    my $record = create_record({ entry=>$buf });
    my $fbgn = $record->{fbgn};
    my $md5 = $record->{md5};
    my $key = $fbgn . '_' . $md5;

    if (!defined $seqs->{$key}) {
      $seqs->{$key} = $record;
    }
    elsif ($record->{md5} ne $seqs->{$key}{md5}) {
      $seqs->{$key} = $record;
    }
    else {
      say STDERR "Skipping $fbgn..." if $debug;
    }
  }
  elsif ($byid) {
    my $record = create_record({ entry=>$buf });
    my $id = $record->{id};
    if (exists $ids->{$id}) {
      $seqs->{$id} = $record;
    }
  }
}

for my $key (keys %{$seqs}) {
  print $seqs->{$key}{fasta};
}


sub read_ids {
  my $file = shift;
  my $ids;
  open (my $fh, '<', $file);
  while (<$fh>) {
    chomp;
    $ids->{$_} = undef;
  }
  close($fh);
  return $ids;
}

sub create_record {
  my ($args) = @_;
  my $entry = $args->{entry};

  my $id     = get_header_value('ID',$entry);
  my $fbgn   = get_fbgn($entry);
  my $len    = get_header_value('length',$entry);
  my $md5    = get_header_value('MD5', $entry);

  return {
    id     => $id,
    fbgn   => $fbgn,
    length => $len,
    md5    => $md5,
    fasta  => '>' . $entry,
  };
}

sub get_fbgn {
  my $header = shift;
  my $id      = get_header_value('ID',$header);
  my $parents = get_header_value('parent',$header);
  my @fbgn = grep { /^FBgn\d+$/ } split(/,/,$parents) if ($parents);

  return $id if $id =~ /^FBgn\d+$/;
  return $fbgn[0] if (scalar @fbgn > 0 && $fbgn[0] =~ /^FBgn\d+$/);
  return undef;
}

sub get_header_value {
  my ($key,$header) = @_;
  $header =~ m/$key=(.*?);/;
  return $1;
}


