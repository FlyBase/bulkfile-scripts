#!/usr/bin/env perl

use strict;
use warnings;

=pod

=head1 NAME

problem_case_filter.pl

=head1 SYNOPSIS

Removes features that can cause problems with certain tools.

cat dmel-all-r5.56.gff | problem_case_filter.pl > problem_free_dmel.gff

=head1 DESCRIPTION

This script was initially written in response to requests from RNA-Seq 
data analysis partners for their use of our data with common RNA-Seq
analysis tools.  These tools had various difficulties with certain
features, requiring that they be removed before further analysis.

When a feature is filtered out, it and all its children are removed
from the file.

See L<Filtered features> for a current list of features that are
removed.

=head2 Filtered features

The current list of problematic features include the following:

=over 5

=item * Genes with these SO terms

=over 5

=item * gene_with_dicistronic_mRNA (SO:0000722)

=item * gene_with_trans_spliced_transcript (SO:0000459)

=back

=back

=head2 Disclaimers & Assumptions

=over 5

=item 1.
Only tested with FlyBase GFF files.  Use with other GFF files may result in incorrect results.

=item 2.
Assumes that parent features come before their children.

=back


=head1 AUTHOR

Josh Goodman, FlyBase

=head1 TODO

=over 5

=item 1. Deal with children coming before parents by implementing two pass processing or pre-sorting.

=back

=head1 LICENSE

 Copyright (c) 2014, Indiana University 
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

my %filter = ();

while (<>) {
    chomp;

    #If we encounter a problem case do this.
    if (/(SO:0000459|SO:0000722)/o) {
        #Get the ID and store it in the %filter hash.
        my ($id) = $_ =~ /\tID=(\S+?);/o;
        print STDERR "Line $. does not have an ID as expected.\n" unless defined $id;
        $filter{$id} = undef if defined $id;
    }
    #Check child features.
    elsif (my ($parent_ids) = $_ =~ /(?:Derives_from|Parent)=(\S+?);/o) {
        #Get the child's ID.
        my ($id) = $_ =~ /\tID=(\S+?);/o;

        #If this child's parent is a problem case, ignore
        #it and store the child ID in the filter hash.
        if (grep {exists $filter{$_} } split(',',$parent_ids)) {
            $filter{$id} = undef if defined $id; #Some children don't have IDs.
        }
        #Otherwise, print it.
        else {
            print "$_\n";
        }

    }
    #Print the GFF line if all else fails.
    else {
        print "$_\n";
    }
}
