# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of RandProt.
# RandProt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# RandProt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with RandProt.  If not, see <http://www.gnu.org/licenses/>.

package ProtKmer;
our $VERSION = 1.00;

use lib '.';
use ParseFasta;
use strict;

# simple code to extract a kmer and sequence length distribution from protein sequences

# used by multiple scripts
sub normalizeSeq {
    # takes a single sequence and changes two rare amino acids into standard equivalents
    # selenocysteine is scored as cysteine, pyrrolysine is scored as lysine
    $_[0] =~ tr/UO/CK/; # this modifies input by reference, so nothing to return
}

# only used internally by countKmersLengthProteome
sub countKmers {
    # a basic tool for getting kmer counts off a sequence
    # sequence is not modified
    # %kmer2c is a pre-existing hash we want to add counts to
    # $k is size of kmer, has to be supplied
    # hint: k=1 for amino acid BG frequencies
    # also works on standard DNA sequences, but degeneracies won't be treated correctly, so be careful!
    # 2012-05-11 10:46:42 EDT
    # update: added initial Met as another case to ignore
    my ($seq, $kmer2c, $k) = @_;
    foreach my $seq2 (split /^M|[BJXZ]+/, $seq) { # degenerate codes are treated as boundaries for kmers, we only count between them...
	my $iMax = length($seq2)-$k;
	for (my $i = 0; $i <= $iMax; $i++) {
	    $kmer2c->{substr $seq2, $i, $k}++; # get kmer, add to count hash
	}
    }
}

sub countKmersLengthProteome {
    # counts kmers in an entire protein FASTA file (if input is scalar, or combined counts on multiple files if input is an array reference), and also returns maps of protein IDs to sequence lengths (for all files, will not work if files have proteins with the same name!), both for the purpose of generating a random genome
    my ($fis, $k) = @_;
    my %kmer2c; # kmer counts
    my %prot2length; # length of every protein
    # get proteome...
    foreach my $fi (ref $fis ? @$fis : $fis) { # handle plural case
	my ($fhi, $protNext) = ParseFasta::parseInit($fi); # get data
	# browse prots
	while ($protNext) {
	    (my $prot, my $seq, $protNext) = ParseFasta::parseNext($fhi, $protNext);
	    $prot2length{$prot} = length $seq; # store protein length
	    # magic subs
	    normalizeSeq($seq);
	    countKmers($seq, \%kmer2c, $k);
	}
	close $fhi;
    }
    return (\%kmer2c, \%prot2length); # return the desired data
}

sub kMax {
    # a utility that may aid in choosing k
    # uses a very loose upper bound of the max k that a dataset this size allows (assuming that you'd like to see every kmer at least once)
    my ($fis, $verbose) = @_;
    # get number
    print "Input: ".(ref $fis ? join(', ', @$fis) : $fis)."\n" if $verbose;
    my $c = ParseFasta::countLetters($fis);
    print "Number of letters: $c\n" if $verbose;
    my $kMax = int log($c)/log(20); # floor of this number is our estimate
    print "kMax: $kMax\n" if $verbose;
    return $kMax;
}

1;
