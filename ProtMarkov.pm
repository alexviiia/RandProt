# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of RandProt.
# RandProt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# RandProt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with RandProt.  If not, see <http://www.gnu.org/licenses/>.

package ProtMarkov;
our $VERSION = 1.01;

# 2015-03-02 v1.01
# - for $n==1, now suffix .RAND\d+ isn't added (because it's trivial, makes sense to keep it simple)
# 2015-06-04 (still v1.01, hadn't been posted yet)
# - debugged $n==1 case (just now tested code, now it definitely works)

use lib '.';
use FileGz;
use strict;

# simple and efficient code to generate random sequences given a precomputed kmer and sequence length distribution.

# constants
# process only standard amino acids, others are a minority and treating them specially doesn't really change distributions overall, but it does make everything a lot slower
our @aas = sort qw(A C D E F G H I K L M N P Q R S T V W Y); # make sure they're sorted for binary searches

sub randomizeGenomeMarkov {
    # this sub takes 
    # - a kmer distribution (precomputed, usually specific to organism of interest)
    # - a hash of IDs mapping to their original protein length
    # - number of randomizations
    # - an output file path, and whether we want it compressed or not
    # and prints a randomized genome with the exact same length and kmer distributions.  It does so via the implied k-1 Markov chain.  Each protein is randomized $n times.
    # the new proteins use the old protein names, suffixed by .RAND$i to list the randomizations
    my ($kmer2p, $id2length, $n, $fo, $comp) = @_;
    # this hash memoizes cumulative distributions per k-1-mer, which are a slow step in our benchmarks
    my %km1mer2cum;
    # make a structure amenable to binary searches, for much faster sampling!
    my ($pCum, $kmersCum) = makeCumFromHash($kmer2p);
    # make format string
    my $spFmt = makeFmtStr($n);
    # to keep memory usage low, we'll print the new sequences directly to the output files, instead of keeping a big hash as we used to
    # obviously this scales optimally with $n (in memory)
    my $fho = FileGz::getOutFh($fo, $comp);
    # navigate ID and length pairs
    while (my ($idOld, $length) = each %$id2length) {
	# number of times we'll do this length
	for (my $i = 0; $i < $n; $i++) {
	    my $idNew = makeRandId($idOld, $spFmt, $i); # name the new random protein sequence
	    my $seqNew = makeRandSeq($kmer2p, $length, $pCum, $kmersCum, \%km1mer2cum); # this sub produces the k-1 Markov random sequence of desired length
	    print $fho ">$idNew\n$seqNew\n"; # print new sequence to output
	}
    }
    close $fho; # done, close output file and return
}

### SEQUENCE CONSTRUCTION

sub makeRandSeq {
    # given a kmer distribution and a target length, this function generates a random sequence using the k-1 Markov model implied by the kmer distribution
    my ($kmer2p, $length, $pCum, $kmersCum, $km1mer2cum) = @_;
    my $seq = 'M'; # always start sequence with Met, only part modelled separately from the Kmer distributions...
    # choosing the first kmer uses this general function, drawing from the overal distribution of kmers
    my $kmer = sampleCumBinarySearch($pCum, $kmersCum);
    $seq .= $kmer; # add to growing sequence, which is now k+1 long!
    # start producing new letters, until sequence is filled
    while ($length > length $seq) {
	# this function will choose the next amino acid (and update the current kmer) using the right k-1 Markov transition probabilities
	($kmer, my $aa) = drawNextAa($kmer2p, $kmer, $km1mer2cum);
	$seq .= $aa; # add new amino acid to growing sequence
    }
    return $seq; # done, return random sequence
}

sub drawNextAa {
    # given a kmer distribution, and a current kmer (turned k-1-mer here, edited by reference!!!), this returns the next kmer and next amino acid with the appropriate distribution
    # assumes all standard kmers have defined probabilities (as they should, a markov model with gaps can lead to forced paths and a higher effective K)
    # NOTE: this new version instead uses a binary search strategy, hoping it will be much faster!
    my ($kmer2p, $kmer, $km1mer2cum) = @_;
    substr $kmer, 0, 1, ''; # this removes the first character, we get a k-1 string
    my $cum = $km1mer2cum->{$kmer}; # get cumulative from memoized hash, but it may not be precomputed yet...
    unless (defined $cum) {
	# we have to calculate this thing the first time we encounter each k-1-mer
	my $pNorm = 0; # we need to get the total probability of all the output possibilities, essentially normalizing them to a conditional
	my @cum; # the cumulative array used in the binary search
	foreach (@aas) { # loop through following amino acid possibilities
	    my $p = $kmer2p->{$kmer.$_}; # get probability, but it'll be undefined if kmer was not observed
	    $pNorm += $p if $p; # add probability, if non-zero
	    push @cum, $pNorm; # put in cumulative array
	}
	$cum = \@cum; # use this array outside this loop
	$km1mer2cum->{$kmer} = $cum; # remember for next time!
    }
    # use handy function to sample using this distribution, and binary search!!!
    # note amino acids are in the same order as cumulative was built!
    my $aa = sampleCumBinarySearch($cum, \@aas);
    return ($kmer.$aa, $aa); # this is the data we want back!
}

### FUNCTIONS FOR SEQ ID FORMAT MAKING/PARSING

sub makeFmtStr {
    # makes format string for random numbers (with zero padding for niceness)
    my ($n) = @_;
    my $spFmt; # leave undefined for n=1 case...
    if ($n > 1) {
	# for 0-padding in protein names, but count will be zero-based too so even $n = 100 should only use two digits
	my $numSigDigits = 1 + int( log($n-1)/log(10) );
	$spFmt = '%0'.$numSigDigits.'d'; # the sprintf format string, calculate it once for the rest of the subroutine
    }
    # done, return
    return $spFmt;
}

sub makeRandId {
    # this makes the randomized ID corresponding to a given protein and rand number (and the format that decides how many zeroes to use for padding, see makeFmtStr() to get that string)
    # chose to make a function so the format is uniform across applications
    my ($id, $spFmt, $i) = @_;
    # if n==1, ID stays the same, otherwise add this numbering...
    # we know whether $spFmt was defined or not
    if (defined $spFmt) {
	return $id . '.RAND' . sprintf $spFmt, $i;
    }
    else {
	return $id; # return same ID as input
    }
}

sub getRealProtFromRandId {
    # given a randomized id, returns the corresponding real protein.  The idea is to supply a function that we can use without having to think about what the formats actually are
    my $protRand = shift;
    if ($protRand =~ /^(.*)\.RAND\d+$/) {
	return $1; # correctly-formatted suffix was found and removed
    }
    else { return $protRand; } # no suffix, return the same as input (assume it wasn't randomized)  Is this the correct behavior?
}

sub getRandNumFromRandId {
    # given a randomized id, returns the corresponding randomization number.  The idea is to supply a function that we can use without having to think about what the formats actually are
    my $protRand = shift;
    if ($protRand =~ /^.*\.RAND(\d+)$/) {
	return $1; # correctly-formatted suffix was found and randomization number extracted
    }
    else { return undef; } # correct format not found, answer is obviously undefined
}

### FUNCTIONS FOR RANDOM KMERS DRAWS

sub makeCumFromHash {
    # this function takes a distribution defined over the keys of a hash, and returns a cumulative distribution (with keys in a random order, but an order that is defined and returned), which we can use to sample from quickly using binary searches!
    my $e2c = shift; # pass hash as reference, elements to weights
    my $cum = 0; # store cumulative as we go
    my @cum; # the cumulative array we want
    my @es; # elements parallel to cumulative array
    while (my ($e, $c) = each %$e2c) {
	$cum += $c; # add to cumulative
	push @cum, $cum; # add to array
	push @es, $e; # add to parallel array
    }
    return (\@cum, \@es);
}

sub sampleCumBinarySearch {
    # assumes we have a very large cumulative distribution in an array, and performs a binary search to "draw" a random element from this distribution
    
    # the following describes Binary Search in general (with least upper bound behavior)
    # search array of numbers @$a for $x.
    # return index of match, or of the least upper bound otherwise
    # IMPORTANT: assumes list @$a is already sorted ascending, it doesn't work otherwise!!!
    # we don't sort here cause it would run the sort at every single call, but we can probably do it much less often than that...
    # based on code stolen from
    # http://staff.washington.edu/jon/dsa-perl/bsearch
    # actually this emulates C version with same name in Dpuc.c, look at that too for more notes...
    my ($cCum, $es) = @_;
    # initialize these appropriately
    my ($l, $u) = (0, @$cCum - 1); # lower, upper end of search interval
    my $x = rand $cCum->[-1]; # the random draw, in relevant range
    my $i; # index of probe
    while ($l < $u) {
	$i = int(($u+$l)/2); # choose the index right at the middle of our range.
	# since we want the least upper bound, upper has to include the probe i always when it needs to be updated!
	if ($cCum->[$i] < $x) { $l = $i+1; } # update bounds accordingly (next line too), get number in array only once!
	else { $u = $i; }
    }
    return $es->[$u]; # returns the element we wanted!!!
}

1;
