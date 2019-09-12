# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of RandProt.
# RandProt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# RandProt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with RandProt.  If not, see <http://www.gnu.org/licenses/>.

# core Perl modules
use FindBin ();
# local modules
use lib '.';
use ProtKmer;
use ProtMarkov;
use strict;

# report percentage of k-mers observed in a given proteome

# get inputs from command line
my ($fiFasta, $k) = @ARGV;

unless ($k) {
    print "# $FindBin::Script: Compute percentage of k-mers observed in a proteome
# " . ProtMarkov::version_string() . "
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $FindBin::Script <input FASTA> <k>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.
    <k>             The size of k-mers to analyze.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Analysis is printed to STDOUT.

This script obtains the number of unique k-mers in the input file, and compares it to the 
theoretical number of possible k-mers, 20^k.  A k with high coverage should be chosen, but 
the largest such k is merely an upper bound, and you may want to choose a smaller k using 
other criteria.
";
    exit 0;
}

# quick validations
$k = int $k; # force to be integers
die "Error: k must be at least 1 (input was k=$k)!\n" unless $k > 0;
print "k: $k\n";

my $kmersTot = 20**$k; # number of unique k-mers from the standard amino acid alphabet
print "Number of theoretical amino acid k-mers: $kmersTot\n";

print "Scanning proteome...\n";
my ($kmer2c) = ProtKmer::countKmersLengthProteome($fiFasta, $k);
my $kmersObs = keys %$kmer2c; # just count the keys, meh
my $kmersPc = sprintf '%0.2f', $kmersObs/$kmersTot*100;

print "Number of observed amino acid k-mers: $kmersObs\n";
print "Proportion of k-mers observed: $kmersPc %\n";
