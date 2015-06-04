# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of RandProt.
# RandProt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# RandProt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with RandProt.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.00';
use lib '.';
use ProtKmer;
use strict;

# report percentage of k-mers observed in a given proteome

# get inputs from command line
my ($fiFasta, $k) = @ARGV;

die "# $0 $VERSION - Compute percentage of k-mers observed in a proteome
# RandProt ".(sprintf '%0.2f', $ProtKmer::VERSION).", viiia.org/randProt
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <input FASTA> <k>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.
    <k>             The size of k-mers to analyze.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Analysis is printed to STDOUT.

See the online manual for more info.
" unless $k;

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
