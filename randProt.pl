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

# randomize genomes by preserving length and kmer distributions (but details are complicated)...

# get inputs from command line
my ($fiFasta, $foRandFasta, $k, $n) = @ARGV;

unless ($n) {
    print "# $FindBin::Script: Make random protein sequences from a high-order Markov model
# " . ProtMarkov::version_string() . "
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $FindBin::Script <input FASTA> <output FASTA> <k> <n>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.
    <output FASTA>  Output protein sequence file in FASTA format.
    <k>             The size of k-mers whose distribution will be preserved.  This is 
                    equivalent to specifying a (k-1)-order Markov model.
    <n>             The number of random sequences per real sequence.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

Output is the desired database of random protein sequences, with n replicates per sequence.
Each output sequence has the same length as the corresponding input sequence.  First the 
(k-1)-order Markov model is learned from the input data, then each sequence is drawn from 
this model.  Output has IDs in a random order.
";
    exit 0;
}

# constant
my $comp = 'gzip';

# quick validations
$k = int $k; # force both of these to be integers
$n = int $n;
die "Error: k must be at least 1 (input was k=$k)!\n" unless $k > 0;
die "Error: n must be at least 1 (input was n=$n)!\n" unless $n > 0;

print "Scanning proteome...\n";
my ($kmer2c, $prot2length) = ProtKmer::countKmersLengthProteome($fiFasta, $k);

print "Generating random proteome...\n";
ProtMarkov::randomizeGenomeMarkov($kmer2c, $prot2length, $n, $foRandFasta, $comp);
