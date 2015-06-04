# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of RandProt.
# RandProt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# RandProt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with RandProt.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.01';
use lib '.';
use ProtKmer;
use ProtMarkov;
use strict;

# randomize genomes by preserving length and kmer distributions (but details are complicated)...

# get inputs from command line
my ($fiFasta, $foRandFasta, $k, $n) = @ARGV;

die "# $0 $VERSION - Make random protein sequences from a high-order Markov model
# RandProt ".(sprintf '%0.2f', $ProtMarkov::VERSION).", viiia.org/randProt
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <input FASTA> <output FASTA> <k> <n>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.
    <output FASTA>  Output protein sequence file in FASTA format.
    <k>             The size of k-mers whose distribution will be preserved.  This is 
                    equivalent to specifying a (k-1)-order Markov model.
    <n>             The number of random sequences per real sequence.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

See the online manual for more info.
" unless $n;

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
