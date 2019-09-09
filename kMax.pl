# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of RandProt.
# RandProt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# RandProt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with RandProt.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.02';
use lib '.';
use ProtKmer;
use ProtMarkov;
use strict;

# output kMax estimate, verbose by default

# get inputs from command line
my ($fiFasta) = @ARGV;

unless ($fiFasta) {
    print "# $0  $VERSION - Compute a weak upper bound on k for k-mer analysis
# RandProt ".(sprintf '%0.2f', $ProtMarkov::VERSION)." - https://github.com/alexviiia/RandProt
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <input FASTA>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Analysis is printed to STDOUT.

See the online manual for more info.
    ";
    exit 0;
}

# constant
my $verbose = 1; # will output analysis to STDOUT

# this does the magic
# also returns kMax (here ignored)
ProtKmer::kMax($fiFasta, $verbose);

