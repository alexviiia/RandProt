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

# output kMax estimate, verbose by default

# get inputs from command line
my ($fiFasta) = @ARGV;

unless ($fiFasta) {
    print "# $FindBin::Script: Compute a weak upper bound on k for k-mer analysis
# " . ProtMarkov::version_string() . "
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $FindBin::Script <input FASTA>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Analysis is printed to STDOUT.

For a given k, there are l = 20^k unique k-mers overall.  For a small k and a large 
proteome, the number of amino acids l is only slightly more than the number of k-mers.
Therefore, the maximum number of k-mers the proteome could possibly contain is

                       k <= k_max = floor( log(l) / log(20) ).

However, usually k must be much smaller than this upper bound. For the input FASTA file, 
this script ouputs l and k_max.  This is a very quick rough analysis.  For a more precise 
assessment on k-mer coverage, use kCov.pl.
";
    exit 0;
}

# constant
my $verbose = 1; # will output analysis to STDOUT

# this does the magic
# also returns kMax (here ignored)
ProtKmer::kMax($fiFasta, $verbose);

