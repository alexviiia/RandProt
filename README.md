<img src="logo.png" alt="RandProt" align="right" />

# RandProt

Generate high-order Markov random protein sequences.
Fast implementation that captures subtleties of protein sequences.

Here you can download our code, and learn how to use it.
This page is the software's manual.

# About RandProt

The goal of this software is to efficiently generate random protein sequences from a high-order Markov model.
Roughtly-speaking, a `(k - 1)`-order Markov model preserves the k-mer distribution of the original data, which is how the model is encoded internally.
The code emphasizes simplicity and computational efficiency.
Note, however, that most of the code is specific to protein sequences, it may not be reused for generating other sequences without significant changes.

The formal specification of the random sequences follows.
The input is a collection of protein sequences.
The output will have, for every input sequence, `n` random "replicate" sequences with the same length (with the same IDs plus, for `n > 1`, the suffix `.RANDi` where `i` goes from `0` to `n - 1` with padded zeroes).
The k-mer distribution is obtained from the input data excluding initial methionine (M) residues, when present.
The random sequences always begin with M, followed by a random k-mer drawn from the global distribution, followed by individual amino acids drawn from the Markov model (that is, conditional on the previous `k - 1` amino acids).
The sequence is terminated when the desired length is attained.

The input is always "normalized" so that only the 20 standard amino acids are considered.
To achieve this, selenocysteines (U) are replaced by cysteines (C), and pyrrolysines (O) by lysines (K). 
Additionally, k-mers with ambiguous amino acid codes (B, J, X, Z) are always ignored.

# Installation

You'll need 

- Perl 5
- `gzip`
- Git versioning control software (optional, to clone repository)
- [`bats`](https://github.com/bats-core/bats-core) (optional, for testing)


## Cloning repository

You can download a ZIP copy of this repository on GitHub.
However, cloning the repository has many advantages.
You can clone it with this command:
```bash
git clone https://github.com/alexviiia/RandProt.git
```
If there are any updates in the future, you can simply update your copy of the repository by typing on a terminal, inside the `RandProt` directory:
```bash
git pull
```

## Running scripts from arbitrary directories

Each script can be run easily from the directory that contains it (and its `*.pm` perl module dependencies).
For example, if you cloned this repository onto the local directory `RandProt/`, on a terminal you can run one of the help messages by typing:
```bash
cd RandProt/ 
perl -w randProt.pl # <ARGS>... 
```

To run the scripts from other directories, you have to specify the directory containing the `*.pm` Perl modules using the `-I` option for `perl`, like so:
```bash
cd ..
perl -IRandProt/ -w RandProt/randProt.pl # <ARGS>...
```
Note that the home directory shortcut `~` doesn't work with `-I`, but you can use `$HOME` instead.
So if the code is in `~/RandProt/`, then this will work:
```bash
perl -I$HOME/RandProt/ -w ~/RandProt/randProt.pl # <ARGS>...
```

## Running tests

You can run automatic tests to make sure the code is running as expected on your system.
The specific errors may help me troubleshooting.

After installing `bats` (available on most standard Linux repositories), run this command:
```bash
bats tests.bats
```
If successful, the output will look like this:
```
 ✓ kMax.pl help
 ✓ kMax.pl main
 ✓ kCov.pl help
 ✓ kCov.pl main 5
 ✓ kCov.pl main 4
 ✓ kCov.pl main 3
 ✓ randProt.pl help
 ✓ randProt.pl main n=1
 ✓ randProt.pl main n=2

9 tests, 0 failures
```

# Running the scripts

## Recommendations

Besides the input sequence file, the main analysis requires positive integer parameters `k` and `n`.
The value of `n` provides the number of random samples per real sequence.
In theory, the larger the better, but the complexity of downsteam analyses may set upper limits on `n`.
I've personally used `n = 100` and `n = 20` for small proteomes, and `n = 1` for a much larger protein database.

The value of `k` controls how much the random sequences will look like real sequences, with small values of `k` yielding less realistic random proteins.
However, a value of `k` that is too large can be bad, as it may recapitulate properties of your real proteins that you do not wish to have in your random sequences (for example, in my unpublished research I found that some protein domains can be well predicted by the presence of 5-mers that represent highly-conserved motifs, so even `k = 5` preserves some biologically-relevant signal that perhaps should be excluded from random sequences).
In my work I favored `k = 3`.

## Sample files

This repository contains the following sample files:

Sample input:

- [sample/PLAF7.fa](sample/PLAF7.fa): the *Plasmodium falciparum* proteome from [PlasmoDB](http://plasmodb.org/plasmo/) version 11.0, with pseudogenes removed, containing 5441 proteins.

Sample outputs:

- [sample/PLAF7.k3.n1.fa](sample/PLAF7.k3.n1.fa): Output from `randProt.pl` with `k = 3` and `n = 1` (see below).
  Note IDs are the same as input (no extra suffixes).
- [sample/PLAF7.k3.n2.fa](sample/PLAF7.k3.n2.fa): Output from `randProt.pl` with `k = 3` and `n = 2` (see below).
  Note IDs have extra suffixes (counting replicates).


## Synopsis of scripts

All scripts give detailed usage instructions when executed without arguments.
The following commands can be run on the sample file placed in the same directory as the code and called from that location.
The input file may be compressed with `gzip` and may be specified with or without the GZ extension.
The output file is compressed with `gzip`, whether the outputs indicate a GZ extension or not.
So the commands as they are below produce and will work entirely with compressed files, without a single GZ extension indicated.

Get a quick estimate of the maximum k to consider 
```bash
perl -w kMax.pl PLAF7.fa 
```

Try a few values of k until a decent coverage is achieved 
```bash
perl -w kCov.pl PLAF7.fa 5 
perl -w kCov.pl PLAF7.fa 4 
perl -w kCov.pl PLAF7.fa 3 
```

Generate your desired random proteome! 
```bash
perl -w randProt.pl PLAF7.fa PLAF7.k3.n2.fa 3 2
```

## Using `kMax.pl`: Compute a weak upper bound on k for k-mer analysis

This is the help message you get by running the script without arguments:
```bash
perl -w kMax.pl
```
```
# kMax.pl: Compute a weak upper bound on k for k-mer analysis
# RandProt 1.04 - https://github.com/alexviiia/RandProt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w kMax.pl <input FASTA>

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
```

This is the output produced on the sample file:
```bash
perl -w kMax.pl sample/PLAF7.fa 
```
```
Input: sample/PLAF7.fa
Number of letters: 4139904
kMax: 5
```

## Using `kCov.pl`: Compute percentage of k-mers observed in a proteome

This is the help message you get by running the script without arguments:
```bash
perl -w kCov.pl 
```
```
# kCov.pl: Compute percentage of k-mers observed in a proteome
# RandProt 1.04 - https://github.com/alexviiia/RandProt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w kCov.pl <input FASTA> <k>

The required inputs are
    <input FASTA>   Input protein sequence file in FASTA format.
    <k>             The size of k-mers to analyze.

Input file may be compressed with gzip, and may be specified with or without the .gz 
extension.  Analysis is printed to STDOUT.

This script obtains the number of unique k-mers in the input file, and compares it to the 
theoretical number of possible k-mers, 20^k.  A k with high coverage should be chosen, but 
the largest such k is merely an upper bound, and you may want to choose a smaller k using 
other criteria.
```

For the sample file, the outputs for `k = 5, 4, 3` are shown below, which shows that `k = 3` should be the maximum value of `k` to consider for this data.
```bash
perl -w kCov.pl sample/PLAF7.fa 5
```
```
k: 5
Number of theoretical amino acid k-mers: 3200000
Scanning proteome...
Number of observed amino acid k-mers: 1226974
Proportion of k-mers observed: 38.34 %
```
```bash
perl -w kCov.pl sample/PLAF7.fa 4
```
```
k: 4
Number of theoretical amino acid k-mers: 160000
Scanning proteome...
Number of observed amino acid k-mers: 148284
Proportion of k-mers observed: 92.68 %
```
```bash
perl -w kCov.pl sample/PLAF7.fa 3
```
```
k: 3
Number of theoretical amino acid k-mers: 8000
Scanning proteome...
Number of observed amino acid k-mers: 7999
Proportion of k-mers observed: 99.99 %
```

This script `kCov.pl` is the more refined way of choosing an upper bound for `k`.
However, `kMax.pl` will be much faster on very large inputs, so it should be run it first to get a rough idea of the range of `k` to test.

## Using `randProt.pl`: Make random protein sequences from a high-order Markov model

This is the help message you get by running the script without arguments:
```bash
perl -w randProt.pl
```
```
# randProt.pl: Make random protein sequences from a high-order Markov model
# RandProt 1.04 - https://github.com/alexviiia/RandProt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w randProt.pl <input FASTA> <output FASTA> <k> <n>

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
```

For the sample file and `k = 3`, I provide here two random outputs, with `n = 1` and `n = 2` replicates per protein, each of which takes a few seconds to run.
```bash
perl -w randProt.pl sample/PLAF7.fa sample/PLAF7.k3.n1.fa 3 1
perl -w randProt.pl sample/PLAF7.fa sample/PLAF7.k3.n2.fa 3 2
```
The output message in all cases is:
```
Scanning proteome...
Generating random proteome...
```
The first two random sequences for the `n = 1` output have this format:
```bash
# assuming default compressed output (repo has uncompressed copies)
zcat sample/PLAF7.k3.n1.fa.gz |head -n 4
# # for uncompressed files, use:
# cat sample/PLAF7.k3.n1.fa |head -n 4
```
```
>PF3D7_0828500
MNCITKNNMVDFFFSTILEKIHNNNNFNYNNKIIDNYILFWNYYMFSFFLKRVMKNDDRKVNEKYAFIKSFSKSDLERNKDVSNDKKCIEGEKKNICYSVPQCTGGDKEKMNINELSITQEEEELTNDYKNIVMKSLINKTFLLCKDNDKSITNQHNNNYHTDADNNLCIYSMLVEYYNSISTLTYIQLVFEHDKQNINIPRQMNNIKLIRQFLGRDLSADHMEKYCKNKKYIGLKVEVYFSHNNHTTDVLFFLDKVEIKHYNKENIIKMKKIVIIWIHQICEKEILLNLYIGVKFNNIKINDYVNGYNDGDYLHLSNYIPDKNKQVAKALNSDINIYIRTI
>PF3D7_1334300
MQLRKNKQSFTCFEVIQDIEIQLGYDEGELVYIKPISSIKNKIRRIRNKYLSFKNDLQKEINEKGRDDDTNVYEKYRTISNIKNINMENGKKNLISICKKSFNVARTKWFERLFLTNLKKNYNKKTYINEQKYQNYKVLDYTQELLSKDQIIFGGGKDGGKKNRKELIFLQPQIRKGYQDKNKKYNGDIKKDTHNIADFIKIVKCSPIKTRGSILHQQKYNIIYCKEEKESDLTISSPFQIQKIQSVKLVQCGYFPFFYRKKELEKNYLSLYKGMNYYKKEILCKSFKSMSSPSDSSKNENYMYRRKSTENTILKSEIYSIKCPYSFTNLFLRKKEEIKQTSAKECDEHVHNECVQQQQPYSDYEIQNVSYDDNDMDFLNKQKYFYMNIKVNDNIHGFIKKVCFFFIQEGANDPNIDENKINMKDKVNYDNNINNYSSDQVTEFSVHSINSNSSINLGR
```
For the `n = 2` output, each ID has a suffix `.RANDi` that distinguishes different replicates `i`:
```basg
zcat sample/PLAF7.k3.n2.fa.gz |head -n 4
```
```
>PF3D7_1239400.RAND0
MKRNNPLSSSLLKIWNWDAKDENACANILINKDKVLNISVHFSTFFYLASLNHIYCLLKRKSISNYDDNIDLNKGNTEQFFFKYRHIIEKYNNNGLNTHVRSLSCLIKYIEAYTQYNLLFHFKSYILLCTNNPNIDDDKVNNNNNEEVVKSSSTDETEKDDDIKRIRKEGVLYICNYNN
>PF3D7_1239400.RAND1
MNDKPKDENNNNYMNPHFHITPASMKNIIDHKDKEEKMKFALNSSENDKNKVTVPIIYQIQSDNISSVWPYMNNNNNKKNDSSLKKDTNNNNEKSKKRNNMDSLYNIMNEELARQILYFYKAQIRASNEYIQENIFIALGMCEQETYHMHMKRRHFITLNYDKRGQKEDIRKRKNYLKK
```
Note that the order in the output is random (and, of course, each sequence is random, so outputs are never the same).

# Original Perl module descriptions

Since my source code is not self-documented, here's a brief description of what each module does.

| File          | Description                                                                                                                                                                                                                                                                                                                               |
|---------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FileGz.pm     | Handles normal and compressed files transparently.                                                                                                                                                                                                                                                                                        |
| ParseFasta.pm | Light-weight tools to handle FASTA files.                                                                                                                                                                                                                                                                                                 |
| ProtKmer.pm   | Functions to normalize protein sequences, and to count k-mers efficiently.                                                                                                                                                                                                                                                                |
| ProtMarkov.pm | Functions that generate random sequences from the k-mer data (and protein lengths). The hardest part is drawing the initial k-mer of a sequence, which entails drawing from a huge categorical distribution with non-uniform probabilies. I encoded a binary-search that draws in O(log(m)), where m =~ 20^k is the number of categories. |


# Software license

This code is released under the GNU GPLv3 (GNU General Public License version 3).
See [LICENSE](LICENSE).

# Compatibility

This code should work for any Perl version >= 5 (tested 5.18, 5.20, 5.22, 5.28).
The code should work on any Linux and MacOS, but let me know otherwise.

## Want a Windows version?

The code will not work on Windows machines because it uses the gzip executable and it also uses Unix pipes. However, if you install the PerlIO::gzip package, or forgo working with compressed files, the code could run (with some adjustments, contact me for more info).

# Citation

2015-11-17. 
Alejandro Ochoa, John D Storey, Manuel Llinás, and Mona Singh. 
Beyond the E-value: stratified statistics for protein domain prediction. 
PLoS Comput Biol. 11 e1004509. 
[Article](http://dx.doi.org/10.1371/journal.pcbi.1004509), 
[arXiv](http://arxiv.org/abs/1409.6384) 2014-09-23.
