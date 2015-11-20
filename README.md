RandProt
===

Generate high-order Markov random protein sequences.
Fast implementation that captures subtleties of protein sequences.

About RandProt
===

The goal of this software is to efficiently generate random protein sequences from a high-order Markov model. Roughtly-speaking, a (k-1)-order Markov model preserves the k-mer distribution of the original data, which is how the model is encoded internally. The code emphasizes simplicity and computational efficiency. Note, however, that most of the code is specific to protein sequences, it may not be reused for generating other sequences without significant changes.

The formal specification of the random sequences follows. The input is a collection of protein sequences. The output will have, for every input sequence, n random sequences with the same length (with the same IDs plus, for n > 1, the suffix .RANDi where i goes from 0 to n−1 with padded zeroes). The k-mer distribution is obtained from the input data excluding initial methionine (M) residues, when present. The random sequences always begin with M, followed by a random k-mer drawn from the global distribution, followed by individual amino acids drawn from the Markov model (that is, conditional on the previous k−1 amino acids). The sequence is terminated when the desired length is attained.

There's more details of lesser importance. The input is always "normalized" so that only the 20 standard amino acids are considered. To achieve this, selenocysteines (U) are replaced by cysteines (C), and pyrrolysines (O) by lysines (K). Additionally, k-mers with ambiguous amino acid codes (B,J,X,Z) are always ignored.

Installation
===

You'll need Perl 5 (without additional Perl packages) and gzip installed.

Synopsis of scripts
===

Get a quick estimate of the maximum k to consider 
```
perl -w kMax.pl PLAF7.fa 
```

Try a few values of k until a decent coverage is achieved 
```
perl -w kCov.pl PLAF7.fa 5 
perl -w kCov.pl PLAF7.fa 4 
perl -w kCov.pl PLAF7.fa 3 
```

Generate your desired random proteome! 
```
perl -w randProt.pl PLAF7.fa PLAF7.k3.n100.fa 3 100
```

More details
===

All scripts give detailed usage instructions when executed without arguments.  I have more detailed recommendations, usage examples (code snippets and test sample inputs and outputs) at my personal website, viiia.org, in [English](http://viiia.org/randProt/?l=en-us) and [Spanish](http://viiia.org/randProt/).


Citation
===

2015-11-17. Alejandro Ochoa, John D Storey, Manuel Llinás, and Mona Singh. Beyond the E-value: stratified statistics for protein domain prediction. PLoS Comput Biol. 11 e1004509. [Article](http://dx.doi.org/10.1371/journal.pcbi.1004509), [arXiv](http://arxiv.org/abs/1409.6384) 2014-09-23.
