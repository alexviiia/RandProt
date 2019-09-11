#!/usr/bin/env bats

# test scripts as run from local directory

@test "kMax.pl help" {
    # missing arguments
    run perl -w kMax.pl
    [ $status -eq 0 ]
}

@test "kMax.pl main" {
    run perl -w kMax.pl sample/PLAF7.fa
    [ $status -eq 0 ]
    # exactly three lines in output
    [ "${#lines[@]}" -eq 3 ]
    # match each of the lines to its correct value (from previous testing)
    [ "${lines[0]}" = "Input: sample/PLAF7.fa" ]
    [ "${lines[1]}" = "Number of letters: 4139904" ]
    [ "${lines[2]}" = "kMax: 5" ]
}

@test "kCov.pl help" {
    # missing arguments
    run perl -w kCov.pl
    [ $status -eq 0 ]
}

@test "kCov.pl main 5" {
    run perl -w kCov.pl sample/PLAF7.fa 5
    [ $status -eq 0 ]
    # exactly 5 lines in output
    [ "${#lines[@]}" -eq 5 ]
    # match each of the lines to its correct value (from previous testing)
    [ "${lines[0]}" = "k: 5" ]
    [ "${lines[1]}" = "Number of theoretical amino acid k-mers: 3200000" ]
    [ "${lines[2]}" = "Scanning proteome..." ]
    [ "${lines[3]}" = "Number of observed amino acid k-mers: 1226974" ]
    [ "${lines[4]}" = "Proportion of k-mers observed: 38.34 %" ]
}

@test "kCov.pl main 4" {
    run perl -w kCov.pl sample/PLAF7.fa 4
    [ $status -eq 0 ]
    # exactly 5 lines in output
    [ "${#lines[@]}" -eq 5 ]
    # match each of the lines to its correct value (from previous testing)
    [ "${lines[0]}" = "k: 4" ]
    [ "${lines[1]}" = "Number of theoretical amino acid k-mers: 160000" ]
    [ "${lines[2]}" = "Scanning proteome..." ]
    [ "${lines[3]}" = "Number of observed amino acid k-mers: 148284" ]
    [ "${lines[4]}" = "Proportion of k-mers observed: 92.68 %" ]
}

@test "kCov.pl main 3" {
    run perl -w kCov.pl sample/PLAF7.fa 3
    [ $status -eq 0 ]
    # exactly 5 lines in output
    [ "${#lines[@]}" -eq 5 ]
    # match each of the lines to its correct value (from previous testing)
    [ "${lines[0]}" = "k: 3" ]
    [ "${lines[1]}" = "Number of theoretical amino acid k-mers: 8000" ]
    [ "${lines[2]}" = "Scanning proteome..." ]
    [ "${lines[3]}" = "Number of observed amino acid k-mers: 7999" ]
    [ "${lines[4]}" = "Proportion of k-mers observed: 99.99 %" ]
}

@test "randProt.pl help" {
    # missing arguments
    run perl -w randProt.pl
    [ $status -eq 0 ]
}

@test "randProt.pl main n=1" {
    file_out='sample/NEW-PLAF7.k3.n1.fa.gz'
    run perl -w randProt.pl sample/PLAF7.fa "$file_out" 3 1 # test n=1 case
    [ $status -eq 0 ]
    # expect these many lines (same as input FASTA file)
    lines=`zcat "$file_out" | wc -l` # capture output of this command into variable
    [ "$lines" = "10882" ]
    # cleanup
    rm "$file_out"
}

@test "randProt.pl main n=2" {
    file_out='sample/NEW-PLAF7.k3.n2.fa.gz'
    run perl -w randProt.pl sample/PLAF7.fa "$file_out" 3 2
    [ $status -eq 0 ]
    # expect these many lines (same as input FASTA file *times 2*)
    lines=`zcat "$file_out" | wc -l` # capture output of this command into variable
    [ "$lines" = "21764" ]
    # cleanup
    rm "$file_out"
}
