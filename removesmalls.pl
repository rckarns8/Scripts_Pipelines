#Date Initiated: 3/25/2020
#Created by: Rachael Storo
#Goal: Filter assembled contigs to only keep those above a certain length, as well as report the number of contigs > a given size (2kb here)


## removesmalls.pl
#!/usr/bin/perl

module load Perl/5.26.1-GCCcore-6.4.0

use strict;
use warnings;

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}


perl /scratch/rck80079/Baker_Tutorial/scripts/removesmalls.pl 2000 final.contigs.fa > contigs-2000.fasta


grep -o '>' contigs-2000.fasta | wc -l > numb_2kb.txt
