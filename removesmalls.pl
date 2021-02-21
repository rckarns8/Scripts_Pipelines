#Date Initiated: 3/25/2020
#Created by: Rachael Storo
#Goal: Filter assembled contigs to only keep those above a certain length, as well as report the number of contigs > a given size (2kb here)


## removesmalls.pl
#!/usr/bin/perl

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


\
