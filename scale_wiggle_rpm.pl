#!opt/bin/perl -w
use strict;

## JBT last updated 10/03/2025
## JMC last updated 03/01/2022 

my($in, $reads, $out)=@ARGV;

open (IN, "$in") or die "Cant open infile $in\n";
open (OUT, ">$out") or die "cant open $out\n";

while (my $line=<IN>) {
    chomp $line;
    if ($line=~/^\D/) {
        print OUT "$line\n";
    } else {
        my @array=split(/\s/, $line);
        my $pos=$array[0];
        my $val=$array[1];
        my $scaled=($val*1000000/$reads);
        print OUT "$pos\t$scaled\n";
    }
}
