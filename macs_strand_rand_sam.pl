#!opt/bin/perl -w
use strict;

# Same version at https://github.com/CalabreseLab/Airn_Xist_manuscript/blob/main/RIPanalysis/macs_strand_rand_sam.pl

my ($in, $out)=@ARGV;

open (IN, "$in") or die "cant open $in in\n";
open (OUT, ">$out") or die "cant open $out in\n";

my ($pos, $neg);

while (my $line=<IN>) {
    chomp $line; #  remove the last trailing newline from the input string

    if ($line=~/^@/) {
	print OUT "$line\n"; # print header
    } else {
	my ($id, $strand, @array)=split(/\t/, $line);
	
	my $rand=rand(2);
	
	my $nstrand;
	if ($rand>1) {
	    $nstrand="16";
	    $neg++;
	} else {
	    $nstrand="0";
	    $pos++;
	}
	
	my $nl="$id\t$nstrand\t";
	foreach my $val (@array) {
	    $nl.="$val\t";
	}
	chop $nl;
	print OUT "$nl\n";
    }
}

print "$pos reads assigned to positive strand, $neg to negative\n";