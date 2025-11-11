#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Time::Local;

#Last updated 6/16/20

#make wig file from bed

###   MAIN   ###

my ($in_sam, $name, $color, $peflag, $log, $binsize)=@ARGV;

if (!defined $peflag) {
    die "paired end flag is not defined; must be y or n ; \n";
}

if (!defined $binsize) {
    die "please specify a binsize, eg 50\n";
}

if (!defined $log) {
    die "please specify y or n for whether data should be log10 normalized\n";
}

my %colors = (
	'pink' => '255,181,197',
	'babyblue' => '135,206,250',
	'grey' => '136,136,136',
	'black' =>'0,0,0',				
	'red' =>'255,0,0',			
	'green' =>'0,51,0',	
	'blue' =>'0,0,255',			
	'yellow' =>'255,204,0',
	'lightgreen' =>'0,112,0',
	'purple' =>'51,0,51',		
	'navyblue'=>'0,0,51',
	'maroon'=> '51,0,0',
	'brown'=>'51,26,0',			
	'orange'=>'255,77,0',
	'magenta'=>'255,0,255',
	);
unless (exists $colors{$color}) 
{
	die "$color is not an accepted color.\nPlease use pink, babyblue,"
	    ."grey, black, red, green, blue, yellow, lightgreen, purple, navyblue maroon, brown,"
	    ." orange, or magenta.\n";
}

my $headerName = $name;
$headerName =~ s/\A.*\///g;
my $outName = $name.".wig";

my ($second,$minute,$houre,@reste) =localtime(time);
print $houre.":".$minute.":".$second." start making wigs...\n";


open (OUT, ">".$outName) or die "Cannot open ".$outName."!\n";
print OUT "track type=wiggle_0 visibility=full" .
	" name=\"".$headerName."\"  color=".$colors{$color}.
	" maxHeightPixels=128:40:11 ".
	" group=\"user\" graphType=bar".
	" priority=100 viewLimits=0:2.8 autoscale=on\n";

my %lengths = (					#based on mm10
	'chr1' => '195471971',
	'chr2' => '182113224',
	'chr3' => '160039680',
	'chr4' => '156508116',
	'chr5' => '151834684',
	'chr6' => '149736546',
	'chr7' => '145441459',
	'chr8' => '129401213',
	'chr9' => '124595110',
	'chr10' => '130694993',
	'chr11' => '122082543',
	'chr12' => '120129022',
	'chr13' => '120421639',
	'chr14' => '124902244',
	'chr15' => '104043685',
	'chr16' => '98207768',
	'chr17' => '94987271',
	'chr18' => '90702639',
	'chr19' => '61431566',
	'chrX' => '171031299',
	#'chrY' => '91744698',
    );

foreach my $chr (sort keys %lengths) {

	my $wighoa = &readsam($in_sam, $chr, $peflag, $binsize);

	if ($wighoa !~ /nochr/) {
		my ($sec,$min,$hour,@rest) = localtime(time);
		print $hour.":".$min.":".$sec." print wig for ".$chr."...\n";
		printer($wighoa, $binsize);
	}
}

### SUBS ###
#store bed hits in an array that counts 1 for every 50bp covered by 
# hit.

sub readsam {
    my ($in, $refchr,  $pe, $bs)=@_;
    my %coords;

      my $togrep = $refchr."[[:space:]]";
    my $temp = $name."_".$refchr."_temp.txt";
    
    print "grep ".$togrep." ".$in." ".$temp."\n";
    `grep $togrep $in > $temp`;
    #my $wc=`wc temp.txt`;
    #print "$wc lines in temp.txt\n";
    
    open (IN, $temp) or die "Cannot open the temporary file ".$temp."!\n";
    
    my $j;
    my $k;
    while (my $line=<IN>) {
	chomp $line;
	my @info=split(/\t/, $line);
	my $lchr=$info[2];
	
	if (($lchr!~/M/) && ($lchr!~/_/)) {
	    
	    if ($refchr !~ /^$lchr$/) {
		print $refchr." does not match:\n".$line."\n";
	    }
	    
	    $k++;
	    
	    #print "$line\n";
	    
	    my ($subreads, $pass, $tot)=splitsam($line, $name, $pe);
	    
	    foreach my $sr (@{$subreads}) {
		my @array=split(/\t/, $sr);
		
		my $chr=$array[1];
		my $start=$array[2];
		my $end=$array[3];
		
		#bins to add seqs to
		my $fbin = int ($start/$bs);
		my $ebin = int ($end/$bs);
		
		#counter to add seq to bins
		my $i = $fbin;
		
		#add a count to each bin of the chr@ between start and
		# end
		while ($i<=$ebin) {
		    @{$coords{$chr}}[$i]+=1;
		    $i++;
		}
	    }
	}
	
	$j++;
    }
    
    close IN;
    `rm -f $temp`;
    
    if (defined $j) {
	print $j." lines read, ".$k." total hits from ".$refchr."\n";
	return \%coords;
    } else {
	return 'nochr';
    }
}

#sub to print out wiggle format
sub printer {
	my ($counts, $bs)=@_;

	foreach my $chr (sort keys %{$counts}) {
		print OUT "variableStep chrom=".$chr." span=50\n";
		my @hits = @{$counts->{$chr}};
		my $i;
		foreach my $hit (@hits) {
			if (defined $hit) {

			    if ($log eq 'y') {
				my $longlog=(log $hit)/(log 10);
				my $log=sprintf("%.2f", $longlog);
				$hit=$log;
			    } elsif ($log eq 'n') {
				$hit=$hit;
			    } else {
				die "please specify y or n for whether data"
				    ."should be log10 normalized -- see usage at top of script\n";
			    }
			    my $pos=$bs*($i);
			    print OUT $pos." ".$hit."\n";
			    
			}
			$i++;
		}
	}
}


#split each read into multiple read equivalents based on cigar string
sub splitsam {
    my ($line, $out, $pe)=@_;

    my $cig="$out"."_CIG";
    open (CIG, ">$cig") or die "cant open $out cig\n";

    my $i="0";
    my $j="0";
    
    chomp $line;
    my @array=split(/\t/, $line);
    my $id=$array[0];
    my $flag=$array[1];  
    my $chr=$array[2];  
    my $start=$array[3];  
    my $mq=$array[4];  
    my $cigar=$array[5]; 
    my $seq=$array[9];
    my $qual=$array[10];

    my $ocigar=$cigar;
    
    #info for each subread
    my @subreads;

    #forward or reverse read
    my $mate="f";

    $j++;
    
    #not yet sure what to do with cigars with these operators
    if ($cigar!~/(P|S|H|X)/) {
	
	#find strand
	#flag webpage:
        #https://broadinstitute.github.io/picard/explain-flags.html

	my $strand;
	if ($pe=~/y/) {
	    if ($flag eq 99) { 
		$strand='+';
	    } elsif ($flag eq 147) {
		$strand='+';
		$mate="r";
	    } elsif ($flag eq 83) {
		$strand='-';
	    } elsif ($flag eq 163) {
		$strand='-';
		$mate="r";
	    } #elsif ($flag eq 73) { 
	    # 	$strand='+';
	    # } elsif ($flag eq 137) {
	    # 	$strand='+';
	    # 	$mate="r";
	    # } elsif ($flag eq 89) {
	    # 	$strand='-';
	    # } elsif ($flag eq 153) {
	    # 	$strand='-';
	    # 	$mate="r";
	    # }
	} else {
	    if ($flag & 16) { 
		$strand='-';
	    } else {
		$strand='+';
	    }
	}
	
	#parsing cigar string

	if (defined $strand) {
	    $i++;
	    #string for cigar length
	    my $cl=length($cigar);
	    #array info for cigar
	    my @cinfo;
	    while ($cl>0) {
		$cigar=~s/^(\d+\D)(.*)/$2/;
		push (@cinfo, $1);
		$cl=length($cigar)
	    }
	    
	    #go through cigar info and make new start end coords, seq for the read
	    my ($end);
	    my $sr;
	    foreach my $action (@cinfo) {
		my $srinfo;
		if ($action=~/(\d+)M/) {
		    $sr++;
		    $end=$start+$1-1;
		    my $nid="$id"."__"."$sr"."$mate";
		    
		    #extract from beginning of seq
		    my $toextract=$end-$start+1;
		    $seq=~s/^([a-zA-Z]{$toextract})(.*)/$2/;
		    my $sseq=$1;
		    $qual=~s/^([\S+]{$toextract})(.*)/$2/;
		    my $squal=$1;
		    $srinfo="$nid\t$chr\t$start\t$end\t$strand\t$sseq\t$squal\t".
			"$action\t$ocigar";
		    push (@subreads, $srinfo);
		    $start=$end+1;
		} elsif  ($action=~/(\d+)N/) {
		    $start=$start+$1;
		} elsif  ($action=~/(\d+)D/) {
		    #read skips \d base in reference
		    $start=$start+$1;
		} elsif  ($action=~/(\d+)I/) {
		    #read has base inserted relative to reference
		    $seq=~s/^([a-zA-Z]{$1})(.*)/$2/;
		    $qual=~s/^([\S+]{$1})(.*)/$2/;
		} 
	    }
	} else {
	    #print out to file with non-canonical cigar strings
	    print CIG "$line\n";
	}
    }
    return (\@subreads, $i, $j);
}