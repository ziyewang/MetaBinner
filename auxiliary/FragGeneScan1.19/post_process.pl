#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($genome_file, $gff_file, $new_gff_file);
my ($freq_start_file, $freq_start1_file);
my ($seq, $head);

GetOptions(
           'genome=s' => \$genome_file,
           'pre=s' => \$gff_file,
           'post=s' => \$new_gff_file,
           );



my %freq_unit =('AAA',0,'AAC',1,'AAG',2,'AAT',3,'ACA',4,'ACC',5,'ACG',6,'ACT',7,
		'AGA',8,'AGC',9,'AGG',10,'AGT',11,'ATA',12,'ATC',13,'ATG',14,'ATT',15,
		'CAA',16,'CAC',17,'CAG',18,'CAT',19,'CCA',20,'CCC',21,'CCG',22,'CCT',23,
		'CGA',24,'CGC',25,'CGG',26,'CGT',27,'CTA',28,'CTC',29,'CTG',30,'CTT',31,
		'GAA',32,'GAC',33,'GAG',34,'GAT',35,'GCA',36,'GCC',37,'GCG',38,'GCT',39,
		'GGA',40,'GGC',41,'GGG',42,'GGT',43,'GTA',44,'GTC',45,'GTG',46,'GTT',47,
		'TAA',48,'TAC',49,'TAG',50,'TAT',51,'TCA',52,'TCC',53,'TCG',54,'TCT',55,
		'TGA',56,'TGC',57,'TGG',58,'TGT',59,'TTA',60,'TTC',61,'TTG',62,'TTT',63);

my $program = $0;
my $dir = substr($0, 0, length($0)-15);
my (@lines, @lines1);


# get frequency model  
$freq_start_file = $dir."train/start";
$freq_start1_file = $dir."train/start1";

open(IN, $freq_start_file);
@lines = <IN>;
chomp(@lines);
close(IN);

open(IN, $freq_start1_file);
@lines1 = <IN>;
chomp(@lines1);
close(IN);

# get out file info
my ($h, %s);

open(IN, $gff_file);
while(my $each_line=<IN>){
    
    if ($each_line=~/^\>(\S+)/){
	$h = $1;
    }else{
	$s{$h} .= $each_line;
    }
}
close(IN);

open OUT, ">$new_gff_file";

$head = "";
$seq = "";
open(SEQ, $genome_file);
while(my $each_line=<SEQ>){
    
    chomp($each_line);
    if ($each_line=~/^\>(\S+)/){
	
	if (length($seq)>0){
	    if (exists $s{$head}){
		print OUT ">".$head."\n";
		call_post(\$seq, \*OUT, \$s{$head},\@lines, \@lines1);
	    }
	}
	$head = $1;
	$seq = "";
    }else{
	$seq .= $each_line;
    }
}
if (length($seq)>0){ #added by YY, Aug, 2013
    if (exists $s{$head}){
	print OUT ">".$head."\n";
	call_post(\$seq, \*OUT, \$s{$head},\@lines, \@lines1);
    }
}

close(OUT);
close(SEQ);

sub call_post{ #\genome_seq, OUT, \$s{$head}, \@lines, \@lines1

    my ($len, $count, $cg);
    my ($genome_seq, $temp_genome_seq);
    my ($id, $i, $j);
    my ($utr, $out, $e_save, $i_save, $codon, $command, $return_seq);
    my $sff;
    my (@freq, @freq1, @lines, @lines1);


    $genome_seq = ${$_[0]};
    $out = $_[1];
    $sff = ${$_[2]};
    @lines = @{$_[3]};
    @lines1 = @{$_[4]};


    #count cg contents
    $temp_genome_seq = $genome_seq;
    $len = length($genome_seq);
    $count = ($temp_genome_seq =~ s/C|G|c|g//g);
    $temp_genome_seq="";
    $cg = sprintf("%.0f", $count*100/$len);

    if ($cg < 26){
	$cg=26;
    }elsif($cg>69){
	$cg=69;
    }

    $id=0;
    for (my $j=($cg-26)*62+1; $j<($cg-26+1)*62; $j += 1){
	my @temp = split(/\s+/, $lines[$j]);
	for (my $i=0; $i<=63; $i++){
	    
	    $freq[$id][$i] = $temp[$i];
	}
	$id+=1;
    }
    
    
    $id=0;
    for (my $j=($cg-26)*62+1; $j<($cg-26+1)*62; $j += 1){
	my @temp = split(/\s+/, $lines1[$j]);
	for (my $i=0; $i<=63; $i++){
	    
	    $freq1[$id][$i] = $temp[$i];
	}
	$id+=1;
    }


    # find the optimal start codon with 30bp up- and downstream of start codon
    my @sffs = split(/\n/, $sff);
    for(my $i=0; $i<=$#sffs; $i+=1){

	my $each_line = $sffs[$i];
	my @temp = split(/\s+/, $each_line);

	$i_save=0;
	if ($temp[2] eq "+"){
	    my $i=0;
	    $codon = substr($genome_seq, $temp[0]-1, 3);
	    
	    while ($codon !~ /TAA|TAG|TGA/ &&$temp[0]-1-$i-35>=0  ){
		
		if ($codon =~ /ATG|GTG|TTG/){
		    $utr = substr($genome_seq, $temp[0]-1-$i-30, 63);
		    
		    my @temp1 = split(//, $utr);
		    my $freq_sum = 0;
		    if ($#temp1==62){
			for($j=0; $j<=60; $j++){
			    my $key = $temp1[$j].$temp1[$j+1].$temp1[$j+2];
			    if (exists $freq_unit{$key}){
				$freq_sum -= log($freq[$j][$freq_unit{$key}]);
			    }else{
				$freq_sum -= log(1/64);
			    }
			}
		    }
		    if ($i==0){
			$e_save = $freq_sum;
			$i_save = 0;
		    }elsif ($freq_sum < $e_save){
			$e_save = $freq_sum;
			$i_save = -1*$i;
		    }
		}
		$i += 3;
		$codon = substr($genome_seq, $temp[0]-1-$i, 3);
	    }
	    print $out eval($temp[0]+$i_save)."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\n";
	    
	}elsif ($temp[2] eq "-"){
	    my $i=0;
	    $codon = substr($genome_seq, $temp[1]-1-2, 3);
	    while ($codon !~ /TTA|CTA|TCA/ && $temp[1]-2+$i+35<length($genome_seq) ){
		
		if ($codon =~ /CAT|CAC|CAA/){
		    $utr = substr($genome_seq, $temp[1]-1-2+$i-30, 63);
		    my @temp1 = split(//, $utr);
		    my $freq_sum = 0;
		    if ($#temp1==62){
			for($j=0; $j<=60; $j++){
			    my $key = $temp1[$j].$temp1[$j+1].$temp1[$j+2];
			    if (exists $freq_unit{$key}){
				$freq_sum -= log($freq1[$j][$freq_unit{$key}]);
			    }else{
				$freq_sum -= log(1/64);
			    }
			}
		    }
		    if ($i==0){
			$e_save = $freq_sum;
			$i_save = 0;
		    }elsif ($freq_sum < $e_save){
			$e_save = $freq_sum;
			$i_save = $i;
		    }
		}
		$i += 3;
		$codon = substr($genome_seq, $temp[1]-1-2+$i, 3);
	    }
	    print $out $temp[0]."\t".eval($temp[1]+$i_save)."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\n";
	}
    }
}
