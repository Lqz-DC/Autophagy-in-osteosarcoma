use strict;
use warnings;

my $pvalueOpt=0.05;
my $normalCount=0;
my $tumorCount=0;

open(RF,"CIBERSORT-Results.txt") or die $!;
open(WF,">CIBERSORT.filter.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	my $RMSE=pop(@arr);
	my $Correlation=pop(@arr);
	my $Pvalue=pop(@arr);
	if($.==1){
		print WF "id\t" . join("\t",@arr) . "\n";
	}
  else{
  	if($Pvalue<$pvalueOpt){
  	  print WF "$sample\t" . join("\t",@arr) . "\n";
  	  my @sampleArr=split(/\-/,$sample);
  	  if($sampleArr[3]=~/^0/){
  	  	$tumorCount++;
  	  }
  	  else{
  	  	$normalCount++;
  	  }
  	}
  }
}
close(WF);
close(RF);

print "normal count: $normalCount\n";
print "tumor count: $tumorCount\n";
