#!/usr/bin/perl -w
use strict;

#########this script is used join the GeneName with topgo output#################
#

die "the first is the ATH_tair10_go2geneid.ALL.map the second is the GO_analysis output,the third is the DE file\n" if @ARGV !=4;

open A,$ARGV[0] or die "the $ARGV[0] cann't be opened";
open B,$ARGV[1] or die "the $ARGV[1] cann't be opened";
open C,$ARGV[2] or die "the $ARGV[2] cann't be opened";
open D,$ARGV[3] or die "the $ARGV[3] cann't be opened";

my %GO;
while (<D>) {
	chomp;
	my @array=split"\t",$_;
    $GO{$array[0]}=$array[1];
}
my %info;
while (<A>) {
	chomp;
	my @array=split"\t",$_;
	my @brray=split (/,/,$array[1]);
	foreach my $gene (@brray) {
		push @{$info{$array[0]}},$gene;
	}
}
my @DE;
while (<C>)  {
	chomp;
	my @array=split "\t",$_;
	push @DE,$array[0];
}

while (<B>) {
	chomp;
	my @array=split "\t",$_;
    if (exists $GO{$array[0]}) {
        $array[1]=$GO{$array[0]};
    }
    $_=join "\t",@array;
	my $str="";
	if ($_!~/^GO\.ID/) {
		if (exists $info{$array[0]}) {
			foreach my $ge (sort @{$info{$array[0]}}) {
				if (grep /^$ge$/, @DE) {
					if ($str eq "") {
						$str=$ge;
					}else{
					$str=$str.",$ge";
					}
				}
			}
			if ($str eq "") {$str="No Gene"};
			print "$_\t$str\n";
			$str="";
		}
	}else{
		print "$_\tDEG_list\n"
	}
}




