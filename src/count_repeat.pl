#!/usr/bin/perl -w
use strict;

#die "input file" if @ARGV !=1;
#open A,$ARGV[0] or die "$ARGV[0] can't be opened";

my %EX;
while  (<>) {
	chomp;
	($EX{$_}) ? ($EX{$_}++) : ($EX{$_}=1);
}
#close A;

foreach my $line (sort keys %EX) {
	print "$line\t$EX{$line}\n";
}
