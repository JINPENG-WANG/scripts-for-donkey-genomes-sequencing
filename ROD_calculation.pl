=pod

1. the perl is to calculate the ROD for each SNPs site.

2. Usage: perl ./ROD_calculation.pl input output


input FILE     the file containing the pi of dun donkeys and non_dun donkeys.

output FILE    the value of ROD per SNPs position

note: the output files contant three columns: chromosome, SNPs position, ROD

=cut




#!/usr/bin/perl -w
use strict;
use IO::File;
my @ARGV;
my $pi_file=shift @ARGV;
my $ROD_file=shift @ARGV;

my $fh_in = IO::File->new("$pi_file",'r');
my $fh_out = IO::File->new(">$ROD_file\n");
my $line = 0;
while(<$fh_in>){
	chomp;
	$line ++;
	if($line>1){
		my ($chr, $pos, $pi_dun, $pi_non_dun)=split /\s+/, $line;
		my $ROD=1-$pi_non_dun/$pi_dun;
		$fh_out->print("$chr\t$pos\t$ROD\n");
	}
}

		


	

