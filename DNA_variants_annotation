=pod
1. the perl is to annotate the VCF-formated SNPs or Indels
2. the perl is only fit for bi allele variant
3. Usage: perl ./annotation.pl step_size gff_file vcf_file output

step_size INT upstream and downstream +- Kb of genes   eg. 2000
gff_file  FILE the gene set file. The gff_file must be convert to 9 columns formats as follows.
1       Hv_IBSC_PGSB_r1 mRNA    41961   45310   1       +       .       ID=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 5UTR    41961   41962   .       +       .       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     41963   42213   .       +       0       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     42300   42338   .       +       1       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     42427   42567   .       +       1       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     42862   43000   .       +       1       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     43072   43164   .       +       0       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     43649   43714   .       +       0       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     43970   44086   .       +       0       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 CDS     44196   44711   .       +       0       Parent=HORVU1Hr1G000010.3;
1       Hv_IBSC_PGSB_r1 3UTR    44712   45310   .       +       .       Parent=HORVU1Hr1G000010.3;
VCF FILE the compressed input VCF file
output FILE the output annotation file
=cut

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 4;

my $step=shift; # eg. gene up down +- 1000bp
my $gff=shift;
my $vcf=shift;
my $output=shift;

my (%x,$total);
open IN,"gzip -dc $vcf | grep -v \"^#\"|" or die $!;
while (<IN>)
{   
    chomp;
    my @a=split;
    $x{$a[0]}{$a[1]}=1;
    $total++;
}
close IN;

my ($mrna,$cds,$utr5,$utr3,$stream);
open IN,"gzip -dc $gff|" or die $!;
while (<IN>)
{
    chomp;
    my @a=split;
    if ($a[2] eq 'mRNA')
    {
        for my $key2($a[3]..$a[4]){$mrna++ if exists $x{$a[0]}{$key2}}

        my $up_st=$a[3]-$step;
        my $up_end=$a[3]-1;
        for my $key2($up_st..$up_end){$stream++ if exists $x{$a[0]}{$key2}}
        my $down_st=$a[4]+1;
        my $down_end=$a[4]+$step;
        for my $key2($down_st..$down_end){$stream++ if exists $x{$a[0]}{$key2}}
    }
    elsif ($a[2] eq 'CDS'){for my $key2($a[3]..$a[4]){$cds++ if exists $x{$a[0]}{$key2}}}
    elsif ($a[2] eq '5UTR'){for my $key2($a[3]..$a[4]){$utr5++ if exists $x{$a[0]}{$key2}}}
    elsif ($a[2] eq '3UTR'){for my $key2($a[3]..$a[4]){$utr3++ if exists $x{$a[0]}{$key2}}}
}
close IN;

my $intergenic=$total-$mrna-$stream;
my $intron=$mrna-$cds-$utr5-$utr3;

open OUT,">$output";
print OUT "cds\t$cds\n5UTR\t$utr5\n3UTR\t$utr3\nintron\t$intron\nstream\t$stream\nintergenic\t$intergenic\n";
close OUT;
