=pod

1. the perl is to calculate the nucleotide diversity (pi) for each SNPs site (Tajima F. Evolutionary relationship of DNA sequences in finite populations. Genetics 1983)
2. the diploid SNPs data will be convert to haploid SNPs data in this perl, by coding all heterozygous genotype calls as missing data.
3. the converted haploid SNPs data will be filtered by an missing rate
4. non-polymorphic variants will be pass 
5. the perl is only fit for : (i) VCFv4.1 or other similar format (GT:AD:DP:GQ:PGT:PID:PL or GT:AD:DP:GQ:PL) (ii) bi allele variant
6. Usage: perl ./vcf_pi.hap.pl maximum_missing_rate input output

maximum_missing_rate FLOAT the maximum proportion of the missing genotype calls after conversion

input FILE     the compressed input VCF file   eg. in.vcf.gz

output FILE    the value of nucleotide diversity per SNPs position

note: the output files contant three columns: chromosome, SNPs position, pi

=cut

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 3;

my $max_miss_rate=shift;
my $vcf=shift;
my $output=shift;

open IN,"gzip -dc $vcf|grep -v \"^#\"|" or die $!;
open OUT,">$output" or die $!;
print OUT "#chr\tpos\tvalue\n";
while (<IN>)
{
    chomp;
    my @a=split;
    my $miss=0;
    my (@geno,%geno);
    for (my $i=9;$i<@a;$i++)
    {        
        if ($a[$i]=~/^([^:])\/([^:]):/)
        {
            if ($1 eq '0' && $2 eq '0'){push @geno,'a';$geno{"0"}=1}
            elsif ($1 eq '1' && $2 eq '1'){push @geno,'b';$geno{"1"}=1}
            else {$miss++}
        }
    }
    next if $miss/(@a-9) > $max_miss_rate;
    if (not exists $geno{"0"} or not exists $geno{"1"}){next}
    my $pi=&pi(@geno);
    print OUT "$a[0]\t$a[1]\t$pi\n";
}
close IN;
close OUT;

sub pi
{
    my @a=@_;
    my $diff=0;
    my $time=0;

    for (my $i=0;$i<@a;$i++)
    {
        for (my $m=0;$m<@a;$m++)
        {
            if ($m>$i)
            {
                $diff++ if $a[$i] ne $a[$m];
                $time++;
            }
        }
    }
    my $value=$diff/$time;
    return $value;
}    
