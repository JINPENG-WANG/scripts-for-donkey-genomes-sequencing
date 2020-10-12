#!/usr/bin/perl -w
use strict;
use warnings;


die "\@ARGV is required" if @ARGV != 6;

my $pi_file_dir=shift;  #high_low.chr4.div
my $zh_file_dir=shift;  
my $style1=shift;
my $style2=shift;
my $percentage=shift; #0.05  
my $out_dir=shift; 


my (%x,%y,%z,%hash1,%hash2,$window);
my @ar=glob "$pi_file_dir/*.div";
my @br=glob "$zh_file_dir/$style1.*.zh";  #high.chr4.zh
my @cr=glob "$zh_file_dir/$style2.*.zh";  #low.chr4.zh

for (@ar)
{
    open IN,"$_" or die $!;
    <IN>;
    while (my $ln=<IN>) #chr10   1       50000   0.00113021953634932     0.00119954696429286     252     -0.909104637276728      0.00105449210087993     0.009898089
    {
        chomp ($ln);
        my @a=split /\t/,$ln;
        next if $a[11] > 1;
        $x{$a[0]}{$a[1]}=$ln; #$x{chr10}{1}=chr10   1       50000   0.00113021953634932     0.00119954696429286     252 
    }
    close IN;
}

for (@br)
{
    open IN,"$_" or die $!; #print "$_\n";
    <IN>;
    while (my $ln=<IN>) #chr6    1       50000   0.107669951146834
    {
        chomp ($ln);
        my @a=split /\t/,$ln;
        $y{$a[0]}{$a[1]}=$ln;   #$y{chr6}{1}=chr6    1       50000   0.107669951146834
    }
    close IN;
}

for (@cr)
{
    open IN,"$_" or die $!;
    <IN>;
    while (my $ln=<IN>)
    {
        chomp ($ln);
        my @a=split /\t/,$ln;
        $z{$a[0]}{$a[1]}=$ln; #$z{chr6}{1}=chr6    1       50000   0.107669951146834
    }
    close IN;
}

for my $key1(sort keys %x)
{
    for my $key2(sort {$a<=>$b} keys %{$x{$key1}})
    {
        $window++;
        my ($zh1,$zh2,$value);
        my @a=split /\t/,$x{$key1}{$key2};
        if (exists $y{$key1}{$key2}){$zh1=(split /\t/,$y{$key1}{$key2})[3]}
        else {$zh1="-"}
        if (exists $z{$key1}{$key2}){$zh2=(split /\t/,$z{$key1}{$key2})[3]}
        else {$zh2="-"}
        $value=join "\t",$a[0],$a[1],$a[2],$a[11],$a[3],$a[6],$zh1,$a[7],$a[10],$zh2;
        $hash1{$window}=$value;
        $hash2{$window}=$a[11];
    }
}

open OUT1,">$out_dir/tend_to_$style1\_$style2.selective" or die $!;
open OUT2,">$out_dir/tend_to_$style1.selective" or die $!;
open OUT3,">$out_dir/tend_to_$style2.selective" or die $!;

print OUT1 "#Chr Start   End     Fst     pi      TjmD    zh1     pi2     TimD2   zh2\n";
print OUT2 "#Chr Start   End     Fst     pi      TjmD    zh1     pi2     TimD2   zh2\n";
print OUT3 "#Chr Start   End     Fst     pi      TjmD    zh1     pi2     TimD2   zh2\n";

my $border=$window*$percentage;
my $n;
my (%chr,%hash);
for (sort {$hash2{$b} <=> $hash2{$a}} keys %hash2)
{
    $n++;
    if ($n <= $border)
    {
        my @a=split /\t/,$hash1{$_};
        if ($a[5] eq "-" or $a[6] eq "-" or $a[8] eq "-" or $a[9] eq "-"){next}        
        if ($a[5]< -1 && $a[6]< -1 && $a[8]<-1 && $a[9]< -1){print OUT1 "$hash1{$_}\n"}
        if ($a[5]< -1 && $a[6]< -1){print OUT2 "$hash1{$_}\n"}
        if ($a[8]<-1 && $a[9]< -1){print OUT3 "$hash1{$_}\n"}
    }
    else {last}
}
close OUT1;
close OUT2;
close OUT3; 




__END__



open OUT,">$out_dir/chr_statistic1" or die $!;
print OUT "$_\t$chr{$_}\n" for sort keys %chr;
close OUT;

open OUT,">$out_dir/chr_statistic2";
print OUT "#Chr Start   End     Fst     pi      TjmD    zh1     pi2     TimD2   zh2\n";
for my $key1(sort keys %hash)
{
    print OUT "$hash{$key1}{$_}\n" for sort {$a<=>$b} keys %{$hash{$key1}};
}
close OUT;

 
