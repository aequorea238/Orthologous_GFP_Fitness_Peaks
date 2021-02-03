#!/usr/bin/perl -w

use strict;

my $line1; my $gi; my $seq; my $i; my $j; my $gi2;
my $codon; my %positions; my $new; my $old; my $line2;
my $line; my $temp; my $file; my $from; my $to; my $pos;

$file = $ARGV[0];
open (INFILE, "$file") or die "Kozyol!";
<INFILE>;
$seq=<INFILE>;

$old=0;
while($seq)
{
chomp($seq);
$codon = substr($seq,0,3);
$seq   = substr($seq,3);

if ($codon eq "---"){$new++;}
else{$old++;$new++;}

$positions{$old}=$new;
}

close(INFILE);

$file = $ARGV[1];
open (INFILE, "$file") or die "Kozyol!";

while(<INFILE>)
{$line=$_;

($pos)=($line =~ m/\w(\d+)/);
($from) =($line =~ m/(\w)\d+\w/);
($to) =($line =~ m/\w\d+(\w)/);

print "$from$positions{$pos}$to\n";

}





