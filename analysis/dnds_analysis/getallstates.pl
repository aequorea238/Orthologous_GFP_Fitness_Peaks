#!/usr/bin/perl -w

use strict;

my ($line, $seq, $codon, $node, $count, $branch, $nodes, %branches, %matrix, $start, $stop, %ds, $ds);


open (INFILE, "GFP.rst.branches") or die "Kozyol!";
while(<INFILE>)
{$line=$_;

if ($line =~ m/Branch/)
{($branch)=($line =~ m/(\d+)/);
 ($nodes)=($line =~ m/(\d+\.\.\d+)/);

$branches{$branch}=$nodes;
($ds)=($line =~ m/s=(.*)\)/);

$ds =~ s/ //g;

$ds{$branch}=$ds;

}

}

close(INFILE);

open (INFILE, "GFP.rst.nodes") or die "Kozyol!";

while(<INFILE>)
{$line=$_;

$node++;
$count=0;

if ($line =~ m/node/){($seq)=($line =~ m/node\s+\S+\s+(.*)/);}
else{($seq)=($line =~ m/\S+\s+(.*)/);}

while($seq)
{
$count++;
($codon)=($seq =~ m/(\S+)/);
($seq)  =($seq =~ m/\S+\s+(.*)/);

$matrix{$node}{$count}=$codon;
}
}

close(INFILE);

foreach $branch (keys %branches)
{
$nodes=$branches{$branch};

($start)=($nodes =~m/(\d+)/);
($stop)=($nodes =~m/\d+\.\.(\d+)/);

$count=0;

print "$branch\t$nodes\t$ds{$branch}\n";


while($count < 306)
{$count++;

print "$count\t$matrix{$start}{$count}\t$matrix{$stop}{$count}\n";

}

}


