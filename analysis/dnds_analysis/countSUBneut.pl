#!/usr/bin/perl -w

use strict;

my ($line, $seq, $codon, $node, $count, $branch, $nodes, %countneu, %countsub, $start, $stop);
my (%code, %count, %substitutions, %neutral, %reverse, $temp, $starts, $stops, $pos);


open (INFILE, "ok_mutations_in_predictions_av_pos_byalignment_substitutions.txt") or die "Kozyol!";

while(<INFILE>)
{$line=$_;

chomp($line);

$substitutions{$line}=1;
}

close(INFILE);


open (INFILE, "ok_mutations_in_predictions_av_pos_byalignment_neutral.txt") or die "Kozyol!";

while(<INFILE>)
{$line=$_;

chomp($line);

$neutral{$line}=1;
}

close(INFILE);

open (INFILE, "GFP.rst.branches.allstates") or die "Kozyol!";

while(<INFILE>)
{$line=$_;

if ($line =~ m/\.\./){$branch=$line; chomp($branch);}
else
{
($pos)=($line =~ m/(\d+)/);
($start)=($line =~ m/\d+\s+(\S+)/);
($stop)=($line =~ m/\d+\s+\S+\s+(\S+)/);
$codon = $start . $pos . $stop;

if ($substitutions{$codon}){$countsub{$branch}++;$count{$branch}++;}
if ($neutral{$codon}){$countneu{$branch}++;$count{$branch}++;}

}

}

open (OUTFILE, ">GFP.rst.branches.allstates.counts") or die "Kozyol!";


foreach $branch (keys %count)
{

if ($countsub{$branch}){}else{$countsub{$branch}=0;}
if ($countneu{$branch}){}else{$countneu{$branch}=0;}

print OUTFILE "$branch\t$countsub{$branch}\t$countneu{$branch}\n";

}
close(OUTFILE);


