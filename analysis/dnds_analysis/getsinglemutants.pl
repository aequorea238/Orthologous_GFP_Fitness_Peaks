#!/usr/bin/perl -w

use strict;

my ($line, $seq, $codon, $node, $f, $s, $t, %branches, %matrix, $new, $stop);
my (%code);


open (INFILE, "codons.txt") or die "Kozyol!";
while(<INFILE>)
{$line=$_;

if ($line =~ m/\w\w\w/)
{($codon)=($line =~ m/(\w\w\w)/);

($f)=($codon =~ m/(\w)/);
($s)=($codon =~ m/\w(\w)/);
($t)=($codon =~ m/\w\w(\w)/);

   if ($f eq "A"){$new = "T" . $s . $t . "\t" . "G" . $s . $t . "\t" . "C" . $s . $t . "\t"; print $new;}
elsif ($f eq "T"){$new = "A" . $s . $t . "\t" . "G" . $s . $t . "\t" . "C" . $s . $t . "\t"; print $new;}
elsif ($f eq "G"){$new = "A" . $s . $t . "\t" . "T" . $s . $t . "\t" . "C" . $s . $t . "\t"; print $new;}
elsif ($f eq "C"){$new = "A" . $s . $t . "\t" . "T" . $s . $t . "\t" . "G" . $s . $t . "\t"; print $new;}

   if ($s eq "A"){$new = $f . "T" . $t . "\t" . $f . "G" . $t . "\t" . $f . "C" . $t . "\t"; print $new;}
elsif ($s eq "T"){$new = $f . "A" . $t . "\t" . $f . "G" . $t . "\t" . $f . "C" . $t . "\t"; print $new;}
elsif ($s eq "G"){$new = $f . "A" . $t . "\t" . $f . "T" . $t . "\t" . $f . "C" . $t . "\t"; print $new;}
elsif ($s eq "C"){$new = $f . "A" . $t . "\t" . $f . "T" . $t . "\t" . $f . "G" . $t . "\t"; print $new;}

   if ($t eq "A"){$new = $f . $s . "T" . "\t" . $f . $s . "G" . "\t" . $f . $s . "C" . "\n"; print $new;}
elsif ($t eq "T"){$new = $f . $s . "A" . "\t" . $f . $s . "G" . "\t" . $f . $s . "C" . "\n"; print $new;}
elsif ($t eq "G"){$new = $f . $s . "A" . "\t" . $f . $s . "T" . "\t" . $f . $s . "C" . "\n"; print $new;}
elsif ($t eq "C"){$new = $f . $s . "A" . "\t" . $f . $s . "T" . "\t" . $f . $s . "G" . "\n"; print $new;}

}
else{print "\n";}

}

close(INFILE);

