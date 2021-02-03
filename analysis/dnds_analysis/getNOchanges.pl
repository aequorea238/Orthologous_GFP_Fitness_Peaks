#!/usr/bin/perl -w

use strict;

my ($line, $seq, $codon, $node, $count, $branch, $nodes, %branches, %matrix, $start, $stop);
my (%code, %mutations, %substitutions, %neutral, %reverse, $temp, $starts, $stops, $pos);

$code{"AAA"}="K";
$code{"AAT"}="N";
$code{"AAG"}="K";
$code{"AAC"}="N";

$code{"ATA"}="I";
$code{"ATT"}="I";
$code{"ATG"}="M";
$code{"ATC"}="I";

$code{"AGA"}="R";
$code{"AGT"}="S";
$code{"AGG"}="R";
$code{"AGC"}="S";

$code{"ACA"}="T";
$code{"ACT"}="T";
$code{"ACG"}="T";
$code{"ACC"}="T";

$code{"TAA"}="X";
$code{"TAT"}="Y";
$code{"TAG"}="X";
$code{"TAC"}="Y";

$code{"TTA"}="L";
$code{"TTT"}="F";
$code{"TTG"}="L";
$code{"TTC"}="F";

$code{"TGA"}="X";
$code{"TGT"}="C";
$code{"TGG"}="W";
$code{"TGC"}="C";

$code{"TCA"}="S";
$code{"TCT"}="S";
$code{"TCG"}="S";
$code{"TCC"}="S";

$code{"GAA"}="E";
$code{"GAT"}="D";
$code{"GAG"}="E";
$code{"GAC"}="D";

$code{"GTA"}="V";
$code{"GTT"}="V";
$code{"GTG"}="V";
$code{"GTC"}="V";

$code{"GGA"}="G";
$code{"GGT"}="G";
$code{"GGG"}="G";
$code{"GGC"}="G";

$code{"GCA"}="A";
$code{"GCT"}="A";
$code{"GCG"}="A";
$code{"GCC"}="A";

$code{"CAA"}="Q";
$code{"CAT"}="H";
$code{"CAG"}="Q";
$code{"CAC"}="H";

$code{"CTA"}="L";
$code{"CTT"}="L";
$code{"CTG"}="L";
$code{"CTC"}="L";

$code{"CGA"}="R";
$code{"CGT"}="R";
$code{"CGG"}="R";
$code{"CGC"}="R";

$code{"CCA"}="P";
$code{"CCT"}="P";
$code{"CCG"}="P";
$code{"CCC"}="P";



$mutations{"AAA"}="TAA	GAA	CAA	ATA	AGA	ACA	AAT	AAG	AAC";
$mutations{"AAT"}="TAT	GAT	CAT	ATT	AGT	ACT	AAA	AAG	AAC";
$mutations{"AAG"}="TAG	GAG	CAG	ATG	AGG	ACG	AAA	AAT	AAC";
$mutations{"AAC"}="TAC	GAC	CAC	ATC	AGC	ACC	AAA	AAT	AAG";

$mutations{"ATA"}="TTA	GTA	CTA	AAA	AGA	ACA	ATT	ATG	ATC";
$mutations{"ATT"}="TTT	GTT	CTT	AAT	AGT	ACT	ATA	ATG	ATC";
$mutations{"ATG"}="TTG	GTG	CTG	AAG	AGG	ACG	ATA	ATT	ATC";
$mutations{"ATC"}="TTC	GTC	CTC	AAC	AGC	ACC	ATA	ATT	ATG";

$mutations{"AGA"}="TGA	GGA	CGA	AAA	ATA	ACA	AGT	AGG	AGC";
$mutations{"AGT"}="TGT	GGT	CGT	AAT	ATT	ACT	AGA	AGG	AGC";
$mutations{"AGG"}="TGG	GGG	CGG	AAG	ATG	ACG	AGA	AGT	AGC";
$mutations{"AGC"}="TGC	GGC	CGC	AAC	ATC	ACC	AGA	AGT	AGG";

$mutations{"ACA"}="TCA	GCA	CCA	AAA	ATA	AGA	ACT	ACG	ACC";
$mutations{"ACT"}="TCT	GCT	CCT	AAT	ATT	AGT	ACA	ACG	ACC";
$mutations{"ACG"}="TCG	GCG	CCG	AAG	ATG	AGG	ACA	ACT	ACC";
$mutations{"ACC"}="TCC	GCC	CCC	AAC	ATC	AGC	ACA	ACT	ACG";

$mutations{"TAA"}="AAA	GAA	CAA	TTA	TGA	TCA	TAT	TAG	TAC";
$mutations{"TAT"}="AAT	GAT	CAT	TTT	TGT	TCT	TAA	TAG	TAC";
$mutations{"TAG"}="AAG	GAG	CAG	TTG	TGG	TCG	TAA	TAT	TAC";
$mutations{"TAC"}="AAC	GAC	CAC	TTC	TGC	TCC	TAA	TAT	TAG";

$mutations{"TTA"}="ATA	GTA	CTA	TAA	TGA	TCA	TTT	TTG	TTC";
$mutations{"TTT"}="ATT	GTT	CTT	TAT	TGT	TCT	TTA	TTG	TTC";
$mutations{"TTG"}="ATG	GTG	CTG	TAG	TGG	TCG	TTA	TTT	TTC";
$mutations{"TTC"}="ATC	GTC	CTC	TAC	TGC	TCC	TTA	TTT	TTG";

$mutations{"TGA"}="AGA	GGA	CGA	TAA	TTA	TCA	TGT	TGG	TGC";
$mutations{"TGT"}="AGT	GGT	CGT	TAT	TTT	TCT	TGA	TGG	TGC";
$mutations{"TGG"}="AGG	GGG	CGG	TAG	TTG	TCG	TGA	TGT	TGC";
$mutations{"TGC"}="AGC	GGC	CGC	TAC	TTC	TCC	TGA	TGT	TGG";

$mutations{"TCA"}="ACA	GCA	CCA	TAA	TTA	TGA	TCT	TCG	TCC";
$mutations{"TCT"}="ACT	GCT	CCT	TAT	TTT	TGT	TCA	TCG	TCC";
$mutations{"TCG"}="ACG	GCG	CCG	TAG	TTG	TGG	TCA	TCT	TCC";
$mutations{"TCC"}="ACC	GCC	CCC	TAC	TTC	TGC	TCA	TCT	TCG";

$mutations{"GAA"}="AAA	TAA	CAA	GTA	GGA	GCA	GAT	GAG	GAC";
$mutations{"GAT"}="AAT	TAT	CAT	GTT	GGT	GCT	GAA	GAG	GAC";
$mutations{"GAG"}="AAG	TAG	CAG	GTG	GGG	GCG	GAA	GAT	GAC";
$mutations{"GAC"}="AAC	TAC	CAC	GTC	GGC	GCC	GAA	GAT	GAG";

$mutations{"GTA"}="ATA	TTA	CTA	GAA	GGA	GCA	GTT	GTG	GTC";
$mutations{"GTT"}="ATT	TTT	CTT	GAT	GGT	GCT	GTA	GTG	GTC";
$mutations{"GTG"}="ATG	TTG	CTG	GAG	GGG	GCG	GTA	GTT	GTC";
$mutations{"GTC"}="ATC	TTC	CTC	GAC	GGC	GCC	GTA	GTT	GTG";

$mutations{"GGA"}="AGA	TGA	CGA	GAA	GTA	GCA	GGT	GGG	GGC";
$mutations{"GGT"}="AGT	TGT	CGT	GAT	GTT	GCT	GGA	GGG	GGC";
$mutations{"GGG"}="AGG	TGG	CGG	GAG	GTG	GCG	GGA	GGT	GGC";
$mutations{"GGC"}="AGC	TGC	CGC	GAC	GTC	GCC	GGA	GGT	GGG";

$mutations{"GCA"}="ACA	TCA	CCA	GAA	GTA	GGA	GCT	GCG	GCC";
$mutations{"GCT"}="ACT	TCT	CCT	GAT	GTT	GGT	GCA	GCG	GCC";
$mutations{"GCG"}="ACG	TCG	CCG	GAG	GTG	GGG	GCA	GCT	GCC";
$mutations{"GCC"}="ACC	TCC	CCC	GAC	GTC	GGC	GCA	GCT	GCG";

$mutations{"CAA"}="AAA	TAA	GAA	CTA	CGA	CCA	CAT	CAG	CAC";
$mutations{"CAT"}="AAT	TAT	GAT	CTT	CGT	CCT	CAA	CAG	CAC";
$mutations{"CAG"}="AAG	TAG	GAG	CTG	CGG	CCG	CAA	CAT	CAC";
$mutations{"CAC"}="AAC	TAC	GAC	CTC	CGC	CCC	CAA	CAT	CAG";

$mutations{"CTA"}="ATA	TTA	GTA	CAA	CGA	CCA	CTT	CTG	CTC";
$mutations{"CTT"}="ATT	TTT	GTT	CAT	CGT	CCT	CTA	CTG	CTC";
$mutations{"CTG"}="ATG	TTG	GTG	CAG	CGG	CCG	CTA	CTT	CTC";
$mutations{"CTC"}="ATC	TTC	GTC	CAC	CGC	CCC	CTA	CTT	CTG";

$mutations{"CGA"}="AGA	TGA	GGA	CAA	CTA	CCA	CGT	CGG	CGC";
$mutations{"CGT"}="AGT	TGT	GGT	CAT	CTT	CCT	CGA	CGG	CGC";
$mutations{"CGG"}="AGG	TGG	GGG	CAG	CTG	CCG	CGA	CGT	CGC";
$mutations{"CGC"}="AGC	TGC	GGC	CAC	CTC	CCC	CGA	CGT	CGG";

$mutations{"CCA"}="ACA	TCA	GCA	CAA	CTA	CGA	CCT	CCG	CCC";
$mutations{"CCT"}="ACT	TCT	GCT	CAT	CTT	CGT	CCA	CCG	CCC";
$mutations{"CCG"}="ACG	TCG	GCG	CAG	CTG	CGG	CCA	CCT	CCC";
$mutations{"CCC"}="ACC	TCC	GCC	CAC	CTC	CGC	CCA	CCT	CCG";

$reverse{"K"}="AAA	AAG";
$reverse{"N"}="AAT	AAC";
$reverse{"M"}="ATG";
$reverse{"I"}="ATA	ATT	ATC";
$reverse{"R"}="AGA	AGG	CGA	CGT	CGG	CGC";
$reverse{"S"}="AGT	AGC	TCA	TCT	TCG	TCC";
$reverse{"T"}="ACA	ACT	ACG	ACC";
$reverse{"Y"}="TAT	TAC";
$reverse{"L"}="TTA	TTG";
$reverse{"F"}="TTT	TTC";
$reverse{"C"}="TGT	TGC";
$reverse{"W"}="TGG";
$reverse{"E"}="GAA	GAG";
$reverse{"D"}="GAT	GAC";
$reverse{"V"}="GTA	GTT	GTG	GTC";
$reverse{"G"}="GGA	GGT	GGG	GGC";
$reverse{"A"}="GCA	GCT	GCG	GCC";
$reverse{"Q"}="CAA	CAG";
$reverse{"H"}="CAT	CAC";
$reverse{"P"}="CCA	CCT	CCG	CCC";


close(INFILE);

open (INFILE, "ok_mutations_in_predictions_av_pos_byalignment.txt") or die "Kozyol!";

while(<INFILE>)
{$line=$_;

($start)=($line =~ m/(\w)/);
($pos)  =($line =~ m/\w(\d+)/);
($stop) =($line =~ m/\w\d+(\w)/);

$starts=$reverse{$start};
$stops =$reverse{$stop};

while($starts)
{
($start)=($starts =~ m/(\S+)/);
($starts)=($starts =~ m/\S+\s+(.*)/);

$temp=$stops;
  while($temp)
  {($codon)=($temp =~ m/(\S+)/);
   ($temp) =($temp =~ m/\S+\s+(.*)/);
   
   if ($mutations{$start} =~ m/$codon/)
    {$substitutions{"$start$pos$codon"}=1;$substitutions{"$codon$pos$start"}=1;
     $neutral{"$start$pos$start"}=1;$neutral{"$codon$pos$codon"}=1;}
  }
}

}

close(INFILE);


open (OUTFILE, ">ok_mutations_in_predictions_av_pos_byalignment_substitutions.txt") or die "Kozyol!";

foreach $codon (keys %substitutions)
{

print OUTFILE "$codon\n";

}
close(OUTFILE);


open (OUTFILE, ">ok_mutations_in_predictions_av_pos_byalignment_neutral.txt") or die "Kozyol!";

foreach $codon (keys %neutral)
{

print OUTFILE "$codon\n";

}
close(OUTFILE);


