#!/usr/bin/perl -w
use strict;

# die "perl $0 <vcf> <sv> <cnv> <name>"unless @ARGV==2;

my $varianttable = $ARGV[0];
my $fsv = $ARGV[1];
my $fcnv = $ARGV[2];
my $id = $ARGV[3];
my $flag = 0;

open IN,"$varianttable";
open OT1,">$id.var.xls";
while(<IN>){
  chomp;
  my @a = split /\t/;
  next if /^CNV|^SV/;
  if(/^Sample/){
    for(my $i = 1; $i < @a; $i++){
      $flag = $i if $a[$i] eq $id;
      print OT1 "$a[0]\t$a[$flag]\n";
      last;
    }
    next;
  }
  print OT1 "$a[0]\t$a[$flag]\n";
}
close IN;

my @cnv = (0) x 2;
if(-e "$fcnv"){
  open IN,"$fcnv";
  while(<IN>){
    chomp;
    next if !/^>/;
    my @a = split;
    $a[3] eq "DEL" ? $cnv[0]++ : $cnv[1]++;
  }
  close IN;
}
print OT1 "CNV deletion\t$cnv[0]\n";
print OT1 "CNV duplication\t$cnv[1]\n";

my @sv = (0) x 5;
if(-e "$fsv"){
  open IN,"$fsv";
  while(<IN>){
    chomp;
    next if /^#/;
    my @a = split;
    my $type;
    foreach my $info (split /\;/, $a[7]){
      my @info = split /\=/, $info;
      $type = $info[1] if $info[0] eq "SVTYPE";
    }
    $type eq "DEL" ? $sv[0]++ :
    $type eq "DUP" ? $sv[1]++ :
    $type =~ /INV/ ? $sv[2]++ : 
    $type =~ /TRA|BND/ ? $sv[3]++ : 
    $sv[4]++;
  }
  close IN;
}
print OT1 "SV DEL\t$sv[0]\n";
print OT1 "SV DUP\t$sv[1]\n";
print OT1 "SV INV\t$sv[2]\n";
print OT1 "SV TRA|BND\t$sv[3]\n";

close OT1;

