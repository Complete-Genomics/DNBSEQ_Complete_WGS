#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

die "perl $0 <name> <hapblock> <result.dir> <ref> <R>"unless @ARGV==5;

my $name = $ARGV[0];
my $ref = $ARGV[3];
my $pd = $ARGV[2];
my $Rscript = "Rscript";
my (%flag, %band);

# build genome txt
open IN,"$ref.fai";
open OT,">$pd/karyotype.$name.genome.txt";
print OT "chr\tstart\tend\n";
while(<IN>){
  chomp;
  my @a = split;
  next if $a[0] =~ /^GL|NC|hs37d5|\_|MT|chrM|\-|EBV/;
  print OT "$a[0]\t1\t$a[1]\n";
}
close IN;
close OT;

# build band txt
  my (%block, $chr, $start, $end) = () x 4;
  $chr = "NA";
  open BLC,"$ARGV[1]";
  while(<BLC>){
    chomp;
    next if /^\*/;
    my @a = split;

    if(/^BLOCK/){
      if($chr ne "NA"){
        push @{$block{$chr}}, "$start:$end";
      }

      my $line = <BLC>;chomp $line;
      next if $line =~ /^\*/;
      my @b = split /\t/, $line;
      $chr = $b[3];
      $start = $b[4];
      $end = $b[4];
    }else{
      $end = $a[4];
    }
  }
  close BLC;
  push @{$block{$chr}}, "$start:$end";

  open FAI,"$pd/karyotype.$name.genome.txt";<FAI>;
  open OT,">$pd/karyotype.$name.band.txt";
  print OT "chr\tstart\tend\tname\tgieStain\n";
  while(<FAI>){
    chomp;
    my @a = split;
    my $pos = 1;
    my ($bn, $bd) = (1) x 2;
    foreach (@{$block{$a[0]}}){
      my ($s, $e) = (split /\:/)[0, 1];
      if($pos < $s){
        print OT "$a[0]\t$pos\t$s\t$a[0]\_$bn\tgneg\n";
        $bn++;
        $pos = $e;
      }elsif($pos == $s){
        $pos = $e;
      }elsif($pos > $e){
        next;
      }
      my $col = ($bd % 2 == 0) ? "stalk" : "gpos25";
      print OT "$a[0]\t$s\t$e\tband\_$bd\t$col\n";
      $bd++;
    }
    if($pos < $a[1]){
      print OT "$a[0]\t$pos\t$a[1]\t$a[0]\_$bn\tgneg\n";
      $bn++;
    }
  }
  close FAI;
  close OT;

# run R
`$Rscript $Bin/haplotype.karyotype.R \\
  $pd/karyotype.$name.genome.txt \\
  $pd/karyotype.$name.band.txt \\
  $pd/$name.haplotype.pdf \\
  $Bin/../../lib`;


