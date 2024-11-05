#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

# die "perl $0 <name> <dir> <ref> <circos> <bcftools> "unless @ARGV==5;

my ($name, $vcf, $fsv, $fcnv, $fbamfinal, $ref, $circos, $bcftools) = @ARGV;

my $cd = ".";
my %flag;
my $win = 1000000;

my $header;
if($ref =~ /hg19.fa$|hs37d5.fa$/){
  $header = "
  karyotype = $circos/data/karyotype/karyotype.human.hg19.txt
  chromosomes = -hsY
  ";
}elsif($ref =~ /hg38.fa$/){
  $header = "
  karyotype = $circos/data/karyotype/karyotype.human.hg38.txt
  chromosomes = -hsY
  ";
}else{
  open IN,"$ref.fai";
  open OT,">$cd/karyotype.txt";
  while(<IN>){
    chomp;
    my @a = split;
    print OT "chr - $a[0] $a[0] 0 $a[1] $a[0]\n";
  }
  close IN;
  close OT;
  $header = "
  karyotype = karyotype.txt
  ";
}

  open MAIN,">circos.$name.conf";
  print MAIN "
  $header
  <<include $circos/etc/housekeeping.conf>>
  <<include $circos/etc/colors_fonts_patterns.conf>>
  <<include $Bin/circos.image.conf>>
  <<include $Bin/circos.ideogram.conf>>
  <<include $Bin/circos.ticks.conf>>
  chromosomes_units = 1000000
  ";

  # cnv and sv
  if($ref =~ /hg38.fa$|hg19.fa$|hs37d5.fa$/){
    if(-e "$fcnv"){
      open CNV,"$fcnv";
      open TXT,">cnv.$name.txt";
      while(<CNV>){
        chomp;
        next if /Start|^>/;
        my @a = split;
        if($a[0] =~ /^chr/){
          $a[0] =~ s/^chr/hs/;
        }else{
          $a[0] = "hs$a[0]";
        }
        print TXT "$a[0]\t$a[1]\t$a[2]\t$a[6]\n";
      }
      close CNV;
      close TXT;
    }else{
      `touch cnv.$name.txt`;
    }
    open CONF,">cnv.$name.conf";
    print CONF "
    <highlight>
    file = cnv.$name.txt
    <<include $Bin/circos.cnv.conf>>
    </highlight>
    ";
    close CONF;

    if(-e "$fsv"){
      #open SV,"$ARGV[1]/../../file/$name/SV/$name.SV.vcf";
      open SV,"$fbamfinal";
      open TXT,">sv.$name.txt";
      while(<SV>){
        chomp;
        my @a = split;
        next if /^#|^EventID/;
        if($a[4] =~ /^chr/){
          $a[4] =~ s/chr/hs/;
          $a[6] =~ s/chr/hs/;
        }else{
          $a[4] = "hs$a[4]";
          $a[6] = "hs$a[6]";
        }
        my $col = ($a[10] eq "DEL") ? 1 :
                  ($a[10] eq "DUP") ? 2 :
                  ($a[10] =~ /INV/) ? 3 :
                  ($a[10] =~ /TRA|BND/) ? 4 :
                  5;
        print TXT "$a[4]\t$a[5]\t$a[5]\t$a[6]\t$a[7]\t$a[7]\tcolor=set2-5-qual-$col\n";
      }
      close SV;

      open SV,"$fsv";
      while(<SV>){
        chomp;
        my @a = split;
        next if /^#/;
        next unless /SOURCE=smoove/;
        my ($type, $end);
        if($a[0] =~ /^chr/){
          $a[0] =~ s/chr/hs/;
        }else{
          $a[0] = "hs$a[0]";
        }

        foreach my $info (split /\;/, $a[7]){
          my @info = split /\=/, $info;
          $type = $info[1] if $info[0] eq "SVTYPE";
          $end = $info[1] if $info[0] eq "END";
        }
        my $col = ($type eq "DEL") ? 1 :
                  ($type eq "DUP") ? 2 :
                  ($type =~ /INV/) ? 3 :
                  ($type =~ /TRA|BND/) ? 4 :
                  5;
        if($type eq "BND"){
          my ($chr2, $pos2);
          my @splitpos = split /\[|\]|\:/, $a[4];
          $chr2 = $splitpos[1];
          $pos2 = $splitpos[2];
          if($chr2 =~ /^chr/){ $chr2 =~ s/chr/hs/; }else{ $chr2 = "hs$chr2"; }
          print TXT "$a[0]\t$a[1]\t$a[1]\t$chr2\t$pos2\t$pos2\tcolor=set2-5-qual-$col\n";
        }else{
          print TXT "$a[0]\t$a[1]\t$a[1]\t$a[0]\t$a[1]\t$a[1]\tcolor=set2-5-qual-$col\n";
        }
      }
      close SV;
      close TXT;
    }else{
      `touch sv.$name.txt`;
    }
    open CONF,">sv.$name.conf";
    print CONF "
    <link>
    file = sv.$name.txt
    <<include $Bin/circos.sv.conf>>
    </link>
    ";
    close CONF;

    print MAIN "
    <highlights>
    <<include cnv.$name.conf>>
    </highlights>
    <links>
    <<include sv.$name.conf>>
    </links>
    ";

  }

  # snp and indel
  treatvcf( "$vcf", "snp", $name );
  treatvcf( "$vcf", "indel", $name );
  print MAIN "
  <plots>
  <<include snp.$name.conf>>
  <<include indel.$name.conf>>
  </plots>
  ";

  # circos
  #`$circos/bin/circos -conf $cd/circos.$name.conf -outputfile $name.circos.png -outputdir $cd`;

sub treatvcf{
  my ($vcf, $type, $name) = (@_);
  my %data;
  my $max = 0;
  open F,"$bcftools view -v ${type}s $vcf | grep -v ^# |";
  open TXT,">$cd/$type.$name.txt";
  while(<F>){
    chomp;
    my @a = split;
    my $key = int($a[1] / $win);
    my $chr = $a[0];
    $chr =~ s/^chr/hs/ if ($ref =~ /hg38.fa$|hg19.fa$/);
    $chr = "hs$chr" if ($ref =~ /hs37d5.fa$/);
    $data{$chr}{$key}++;
  }
  close F;
  foreach my $chr (sort keys %data){
    foreach my $pos (sort {$a <=> $b} keys %{$data{$chr}}){
      my $s = $pos * $win;
      my $e = $s + $win;
      my $r = $data{$chr}{$pos} / $win;
      $max = ($max > $r) ? $max : $r;
      print TXT "$chr\t$s\t$e\t$r\n";
    }
  }
  close TXT;
  open CONF,">$type.$name.conf";
  print CONF "
  <plot>
  file = $type.$name.txt
  max = $max
  <<include $Bin/circos.$type.conf>>
  <<include $Bin/circos.bg.conf>>
  </plot>
  ";
  close CONF;
};
