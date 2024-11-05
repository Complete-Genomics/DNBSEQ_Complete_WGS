#!/usr/bin/perl -w
use strict;

die "perl $0 <sample> <hapdir> <vcfdir> <statdir>"unless @ARGV==4;

my (%name, %name_flag);
`mkdir -p $ARGV[1]`;
my $list = "chr\tswitch rate\tmismatch rate\tflat rate\tmissing rate\tphased count\tAN50\tN50\tmax block snp frac\tphasing rate";

process($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

sub process{
  my ($name, $inputdir, $vcfdir, $outputdir) = (@_);

  my %stat = ();
  my @chr = ();
  my $chr;
  open RESULT,"$inputdir/$name.hapcut_stat.txt";
  while(<RESULT>){
    chomp;
    next unless /\S+/;
    next if /compare L0 with giab/;

    if(!/\:/ && !/combine all chrs/){
      $chr = $_;
      push @chr, $chr;
      $stat{$chr}{"chr"} = $chr;
      next;
    }
    elsif(/combine all chrs/){
      $chr = "chrAll";
      push @chr, $chr;
      $stat{$chr}{"chr"} = $chr;
      next;
    }
    my @b = split /\:/;
    $b[1] =~ s/^\s+//;
    $stat{$chr}{$b[0]} = $b[1];
  }
  close RESULT;

  foreach my $chr (@chr){
    next if $chr eq "chrAll";
    open FILE,"$vcfdir/$name.$chr.hetsnp.vcf";
    while(<FILE>){
      chomp;
      next if /^#/;
      my @b = split /\t/;
      
      $stat{$chr}{"het"} += 1;
      $stat{"chrAll"}{"het"} += 1;
    }
    close FILE;
    $stat{$chr}{"phasing rate"} = (defined $stat{$chr}{"het"} > 0) ? $stat{$chr}{"phased count"} / $stat{$chr}{"het"} : 0;
  }

  $chr = "chrAll";
  $stat{$chr}{"phasing rate"} = (defined $stat{$chr}{"het"} > 0) ? $stat{$chr}{"phased count"} / $stat{$chr}{"het"} : 0;

  open OT,">$outputdir/$name.haplotype.xls";
  print OT "$list\n";
  foreach my $chr (@chr){
    print OT "$chr";
    foreach my $key(split /\t/, $list){
      next if $key eq "chr";
      $stat{$chr}{$key} = (defined $stat{$chr}{$key}) ? $stat{$chr}{$key} : 0;
      print OT "\t$stat{$chr}{$key}";
    }
    print OT "\n";
  }
  close OT;

};
