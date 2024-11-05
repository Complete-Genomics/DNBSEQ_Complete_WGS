#!/usr/bin/perl -w
use strict;

# die "perl $0 <indir> <out> <sample>"unless @ARGV==3;

my ($flgstat, $stats, $info1, $info2, $insertsize, $fout, $sample) = @ARGV;
my %stat;
#my $aligntable = "Sample name,Mapped reads,Mapped bases(bp),Mapping rate,Paired mapping rate,Mismatch bases(bp),Mismatch rate,Duplicate reads,Duplicate rate,Dup depth,Average sequencing depth,Coverage,Coverage at least 4X,Coverage at least 10X,Coverage at least 20X,Mean insert size";
my $aligntable = "Sample name,Mapped reads,Mapped bases(bp),Mapping rate,Paired mapping rate,Mismatch bases(bp),Mismatch rate,Duplicate reads,Duplicate rate,Average sequencing depth,Coverage,Coverage at least 4X,Coverage at least 10X,Coverage at least 20X,Mean insert size";

$stat{$sample}{"Sample name"} = $sample;
my @chimerical = ();
open IN3,"$flgstat";
while(<IN3>){
	chomp;
	my @a = split /\s+/;
	$chimerical[0] = $a[0] if /total/;
	$chimerical[1] = $a[0] if /mapped \(/;
	$chimerical[2] = $a[0] if /properly paired/;
	$chimerical[3] = $a[0] if /singletons/;
  if(/properly paired \((\S+) : /){
    $stat{$sample}{"Paired mapping rate"} = $1;
  }
}
close IN3;
open IN4,"$stats";
while(<IN4>){
  chomp;
  next if /^#/;
  if(/^SN/){ 
    my @a = split /\t/;
    $stat{$sample}{"Clean reads"} = $a[2] * 2 if $a[1] eq "1st fragments:";
    $stat{$sample}{"Mapped reads"} = $a[2] if $a[1] eq "reads mapped:";
    $stat{$sample}{"Duplicate reads"} = $a[2] if $a[1] eq "reads duplicated:";
  }
}
close IN4;

if($chimerical[0] > 0){
	$stat{$sample}{"Mapping rate"} = int( $stat{$sample}{"Mapped reads"} / $stat{$sample}{"Clean reads"} * 10000 + 0.5) / 100;$stat{$sample}{"Mapping rate"} .= "%";
	$stat{$sample}{"Duplicate rate"} = int($stat{$sample}{"Duplicate reads"} / $stat{$sample}{"Clean reads"} * 10000 + 0.5) / 100;$stat{$sample}{"Duplicate rate"} .= "%";
	$stat{$sample}{"Chimerical rate"} = int(($chimerical[1] - $chimerical[2] - $chimerical[3]) / $chimerical[0] * 10000 + 0.5) / 100;$stat{$sample}{"Chimerical rate"} .= "%";
}
else{
  $stat{$sample}{"Mapping rate"} = 0;
  $stat{$sample}{"Duplicate rate"} = 0;
  $stat{$sample}{"Chimerical rate"} = 0;
}

open IN1,"$info1";
while(<IN1>){
	chomp;
	my @a = split /\t/;
	$stat{$sample}{$a[0]} = $a[1];
}
close IN1;
$stat{$sample}{"Clean bases(bp)"} = $stat{$sample}{"Clean reads"} * $stat{$sample}{"Readlength"};
$stat{$sample}{"Mapped bases(bp)"} = $stat{$sample}{"Mapped reads"} * $stat{$sample}{"Readlength"};

open IN2,"$info2";
while(<IN2>){
	chomp;
	my @a = split /\t/;
	$stat{$sample}{$a[0]} = $a[1];
}
close IN2;

# insert size
print STDERR "processing $sample insert size ......\n";
open INSERT,"$insertsize";
while(<INSERT>){
  chomp;
  if(/^MEDIAN_INSERT_SIZE/){
    my @id = split /\s+/;
    my $id;
    for(my $k = 0; $k < @id; $k++){
      $id = $k if $id[$k] eq "MEAN_INSERT_SIZE";
    }
    my $value = <INSERT>;chomp $value;
    $stat{$sample}{"Mean insert size"} = int((split /\s+/, $value)[$id] * 100 + 0.5) / 100;
  }
}
close INSERT;

my $rate = $stat{$sample}{"Duplicate rate"};
$rate =~ s/%$//;
$rate /= 100;
#$stat{$sample}{"Dup depth"} = $stat{$sample}{"Average sequencing depth"} * (1 - $rate);
#remove dup when calculating ave depth in bam2depth, so dup depth is ave depth,and ave depth is all data * mapping rate, 20200716_rjh

open OT1,">$fout"; 
foreach my $key (split /\,/, $aligntable){
	chomp $key;
	print OT1 "$key\t$stat{$sample}{$key}\n";	
}
close OT1;

