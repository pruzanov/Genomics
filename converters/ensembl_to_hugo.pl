#!/usr/bin/perl -w

# Script finds ensembl ids and change them to HUGO ids
use FindBin qw($Bin);
use strict;

my $file = shift @ARGV;
my $mapfile = shift @ARGV;

my $USAGE = "./ensembl_to_hugo.pl [FILE] [MAPPINGS] > [OUT]";
$mapfile ||="$Bin/id_mapping.GRCh38.92";

die "Need file to work with, $USAGE" if (!$file || ! -e $file);

my %ids = ();

open(MAP,"<$mapfile") or die "Cannot open mappings file for reading";

while (<MAP>) {
  chomp;
  my @tmp = /\t/ ? split("\t") : split(" ");
  if ($ids{$tmp[0]}) { print STDERR "Duplicate HUGO id for gene [$tmp[0]]!\n"; } else { $ids{$tmp[0]} = $tmp[1]; }
}
close MAP;


open(FILE, "<$file") or die "Couldn't open file for reading";
while (<FILE>) {
  chomp;
  my @tmp  = split("\t");
  if ($ids{$tmp[0]}) {
     $tmp[0] = $ids{$tmp[0]};
  }

  print join ("\t", @tmp)."\n";
}
close FILE;


