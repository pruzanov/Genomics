#!/usr/bin/perl -w

#===========================================================
# Script for changing refs in sam (bam) file header and body
#===========================================================

=head2 convert.contigs.pl

 convert.contigs.pl --input [bam or sam] --output [out bam] --ref [reference assembly, need to have fai] --samtools [default is system-wide samtools]
 
 script will fix header and field 3 in all lines
  
=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use constant DEBUG=>0;

my $USAGE = "convert.contigs.pl --input [bam or sam] --ref [reference assembly file] --output [out bam] --samtools [default is system-wide samtools]\n";

# refs point to correct chromosome ids
my %refs = ("23"    => "X",
            "24"    => "Y",
            "25"    => "M",
            "chr23" => "chrX",
            "chr24" => "chrY",
            "chr25" => "chrM");

my($in, $out, $samtools, $ref, $strict);
my $results = GetOptions('samtools=s'  => \$samtools,
                         'input=s'     => \$in,
                         'output=s'    => \$out,
                         'ref|fasta=s' => \$ref,
                         'strict'      => \$strict,
                         'out=s'       => \$out);

if(!$in || ! -e $in || !$ref) {die "We need an input file and reference, $USAGE";}
$samtools ||= "samtools";
$out      ||= $in.".conv.bam";
$out =~s/\.bam\.conv\./\.conv\./; # just in case
my $temp_sam   = $out;
my $chr_prefix = 1; # assume that we must have chr prefix 
$temp_sam.="_temp$$";

my $contigs = &read_contigs($ref);

my $file_fork = $in =~ /bam$/ ? "$samtools view -h $in | " : "$samtools view -hS $in | ";
open(IN, $file_fork) or die "Couldn't open file for reading";
open(OUT,">$temp_sam") or die "Couldn't write to file [$temp_sam]";


print STDERR "Converting starts\n" if DEBUG;
my $count = 0; # modified line's counter
my @headerlines = ();

LINE:
while(my $line = <IN>) {
 if ($line=~/^\@/) {
   if ($line=~/^\@SQ\tSN:(\S+)/ && $refs{$1}) {
    print STDERR "BEFORE: $line" if DEBUG;
    my $r = $refs{$1};
    $line=~s/SN:\S+/SN:$r/;
    print STDERR "AFTER: $line" if DEBUG;
   }
   if ($chr_prefix) {
     if ($line!~/SN:chr/) {$line=~s/SN:/SN:chr/;$count++;}
   } elsif ($line=~/SN:chr/) {
     $line=~s/SN:chr/SN:/;
   }
   
   # Check if we have this contig
   if ($line=~/SN:(\S+)/ && !$contigs->{$1}) {
     die "Unknown contig $1 detected, aborting";
   }
   push @headerlines, $line;
 } else {
   if (@headerlines > 1) {
     map {print OUT $_} @headerlines;
     @headerlines = ();
   }
   my @temp = split("\t",$line);

   # Append or remove chr if needed
   if ($chr_prefix) {
     if ($temp[2]!~/^chr/) {
       $temp[2] = "chr$temp[2]";
       $count++;
     }
   } else {
     $count++ if $temp[2]=~/^chr/;
     $temp[2]=~s/^chr//;
   }
 
   if ($refs{$temp[2]}) {
     $temp[2] = $refs{$temp[2]};
     $count++;
   }
   
   # Skip if we have this contig
   if (!$contigs->{$temp[2]}) {
       print STDERR "Skipping a line since contig [$temp[2]] is unknown\n";
       next LINE;
   }
   
   print OUT join("\t", @temp);
 }
}
close IN;
close OUT;

print STDERR "Converting to bam format...\n" if DEBUG;
`$samtools view -Sh $temp_sam -b > $out && rm $temp_sam`;
print STDERR "Converting finished, $count lines modified\n" if DEBUG;


=head2 read_contigs

 Function for reading contig information from a fai file
 
 User may supply a fata file, this function will try to find fai
 by using the path to fasta and assume that .fai is in the same 
 directory. 

=cut

sub read_contigs {

 my $file = shift @_;
 my $path = $file;
 my $refs = {};

 if ($file!~/.fai/) {
   $path = $file.".fai";
   if (! -e $path) {
     # make another attempt
     $path = $file;
     $path=~s/(.*)\./$1/;
     $path.=".fai";
     -e $path or die "Index file for reference File [$file] could not be found, make sure the reference is indexed with samtools faidx";
   }
 }

 open(REF, "<$path") or die "Couldn't read from reference file [$path]";
 while(<REF>) {
   my @temp = split("\t");
   $refs->{$temp[0]} = $temp[1];
   if ($temp[0]!~/^chr/) {$chr_prefix = 0;}
 }
 close REF;

 return $refs;
}

