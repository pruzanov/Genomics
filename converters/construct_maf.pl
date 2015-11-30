#!/usr/bin/perl -w

=head1 Construct maf files for D-ToxoG tool

 This script will take in a bam and vcf file (for the same sample)
 And construct a .maf file with the following fields:

 vcf processing (output of IntersectBed, optional filtering for PASS calls)
 1. chromosome    -
 2. start of SNV  - 
 3. end           - the same as start
 4. ref allele    -
 5. tumor allele 1 -
 6. tumor allele 2 -

 7. tumor name  (barcode) - from file metadata or TUMOR
 8. normal name (barcode) - from file metadata or NORMAL
 9. ref context - 3-21 base (odd number) - samtools faidx

 Processing bam file with samtools (having hash with all calls)
 10.  F1R2 alt allele #
 11. F2R1 alt allele #
 12. F1R2 ref allele #
 13. F2R1 ref allele #
 14. Foxog (depends on mutation type)
 15. Variant type
 16. i_picard_oxoQ value (OXIDATION_Q produced with Picard's CollectOxoGMetrics tool)

 MARGIN defines how many bases we want on each side (with 1 we get 3, with 2 we get 5 etc.)
 when extracting context for a SNP - having margin set to anything other than 1 crashes D-ToxoG 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;
use lib "$Bin/lib/";
use SNVCall;
use constant DEBUG=>0;
use constant MARGIN=>1;

=head1 USAGE

 construct_maf.pl --bam [bam file] 
                  --vcf [vcf file]
                  --picard-oxog-file
                  --samtools [samtools] 
                  --flter [i.e. PASS, to limit the number of observed SNPs]
                  --isect-bed [intersectBed] 
                  --tumor-id [optional, tumor id] 
                  --normal-id [optional, normal id] 
                  --out [optional, filename]

=cut

my @colnames = qw(Chromosome Start_position End_position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode Matched_Norm_Sample_Barcode ref_context i_t_ALT_F1R2 i_t_ALT_F2R1 i_t_REF_F1R2 i_t_REF_F2R1 i_t_Foxog Variant_Type i_picard_oxoQ);
my ($bam,$vcf,$poxog,$samtools,$isectbed,$tumor_id,$normal_id,$ref_fasta,$outfile,$filter);
my $result = GetOptions ('bam=s'      => \$bam,
                         'vcf=s'      => \$vcf,
                         'picard-oxog-file=s' => \$poxog,
                         'filter=s'   => \$filter,
                         'samtools=s' => \$samtools,
                         'ref-fasta=s'=> \$ref_fasta,
                         'isect-bed=s'=> \$isectbed,
                         'tumor-id=s' => \$tumor_id,
                         'normal-id=s'=> \$normal_id,
                         'out=s'      => \$outfile);

my $USAGE = "construct_maf.pl --bam [bam file] --vcf [vcf file] --picard-oxog-file [] --samtools [samtools] --ref-fasta [i.e. hg19.fa, for context extraction] --isect-bed [intersectBed] --tumor-id [optional, tumor id] --normal-id [optional, normal id] --filter [optional, i.e. PASS to limit the list of calls] --out [optional, filename]\n";
if (!$bam || !$vcf || !$poxog || !$samtools || !$isectbed || !$ref_fasta) {die $USAGE;}

#===============MAIN BLOCK===============//
my %read2snp = ();
my %poxo_q   = %{&read_poxo($poxog)};

# read information from genotype file
my %vars = %{&init_varhash};

# read information from alignment file
&read_alignments();

open(OUT,">$outfile") or die "Cannot write to output file [$outfile]";
print OUT join("\t",@colnames)."\n";
map{print OUT $vars{$_}->to_string()."\n"} (keys %vars);

close OUT;
#===========MAIN BLOCK ENDS=============//

=head2

 Read oxoG Q score from the file produced with picard-tools

 Scores may not be available for all context sequences (triplets)
 in such cases 0.00 will be used 

=cut

sub read_poxo {

 my $file = shift;
 my %scores = ();
 open(P,"<$file") or die "Couldn't read Picard metric file from [$poxog]";
 while(<P>) {
   next if (/CONTEXT/ || /^\s*$/ || /^#/); # skip header, we need fields 3 and 12 only
   chomp;
   my @temp = split("\t");
   if ($temp[2] =~/[A-T]{3}/i) {
     $scores{uc($temp[2])} = sprintf "%.2f", $temp[11];
   }
 }
 close P;
 return \%scores;
}

=head2 Subroutine for getting allele info

 Here, we use intersectBed from bedtools to get all alleles
 that are overlapping our bam reads 

 if we don't have a vcf entry yet, we initialize it first
 calling SNVCall->new(chr, start, ref, alt)

 and/or adding read id to the hash that contains only unique ids 
 (read + mate id), as in 
 HWI-ST801:93:C12P2ACXX:2:2102:15104:40299/2 (note /2 at the end)

=cut

sub init_varhash {

 my %varhash = ();
 my $vcf_pipe = "$isectbed -abam $bam -b $vcf -wb -bed | ";
 open(VCF, "$vcf_pipe") or die "Couldn't create a vcf pipe [$vcf_pipe]";
 while(<VCF>) {
  if (/^(chr.+)(chr.*)$/) {
   my @baminfo = split("\t", $1);
   my @vcfinfo = split("\t", $2);
   my $readid = $baminfo[3];
   my $snp_id = join(":",($vcfinfo[0], $vcfinfo[1]));
   
   if ($filter && $vcfinfo[6] ne $filter) {next;}
   $read2snp{$readid} = $snp_id;

   if (!$varhash{$snp_id} && $vcfinfo[3]=~/^[A-T]{1}$/ && $vcfinfo[4]=~/^[A-T]{1}(\,.*)*$/)  { 
       $varhash{$snp_id} = SNVCall->new(@vcfinfo[0..1],@vcfinfo[3..4]);
       $varhash{$snp_id}->set_tumor_id($tumor_id);
       $varhash{$snp_id}->set_normal_id($normal_id);
       $varhash{$snp_id}->set_context(&get_context(@vcfinfo[0..1],MARGIN));

       # This will work only if MARGIN is set to 1
       if (my $poxo = $poxo_q{$varhash{$snp_id}->get_context()}) {
          $varhash{$snp_id}->set_poxoq($poxo);
       }
       
   } else {
     next;
   }
  
   $varhash{$snp_id}->add_read($readid);
  }
 }
 print STDERR "Got ".scalar(keys(%varhash))." variants from intersect analysis\n";

 close VCF;

 if (DEBUG) {
  print STDERR "Press any key to continue...\n";
  my $answer = <STDIN>;
 }
 return \%varhash;
}

=head2 Subroutine for bam parsing

 This part deals with alignments and the goal here is
 to get all info for bam reads overlapping our SNPs

=cut

sub read_alignments {

 my $bam_pipe = "$isectbed -abam $bam -b $vcf -wa | $samtools view - | ";
 open(BAM, "$bam_pipe") or die "Couldn't create a bam pipe";
 my $count = 0;
 while(<BAM>) {
  my @baminfo = split("\t");
  my($read,$tag,$chr,$start,$cig,$seq) = (@baminfo[0..3],$baminfo[5],$baminfo[9]);
  my $readid = $read;

  if ($tag & 64) {
    $readid.="/1";
  } elsif($tag & 128) {
    $readid.="/2";
  }

  if (!$read2snp{$readid}) {
    print STDERR "Read [$readid] could not be identified!\n" if DEBUG;
    next;
  }

  $vars{$read2snp{$readid}}->update_read($readid,$tag,$chr,$start,$cig,$seq);
  $count++;
 }
 
 close BAM;

 print STDERR $count." reads updated\n" if DEBUG;
 
 if (DEBUG) {
  print STDERR "Press any key to continue...\n";
  my $answer = <STDIN>;
 }

}

=head2 Getting nucleotide context for a SNP
 
 Borrowed this from another script
 
 We get surronding n reference bases on each side of our SNP
 using samtools

 samtools faidx hg19 chr22:22000-22002
 
=cut

sub get_context {

 my($chr,$start,$margin) = @_;
 my $site = join(":", ($chr, join("-", ($start - $margin, $start + $margin))));
 print STDERR "Checking $site\n" if DEBUG;

 my $bases =  `$samtools faidx $ref_fasta $site | grep -v $chr`;
 chomp($bases);
 my $contsize = 1+MARGIN*2;
 print STDERR "Problem with context [$bases] requested for $site\n" if $bases !~/[A-T]{$contsize}/i;
 return uc($bases);
}




