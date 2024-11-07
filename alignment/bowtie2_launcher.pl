#!/usr/bin/perl -w

=head2 Bowtie2 launcher
 
 align RNA-Seq data using Bowtie2 (we expect adaptor-trimmed bams with low complexity reads removed

 sampleR_file.R1.fastq.gz	sampleR_file.R2.fastq.gz	rg_line		sample_Alignments.bam
 sampleT_file.R1.fastq.gz	sampleT_file.R2.fastq.gz

 For single-ended reads put the same file in both columns,

 sampleR_file.R1.fastq.gz       sampleR_file.R1.fastq.gz        rg_line         sample_Alignments.bam

 bowtie2_launcher.pl --list [mylist.txt] --bowtie2 [path to Bowtie2 executable] --threads 6 --k 10
                     --index-base /mt/dir/bowtie2-indexes/hg38 --adaptors [String(s) of adaptors to clip]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my %strand = ("FR" => "--fr ","RF" => "--rf ","FF" => "--ff ");

my($list, $bowtie, $index_base, $threads, $adaptors, $samtools, $k, $ref_fasta, $stranded);
my $USAGE = "bowtie2_launcher.pl --list [mylist.txt] --bowtie [path to BOWTIE2 executable] --ref-fasta[Reference fasta] --threads [threads to use] --stranded [Optional, not used by default but may be FR RF or FF] --index-base [basename Bowtie2 indexes] --samtools [path to samtools binary] --k [Optional, if passed a non-zero value multimappers will be generated]\n";
my $result = GetOptions("list=s"          => \$list,
                        "index-base=s"    => \$index_base,
                        "ref-fasta=s"     => \$ref_fasta,
                        "samtools=s"      => \$samtools,
                        "k=i"             => \$k,
                        "stranded=s"      => \$stranded,
                        "threads=i"       => \$threads,
                        "bowtie2|bowtie=s"=> \$bowtie);

if (!$list || !$bowtie || !$samtools || !$ref_fasta) { die $USAGE; }
$k ||=0;
$threads ||=6;
if (!$index_base) {
  print STDERR "It is recommended that user specifies the index directory, reverting to default (hg38)\n";
  $index_base = "/.mounts/labs/gsiprojects/gsi/reference/indexes/bowtie2/hg38/hg38_human";
}

if ($stranded && $strand{$stranded}) {
  $stranded = $strand{$stranded};
} else {
  $stranded = undef;
}

my $SGEscript = <<'BOWTIE_SCRIPT';
#!/usr/bin/perl
# : Align RNAseq data with BOWTIE2
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=36G
#$ -pe smp THREADS_TAG
#$ -o /dev/null
#$ -l h_rt=36:00:00
#$ -e BowtieAlign.e
#$ -S /usr/bin/perl
#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $bowtie          = "BOWTIE_TAG";
my $threads         = "THREADS_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $indexbase       = "IDXBASE_TAG";
my $fasta           = "FASTA_TAG";
my $k               = "K_TAG";
my $stranded        = "STRANDED_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $fileR2, $outbam, $rg);

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        # Below only lines that we want get processed: 
        if (/^#/) {
          print STDERR "Skipping Commented Line\n";
          next;
        }

        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }


        ($fileR1, $fileR2, $rg, $outbam) = split("\t");
        if (!-e $fileR1 || !-e $fileR2 ) { die "Incorrect input parameters! Check if files exist"; }
}
close INPUTFILE;

my $bamdir = dirname($outbam);
my $tmpdir = $bamdir."/tmp_$$";
if (!-d $bamdir) { `mkdir -p $bamdir`; }
if (!-d $tmpdir) { `mkdir -p $tmpdir`; }

# Construct prefix
my($rg_id,$rg_barcode,$rg_library,$rg_sample) = ("Bowtie2","NoIndex","LIBRARY","SAMPLE");
$rg_library = $1 if $rg=~/LB:(\S+)/;
$rg_barcode = $1 if $rg=~/PU:(\S+)/;
$rg_sample  = $1 if $rg=~/SM:(\S+)/;
$rg_id      = $1 if $rg=~/ID:(\S+)/;

my $command = "$bowtie";
   #$command.= " --very-sensitive-local";
   $command.= " --sensitive-local";
   $command.= " -x $indexbase ";
   $command.= " -p $threads ";
   $command.= $fileR1 eq $fileR2 ? " -U $fileR1 " : " -1 $fileR1 -2 $fileR2 ";
   if ($k > 0) { $command.=" -k $k "; }
   $command.= $stranded if $stranded ne "STRANDED_TAG";
   $command.= " --fr"; # This is for stranded ERNA samples
   $command.= " --rg-id $rg_id";
   $command.= " --rg PL:Illumina";
   $command.= " --rg PU:$rg_barcode";
   $command.= " --rg LB:$rg_library";
   $command.= " --rg SM:$rg_sample -S /dev/stdout ";
   $command.= "| samtools sort -O bam --reference $fasta -T $tmpdir -o $outbam -";
print STDERR $command."\n\n";

`$command`;

BOWTIE_SCRIPT

my $JobName = "BOWTIE2_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/THREADS_TAG/$threads/g;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/BOWTIE_TAG/$bowtie/;
$SGEscript =~s/IDXBASE_TAG/$index_base/;
$SGEscript =~s/FASTA_TAG/$ref_fasta/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/K_TAG/$k/;
$SGEscript =~s/STRANDED_TAG/$stranded/ if $stranded;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/bowtie_script$$.pl") or die "Couldn't write to [$Bin/bowtie_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V  $Bin/bowtie_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
