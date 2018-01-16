#!/usr/bin/perl -w

=head2 Bowtie2 launcher
 
 align RNA-Seq data using Bowtie2 (we expect adaptor-trimmed bams with low complexity reads removed

 sampleR_file.R1.fastq.gz	sampleR_file.R2.fastq.gz	sampleR_Alignments.bam	RGline[Optional]
 sampleT_file.R1.fastq.gz	sampleT_file.R2.fastq.gz	sampleT_Alignments.bam  ...

 bowtie2_launcher.pl --list [mylist.txt] --bowtie2 [path to Bowtie2 executable] --threads 6 
                     --index-base /mt/dir/bowtie2-indexes/hg38 --adaptors [String(s) of adaptors to clip]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $bowtie, $index_base, $threads, $adaptors, $samtools);
my $USAGE = "bowtie2_launcher.pl --list [mylist.txt] --bowtie [path to BOWTIE2 executable] --threads [threads to use] --index-base [basename Bowtie2 indexes] --samtools [path to samtools binary]\n";
my $result = GetOptions("list=s"          => \$list,
                        "index-base=s"    => \$index_base,
                        "samtools=s"      => \$samtools,
                        "threads=i"       => \$threads,
                        "bowtie2|bowtie=s"=> \$bowtie);

if (!$list || !$bowtie || !$samtools) { die $USAGE; }
$threads ||=6;
if (!$index_base) {
  print STDERR "It is recommended that user specifies the index directory, reverting to default\n";
  $index_base = "/.mounts/labs/prod/scratch/pruzanov/hg38_indexes/bowtie2/hg38_genome";
}

my $SGEscript = <<'BOWTIE_SCRIPT';
#!/usr/bin/perl
# : Align RNAseq data with BOWTIE2
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=14G
#$ -pe smp THREADS_TAG
#$ -o /dev/null
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

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $fileR2, $outbam);

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        # Below only lines that we want get processed: 
        if (/^#/) {
          print STDERR "Skipping Line $Line\n";
          exit;
        }

        ($fileR1, $fileR2, $outbam, $rg) = split("\t");
        if (!-e $fileR1 || !-e $fileR2 ) { die "Incorrect input parameters! Check if files exist"; }
}
close INPUTFILE;

my $command = "set -o pipefail && $bowtie";
   $command.= " --fast-local";
   $command.= " -x $indexbase ";
   $command.= " -p $threads ";
   $command.= " -1 $fileR1 -2 $fileR2 ";
   if ($rg) {
    my @rgs = split(";",$rg);
    foreach my $rg (@rgs) {
      $rg =~ /^ID\:/ ? $command.="--rg-id $rg " : $command.="--rg $rg ";
    }
   }
   $command.= "-S /dev/stdout ";
   $command.= "| $samtools view -Sbh - > ".$outbam;
print STDERR $command."\n\n";

`$command`;

BOWTIE_SCRIPT

my $JobName = "BOWTIE2_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/THREADS_TAG/$threads/g;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/BOWTIE_TAG/$bowtie/;
$SGEscript =~s/IDXBASE_TAG/$index_base/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/bowtie_script$$.pl") or die "Couldn't write to [$Bin/bowtie_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V  $Bin/bowtie_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
