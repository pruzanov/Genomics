#!/usr/bin/perl -w

=head2 BWA launcher
 
 align Single-end data using BWA (Paired is handled with workflow)

 sampleR_file.R1.fasq.gz	sampleR_file.R2.fasq.gz		@RG\tID:210..-GATCAG_1\tLB:LIB01\tPL:ILLUMINA\tPU:210..-GATCAG_1\tSM:LIB01	outdir
 sampleT_file.R1.fastq.gz	sampleT_file.R2.fastq.gz	@RG\tID:210..-GATCAG_1\tLB:LIB02\tPL:ILLUMINA\tPU:210..-GATCAG_1\tSM:LIB02	outdir
 
 bwa_launcher.pl --list [mylist.txt] --bwa [path to bwa-mem2 executable] --bwa-param [optional parameters] --threads 8 --ref-base rererence_base

 THIS HAS BEEN SWITCHED TO bwa-mem2!

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $bwa, $ref, $threads, $samtools, $params, $prefix);
my $USAGE = "bwa_launcher.pl --list [mylist.txt] --bwa [path to BWA-MEM2 executable] --threads [threads to use] --ref-base [/some/path/hg19_indecbase] --samtools [samtools] --prefix [prefix]\n";
my $result = GetOptions("list=s"          => \$list,
                        "ref-base=s"      => \$ref,
                        "threads=i"       => \$threads,
                        "samtools=s"      => \$samtools,
                        "bwa-param=s"     => \$params,
                        "prefix=s"        => \$prefix,
                        "bwa=s"           => \$bwa);

$threads ||=6;
$params  ||="";
if (!$list) { die $USAGE; }
if (!$ref) {
  print STDERR "It is recommended that user specifies the reference file with bwa indexes in the same dir, default is hg19_random\n";
  $ref = "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.12/hg19_random.fa";
}
if (!$bwa) {
  print STDERR "Default bwa will be used\n";
  $bwa = "/u/pruzanov/Downloads/bwa-0.7.13/bwa";
}
if (!$samtools) {
 print STDERR "Default samtools will be used\n";
 $samtools = "/u/pruzanov/bin/samtools";
}


my $SGEscript = <<'BWA_SCRIPT';
#!/usr/bin/perl
# : Align data with BWA
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=48G
#$ -l h_rt=48:00:00
#$ -pe smp THREADS_TAG
#$ -o /dev/null
#$ -e BwaAlign.e
#$ -S /usr/bin/perl
#$ -hold_jid FILT_3885390
#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $bwa             = "BWA_TAG";
my $threads         = "THREADS_TAG";
my $refbase         = "REF_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $prefix          = "PREFIX_TAG";
my $params          = "PARAMS_TAG";

$prefix = "" if $prefix=~/_TAG/;
open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $fileR2, $rg, $outdir);

my $Count=0;
while (<INPUTFILE>) {
        next if /^#/;
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        # Below only lines that we want get processed: 
        if (/^#/) {
          print STDERR "Skipping Line $Line\n";
          exit;
        }

        ($fileR1, $fileR2, $rg, $outdir) = split("\t");
        
        if (!-e $fileR1 || !$fileR2 || !$rg || !$outdir) { die "Incorrect input parameters! Check if files exist"; }
        if (!-d $outdir) { `mkdir -p $outdir`; }
}
close INPUTFILE;

# Construct prefix
$prefix ||= basename($fileR1);
$prefix=~s/_L000(\d)_.*/_$1/;
$prefix=~s/\.fq$//;
$prefix=~s/\..*//;


$outdir =~s!//!/!g;
if ($outdir !~ m!/$!) {$outdir.="/";}
$prefix = $outdir.$prefix;
my $sorted = $prefix.".sorted.bam";
my $sort_prefix = $outdir."SORT_".$Line;

# Support for single-end reads
if ($fileR2 eq $fileR1) {
  $fileR2 = "";
}

my $command;
 
$command = "$bwa mem -t $threads $refbase $fileR1 $fileR2 -R \'$rg\' -M  $params | $samtools view -Su - -b | $samtools sort -T $sort_prefix -o $sorted -"; 

print STDERR $command."\n";

`$command`;

BWA_SCRIPT

my $JobName = "BWA_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/THREADS_TAG/$threads/g;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/BWA_TAG/$bwa/;
$SGEscript =~s/REF_TAG/$ref/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/PREFIX_TAG/$prefix/ if $prefix;
$SGEscript =~s/PARAMS_TAG/$params/;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/bwa_script$$.pl") or die "Couldn't write to [$Bin/bwa_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/bwa_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
