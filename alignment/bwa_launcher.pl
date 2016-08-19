#!/usr/bin/perl -w

=head2 BWA launcher
 
 align Single-end data using BWA (Paired is handled with workflow)

 sampleR_file.R1.fasq.gz	outdir
 sampleT_file.R1.fastq.gz	outdir

 bwa_launcher.pl --list [mylist.txt] --bwa [path to bwa executable] --threads 8 --ref-fasta fererence.fa [the same dir needs to have indexes too]

=cut
 
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my($list, $bwa, $ref, $threads, $samtools);
my $USAGE = "bwa_launcher.pl --list [mylist.txt] --bwa [path to BWA executable] --threads [threads to use] --ref-fasta [hg19_random.fa] --samtools [samtools]\n";
my $result = GetOptions("list=s"          => \$list,
                        "ref-fasta=s"     => \$ref,
                        "threads=i"       => \$threads,
                        "samtools=s"      => \$samtools,
                        "bwa=s"           => \$bwa);

$threads ||=6;
if (!$list) { die $USAGE; }
# Some Defaults, TODO: user should change these
if (!$ref) {
  print STDERR "It is recommended that user specifies the reference file with bwa indexes in the same dir\n";
  $ref = "$HOME/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.12/hg19_random.fa";
}
if (!$bwa) {
  print STDERR "Default bwa will be used\n";
  $bwa = "$HOME/Downloads/bwa-0.7.13/bwa";
}
if (!$samtools) {
 print STDERR "Default samtools will be used\n";
 $samtools = "$HOME/bin/samtools";
}


my $SGEscript = <<'BWA_SCRIPT';
#!/usr/bin/perl
# : Align data with BWA
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -pe smp THREADS_TAG
#$ -o /dev/null
#$ -e BwaAlign.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $bwa             = "BWA_TAG";
my $threads         = "THREADS_TAG";
my $fasta           = "FASTA_TAG";
my $samtools        = "SAMTOOLS_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $prefix, $outdir);

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

        ($fileR1, $prefix, $outdir) = split("\t");
        
        if (!-e $fileR1 || !$outdir) { die "Incorrect input parameters! Check if files exist"; }
        if (!-d $outdir) { `mkdir -p $outdir`; }
}
close INPUTFILE;

# Construct prefix
$prefix ||= basename($fileR1);
# OICR specific
$prefix =~s/SWID_\d+_//;        # remove possible SWID part
my $s = $1 if $prefix =~/([A-Z]+)_\d+/;
$s ||= "SAMPLE";
my $rg ="\'\@RG\\tID:$prefix\\tSM:$s\'";

$outdir =~s!//!/!g;
if ($outdir !~m!/$!) {$outdir.="/";}
$prefix = $outdir.$prefix;
my $unsorted = $prefix.".unsorted.bam";

my $command = "$bwa mem -t $threads -R $rg -a -B 7 $fasta $fileR1 | $samtools view -Su - -b > $unsorted"; # | $samtools sort - $prefix";
 
print STDERR $command."\n";

`$command`;

BWA_SCRIPT

my $JobName = "BWA_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/THREADS_TAG/$threads/g;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/BWA_TAG/$bwa/;
$SGEscript =~s/FASTA_TAG/$ref/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/bwa_script$$.pl") or die "Couldn't write to [$Bin/bwa_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/bwa_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
