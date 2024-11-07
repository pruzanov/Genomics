#!/usr/bin/perl -w

=head2 BWA samse launcher
 
 convert Single-end sai files using BWA

 idS1	sampleS1.sai	sampleS1_file.R1.fasq.gz	outdir
 idS2	sampleS2.sai	sampleS2_file.R1.fastq.gz	outdir

 bwa_samse_launcher.pl --list [mylist.txt] --bwa [path to bwa executable] --threads 8 --ref-fasta fererence.fa [the same dir needs to have indexes too]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $bwa, $ref, $samtools);
my $USAGE = "bwa_launcher.pl --list [mylist.txt] --bwa [path to BWA executable] --ref-fasta [hg19_random.fa] --samtools [samtools]\n";
my $result = GetOptions("list=s"          => \$list,
                        "ref-fasta=s"     => \$ref,
                        "samtools=s"      => \$samtools,
                        "bwa=s"           => \$bwa);

if (!$list) { die $USAGE; }
if (!$ref) {
  print STDERR "It is recommended that user specifies the reference file with bwa indexes in the same dir\n";
  $ref = "/.mounts/labs/PDE/data/RegressionTests/BWA/workflow/input_data/0.6.2/hg19_random.fa";
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
# : Convert sai data with BWA
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -o /dev/null
#$ -e BwaConvert.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $bwa             = "BWA_TAG";
my $fasta           = "FASTA_TAG";
my $samtools        = "SAMTOOLS_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($file, $fileFQ, $prefix, $outdir);

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

        ($prefix, $file, $fileFQ, $outdir) = split("\t");
        
        if (!-e $file || !-e $fileFQ || !$outdir) { die "Incorrect input parameters! Check if files exist"; }
        if (!-d $outdir) { `mkdir -p $outdir`; }
}
close INPUTFILE;

# Construct prefix
$prefix ||= basename($file);
$prefix =~s/SWID_\d+_//;        # remove possible SWID part

$outdir =~s!//!/!g;
if ($outdir !~m!/$!) {$outdir.="/";}
$prefix = $outdir.$prefix;
my $unsorted = $prefix.".unsorted.bam";

my $command = "$bwa samse $fasta $file $fileFQ | $samtools view -Su - -b > $unsorted"; 
 
print STDERR $command."\n";

`$command`;

BWA_SCRIPT

my $JobName = "SAMSE_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/BWA_TAG/$bwa/;
$SGEscript =~s/FASTA_TAG/$ref/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/samse_script$$.pl") or die "Couldn't write to [$Bin/samse_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/samse_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
