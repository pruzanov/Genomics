#!/usr/bin/perl -w

=head2 HISAT2 launcher
 
 align RNA-Seq data using HISAT2

 sampleR_file.R1.fastq.gz	sampleR_file.R2.fastq.gz	sample_Alignments.bam
 sampleT_file.R1.fastq.gz	sampleT_file.R2.fastq.gz

 hisat2_launcher.pl --list [mylist.txt] --hisat [path to HISAT2 executable] --threads 6 --index-base /mt/dir/hisat2-indexes/hg38 --adaptors [String(s) of adaptors to clip]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $hisat, $index_base, $threads, $adaptors, $samtools);
my $USAGE = "hisat2_launcher.pl --list [mylist.txt] --hisat [path to HISAT2 executable] --threads [threads to use] --index-base [basename HISAT2 indexes] --samtools [path to samtools binary]\n";
my $result = GetOptions("list=s"          => \$list,
                        "index-base=s"    => \$index_base,
                        "samtools=s"      => \$samtools,
                        "threads=i"       => \$threads,
                        "hisat=s"         => \$hisat);

if (!$list || !$hisat || !$samtools) { die $USAGE; }
$threads ||=6;
if (!$index_base) {
  print STDERR "It is recommended that user specifies the index directory, reverting to default\n";
  $index_base = "/.mounts/labs/prod/scratch/pruzanov/hg38_indexes/hisat2/hg38_genome";
}

my $SGEscript = <<'HISAT_SCRIPT';
#!/usr/bin/perl
# : Align RNAseq data with HISAT2
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=14G
#$ -pe smp THREADS_TAG
#$ -o /dev/null
#$ -e HisatAlign.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $hisat           = "HISAT_TAG";
my $threads         = "THREADS_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $indexbase       = "IDXBASE_TAG";
my $adaptors        = "ADAPTOR_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $fileR2, $outbam, $rgline);

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

        ($fileR1, $fileR2, $outbam, $rgline) = split("\t");
        if (!-e $fileR1 || !-e $fileR2 ) { die "Incorrect input parameters! Check if files exist"; }
}
close INPUTFILE;

# Construct prefix
my $command = "$hisat -p $threads ";
   $command.= "-x $indexbase ";
   $command.= "--pen-cansplice 4 ";
   $command.= "-k 10 ";
   $command.= "--fr ";
   if ($rgline) {
    my @rgs = split(";",$rgline);
    foreach my $rg (@rgs) {
      $rg =~ /^ID\:/ ? $command.="--rg-id $rg " : $command.="--rg $rg ";
    }
   }    
   $command.= "--secondary ";
   $command.= "-1 $fileR1 -2 $fileR2 ";
   $command.= "| $samtools view -h - -b > ".$outbam;
print STDERR $command."\n\n";

`$command`;

HISAT_SCRIPT

my $JobName = "HISAT_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/THREADS_TAG/$threads/g;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/HISAT_TAG/$hisat/;
$SGEscript =~s/IDXBASE_TAG/$index_base/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/hisat_script$$.pl") or die "Couldn't write to [$Bin/hisat_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -q production $Bin/hisat_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
