#!/usr/bin/perl -w

=head2 STAR launcher
 
 align RNA-Seq data using STAR

 sampleR_file.R1.fastq.gz	sampleR_file.R2.fastq.gz
 sampleT_file.R1.fastq.gz	sampleT_file.R2.fastq.gz

 star_launcher.pl --list [mylist.txt] --star [path to STAR executable] --threads 6 --index-dir /mt/dir/star-indexes/75b --adaptors [String(s) of adaptors to clip]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $star, $index_dir, $threads, $adaptors);
my $USAGE = "star_launcher.pl --list [mylist.txt] --star [path to STAR executable] --threads [threads to use] --index-dir [Directory with STAR indexes]\n";
my $result = GetOptions("list=s"          => \$list,
                        "index-dir=s"     => \$index_dir,
                        "adaptors=s"      => \$adaptors,
                        "threads=i"       => \$threads,
                        "star=s"     => \$star);

if (!$list || !$star) { die $USAGE; }
$threads ||=6;
if (!$index_dir) {
  print STDERR "It is recommended that user specifies the index directory, reverting to default specific to 100b reads\n";
  $index_dir = "$HOME/data/reference/hg19_random/star2.5/100b";
}
$adaptors ||= "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

my $SGEscript = <<'STAR_SCRIPT';
#!/usr/bin/perl
# : Align RNAseq data with STAR
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -pe smp THREADS_TAG
#$ -o /dev/null
#$ -e StarAlign.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $star            = "STAR_TAG";
my $threads         = "THREADS_TAG";
my $indexdir        = "IDXDIR_TAG";
my $adaptors        = "ADAPTOR_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $fileR2, $prefix, $outdir);

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

        ($fileR1, $fileR2, $outdir) = split("\t");
        if (!-e $fileR1 || !-e $fileR2 || !-d $outdir) { die "Incorrect input parameters! Check if files exist"; }
}
close INPUTFILE;

# Construct prefix
$prefix = basename($fileR1);
$prefix =~s/SWID_\d+_//;        # remove possible SWID part
$prefix =~s/^(\w+?_\d+)_.*/$1/; # extract donor id
$outdir =~s!//!/!g;
if ($outdir !~m!/$!) {$outdir.="/";}
$prefix = $outdir.$prefix;
# TODO: Add RG
my $command = "$star --runThreadN $threads ";
   $command.= "--genomeDir $indexdir ";
   $command.= "--readFilesIn $fileR1 $fileR2 ";
   $command.= "--clip3pAdapterSeq $adaptors ";
   $command.= "--readFilesCommand zcat ";
   $command.= "--outFileNamePrefix $prefix ";
   $command.= "--outSAMmultNmax -1 ";
   $command.= "--outSAMtype BAM SortedByCoordinate";

#print STDERR $command."\n";

`$command`;

STAR_SCRIPT

my $JobName = "STAR_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/ADAPTOR_TAG/$adaptors/g;
$SGEscript =~s/THREADS_TAG/$threads/g;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/STAR_TAG/$star/;
$SGEscript =~s/IDXDIR_TAG/$index_dir/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/star_script$$.pl") or die "Couldn't write to [$Bin/star_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/star_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
