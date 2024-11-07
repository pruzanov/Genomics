#!/usr/bin/perl -w

=head2 Coverage Maker 2 (Unstranded)

 Development: make filtering run faster with samtools collate, pipe into modified script
 that will expect only primary alignments
 
 This version will run a script to extract unique alignments for reads,
 removing supplemental, secondary and extra primary alignments
 Coverage with bedtools will run on this filtered file.
 
 ./cvrg_maker2.pl --list [list of bam files] 
                  --ref [reference intervals] 
                  --coverage-bed [path to covergeBed] 
                  --tag [optional] 
                  --memory [Optional, gigabytes]
                  --samtools [samtools needed only when using filter]
                  --outdir [optional]
                  --filter [optional, filtered reads only]
                  --fasta [optional, need ref fasta when working with CRAM files]


 modules required:
 
 hg38/p12 samtools/1.14 xenoclassify/1.0 (for pysam)

 in order this to work also need to ln -s coverageMaker.py

=head2 On Bedtools

 seems that after 2.25 the logic of coverageBed has changed, will need to calculate coverage
 differently, will need to modify this script. Use 2.17!

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $covbed, $refbed, $tag, $mem, $filter, $samtools, $fasta);
my $USAGE = "cvrg_maker2.pl --list [list of bam files] --ref [reference intervals] --coverage-bed [path to coverageBed] --outdir [outdir, optional] --tag [optional suffix for outputs] --memory [Gigabytes allocated, default is 10] --filter [Optional flag, will apply -F 256 -q 30 samtools filter] --samtools [path to samtools]\n"; 
my $result = GetOptions("list=s"         => \$list,
                        "coverage-bed=s" => \$covbed,
                        "ref=s"          => \$refbed,
                        "memory|mem=i"   => \$mem,
                        "filter|f"       => \$filter,
                        "samtools=s"     => \$samtools,
                        "tag=s"          => \$tag,
                        "outdir=s"       => \$outdir);

if (!$list || !$refbed || !$covbed) { die $USAGE; }
if ($filter && !$samtools){ $samtools = "samtools";}
if (!$outdir) {
  print STDERR "Outdir will be the same as input dir for each file\n";
  $outdir = "";
} else {
  if (!-d $outdir) {
    print STDERR "Directory $outdir does not exists, create it (y/n)?\n";
    my $answer = <STDIN>;
    if ($answer !~/^y/i) {die "Aborting...\n";}
    
    `mkdir $outdir`;
  }
}
$tag ||="RM2";

my $SGEscript = <<'COVERAGE_SCRIPT';
#!/usr/bin/perl
# : Calculate coverage for files in a list
#$ -t 1-FILES_TAG
#$ -tc 10
#$ -l h_vmem=MEM_TAG
#$ -l h_rt=96:00:00
#$ -cwd
#$ -o /dev/null
#$ -e Coverage.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;

$ENV{'LD_LIBRARY_PATH'} = "LDPATH_TAG";

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $ResultOutputDir = "OUTDIR_TAG";   #Where we will put the results
my $FileList        = "INLIST_TAG";
my $coverageBed     = "COVBED_TAG";
my $refBed          = "REFBED_TAG";
my $tag             = "SUFFIX_TAG";
my $filter          = "FILTER_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $preprocessScript= "PREPSCRIPT_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $File;

my $Count=0;
while (<INPUTFILE>) {
        if (/^#/) { next;}
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }
        #Below here only processed for wanted lines: 
        $File = $_;
}
close INPUTFILE;

my $out = $File;
$out =~ s!(.*/)!!;
my $filedir = $1;


# Skip directory if there is a txt file
#my @files = `ls $filedir/*txt $filedir/*txt.gz`;
#if (@files > 0) { 
#    print STDERR "Found a coverage file, will not proceed for [$filedir]";
#    exit;
#}


# TODO: break the pipeline into three qsub operations: filtering, de-duplication, coverage
$out =~ s!\.bam$!.coverage!;
$out =~ s!\.cram$!.coverage!;
$ResultOutputDir ||= $filedir;
$out = join("/",($ResultOutputDir,$out));

if ($tag=~/\S+/ && $tag!~/^\./) { $tag = ".".$tag; }
if (!$File || !-e $File || !-s $File) { die "Couldn't read from File\n";}

my $cover = $out.$tag.".txt";
my $prebam = $out;
my $filtered = $out;
$prebam=~s/.coverage.*/_filtered_cleaned.bam/;
$filtered=~s/.coverage.*/_filtered.bam/;
print STDERR "Starting preprocessing [$File]\n";

# 1.Pre-filter, remove all non-primary alignments
if ($filtered && !-e $filtered && !-e $prebam) {
   `samtools collate -f --no-PG $File -o $filtered`;
}

# 2.Uniquify primary alignments
if ($filtered && -e $filtered && !-e $prebam) {
    `python3 $preprocessScript -b $filtered -o $ResultOutputDir -i $Line`;
    # remove filtered bam
    if ($prebam && -s $prebam) { `rm $filtered`; }
}

# need to check this
# 3.Coverage
if ($prebam && -e $prebam) {
  if ($filter) {
    `$samtools view -F 260 -q30 $prebam -b | $coverageBed -abam stdin -b $refBed -split > $cover`;
  } else {
    `$coverageBed -abam $prebam -b $refBed -split > $cover`;
  }
}

if (-f $cover && -s $cover) { `rm $prebam`; }

COVERAGE_SCRIPT

my $JobName  = "CVRG_$$";
my $files = `grep -v ^# $list | wc -l`;
my $prepScript = "$Bin/coverageMakerLite.py";
chomp($files);
$files =~s/\s+.*//;

$mem||=10;
$mem.="G" if $mem!~/G$/;

my $ldpath = $ENV{'LD_LIBRARY_PATH'};

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/LDPATH_TAG/$ldpath/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/COVBED_TAG/$covbed/;
$SGEscript =~s/PREPSCRIPT_TAG/$prepScript/;
$SGEscript =~s/REFBED_TAG/$refbed/;
$SGEscript =~s/FASTA_TAG/$fasta/ if $fasta;
$SGEscript =~s/OUTDIR_TAG/$outdir/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/SUFFIX_TAG/$tag/;
$SGEscript =~s/MEM_TAG/$mem/;
$filter ? $SGEscript =~s/FILTER_TAG/1/ : $SGEscript =~s/FILTER_TAG/0/;
$samtools ? $SGEscript =~s/SAMTOOLS_TAG/$samtools/ : $SGEscript =~s/SAMTOOLS_TAG/samtools/;

# Now, print the coverage script to the working directory and qsub it

open(SCRIPT,">$Bin/coverage_script2_$$.pl") or die "Couldn't write to [$Bin/coverage_script2_$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/coverage_script2_$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
