#!/usr/bin/perl -w

=head2 Coverage Maker (Unstranded)
 
 ./cvrg_maker.s.pl --list [list of bam files] --ref [reference intervals] --coverage-bed [path to covergaeBed] --outdir [outdir, optional] --memory [Optional, gigabytes]

=head2 On Bedtools

 seems that after 2.25 the logic of coverageBed has changed, will need to calculate coverage
 differently, will need to modify this script

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $covbed, $refbed, $tag, $mem, $filter, $samtools);
my $USAGE = "cvrg_maker.s.pl --list [list of bam files] --ref [reference intervals] --coverage-bed [path to coverageBed] --outdir [outdir, optional] --tag [optional suffix for outputs] --memory [Gigabytes allocated, default is 10]\n"; 
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
$tag ||="";

my $SGEscript = <<'COVERAGE_SCRIPT';
#!/usr/bin/perl
# : Calculate coverage for files in a list
#$ -t 1-FILES_TAG
#$ -l h_vmem=MEM_TAG
#$ -l h_rt=56:00:00
#$ -cwd
#$ -o /dev/null
#$ -e Coverage.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $ResultOutputDir = "OUTDIR_TAG";   #Where we will put the results
my $FileList        = "INLIST_TAG";
my $coverageBed     = "COVBED_TAG";
my $refBed          = "REFBED_TAG";
my $tag             = "SUFFIX_TAG";
my $filter          = "FILTER_TAG";
my $samtools        = "SAMTOOLS_TAG";

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
$out =~ s!\.bam!.coverage!;
$ResultOutputDir ||= $filedir;
$out = join("/",($ResultOutputDir,$out));

if ($tag=~/\S+/ && $tag!~/^\./) { $tag = ".".$tag; }
if (!$File || !-e $File || !-s $File) { die "Couldn't read from File\n";}

my $cover = $out.$tag.".txt";

# need to check this
# First read is always synthesized from first cDNA strand ? (antisense)
if ($File && -e $File) {
  if ($filter) {
    `$samtools view -F 260 -q30 $File -b | $coverageBed -abam stdin -b $refBed -split > $cover`; # -split
  } else {
    `$coverageBed -abam $File -b $refBed -split > $cover`; # -split
  }
}

COVERAGE_SCRIPT

my $JobName = "CVRG_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$mem||=10;
$mem.="G" if $mem!~/G$/;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/COVBED_TAG/$covbed/;
$SGEscript =~s/REFBED_TAG/$refbed/;
$SGEscript =~s/OUTDIR_TAG/$outdir/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/SUFFIX_TAG/$tag/;
$SGEscript =~s/MEM_TAG/$mem/;
$filter ? $SGEscript =~s/FILTER_TAG/1/ : $SGEscript =~s/FILTER_TAG/0/;
$samtools ? $SGEscript =~s/SAMTOOLS_TAG/$samtools/ : $SGEscript =~s/SAMTOOLS_TAG/samtools/;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/coverage_script$$.pl") or die "Couldn't write to [$Bin/coverage_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/coverage_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
