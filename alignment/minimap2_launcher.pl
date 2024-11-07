#!/usr/bin/perl -w

=head2 minimap2 launcher
 
 align Single-end data using minimap2 

 sampleR_file.fasq.gz
 sampleT_file.fastq.gz
 ...

 minimap2_launcher.pl --list [mylist.txt] --minimap2 [path to minimap2 executable] --ref-base rererence_base

 By default, alignments in .bam format. Use --paf or -p to output in PAF format

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $mm2, $ref, $samtools, $prefix, $paf);
my $USAGE = "minimap2_launcher.pl --list [mylist.txt] --minimap2 [path to minimap2 executable] --ref [/some/path/mm2_idx.mmi] --samtools [samtools] --prefix [prefix]\n";
my $result = GetOptions("list=s"          => \$list,
                        "ref=s"           => \$ref,
                        "samtools=s"      => \$samtools,
                        "prefix=s"        => \$prefix,
                        "paf|p"           => \$paf,
                        "minimap2=s"      => \$mm2);

if (!$list || !$mm2) { die $USAGE; }
if (!$ref) {
  print STDERR "It is recommended that user specifies the reference file with bwa indexes in the same dir, default is hs1\n";
  $ref = " ";
}
if (!$samtools) {
 print STDERR "Default samtools will be used\n";
 $samtools = "samtools";
}


my $SGEscript = <<'MM2_SCRIPT';
#!/usr/bin/perl
# : Align data with minimap2
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=48:00:00
#$ -o /dev/null
#$ -e MM2Align.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $mm2             = "MM2_TAG";
my $refidx          = "REF_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $prefix          = "PREFIX_TAG";
my $paf             = "PAF_TAG";

$prefix = "" if $prefix=~/_TAG/;
open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($fileR1, $outdir);

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

        ($fileR1, $outdir) = split("\t");
        
        if (!-e $fileR1 || !$outdir) { die "Incorrect input parameters! Check if files exist"; }
        if (!-d $outdir) { `mkdir -p $outdir`; }
}
close INPUTFILE;

# Construct prefix
$prefix ||= basename($fileR1);
$prefix=~s/_L000(\d)_.*/_$1/;
$prefix=~s/\.fq$//;
$prefix=~s/\..*//;

my $format = $paf=~/TAG$/ ? ".bam" : ".paf";

$outdir =~s!//!/!g;
if ($outdir !~ m!/$!) {$outdir.="/";}
$prefix = $outdir.$prefix;
my $sorted = $prefix.".sorted".$format;
my $sort_prefix = $outdir."SORT_".$Line;

my $command;

if ($paf=~/TAG$/) {
  $command = "minimap2 -a $refidx $fileR1 | $samtools view -Su - -b | $samtools sort -T $sort_prefix -o $sorted -";
} else {
  $command = "minimap2 -x map-ont $refidx $fileR1 > $sorted";
}

print STDERR $command."\n";

`$command`;

MM2_SCRIPT

my $JobName = "MM2_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/MM2_TAG/$mm2/;
$SGEscript =~s/REF_TAG/$ref/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/PREFIX_TAG/$prefix/ if $prefix;
$SGEscript =~s/PAF_TAG/paf/ if $paf;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/mm2_script$$.pl") or die "Couldn't write to [$Bin/mm2_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/mm2_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
