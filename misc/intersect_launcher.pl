#!/usr/bin/perl -w

=head2 Intersect Maker
 
 A simple intersect checker for bam files. Will produce
 a bam file of reads intersecting with intervals in supplied .bed

 given a tag, will make a file with this tag, like so:

 myinput.bam
 myinput.parsed.bam (where 'parsed' passed using tag argument)

 ./intersect_maker.s.pl --list [list of bam files] --ref [reference intervals] --tag [for making output names] --intersect-bed [path to intersectBed] --outdir [outdir, optional]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $intbed, $refbed, $ioptions, $tag);
my $USAGE = "intersect_launcher.pl --list [list of bam files] --ref [reference intervals] --intersect-bed [path to IntersectBed] --outdir [outdir, optional] --tag [optional suffix for outputs] --ioptions [optonal, options for intersectBed]\n"; 
my $result = GetOptions("list=s"         => \$list,
                        "intersect-bed=s"=> \$intbed,
                        "ref=s"          => \$refbed,
                        "tag=s"          => \$tag,
                        "ioptions=s"     => \$ioptions,
                        "outdir=s"       => \$outdir);

if (!$list || !$refbed || !$intbed || !$tag) { die $USAGE; }
$ioptions ||=" ";
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

my $SGEscript = <<'INTERSECT_SCRIPT';
#!/usr/bin/perl
# : Calculate intersect for files in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=6G
#$ -o /dev/null
#$ -e /dev/null
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $ResultOutputDir = "OUTDIR_TAG";   #Where we will put the results
my $FileList        = "INLIST_TAG";
my $intersectBed    = "INTBED_TAG";
my $refBed          = "REFBED_TAG";
my $options         = "OPTIONS_TAG";
my $tag             = "SUFFIX_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $File;

my $Count=0;
while (<INPUTFILE>) {
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
$out =~ s!\.bam!.isect!;
$ResultOutputDir ||= $filedir;
$out = join("/",($ResultOutputDir,$out));

if ($tag=~/\S+/ && $tag!~/^\./) { $tag = ".".$tag; }
if (!$File || !-e $File || !-s $File) { die "Couldn't read from File\n";}

if ($File && -e $File) {
  my $outPath = $out.$tag;
  my $type;
  if ($File=~/bam$/) {
      $type =  "-u -abam";
      $options =~/-bed/ ? $outPath.=".bed" : $outPath.=".bam";
  } else {
      $type = "-u -a";
      $outPath.=".bed";
  }
  
  # The line that needs to be customized to suit the needs of the user:
  
  `$intersectBed -f 0.6 $type $File -b $refBed $options > $outPath`;
}

INTERSECT_SCRIPT

my $JobName = "ISECT_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/INTBED_TAG/$intbed/;
$SGEscript =~s/REFBED_TAG/$refbed/;
$SGEscript =~s/OUTDIR_TAG/$outdir/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/OPTIONS_TAG/$ioptions/;
$SGEscript =~s/SUFFIX_TAG/$tag/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/intersect_script$$.pl") or die "Couldn't write to [$Bin/intersect_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/intersect_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
