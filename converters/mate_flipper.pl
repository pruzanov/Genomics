#!/usr/bin/perl -w

=head2 Mate Flipper

 ./mate_flipper.pl --list [list of files] --outdir [output dir, optional]
 
 This script may be useful for those who want to get stranded coverage for paired-end reads
 which come from RNA-seq experiments 
 (let's just say you are interested in antisense transcripts)

 The problem is that SAM format sets binary tags so that mates are registered as sitting on
 opposite strands (which is true) and even if you are using stranded sequencing chemistry
 i.e. stranded TrueSeq

 Sequence of a read is always reported in the same orientation as the reference (no rev-complement)
 but if you try to use BEDtools which works fine with single-end sequencing you may encounter 
 a problem: the proportion of +/- reads will be pretty much balanced since no matter which
 transcript is read - sense or antisense. 

 This script flips the binary tag for the 2nd mate so that two reads in a pair are sitting in the
 same orientation and BEDtools counts them both toward either sense or antisense coverage, intersect etc.

 Note that the script works with BWA alignments generated with --library-type set to" fr-firststrand"
 in case of "fr-secondstrand" script needs to modified (it may be done automatically in a future)

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $samtools);
my $USAGE = "mate_flipper.pl --list [list of bam files] --samtools [path to samtools executable] --outdir [output dir, optional]\n";
my $result = GetOptions("list=s"     => \$list,
                        "samtools=s" => \$samtools,
                        "outdir=s"   => \$outdir);

if (!$list || !$samtools) { die $USAGE; }
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

my $SGEscript = <<'CONVERT_SCRIPT';
#!/usr/bin/perl
# : Calculate coverage for files in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub
my $samtools = "SAMTOOLS_TAG";

my $ResultOutputDir = "OUTDIR_TAG";   #Where we will put the results
my $FileList        = "INLIST_TAG";
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
$out =~ s!\.bam!\.mateflip.bam!;
$ResultOutputDir ||= $filedir;
$out = join("/",($ResultOutputDir,$out));

if (!$File || !-e $File || !-s $File) { die "Couldn't read from File\n";}

open BAM, "$samtools view -h $File |" or die "Cannot fork!"; 
my $temp_sam = join("/",($ResultOutputDir,"temp_$$.sam"));
open(OUT,">$temp_sam") or die "Couldn't write to temp SAM file [$temp_sam]";
while (<BAM>) {
 if (/^\@/) {
    print OUT $_;
    next;
 } # header lines

 my @temp = split("\t");

 if ($temp[1] & 64) { # If the read is read#1, mark as properly paired and 
     $temp[1] |=2;  # Mark as properly paired
     $temp[1] ^=32; # Reverse mate
 } elsif ($temp[1] & 128) {
     $temp[1] |=2;  # Mark as properly paired
     $temp[1] ^=16; # Reverse read
 }

 print OUT join("\t",@temp);

}

close BAM;
close OUT;

# Convert to BAM and rm sam file:

`$samtools view -Sh $temp_sam -b > $out`;
`rm $temp_sam`;

CONVERT_SCRIPT

my $JobName = "CONV_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/OUTDIR_TAG/$outdir/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/convert_script$$.pl") or die "Couldn't write to [$Bin/convert_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/convert_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
