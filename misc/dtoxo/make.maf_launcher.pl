#!/usr/bin/perl -w

=head2 Make MAF launcher
 
 make MAF files for DToxoG using pre-compiled .m files from Picard, .bam and vcf files
 but the list of files (tab-delimited) needs to be formatted properly:

 sampleP_file.bam	sample.mutect.vcf	output.maf
 ...

 (metric files will be in the same dir as input)

 make.maf_launcher.pl --list [mylist.txt] --maf-script [path to maf maker script]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $mafmaker, $isectbed, $samtools);
my $USAGE = "make.maf_launcher.pl --list [mylist.txt] --maf-script [path to maf maker script] --isect-bed [path to IntersectBed] --samtools [path to samtools]\n";
my $result = GetOptions("list=s"         => \$list,
                        "maf-script=s"   => \$mafmaker,
                        "isect-bed=s"    => \$isectbed,
                        "samtools=s"     => \$samtools);

if (!$list || !$mafmaker || !$samtools || !$isectbed) { die $USAGE; }

my $SGEscript = <<'MAF_SCRIPT';
#!/usr/bin/perl
# : Make MAF files for using with DToxoG
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -o /dev/null
#$ -e Mafmaker.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $refFasta  = "/.mounts/labs/PDE/data/MutectStrelka/hg19/hg19_random.fa";
my $FileList        = "INLIST_TAG";
my $mafmaker        = "MAFSCRIPT_TAG";
my $isectBed        = "ISECT_TAG";
my $samtools        = "SAMTOOLS_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $inputLine;

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        #Below here only processed for wanted lines: 
        $inputLine = $_;
        if ($inputLine=~/^#/) {
          print STDERR "Skipping Line $Line\n";
          exit;
        }
}
close INPUTFILE;

chomp($inputLine);

my($bam,$vcf,$output) = split("\t",$inputLine);
if (!-e $bam || !-e $vcf) { die "Couldn't find inputs"; }
my $poxo = $bam.".m";

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = $mafmaker;

$command.= " --bam $bam";
$command.= " --vcf $vcf";
$command.= " --picard-oxog-file $poxo";
$command.= " --samtools $samtools";
$command.= " --ref-fasta $refFasta";
$command.= " --isect-bed $isectBed";
$command.= " --tumor-id TUMOR";
$command.= " --normal-id NORMAL";
$command.= " --filter PASS";
$command.= " --out $output";

#print STDERR $command."\n";

`$command`;

MAF_SCRIPT

my $JobName = "MAKEMAF_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/MAFSCRIPT_TAG/$mafmaker/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/ISECT_TAG/$isectbed/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;



# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/mafmaker_script$$.pl") or die "Couldn't write to [$Bin/mafmaker_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/mafmaker_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
