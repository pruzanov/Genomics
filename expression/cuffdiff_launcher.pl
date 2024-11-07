#!/usr/bin/perl -w

=head2 Cuffdiff Launcher
 
 .bam files will be processed by cuffdiff tool
 but the list of files need to be formatted properly:

 sample1_file.bam	sample2_file.bam	output_dir	label1,label2

 OR

 sample1_file.bam:label1       sample2_file.bam:label2        output_dir

 Will create cuffdiff output in output_dir

 cuffdiff_launcher.pl --list [mylist.txt] --cuffdiff [path to cuffdiff] --ref-gtf [transcripts.gtf]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $cuffdiff, $refgtf);
my $USAGE = "cuffdiff_launcher.pl --list [mylist.txt] --cuffdiff [path to cuffdiff] --ref-gtf [transcripts.gtf]\n";
my $result = GetOptions("list=s"         => \$list,
                        "cuffdiff=s"     => \$cuffdiff,
                        "ref-gtf=s"      => \$refgtf);

if (!$list || !$refgtf || !$cuffdiff) { die $USAGE; }

my $SGEscript = <<'CUFFDIFF_SCRIPT';
#!/usr/bin/perl
# : Analyze Differrential expression for files in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=24G
#$ -l h_rt=48:00:00
#$ -pe smp 1
#$ -o /dev/null
#$ -e Cuffdiff.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $cuffDiff        = "CUFFDIFF_TAG";
my $refGtf          = "REFGTF_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $Comparison;

my $Count=0;
while (<INPUTFILE>) {
        next if /^#/;
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {   next;   }       #Skip until the lines we want:
        if ($Count > $Line )       {   last;   }

        #Below here only processed for wanted lines: 
        $Comparison = $_;
}
close INPUTFILE;

my @out = split("\t", $Comparison);
if (!@out || @out < 3) {die "MALFORMED INSTRUCTION FILE!";}

my $ResultOutputDir ||= $out[2];
if ( ! -d $ResultOutputDir) {`mkdir $ResultOutputDir`;}

# Check if all files exist
foreach(@out[0..1]) {
  my @files = split(",");
  map{if (!$_ || !-e $_ || !-s $_) { die "Couldn't read from File\n";}} @files;
}

my($label1,$label2);


if ($out[0] =~/\:/) {
 $label1 = $`;
} else {
 $label1 = $1 if $out[0]=~m![A-Z]{4,5}_*\d+!;
}

if ($out[1]=~/\:/) {
 $label2 = $`;
} else  {
 $label2 = $1 if $out[1]=~m![A-Z]{4,5}_*\d+!; 
}

if ($out[2] && $out[2]=~/\,/) {($label1,$label2) = split(",",$out[2]);}

my $labels ="";
if ($label1 && $label2) {$labels = "--labels $label1".",".$label2;}

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
# --library-type fr-firststrand
my $c = "$cuffDiff -v $refGtf $out[0] $out[1] -o $ResultOutputDir $labels -p 1 --no-update-check --library-norm-method classic-fpkm --max-mle-iterations 2000";
print STDERR "Executing command:\n";
print STDERR $c."\n";
`$c`;

CUFFDIFF_SCRIPT

my $JobName = "CUFDIF_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/CUFFDIFF_TAG/$cuffdiff/;
$SGEscript =~s/REFGTF_TAG/$refgtf/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/cuffdiff_script$$.pl") or die "Couldn't write to [$Bin/cuffdiff_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/cuffdiff_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
