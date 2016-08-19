#!/usr/bin/perl -w

=head2 Cufflinks Launcher (Multimapped fraction set to 1.0)
 
 .bam files will be processed by cufflinks tool
 but the list of files need to be formatted properly:

 sample1_file.bam	output_dir

 Will create cufflinks output in output_dir

 cufflinks_launcher.pl --list [mylist.txt] --cufflinks [path to cufflinks] --ref-gtf [transcripts.gtf]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $cufflinks, $refgtf);
my $USAGE = "cufflinks_launcher.pl --list [mylist.txt] --cufflinks [path to cufflinks] --ref-gtf [transcripts.gtf]\n";
my $result = GetOptions("list=s"         => \$list,
                        "cufflinks=s"    => \$cufflinks,
                        "ref-gtf=s"      => \$refgtf);

if (!$list || !$refgtf || !$cufflinks) { die $USAGE; }

my $SGEscript = <<'CUFFLINKS_SCRIPT';
#!/usr/bin/perl
# : Analyze expression (transcript abundance) for files in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=24G
#$ -o /dev/null
#$ -e Cufflinks.multi.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $cuffLinks       = "CUFFLINKS_TAG";
my $refGtf          = "REFGTF_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $Comparison;

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

        #Below here only processed for wanted lines: 
        $Comparison = $_;
}
close INPUTFILE;

my @out = split("\t", $Comparison);
if (!@out || @out != 2) {die "MALFORMED INSTRUCTION FILE!";}

my $ResultOutputDir ||= $out[1];
if ( ! -d $ResultOutputDir) {`mkdir $ResultOutputDir`;}

if (!$out[0] || !-e $out[0] || !-s $out[0]) { die "Couldn't read from File\n";}

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
`$cuffLinks $out[0] -G $refGtf -o $ResultOutputDir --library-type fr-firststrand -q -u –no-update-check –max-multiread-fraction 1.0 –min-frags-per-transfrag 100`;

CUFFLINKS_SCRIPT

my $JobName = "CUFLINK_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/CUFFLINKS_TAG/$cufflinks/;
$SGEscript =~s/REFGTF_TAG/$refgtf/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/cufflinks_script$$.pl") or die "Couldn't write to [$Bin/cufflinks_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/cufflinks_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
