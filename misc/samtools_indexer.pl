#!/usr/bin/perl -w

=head2 samtools indexer
 
 index .bam files from a list

 sampleR_file.bam
 sampleT_file.bam

 samtools_indexer.pl --list [mylist.txt] --samtools [path to samtools]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $samtools);
my $USAGE = "samtools_indexer.pl --list [mylist.txt] --samtools [path to samtools]\n";
my $result = GetOptions("list=s"         => \$list,
                        "samtools=s"     => \$samtools);

if (!$list || !$samtools) { die $USAGE; }

my $SGEscript = <<'SAMTOOLS_SCRIPT';
#!/usr/bin/perl
# : Index .bam files with samtools
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=8G
#$ -o /dev/null
#$ -e Indexing.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $samtools        = "SAMTOOLS_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $bam;

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
        $bam = $_;
        $bam =~s/\s+.*//;
        $bam =~s/\t.*//;
}
close INPUTFILE;

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "$samtools index $bam";

#print STDERR $command."\n";

`$command`;

SAMTOOLS_SCRIPT

my $JobName = "BAMIDX_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/samtools_script$$.pl") or die "Couldn't write to [$Bin/samtools_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/samtools_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
