#!/usr/bin/perl

=head2 samtools flagstat launcher
 
 collect stats on bam files in supplied list

 file01.bam
 file02.bam
 file03.bam

 samstat_launcher.pl --list [mylist.txt] --samtools [path to samtools] --report [Optional summary report file]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $samtools, $report);
my $USAGE = "samstat_indexer.pl --list [mylist.txt] --samtools [path to samtools] --report [Optional summary report file]\n";
my $result = GetOptions("list=s"         => \$list,
                        "samtools=s"     => \$samtools,
                        "report=s"       => \$report);

if (!$list || !$samtools ||!$report) { die $USAGE; }

my $SGEscript = <<'SAMSTAT_SCRIPT';
#!/usr/bin/perl
# : Flagstat .bam files with samtools
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=4G
#$ -o /dev/null
#$ -hold_jid PCSI612.mrg
#$ -e Stats.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $samtools        = "SAMTOOLS_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($bam,$report);

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        #Below here only processed for wanted lines: 
        $bam = $_;
        $bam =~s/\s+.*//;
        $bam =~s/\t.*//;
        $report = $bam;
        $report=~s/bam$/stat/;
}
close INPUTFILE;

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "$samtools flagstat $bam > $report";

#print STDERR $command."\n";

`$command`;

SAMSTAT_SCRIPT

my $JobName = "SAMSTAT_$$";
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

# Ratio calculating script
my $SGECollector = <<'REPORT_SCRIPT';
#!/usr/bin/perl
# : Collect reports in a single .stat file
#$ -cwd
#$ -l h_vmem=2G
#$ -hold_jid STATJOB_TAG
#$ -o /dev/null
#$ -e Summary.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $listfile      = "INLIST_TAG";
my $report        = "REPORT_TAG";
my @lines         = ();

open(IN, "<$listfile") or die "Couldn't open file [$listfile] for reading";
while(<IN>) {
 chomp;
 my $repfile = $_;
 $repfile=~s/bam$/stat/;
 `echo $repfile >> $report`;
 
 if (-e $repfile) {
   `cat $repfile >> $report`;
 } else {
   `echo "ERROR collecting stats for this file" >> $report`;
 }

}
close IN;

REPORT_SCRIPT

my $CounterName = "COUNT_$$";

$SGECollector =~s/STATJOB_TAG/$JobName/;
$SGECollector =~s/JOBNAME_TAG/$CounterName/;
$SGECollector =~s/REPORT_TAG/$report/;
$SGECollector =~s/INLIST_TAG/$list/;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/counter_script$$.pl") or die "Couldn't write to [$Bin/counter_script$$.pl]";
print SCRIPT $SGECollector;
close SCRIPT;

# qsub
sleep(10); # wait 10 sec

print STDERR "Submitting Collector script\n";

$SGEResult = `qsub $Bin/counter_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";

