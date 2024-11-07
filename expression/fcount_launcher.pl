#!/usr/bin/perl -w

=head2 launch featureCounter from SUBREAD package

 sample1_file.bam
 sample2_file.bam

 fcount_launcher.pl --list [mylist.txt] -f [count features (exons), not meta-features (genes)] --gtf [path to GTF file] --id [Optional, Default:gene_name] --feature-counter [path to featureCounter]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $fcounter, $features, $id, $gtf);
my $USAGE = "fcount_launcher.pl --list [mylist.txt] -f [count features (exons), not meta-features (genes)] --gtf [path to GTF file] --id [Optional, Default:gene_name] --feature-counter [path to featureCounter]\n";
my $result = GetOptions("list=s"            => \$list,
                        "f"                 => \$features,
                        "id=s"              => \$id,
                        "gtf=s"             => \$gtf,
                        "feature-counter=s" => \$fcounter);

if (!$list || !$gtf) { die $USAGE; }
$id ||="gene_name";
$fcounter ||="/u/pruzanov/Downloads/subread-2.0.6-Linux-x86_64/bin/featureCounts";

my $SGEscript = <<'FCOUNT_SCRIPT';
#!/usr/bin/perl
# : count reads per gene with featureCount
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=6:00:00
#$ -o /dev/null
#$ -e fCounter.e
#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $fcounter        = "FCOUNT_TAG";
my $ftag            = "F_TAG";
my $gtf             = "GTF_TAG";
my $id              = "ID_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $bam;

my $Count=0;
while (<INPUTFILE>) {
        next if /^#/;
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)  {  next;  }       #Skip until the lines we want:

        #Below here only processed for wanted lines: 
        $bam = $_;
        $bam =~s/\s+.*//;
        $bam =~s/\t.*//;
        last;
}
close INPUTFILE;
my $output = $bam;
$output=~s/bam$/counts/;

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "$fcounter -t exon -g $id -p -a $gtf -o $output $bam ";
$command.= "-f " if $ftag!~/TAG$/;

print STDERR $command."\n";

`$command`;

FCOUNT_SCRIPT

my $JobName = "FCNT_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/FCOUNT_TAG/$fcounter/;
if ($features) {$SGEscript =~s/F_TAG/ENABLE/;}
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/GTF_TAG/$gtf/;
$SGEscript =~s/ID_TAG/$id/;
# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/fcount_script$$.pl") or die "Couldn't write to [$Bin/fcount_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/fcount_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
