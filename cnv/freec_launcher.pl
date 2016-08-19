#!/usr/bin/perl -w

=head2 FREEC launcher
 
 process bam data (controls ONLY) using FREEC

 mybam_file1.bam	my_output_dir1
 mybam_file2.bam        my_output_dir2

 freec_launcher.pl --list [mylist.txt] --freec [dir for FREEC] --freec-script freec_script --r-libs [directory with R libraries] --samtools samtools

 need external FREEC launching script

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $freec, $freec_script, $rlibdir, $samtools, $type, $lenfile, $chrfiles);
my $USAGE = "freec_launcher.pl --list [mylist.txt] --freec [freec executable] --freec-script freec_script --samtools samtools --r-libs rlib-dir --chr-files dir_with_fasta";
my $result = GetOptions("list=s"             => \$list,
                        "r-libs=s"           => \$rlibdir,
                        "samtools=s"         => \$samtools,
                        "lenfile=s"          => \$lenfile,
                        "chr-files=s"        => \$chrfiles,
                        "freec=s"            => \$freec,
                        "freec-script=s"     => \$freec_script);
$freec_script ||="$Bin/scripts/launchFREEC.g.pl";
$lenfile  ||="$Bin/scripts/hg19_chr.len";
$chrfiles ||=$NV{HOME}."/data/reference/hg19_random/fasta/UCSC/";
if (!$samtools || !$freec || !$rlibdir || !$list) { die $USAGE; }

my $SGEscript = <<'FREEC_SCRIPT';
#!/usr/bin/perl
# : Run FREEC pipeline for germline CNV analysis
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=8G
#$ -o /dev/null
#$ -e FREEC.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
use constant DEBUG=>1;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $freec           = "FREEC_TAG";
my $rlibdir         = "RLIB_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $launcher_script = "LAUNCHER_TAG";
my $lenfile         = "LENFILE_TAG";
my $chrfiles        = "CHRFILES_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($bam, $outdir);

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

        ($bam, $outdir) = split("\t");
        if (!-e $bam || ! -d $outdir) { die "Incorrect input parameters! Check if files exist"; }
        my $freec_dir=$outdir."/freec";
        $freec_dir=~s!//!/!g;
        if (-d $freec_dir) {
            print STDERR "FREEC analysis has been run for $bam, won't continue\n";
            exit;
        } 

        `mkdir $freec_dir`;
        $outdir = $freec_dir;
}
close INPUTFILE;

# Construct prefix
my $prefix = basename($bam);
$prefix =~s/SWID_\d+_//;                  # remove possible SWID part
$prefix = $1 if $prefix=~/^(\S+?_\d+)/;   # extract donor id
print STDERR "Prefix defined a $prefix\n" if DEBUG;

if ($outdir !~m!/$!) {$outdir.="/";}
# my $data_dir = $outdir."data";
# if (!-d $data_dir) { `mkdir $data_dir`; }

my $command = "$launcher_script --input $bam --freec $freec --samtools $samtools --id $prefix --r-libdir $rlibdir --outdir $outdir --lenfile $lenfile --chr-files $chrfiles";
print STDERR "Running FREEC via external script as [$command] ...\n" if DEBUG;

`$command`;

FREEC_SCRIPT

my $JobName = "FREEC_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/FREEC_TAG/$freec/;
$SGEscript =~s/LAUNCHER_TAG/$freec_script/;
$SGEscript =~s/RLIB_TAG/$rlibdir/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/LENFILE_TAG/$lenfile/;
$SGEscript =~s/CHRFILES_TAG/$chrfiles/;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/freec_script$$.pl") or die "Couldn't write to [$Bin/freec_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/freec_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
