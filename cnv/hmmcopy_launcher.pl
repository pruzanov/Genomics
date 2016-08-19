#!/usr/bin/perl -w

=head2 HMMcopy launcher
 
 process bam data (controls ONLY) using HMMcopy

 mybam_file1.bam	my_output_dir1
 mybam_file2.bam        my_output_dir2

 hmmcopy_launcher.pl --list [mylist.txt] --hmm-dir [dir for HMMCOPY] --script-dir script-dir --r-libs [directory with R libraries] --ref-dir [dir with ref files for HMMcopy]

 need external HMMCOPY launching script

 A configuration file should point to mapping and cg files
 Comments are alllowed, spaces - no
 
 cg_file=/my/dir/cg_file
 map_file=/my/dir/map_file

 In order for this to work we need to:

 1. index all input bam files
 2. run 'module load R/3.2.1' before launching

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $hmmbindir, $scriptdir, $rlibdir, $refconf, $outbase);
my $USAGE = "hmmcopy_launcher.pl --list [mylist.txt] --hmm-bindir [bin dir for HMMCOPY] --script-dir script-dir --r-libs [directory with R libraries] --ref-conf [conf file for HMMcopy]";


my $result = GetOptions("list=s"             => \$list,
                        "r-libs=s"           => \$rlibdir,
                        "hmm-bindir=s"       => \$hmmbindir,
                        "ref-conf=s"         => \$refconf,
                        "out-base=s"         => \$outbase);
# DEFAULTS
$scriptdir  ||= "$Bin/scripts";
$refconf    ||="$scriptdir/HMMCOPY.conf";
if (!$hmmbindir || !$rlibdir || !$list || !-d $hmmbindir) { die $USAGE; }

my $SGEscript = <<'HMMCOPY_SCRIPT';
#!/usr/bin/perl
# : Run HMMCOPY pipeline for germline CNV analysis
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=6G
#$ -o /dev/null
#$ -e HMMCOPY.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
use constant DEBUG=>0;

my $FileList         = "INLIST_TAG";
my $hmmbindir        = "HMMBINDIR_TAG";
my $rlibdir          = "RLIB_TAG";
my $hmm_script       = "LAUNCHER_TAG";
my $hmm_convert      = "CONVERTER_TAG";
my $hmm_config       = "HMMCONFIG_TAG";

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($bam, $outdir);

my $Count=0;
while (<INPUTFILE>) {
        print STDERR "Reading from $FileList\n" if DEBUG;
        chomp;
        $Count++;  #Increment the line counter
        print STDERR "Will see if we have line $Line, current is $Count\n" if DEBUG;
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        # Below only lines that we want get processed: 
        if (/^#/) {
          print STDERR "Skipping Line $Line\n";
          exit;
        }

        ($bam, $outdir) = split("\t");
        print STDERR "Extracted $bam in $outdir\n" if DEBUG;
        if (!-e $bam || ! -d $outdir) { die "Incorrect input parameters! Check if files exist"; }
        
      
        my $hmm_dir=$outdir."/hmmcopy";
        $hmm_dir=~s!//!/!g;
      
        if (-d $hmm_dir) {
            print STDERR "HMMCOPY analysis has been run for $bam, won't continue\n";
            exit;
        } 
        print STDERR "Will launch analysis for $bam in $hmm_dir\n" if DEBUG;

        `mkdir $hmm_dir`;
        $outdir = $hmm_dir;
}
close INPUTFILE;

# Construct prefix
my $prefix = basename($bam);
$prefix =~s/SWID_\d+_//;                  # remove possible SWID part
$prefix = $1 if $prefix=~/^(\S+?_\d+)/;   # extract donor id
print STDERR "Prefix defined as $prefix\n" if DEBUG;

if ($outdir !~m!/$!) {$outdir.="/";}
my $data_dir = $outdir."data/";
my $output_base = $data_dir.$prefix;
if (!-d $data_dir) { `mkdir $data_dir`; }

# 1. HMMCOPY convert script

my $outwig = join("/",($data_dir,basename($bam)));
$outwig =~s/.bam$/_reads.wig/;
my $command = "$hmm_convert --read-counter $hmmbindir/readCounter --input $bam --output $outwig";

print STDERR "Running HMMCOPY converting step via external script as [$command] ...\n" if DEBUG;

`$command`;

# 2. HMMCOPY analysis
# TODO read the files from config
my($cg_file,$map_file);

open(CONF, "<$hmm_config") or die "Couldn't open Config file for hmm copy analysis";
while(<CONF>) {
  chomp;
  my @data = split("=",$_);
  if ($data[0] eq "cg_file")  {$cg_file = $data[1];}
  if ($data[0] eq "map_file") {$map_file = $data[1];}
}
close CONF;

unless (-e $cg_file && -e $map_file) { die "Couldn't find parameters for HMM copy analysis, aborting"; }
$command = "Rscript $hmm_script $outwig $cg_file $map_file $output_base";
unless (-e $hmm_script) { die "Could not finf HMMcopy R script"; }

print STDERR "Running HMMCOPY via external script as [$command]\n" if DEBUG;

`$command`;

HMMCOPY_SCRIPT

my $JobName = "HMMCOPY_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

my $hmm_convert = join("/",($scriptdir,"convertHMMcopy.pl"));
my $hmm_launch  = join("/",($scriptdir,"run_HMMcopy.g.r"));

$SGEscript =~s/HMMBINDIR_TAG/$hmmbindir/;
$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/LAUNCHER_TAG/$hmm_launch/;
$SGEscript =~s/CONVERTER_TAG/$hmm_convert/;
$SGEscript =~s/HMMCONFIG_TAG/$refconf/;
$SGEscript =~s/RLIB_TAG/$rlibdir/;
$SGEscript =~s/INLIST_TAG/$list/;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/hmm_script$$.pl") or die "Couldn't write to [$Bin/hmm_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

print STDERR "Submitting script\n";
my $SGEResult = `qsub -V -cwd $Bin/hmm_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
