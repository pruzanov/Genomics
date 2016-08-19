#!/usr/bin/perl -w

=head2 BICseq2 launcher
 
 process bam data (controls ONLY) using BICseq2

 mybam_file1.bam	my_output_dir1
 mybam_file2.bam        my_output_dir2

 bicseq_launcher.pl --list [mylist.txt] --bicseq-norm [dir for bicseq-norm part] --bicseq-seq [dir for bicseq-seg part] --chrinfo [tmplate config for bicseq] --samtools samtools-U

 special samtools needed

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $bicseq_norm, $bicseq_seg, $samtools, $chrinfo, $tmpdir);
my $USAGE = "bicseq_launcher.pl --list [mylist.txt] --bicseq-norm [dir for bicseq-norm part] --bicseq-seq [dir for bicseq-seg part] --chrinfo [chrom info file] --samtools samtools-U";
my $result = GetOptions("list=s"             => \$list,
                        "chrominfo=s"        => \$chrinfo, # file with information on fasta files and 
                        "samtools=s"         => \$samtools,
                        "temp-dir|temp=s"    => \$tmpdir,
                        "bicseq-norm=s"      => \$bicseq_norm,
                        "bicseq-seg=s"       => \$bicseq_seg);

$chrinfo ||="$Bin/chrominfo.txt";
if (!$list || !$bicseq_norm || !$bicseq_seg || ! -e $chrinfo || !$samtools) { die $USAGE; }
$tmpdir ||="$HOME/Results/BicSeq/tmp";

my $SGEscript = <<'BICSEQ_SCRIPT';
#!/usr/bin/perl
# : Run BiqSeq2 pipeline for germline CNV analysis
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=8G
#$ -o /dev/null
#$ -e BICseq.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
use constant DEBUG=>1;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $bicnorm         = "BICNORM_TAG";
my $bicseg          = "BICSEG_TAG";
my $infofile        = "INFO_TAG";
my $samtools        = "SAMTOOLS_TAG";
my $temp            = "TEMP_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($bam, $outdir, $bicseq_dir);

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
        $bicseq_dir=$outdir."/bicseq";
        $bicseq_dir=~s!//!/!g;
        if (-d $bicseq_dir) {
            print STDERR "BicSeq analysis has been run for $bam, wont't continue\n";
            exit;
        } 

        `mkdir $bicseq_dir`;
        $outdir = $bicseq_dir;
}
close INPUTFILE;

# Construct prefix
my $prefix = basename($bam);
# OICR specific
$prefix =~s/SWID_\d+_//;        # remove possible SWID part
$prefix = $1 if $prefix=~/^(\S+?_\d+)/; # extract donor id
print STDERR "Prefix defined a $prefix\n" if DEBUG;
$prefix.=".";

if ($outdir !~m!/$!) {$outdir.="/";}
my $data_dir = $outdir."data";
if (!-d $data_dir) { `mkdir $data_dir`; }

# 1. Run Uniquely mapped reads analysis with custom samtools
# 2. Run Biq-Seq2-norm
# 3. Run Biq-Seq2-seg

my $command = "$samtools view -U BWA,$data_dir/$prefix,N,N $bam";
print STDERR "Running command [$command]\n" if DEBUG;
`$command`;

# Make the config file for BICseq2-norm
my $bicnorm_config = $outdir."bicnorm.conf";
open(CONF,">$bicnorm_config") or die "Cannot write to config file [$bicnorm_config]";
open(INFO,"<$infofile") or die "Couldn't read from info file [$infofile]";

my $first = <INFO>;
print CONF $first;

while(<INFO>) {
  my @temp = split("\t");
  $temp[3] = "$outdir/data/$prefix$temp[3]";
  $temp[4] = "$outdir/data/$prefix$temp[4]";
  if (-e $temp[3] && -s $temp[3]) {
     print CONF join("\t",@temp);
  } else { 
     print STDERR "File $temp[3] does not exist\n" if DEBUG;
  }
}
close CONF;
close INFO;

my $bicseg_config = $outdir."bicseg.conf";
`cut -f 1,5 $bicnorm_config | grep -v chrM > $bicseg_config`;
print STDERR "Made config files\n" if DEBUG;

# chromName  	faFile 	MapFile 	readPosFile 	binFileNorm
# chr1 	chr1.fa 	hg18.CRC.50mer.chr1.txt 	chr1.seq 	chr1.norm.bin
# chr2 	chr2.fa 	hg18.CRC.50mer.chr2.txt 	chr2.seq 	chr2.norm.bin

my $gam_file = $bicnorm_config;
$gam_file =~s/.conf/.gam/;
$command = "$bicnorm -l=100 -s=20 $bicnorm_config $gam_file --tmp $temp"; 
print STDERR "Running BicNorm as [$command]...\n" if DEBUG;
`$command`;

# Make the config file for BICseq2-seg

# chromName 	binFileNorm
# chr1 	        chr1.norm.bin
# chr2 	        chr2.norm.bin

my $final_output = $outdir.$prefix."final.out";
$command = "$bicseg $bicseg_config $final_output --lambda 0.5 --tmp $temp";
print STDERR "Running BicSeq as [$command]...\n" if DEBUG;
`$command`;

BICSEQ_SCRIPT

my $JobName = "BIC_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/BICNORM_TAG/$bicseq_norm/;
$SGEscript =~s/BICSEG_TAG/$bicseq_seg/;
$SGEscript =~s/INFO_TAG/$chrinfo/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/SAMTOOLS_TAG/$samtools/;
$SGEscript =~s/TEMP_TAG/$tmpdir/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/bicseq_script$$.pl") or die "Couldn't write to [$Bin/bicseq_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/bicseq_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
