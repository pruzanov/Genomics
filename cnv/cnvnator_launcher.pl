#!/usr/bin/perl -w

=head2 CNVnator pipeline launcher

 will work with lists of bam files
 better to nave bam files in separate directories
 (each bam in an indivdual dir)

 mybam_dir1/mybam_file1.bam
 mybam_dir2/mybam_file2.bam

 will need list of files, will use defaults but may
 accept path to samtools, cnvnator

 This is a cluster version, needs to be run on a headnode capable of submitting sge jobs (qsub)

 need to make sure we have root module

 module load root

 after that which root 

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin);
use constant DEBUG=>1;

my($list, $samtools, $rootdir, $cnvnator, $genome, $bin, $refFai, $q);
my $USAGE = "cnvnator_launcher.pl --list [mylist.txt] --samtools [Optional path to samtools] --cnvnator [Optional path to cnvnator binary] --genome [Optional, i.e. hg19] --bin [binsize in bases] --queue [optional queue]\n";
my $result = GetOptions("list=s"         => \$list,
                        "genome=s"       => \$genome,
                        "bin=i"          => \$bin,
                        "ref-fai|fai=s"  => \$refFai,
                        "queue|q=s"      => \$q,
                        "cnvnator=s"     => \$cnvnator,
                        "samtools=s"     => \$samtools);

# Defaults
$refFai   ||= "$HOME/data/reference/hg19_random/fasta/UCSC/hg19_random.fa.fai";
$samtools ||= $ENV{HOME}."/bin/samtools";
$rootdir  = `which root`;
$rootdir  =~s!bin/root!!;
chomp($rootdir);
$cnvnator ||= $ENV{HOME}."/bin/cnvnator";
my $contig_dir=dirname($refFai);
$genome   ||= "hg19";


if (!$rootdir || !-d $rootdir) {die "root needs to be loaded with \"module load root\" to run this analysis\n";}

# ======================================
# Make sure environment variables r set
# ======================================

$ENV{ROOTSYS} = $rootdir;
unless ($ENV{LD_LIBRARY_PATH} =~/$rootdir/) {my $ld = join(":",($rootdir."/lib",$ENV{LD_LIBRARY_PATH})); $ENV{LD_LIBRARY_PATH}=$ld;}
if (!$list || !$bin) { die $USAGE; }

my $SGEscript = <<'CNVNATOR_SCRIPT';
#!/usr/bin/perl
# : Call germline CNVs for bams in a list
#$ -t 1-CHROMS_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -o /dev/null
#$ -e CNVnator.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
use constant DEBUG=>0;

my $Line=$ENV{"SGE_TASK_ID"};     # This is passed as an environment variable by SGE / Qsub

my $ChromList       = "REFFAI_TAG";
my $genome          = "GENOME_TAG";
my $bin             = "BIN_TAG";
my $file            = "INFILE_TAG";
my $cnvnator        = "CNVNATOR_TAG";

my $CHROM;

open INPUTFILE , "$ChromList" or die "Cannot open list of files '$ChromList' containing inputs\n";
my $Comparison;

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }
        
        $CHROM = $_;
        chomp($CHROM);
}
close INPUTFILE;

exit if $CHROM=~/M/ || $CHROM=~/_/;

my @chrinfo = split("\t",$CHROM); #(":",$CHROM);
my $chrom = $chrinfo[0];

my $TMPDIR = dirname($file)."/tmpfiles/";
if (!-d $TMPDIR) { `mkdir $TMPDIR`; }

print STDERR "Starting to work on ".basename($file)."\n";
print STDERR "Working on $chrom...\n" if DEBUG;
chomp($chrom);

my $callfile = $TMPDIR.basename($file);
$callfile=~s/bam$/$chrom.calls/;
my $rootfile = $TMPDIR.basename($file);
$rootfile=~s/bam$/$chrom.root/;
my $contig_dir = dirname($ChromList);

`$cnvnator -root $rootfile -chrom $chrom -tree $file`;
 print STDERR "Produced root file for $chrom...\n" if DEBUG;

`$cnvnator -genome $genome -root $rootfile -chrom $chrom -his $bin -d $contig_dir`;
 print STDERR "Produced hist data for $chrom...\n" if DEBUG;

`$cnvnator -genome $genome -root $rootfile -chrom $chrom -stat $bin`;
 print STDERR "Produced stat data for $chrom...\n" if DEBUG;

`$cnvnator -genome $genome -root $rootfile -chrom $chrom -partition $bin`;
 print STDERR "Partitioned read densities for $chrom...\n" if DEBUG;

`$cnvnator -root $rootfile -chrom $chrom -call $bin > $callfile`;

print STDERR "Completed CNV analysis for $chrom, calls written into [$callfile]...\n" if DEBUG;

CNVNATOR_SCRIPT

# Collector Script
my $SGECollector = <<'COLLECT_SCRIPT';
#!/usr/bin/perl
# : Collect CNVnator results concat them
#$ -cwd
#$ -l h_vmem=2G
#$ -hold_jid CNVNATOR_TAG
#$ -o /dev/null
#$ -e CNVcat.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
use constant DEBUG=>0;

my $final_calls = "FILE_TAG";
my $filedir = dirname($final_calls)."/tmpfiles";
$final_calls =~s/bam$/final.calls/;


opendir(DIR, $filedir) or die "Couldn't read files from $filedir";
my @final_calls = map{$filedir."/".$_} (grep {/.calls$/} readdir(DIR));
closedir(DIR);

my $tempfile = join("/",($filedir,"temp_$$"));
my $to_cat;
map {if(-e $_ && (stat($_))[7] > 0){$to_cat.=" $_"}} @final_calls;

print STDERR "List of files to merge: [$to_cat]\n" if DEBUG;

if ($to_cat=~/\S+/) { # list is not empty
 my $command = "cat $to_cat > $tempfile && mv $tempfile $final_calls";
 print STDERR $command;
 `$command`;
} else {
 print STDERR "Analysis did not produce any calls\n";
}

COLLECT_SCRIPT


my $chromfile = $refFai;
my $chroms  = `wc -l $chromfile`;
chomp($chroms);
$chroms =~s/\s+.*//;

$SGEscript =~s/CHROMS_TAG/$chroms/;
$SGEscript =~s/REFFAI_TAG/$chromfile/;
$SGEscript =~s/GENOME_TAG/$genome/;
$SGEscript =~s/BIN_TAG/$bin/;
$SGEscript =~s/CNVNATOR_TAG/$cnvnator/;

# ==========================================
# For each file in the list run the pipeline
# ==========================================
my $Count = 0;
open(INLIST,"<$list") or die "Cannot read from the list with files";
while(<INLIST>) {
 chomp;
 $Count++;
 next if /^#/; # skip commented out files
 
 if (!-e $_) {
     print STDERR "File $_ does not exist, skipping\n";
     next;
 }

 my $file = $_;

 my $JobName = "CNVNTR_$Count"."_$$";
 my $SGEscript_next = $SGEscript;
 $SGEscript_next =~s/JOBNAME_TAG/$JobName/;
 $SGEscript_next =~s/INFILE_TAG/$file/;

 my $scriptFile = "$Bin/cnvnator_script$$.".$Count.".pl";

 open(SCRIPT,">$scriptFile") or die "Couldn't write to [$scriptFile]";
 print SCRIPT $SGEscript_next;
 close SCRIPT;

 # qsub main script (array job)
 print STDERR "Submitting script for [$file]\n";

 my $SGEResult = $q ? `qsub -V -q $q $scriptFile` : `qsub -V $scriptFile`;

 $SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
 print "# SGE QSub launch result was:'$SGEResult'\n";

 # qsub collector
 my $CJobName = "CNVCOLL_$Count"."_$$";
 #make aggregation script for merging vcf files
 my $SGECollnext = $SGECollector;
 $SGECollnext=~s/FILE_TAG/$file/;
 $SGECollnext=~s/CNVNATOR_TAG/$JobName/;
 $SGECollnext=~s/JOBNAME_TAG/$CJobName/;
 
 my $scriptCFile = "$Bin/cnvnator_collector_script$$.".$Count.".pl";
 open(SCRIPT,">$scriptCFile") or die "Couldn't write to [$scriptCFile]";
 print SCRIPT $SGECollnext;
 close SCRIPT;

 print STDERR "Submitting collector script\n";
 
 my $SGECResult = `qsub -V $scriptCFile`;
 print "# SGE QSub launch result was:'$SGECResult'\n";
}

close INLIST;
