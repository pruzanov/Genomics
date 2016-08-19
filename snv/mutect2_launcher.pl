#!/usr/bin/perl -w

=head2 MuTect2 Launcher
 
 .bam files will be processed by MuTect2 tool
 but the list of files (tab-delimited) needs to be formatted properly:

 sampleR_file.bam	sampleT_file.bam	output_vcf

 Will create mutect-produced vcf given full path as output_vcf

 Note: This needs to be run from a qrsh session, gatk module needs to be loaded with
 
 module load gatk/3.5.0.nightly.2016.05.02

 mutect2_launcher.pl --list [mylist.txt] 

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir);
my $USAGE = "mutect2_launcher.pl --list [mylist.txt]\n";
my $result = GetOptions("list=s" => \$list);

my $refFai = $ENV{HOME}."/data/MutectStrelka/hg19/hg19_random.fa.fai";
if (!$list) { die $USAGE; }

my $SGEscript = <<'MUTECT_SCRIPT';
#!/usr/bin/perl
# : Call SNV for paired normal/tumor in a list
#$ -t 1-CHROMS_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -o /dev/null
#$ -e Mutect2.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     # This is passed as an environment variable by SGE / Qsub

my $referenceFasta  = $ENV{HOME}."/data/MutectStrelka/hg19/hg19_random.fa";
my $cosmicVcf       = $ENV{HOME}."/data/MutectStrelka/cosmic/hg19_cosmic_v54_120711.vcf";
my $dbSNPvcf        = $ENV{HOME}."/data/MutectStrelka/dbSNP/dbsnp_138.hg19.leftAligned.vcf.gz";
my $ChromList       = "REFFAI_TAG";
my $FileTumor       = "TUMORFILE_TAG";
my $FileNormal      = "NORMALFILE_TAG";
my $FileOut         = "OUTFILE_TAG";

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

my @chrinfo = split("\t",$CHROM); #(":",$CHROM);
my $CHROMID = $chrinfo[0];
my $CHR=$chrinfo[0].":1-".$chrinfo[1];

my $TMPDIR = dirname($FileOut)."/tmpfiles";
if (!-d $TMPDIR) { `mkdir $TMPDIR`; }

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar \$GATKROOT/GenomeAnalysisTK.jar";
$command.= " -T MuTect2";
$command.= " -R $referenceFasta";
$command.= " -I:tumor $FileTumor";
$command.= " -I:normal $FileNormal";
$command.= " --dontUseSoftClippedBases";
$command.= " --dbsnp $dbSNPvcf";
$command.= " --cosmic $cosmicVcf";
$command.= " --intervals $CHR";
$command.= " --max_alt_alleles_in_normal_count 6";
$command.= " -U ALLOW_SEQ_DICT_INCOMPATIBILITY";
$command.= " -o $TMPDIR/mutect.$CHROMID.vcf";

print STDERR $command;
`$command`;

MUTECT_SCRIPT

# Collector Script
my $SGECollector = <<'COLLECT_SCRIPT';
#!/usr/bin/perl
# : Collect vcf files and merge them
#$ -cwd
#$ -l h_vmem=2G
#$ -hold_jid MUTECT_TAG
#$ -o /dev/null
#$ -e VcfMerge.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $outfile = "OUTFILE_TAG";
my $TMPDIR  = dirname($outfile)."/tmpfiles";

opendir(DIR, $TMPDIR) or die "Couldn't read files from $TMPDIR";
my @vcfs = map{$TMPDIR."/".$_} (grep {/.vcf$/} readdir(DIR));
closedir(DIR);

my $command = "vcf-concat ".join(" ",@vcfs)." | vcf-sort > $outfile";

print STDERR $command;
`$command`;

COLLECT_SCRIPT


my $chromfile = $refFai;
my $chroms  = `wc -l $chromfile`;
chomp($chroms);
$chroms =~s/\s+.*//;

$SGEscript =~s/CHROMS_TAG/$chroms/;
$SGEscript =~s/REFFAI_TAG/$chromfile/;

my $gatk     = "gatk/3.5.0.nightly.2016.05.02";
my $vcftools = "vcftools";

# Now, print this script to the working directory and qsub it
open INPUTFILE , "<$list" or die "Cannot open list of files '$list' containing inputs\n";
my $Count = 0;
while (<INPUTFILE>) {
 chomp;
 $Count++;
 
 #Below here only processed for wanted lines: 
 if (/^#/) {
     print STDERR "Skipping Line $Count\n";
     next;
 }

 my($inN, $inT, $outVcf) = split("\t");
 next if (!$inN || !$inT || !$outVcf);

 my $JobName = "MUTECT2_$Count"."_$$";
 my $SGEscript_next = $SGEscript;
 $SGEscript_next =~s/JOBNAME_TAG/$JobName/;
 $SGEscript_next =~s/NORMALFILE_TAG/$inN/;
 $SGEscript_next =~s/TUMORFILE_TAG/$inT/;
 $SGEscript_next =~s/OUTFILE_TAG/$outVcf/;
 my $scriptFile = "$Bin/mutect_script$$.".$Count.".pl";

 open(SCRIPT,">$scriptFile") or die "Couldn't write to [$scriptFile]";
 print SCRIPT $SGEscript_next;
 close SCRIPT;

 # qsub
 print STDERR "Submitting script\n";

 my $SGEResult = `qsub -V $scriptFile`;

 $SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
 print "# SGE QSub launch result was:'$SGEResult'\n";

 my $CJobName = "M2COLL_$Count"."_$$";
 #make aggregation script for merging vcf files
 my $SGECollnext = $SGECollector;
 $SGECollnext=~s/MUTECT_TAG/$JobName/;
 $SGECollnext=~s/OUTFILE_TAG/$outVcf/;
 $SGECollnext=~s/JOBNAME_TAG/$CJobName/;
 
 my $scriptCFile = "$Bin/mutect_collector_script$$.".$Count.".pl";
 open(SCRIPT,">$scriptCFile") or die "Couldn't write to [$scriptCFile]";
 print SCRIPT $SGECollnext;
 close SCRIPT;

 my $SGECResult = `qsub -V $scriptCFile`;
 print "# SGE QSub launch result was:'$SGECResult'\n";
}
close INPUTFILE;




