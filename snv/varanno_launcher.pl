#!/usr/bin/perl -w

=head2 VariantAnnotator Launcher
 
 USAGE:

 varanno_launcher.pl --list [mylist.txt] --gatk [path to GATK jar] --tabix [for optional indexing, if tabix indexes need to be produced]

 .bam and .vcf files will be processed by GATK's VariantAnnotator walker
 but the list of files (tab-delimited) needs to be formatted properly:

 sample_file.bam,(...bam)	variants.vcf	output_vcf

 first column can be a of files (ideally, one tumor and one normal bam)

 Will create gatk-produced .annotated.vcf given full path as output_vcf
 if output_vcf is not present, will create annotated vcf in the same 
 directory where the input vcf is, adding suffix va_annotated

 varanno_launcher.pl --list [mylist.txt] --gatk [path to GATK jar] --tabix [for optional indexing, if tabix indexes need to be produced]
 List of metrics is hard-coded, may be modified below

 in order for this to work we need java and Gatk available

 (module load needs to be used with appropriate modules for required tools above)
 
 modules vcan be loaded on a headnode, the script should capture the environment and the same
 modules will be available on SGE node(s) once the jobs get submitted.

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);
use constant DEBUG=>0;

my($list, $outdir, $gatk, $tabix);
my $USAGE = "varanno_launcher.pl --list [mylist.txt] --gatk [path to GATK jar]\n";
my $result = GetOptions("list=s"         => \$list,
                        "tabix=s"        => \$tabix,
                        "gatk=s"         => \$gatk);

$tabix ||=$ENV{HOME}."/bin/tabix";
if (!$list || !$gatk) { die $USAGE; }

my $SGEscript = <<'GATK_SCRIPT';
#!/usr/bin/perl
# : Annotate Variants, add additional metrics
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=12G
#$ -o /dev/null
#$ -e VarAnno.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $referenceFasta  = "/.mounts/labs/PDE/data/MutectStrelka/hg19/hg19_random.fa";
my $cosmicVcf       = "/.mounts/labs/PDE/data/MutectStrelka/cosmic/hg19_cosmic_v54_120711.vcf";
my $dbSNPvcf        = "/.mounts/labs/PDE/data/MutectStrelka/dbSNP/dbsnp_142.CAF.chr.vcf.gz";
my $java            = "java"; # This may be loaded as module, Java 7 recommended
my $FileList        = "INLIST_TAG";
my $gatk            = "GATK_TAG";
my $tabix           = "TABIX_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $Comparison;

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        #Below here only processed for wanted lines: 
        $Comparison = $_;
        if ($Comparison=~/^#/) {
          print STDERR "Skipping Line $Line\n";
          exit;
        }
}
close INPUTFILE;

my($bam, $input_vcf, $vcf) = split("\t",$Comparison);
my @bams = split(",",$bam);
my $vfile = basename($input_vcf);
my $sample = $1 if $vfile=~/SWID_\d+_(\S+?_\d+)_/ || $vfile=~/^(\S+?_\d+)_/;
my $idx = $input_vcf.".tbi";

if (! -e $idx) {`$tabix -p vcf $input_vcf`;}

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
# ============= QC Metrics =======================
my @metrics = ("StrandOddsRatio","AS_FisherStrand","AS_StrandOddsRatio","AS_BaseQualityRankSumTest","AS_MappingQualityRankSumTest");

my $command = "$java -Xmx10000M -jar ";
if (!$vcf) {
 $vcf  = $input_vcf;
 $vcf =~s/\.vcf.*/.va_annotated.vcf/;
}

$command.=    "$gatk -T VariantAnnotator ";
$command.=    "--variant:$sample $input_vcf ";
$command.=    "-R $referenceFasta ";
$command.=    "--dbsnp $dbSNPvcf ";
foreach my $bam(@bams) {
 $command.=    "-I $bam ";
}
foreach my $m(@metrics) {
 $command.=    "-A $m ";
}
$command.=    "-o ".$vcf;

`$command`;

GATK_SCRIPT

my $JobName = "VARA_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/GATK_TAG/$gatk/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/TABIX_TAG/$tabix/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/varanno_script$$.pl") or die "Couldn't write to [$Bin/varanno_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/varanno_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
