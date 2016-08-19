#!/usr/bin/perl -w

=head2 HaplotypeCaller Launcher
 
 .bam files will be processed by iGATK's HaplotypeCaller tool
 but the list of files (tab-delimited) needs to be formatted properly:

 sampleR_file.bam	output_vcf

 Will create gatk-produced vcf given full path as output_vcf

 hcaller_launcher.pl --list [mylist.txt] --gatk [path to GATK jar]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $gatk);
my $USAGE = "hcaller_launcher.pl --list [mylist.txt] --gatk [path to GATK jar]\n";
my $result = GetOptions("list=s"         => \$list,
                        "gatk=s"         => \$gatk);

if (!$list || !$gatk) { die $USAGE; }

my $SGEscript = <<'GATK_SCRIPT';
#!/usr/bin/perl
# : Call SNV for tumor .bam in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=8G
#$ -o /dev/null
#$ -e HCaller.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $referenceFasta  = $ENV{HOME}."/data/MutectStrelka/hg19/hg19_random.fa";
my $cosmicVcf       = $ENV{HOME}."/data/MutectStrelka/cosmic/hg19_cosmic_v54_120711.vcf";
my $dbSNPvcf        = $ENV{HOME}."/data/MutectStrelka/dbSNP/dbsnp_142.CAF.chr.vcf.gz";
my $java            = $ENV{HOME}."/bin/java";
my $FileList        = "INLIST_TAG";
my $gatk            = "GATK_TAG";

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

my $bam = $Comparison;

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "$java -Xmx4000M -jar ";

my $vcf = $bam;
$vcf =~s/\.bam$/.hc.vcf/;

$command.=    "$gatk -T HaplotypeCaller ";
$command.=    "--reference_sequence $referenceFasta ";
$command.=    "--genotyping_mode DISCOVERY ";
$command.=    "--dbsnp $dbSNPvcf ";
$command.=    "-I $bam ";
$command.=    "-stand_call_conf 20 ";
$command.=    "-stand_emit_conf 10 ";
$command.=    "--filter_reads_with_N_cigar ";
$command.=    "-out_mode EMIT_VARIANTS_ONLY ";
$command.=    "-o ".$vcf;

`$command`;

GATK_SCRIPT

my $JobName = "HCALL_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/GATK_TAG/$gatk/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/gatk_script$$.pl") or die "Couldn't write to [$Bin/gatk_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/gatk_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
