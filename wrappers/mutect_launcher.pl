#!/usr/bin/perl -w

=head2 MuTect Launcher
 
 .bam files will be processed by MuTect tool
 but the list of files (tab-delimited) needs to be formatted properly:

 sampleR_file.bam	sampleT_file.bam	output_vcf

 Will create mutect-produced vcf given full path as output_vcf

 mutect_launcher.pl --list [mylist.txt] --mutect [path to mutect jar]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $mutect, $fasta, $dbsnp, $cosmic);
my $USAGE = "mutect_launcher.pl --list [mylist.txt] --mutect [path to MuTect jar]\n";
my $result = GetOptions("list=s"         => \$list,
                        "mutect=s"       => \$mutect,
                        "ref-faste=s"    => \$fasta,
                        "dbsnp=s"        => \$dbsnp,
                        "cosmic=s"       => \$cosmic);

if (!$list || !$mutect) { die $USAGE; }
# Some PDE-specific defaults:
$fasta  ||= "/.mounts/labs/PDE/data/MutectStrelka/hg19/hg19_random.fa";
$dbsnp  ||= "/.mounts/labs/PDE/data/MutectStrelka/dbSNP/dbsnp_138.hg19.leftAligned.vcf.gz";
$cosmic ||= "/.mounts/labs/PDE/data/MutectStrelka/cosmic/hg19_cosmic_v54_120711.vcf";


my $SGEscript = <<'MUTECT_SCRIPT';
#!/usr/bin/perl
# : Call SNV for paired normal/tumor in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=6G
#$ -o /dev/null
#$ -e /dev/null
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $referenceFasta  = "FASTA_TAG";
my $cosmicVcf       = "COSMIC_TAG";
my $dbSNPvcf        = "DBSNP_TAG";
my $FileList        = "INLIST_TAG";
my $muTect          = "MUTECT_TAG";

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
}
close INPUTFILE;

my @out = split("\t", $Comparison);
if (!@out || @out != 3) {die "MALFORMED INSTRUCTION FILE!";}

map{if (!$_ || !-e $_ || !-s $_) { die "Couldn't read from File\n";}} @out[0..1];
my $id = $1 if ($out[0]=~/([A-Z]{3,4,5}_\d+)/);

my $normalName = $id."_R";
my $tumorName  = $id."_T";

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "java -Xmx3000M -jar ";

my $outdir = $out[2]=~m!(.*/)! ? $1 : "./";

# This is the place for customization:
$command.=    "$muTect --analysis_type MuTect ";
$command.=    "--reference_sequence $referenceFasta ";
$command.=    "--cosmic $cosmicVcf ";
$command.=    "--dbsnp $dbSNPvcf ";
$command.=    "--input_file:normal $out[0] ";
$command.=    "--input_file:tumor $out[1] ";
$command.=    "--only_passing_calls ";
$command.=    "--normal_sample_name $normalName ";
$command.=    "--tumor_sample_name  $tumorName ";
$command.=    "--vcf $out[2]";

#print STDERR $command."\n";

`$command`;

MUTECT_SCRIPT

my $JobName = "MUTECT_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/MUTECT_TAG/$mutect/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/FASTA_TAG/$fasta/;
$SGEscript =~s/DBSNP_TAG/$dbsnp/;
$SGEscript =~s/COSMIC_TAG/$cosmic/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/mutect_script$$.pl") or die "Couldn't write to [$Bin/mutect_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

#my $SGEResult = `qsub $Bin/mutect_script$$.pl`;

#$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
#print "# SGE QSub launch result was:'$SGEResult'\n";
