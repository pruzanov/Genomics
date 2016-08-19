#!/usr/bin/perl -w

=head2 Add RG with picard
 
 Will add read groups to bam files

 add_rg.pl --list [mylist.txt] --addrg [path to Picard AddOrReplaceReadGroups.jar]

 mylist.txt whould be a two-column file that looks like this:

 ERNA0051Aligned.bam       ID:150213_D00355_0080_BC5UR0ANXX_CCGTCC_3 SM:ERNA_0051_nn_C LB:ERNA_0051_nn_C_PE_366_WT PU:150213_D00355_0080_BC5UR0ANXX_CCGTCC_3 PI:100 PL:ILLUMINA

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $addrg);
my $USAGE = "add_rg.pl --list [mylist.txt] --addrg [path to Picard AddOrReplaceReadGroups.jar]\n";
my $result = GetOptions("list=s"         => \$list,
                        "addrg=s"        => \$addrg);

if (!$list || !$addrg) { die $USAGE; }

my $SGEscript = <<'PICARD_SCRIPT';
#!/usr/bin/perl
# : Add read groups to bam files in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=6G
#$ -o /dev/null
#$ -e AddRg.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $addRg           = "ADDRG_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my($infile, $rgline);
my $Count=0;
while (<INPUTFILE>) {
        chomp;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }

        if (/^#/) {
          print STDERR "Line $Line skipped, Aborting...\n";
          exit;
        }

        #Below here only processed for wanted lines: 
        ($infile,$rgline) = split("\t",$_);
}
close INPUTFILE;
my($rgid,$rgpl,$rgpi,$rgpu,$rglb,$rgsm);

if (!$infile || !-e $infile || !-s $infile) { die "Couldn't read from File\n";}
$rgline=~s/\:/ /g;
if ($rgline=~/ID (\S+) SM (\S+) LB (\S+) PU (\S+) PI (\S+) PL (\S+)/) {
  ($rgid,$rgsm,$rglb,$rgpu,$rgpi,$rgpl) = ($1,$2,$3,$4,$5,$6);
} else {
  print STDERR "The RG line should have format ID:\S+ SM:\S+ LB:\S+ PU:\S+ PI:\d+ PL:\S+\n";
}
my $libid   = $1 if ($infile=~/([A-Z]{3,5}\d+)/);
my $outfile = $infile;
$outfile=~s/.bam$/.rg.bam/;

# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $command = "java -Xmx3000M -jar $addRg ";

$command.= "INPUT=$infile " ;
$command.= "OUTPUT=$outfile "; 
$command.= "RGID=$rgid ";
$command.= "RGSM=$rgsm ";
$command.= "RGLB=$rglb ";
$command.= "RGPU=$rgpu ";
$command.= "RGCN=OICR ";
$command.= "RGPL=$rgpl ";
$command.= "CREATE_INDEX=true";

`$command`;
 print STDERR $command;

PICARD_SCRIPT

my $JobName = "PICARD_$$";
my $files = `wc -l $list`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/ADDRG_TAG/$addrg/;
$SGEscript =~s/INLIST_TAG/$list/;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/add.readgroups_script$$.pl") or die "Couldn't write to [$Bin/add.readgroups_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/add.readgroups_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
