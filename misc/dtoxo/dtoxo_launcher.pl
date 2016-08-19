#!/usr/bin/perl -w

=head2 Run DToxoG on pre-compiled MAF files
 
 Run DToxoG
 but the list of files (tab-delimited) needs to be formatted properly:

 sample_file.maf	output_dir
 ...

 Note: This needs to be run from a qrsh session, matlab module needs to be loaded with
 
 module load matlab

 dtoxo_launcher.pl --list [mylist.txt] --dtoxo-dir [path to DToxo-G dir containing script]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $dtoxo);
my $USAGE = "dtoxo_launcher.pl --list [mylist.txt] --dtoxo-dir [path to dtoxo-g directory containing script]\n";
my $result = GetOptions("list=s"         => \$list,
                        "dtoxo-dir=s"    => \$dtoxo);

if (!$list || !$dtoxo) { die $USAGE; }

my $SGEscript = <<'DTOXO_SCRIPT';
#!/usr/bin/perl
# : Run DToxoG
#$ -cwd
#$ -l h_vmem=12G
#$ -o /dev/null
#$ -hold_jid MAKEMAF_5466
#$ -e Dtoxo.e
#$ -V
#$ -S /bin/bash

#$ -N JOBNAME_TAG

use strict;
use File::Basename;
use FindBin qw($Bin);

my $refFasta  = "/.mounts/labs/PDE/data/MutectStrelka/hg19/hg19_random.fa";
my $FileList  = "INLIST_TAG";
my $dtoxo_dir = "DTOXO_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $inputLine;
my $Count = 1;

while (<INPUTFILE>) {
        chomp;
        #Below here only processed for wanted lines: 
        $inputLine = $_;
        if ($inputLine=~/^#/) {
          print STDERR "Skipping Line $Count\n";
          $Count++;
          next;
        }
  chomp($inputLine);

  my($inputMaf,$outdir) = split("\t",$inputLine);
  if (!-e $inputMaf || !-d $outdir) { die "Couldn't find inputs"; }
  my($filename, $dirs, $suffix) = fileparse($inputMaf);
  $filename =~s/.maf$/.parsed.maf/;

  # Let's assume that we have valid .bam files at this point (may introduce a check later, if needed

  my $command.= "matlab -nodesktop -nosplash -r ";
  $command.= " \"addpath $dtoxo_dir;";
  $command.= "startFilterMAFFile(\'$inputMaf\', \'$filename\', \'$outdir\',0,0,\'0.96\');";
  $command.= "exit\"";

  print STDERR $command."\n";
  system($command);


  $Count++;
 }
 close INPUTFILE;

DTOXO_SCRIPT

my $JobName = "DTOX_$$";

$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/DTOXO_TAG/$dtoxo/;
$SGEscript =~s/INLIST_TAG/$list/;

# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/dtoxo_script$$.pl") or die "Couldn't write to [$Bin/dtoxo_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub -V $Bin/dtoxo_script$$.pl`;


$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
