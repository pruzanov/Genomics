#!/usr/bin/perl -w

# TODO : this is still WIP

=head2 StringTie Launcher  
 
 .bam files will be processed by StringTie tool
 but the list of files need to be formatted properly:

 sample1_file.bam	output_dir

 Will create strigtie output in output_dir

 stringtie_launcher.pl --list [mylist.txt] --stringtie [path to stringtie] --ref-gtf [transcripts.gtf]

=cut
 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my($list, $outdir, $stringtie, $refgtf, $threads);
my $USAGE = "stringtie_launcher.pl --list [mylist.txt] --stringtie [path to stringtie] --ref-gtf [transcripts.gtf]\n";
my $result = GetOptions("list=s"         => \$list,
                        "stringtie=s"    => \$stringtie,
                        "threads=i"      => \$threads,
                        "ref-gtf=s"      => \$refgtf);

if (!$list || !$refgtf || !$stringtie) { die $USAGE; }
$threads ||=2; # Default number of threads

my $SGEscript = <<'STRINGTIE_SCRIPT';
#!/usr/bin/perl
# : Analyze expression (transcript abundance) for files in a list
#$ -t 1-FILES_TAG
#$ -cwd
#$ -l h_vmem=24G
#$ -o /dev/null
#$ -pe smp THREADS_TAG
#$ -e Stringtie.e
#$ -S /usr/bin/perl

#$ -N JOBNAME_TAG

use strict;
use File::Basename;

my $Line=$ENV{"SGE_TASK_ID"};     #This is passed as an environment variable by SGE / Qsub

my $FileList        = "INLIST_TAG";
my $stringtie       = "STRINGTIE_TAG";
my $refGtf          = "REFGTF_TAG";
my $cores           = "THREADS_TAG";

open INPUTFILE , "$FileList" or die "Cannot open list of files '$FileList' containing inputs\n";
my $Comparison;

my $Count=0;
while (<INPUTFILE>) {
        chomp;
        next if /^#/;
        $Count++;  #Increment the line counter
        if ($Count < $Line)        {       next;   }       #Skip until the lines we want:
        if ($Count > $Line )                  {               last;   }
        
        #Below here only processed for wanted lines: 
        $Comparison = $_;
}
close INPUTFILE;

my @out = split("\t", $Comparison);
if (!@out || @out != 2) {die "MALFORMED INSTRUCTION FILE!";}

my $ResultOutputDir ||= $out[1];
if ( ! -d $ResultOutputDir) {`mkdir -p $ResultOutputDir`;}

if (!$out[0] || !-e $out[0] || !-s $out[0]) { die "Couldn't read from File\n";}

my $label ="";
my $threads = THREADS_TAG;
# Let's assume that we have valid .bam files at this point (may introduce a check later, if needed
my $outFile = basename($out[0]).".gtf";

`$stringtie $out[0] -G $refGtf -o $ResultOutputDir/$outFile -p $cores --fr -l $label -B`;

STRINGTIE_SCRIPT

my $JobName = "STRING_$$";
my $files = `grep -v ^# $list | wc -l`;
chomp($files);
$files =~s/\s+.*//;

$SGEscript =~s/FILES_TAG/$files/;
$SGEscript =~s/JOBNAME_TAG/$JobName/;
$SGEscript =~s/STRINGTIE_TAG/$stringtie/;
$SGEscript =~s/REFGTF_TAG/$refgtf/;
$SGEscript =~s/INLIST_TAG/$list/;
$SGEscript =~s/THREADS_TAG/$threads/g;


# Now, print this script to the working directory and qsub it

open(SCRIPT,">$Bin/stringtie_script$$.pl") or die "Couldn't write to [$Bin/stringtie_script$$.pl]";
print SCRIPT $SGEscript;
close SCRIPT;

# qsub
print STDERR "Submitting script\n";

my $SGEResult = `qsub $Bin/stringtie_script$$.pl`;

$SGEResult =~ s/[\n\s]+$//g;    $SGEResult =~ s/[\r\n]/\n#:  /g;
print "# SGE QSub launch result was:'$SGEResult'\n";
