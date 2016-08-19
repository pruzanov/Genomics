#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use constant DEBUG=>1;

=head2 Running HMMcopy germline (normal only) analysis

 HMMcopy is a tool for analyzing copy-number variation
 This script launches HMMcopy analysis using wig file 
 prepared with readCounter binary

 also needs some pre-generated support data

=cut

my $USAGE = "launchHMMcopy.g.pl --r-libdir [Path to directory with R modules] --wig [input wig] --cg-file [cg file] --map-file [map file] --output-base [basename for output]";

# Required parameters
my($hmm_script,$rlibdir,$wig,$cg_file,$map_file,$output_base);
my $results = GetOptions( "r-libdir=s"      =>  \$rlibdir,
                          "wig=s"           =>  \$wig,
                          "cg-file=s"       =>  \$cg_file,
                          "map-file=s"      =>  \$map_file,
                          "hmm-script=s"    =>  \$hmm_script,
                          "output-base=s"   =>  \$output_base);

if (!$wig || !$rlibdir || !$cg_file || !$map_file || !$output_base) { die $USAGE; }
$hmm_script ||="$Bin/run_HMMcopy.g.r";

unless (-e $hmm_script) { die "Cannot find R script for HMM cpy analysis"; }

# Make output directory if it does not exists                            
my $outdir = dirname($output_base);
`mkdir -p $outdir` if (!-e $outdir || !-d $outdir);

$ENV{R_LIBS} = $rlibdir;

#=====================================
# RUNNING HMMcopy script here
#=====================================
print STDERR "Running Rscript $hmm_script $wig $cg_file $map_file $output_base\n" if DEBUG;
`Rscript $hmm_script $wig $cg_file $map_file $output_base`;
