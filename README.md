# Genomics

Scripts for running various tools, data conversion and visualization
The scripts are simple wrappers for common NGS-processing tools

all of the scripts work with lists of inputs and launch sge jobs using qsub command.
The jobs are submitted as array jobs which simplifies job mainainance

# File directory structure

* converters - file manipulation and producing certain formats
* cnv - Copy Number Variation analysis
* snv - SNV and indel analysis
* alignment - wrappers for alignment tools
* misc - anything else
