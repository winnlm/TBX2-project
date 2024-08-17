# TBX2-project

#UNIX
The first steps with CUT and Run is trimming and FASTQ paired files.

The Amik.slurm script has the capability to conduct trimming on all pairs in parralel.

This results in two of the most important files we will use further down the line: BigWig and Bed files.

However, these files are not normalized to account for the signal to noise ratio so we have to use a control to reduce the noise.

This is done in the Amik2.slurm file.

The bed file converts bed files into BigBed for ease of use on University of Santa Cruz's Genome Browser

Bed2Bam is a script that reverse the changes and converts bed files to bam files. This is useful if you want to use averaged bam files inside of the amik2.slurm script.

#RStudio

The Diff_Bind_CSV.Rmd is used to create a differential binding affinity analysis to get the differential enriched or depleted sites of the specified DNA binding protein.

The SCCcrosscompare.Rmd is an unused script that has not been tested beacause the package system wasn't calibrated. The function is to create an upset plot that shows intersections between groups.

The remmap enrich.Rmd is used to plot the significant binding sites of the depletion and enrichment to a catalog genome to see regions of significant intersections.


Contact
email: eleuthexius@gmail.com


