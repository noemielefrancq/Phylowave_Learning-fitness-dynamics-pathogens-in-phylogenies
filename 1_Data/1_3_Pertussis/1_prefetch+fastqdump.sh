#########################################################
## Download data from sra
## i must be an accession number
#########################################################

#!/bin/sh

~/softwares/sratoolkit.2.11.0-ubuntu64/bin/fastq-dump.2.11.0 --split-files ${i};
~/softwares/DSRC-master/bin/dsrc c -t1 ${i}_1.fastq ${i}_1.fastq.dsrc;
~/softwares/DSRC-master/bin/dsrc c -t1 ${i}_2.fastq ${i}_2.fastq.dsrc;
rm ${i}_1.fastq  ${i}_2.fastq 

