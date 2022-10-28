#!/usr/bin/perl


$pwd=`pwd`;
chomp $pwd;

while(<>){


chomp;

$cmd="qsub -e $pwd/$_\.log -o $pwd/$_\.qsub -b y \"/sysapps/cluster/software/HTSeq/0.9.1-goolf-1.7.20-Python-2.7.9/bin/htseq-count   $_ $pwd/Homo_sapiens.GRCh38.95.gtf  -f bam -s reverse  >  $_\.count\" ";


print "$cmd\n";
system($cmd);

}
