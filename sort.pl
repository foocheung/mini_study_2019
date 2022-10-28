#!/usr/bin/perl


$pwd=`pwd`;
chomp $pwd;


open BAM, '-|', 'ls bam/*.bam';
@bam=(<BAM>);


foreach $b (@bam){
chomp $b;


$cmd="/sysapps/cluster/software/SAMtools/1.9-goolf-1.7.20/bin/samtools  sort $pwd/$b -o $pwd/$b\.2 ";

print "$cmd\n";
system($cmd);


}
