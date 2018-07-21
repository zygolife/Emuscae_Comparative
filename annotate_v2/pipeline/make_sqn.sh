#!/usr/bin/bash
#SBATCH -p batch --time 4:00:00 --mem 4G --nodes 1 --ntasks 1 --out tbl2asn.log

module load tbl2asn

tbl2asn -p . -y "Annotated using funannotate v0.8.0" -N 1 -t /bigdata/bioinfo/pkgadmin/opt/linux/centos/7.x/x86_64/pkgs/funannotate/git-live/lib/test.sbt -M n -Z ../../predict_results/Rhizopus_azygosporus.discrepency.report.txt -j "[organism=Rhizopus azygosporus]" -V b -c fx -T -a r10u -l paired-ends
