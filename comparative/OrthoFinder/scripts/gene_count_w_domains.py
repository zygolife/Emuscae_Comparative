#!/usr/bin/env python3

import os, sys, re, argparse, csv

parser = argparse.ArgumentParser(description='Summarize OG family counts and add Pfam and Domain columns.')

parser.add_argument('--domains', default='../domains', help='Domains folder')
parser.add_argument('-o','--out', default='Orthologs_Summary.tsv', help='Output summary table file')
parser.add_argument('--orthogroups', default='OrthoFinder_Results/Orthogroups/Orthogroups.tsv',
                    help='Orthogroup folder')
args = parser.parse_args()

domains = {}
domain_types = set()

# function for adding to domain dictionary db
def add_domain_obs(domain_db, gene_name, domain_type, domain_name):
    if gene_name not in domain_db:
        # init the essential nested dictionaries if
        # this is the first time
        domains[gene_name] = {domain_type: {domain_name: 1}}
    elif domain_type not in domain_db[gene_name]:
        # domain type (eg Pfam) not seen
        # so create it and also
        domain_db[gene_name][domain_type] = {domain_name: 1}
    elif domain_name not in domains[gene_name][domain_type]:
        domain_db[gene_name][domain_type][domain_name] = 1
    else:
        # add if it already exists
        domain_db[gene_name][domain_type][domain_name] += 1
#    print("updated db so that {} -> {}:{}".format(gene_name,domain_type,domain_name))

for dtype in os.listdir(args.domains):
    dpath=os.path.join(args.domains,dtype)
    if dtype == 'Pfam':
        for dfile in os.listdir(dpath):
            if dfile.endswith('.domtbl'):
                with open(os.path.join(dpath,dfile),"r") as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        row = line.split()
                        gname = row[0]
                        domainname = row[3]
                        domainacc  = row[4]
                        add_domain_obs(domains,gname,dtype,domainname)
    elif dtype == "CAZY":
        for dfile in os.listdir(dpath):
            if dfile.endswith('.run_dbcan'):
                hotpepfile = os.path.join(dpath,dfile,'overview.txt')
                if os.path.exists(hotpepfile):
                    with open(hotpepfile,"r") as fh:
                        for line in fh:
                            if line.startswith('Gene ID'):
                                continue

                            row = line.strip().split("\t")

                            if row[-1] == "3": # only handle hotpep full hits
                                gname = row[0]
                                hotpep = row[2]
                                add_domain_obs(domains,gname,dtype,hotpep)

    elif dtype == "MEROPS":
        for dfile in os.listdir(dpath):
            if dfile.endswith('.blasttab'):
                with open(os.path.join(dpath,dfile),"r") as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        row = line.split()
                        gname = row[0]
                        domainname = row[1]
                        add_domain_obs(domains,gname,dtype,domainname)
    else:
        print("did not process domain dir {}".format(dtype))
        continue
    domain_types.add(dtype)

print(len(domains),"proteins scored")

with open(args.orthogroups,"r") as fh, open(args.out,"wt") as ofh:
    writer = csv.writer(ofh,delimiter="\t")
    rdr    = csv.reader(fh,delimiter="\t")
    header = next(rdr)
    header.append('Total')
    header.extend(sorted(domain_types))
    writer.writerow(header)
    for row in rdr:
        outrow = [row[0]]
        domain_counts = {}
        total_genes = 0
        # process all the genes in a row by group
        for n in range(1,len(row)):
            if len(row[n]) == 0:
                outrow.append(0)
                continue
            genes = row[n].split(", ")
            total_genes += len(genes)
            outrow.append(len(genes))
            for gene in genes:
                if gene in domains:
                    for dtype in domains[gene]:
                        if dtype not in domain_counts:
                            domain_counts[dtype] = {}
                        for domain,dcount in domains[gene][dtype].items():
                            if domain in domain_counts[dtype]:
                                domain_counts[dtype][domain] += dcount
                            else:
                                domain_counts[dtype][domain] = dcount
        outrow.append(total_genes)
        for dtype in sorted(domain_types):
            dvalues = []
            if dtype in domain_counts:
                for dname,dcount in sorted(domain_counts[dtype].items(),
                                           key=lambda item: item[1]):
                    dvalues.append("{}:{}".format(dname,dcount))
            outrow.append(",".join(dvalues))
        writer.writerow(outrow)
