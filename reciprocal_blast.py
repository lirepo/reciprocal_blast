#!/usr/bin/python2.7
import os
import time
import re
from Bio import SeqIO

def populate_dbs(orgs, data_dir, blast_bin_path, db_dir):
    """Given a list of genome accession numbers, construct the databases
    required for BLAST search."""
    file_ext = '.faa'
    for org in orgs:
        print "Constructing database for", org
        os.system(os.path.join(blast_bin_path, 'makeblastdb') + 
                  ' -in %s -title %s -dbtype prot -out %s -logfile %s' %
                  (os.path.join(data_dir, org+file_ext), # faa file
                   org,                                  # accession number
                   os.path.join(db_dir, org),            # db output file
                   os.path.join(db_dir, org+'.log')))    # logfile
    
def reciprocal_blasts(orgs):
    """Compute reciprocal blast hits for orgs"""
    results_contents = os.listdir(BLAST_RESULTS_PATH)
    file_ext = '.faa'
    for org1 in orgs:
        for org2 in orgs:
            print "starting on: ", org1, org2, "at", time.ctime()
            out_file = "results_%s_%s.txt" % (org1, org2)
            if org1 == org2:
                print "skipping", org1, org2
                continue            
            full_fasta_file = os.path.join(DATA_PATH, org1+file_ext)
            full_db_path = os.path.join(DB_PATH, org2)
            if out_file in results_contents:
                print "skipping", org1, org2
                continue
            os.system(os.path.join(BLAST_BIN_PATH, 
                      ('blastp -query %s -task blastp -db %s -out %s -outfmt 5'
                       % (full_fasta_file, full_db_path, 
                          os.path.join("blast_results", out_file)))))
            print "finished", org1, org2, "at", time.ctime()

def head(xs):
    return xs[0]

def choose2(xs):
    """Return list of choose(xs,2) pairs, retaining ordering on xs"""
    return [(x1,x2) for i,x1 in enumerate(xs) for x2 in xs[i+1:]]

def collate_reciprocal_blasts(orgs, data_dir, blast_results_dir, reciprocal_results_dir):
    # Create locus tag dictionary of the form {genome accession: {gi:locus tag}}
    locus_tag_dicts = {}
    for org in orgs:
        genome_file = os.path.join(data_dir, org+'.gbk')
        locus_tag_dicts[org] = make_locus_tag_dict(genome_file)
    annotation_dicts = make_annotation_dicts(data_dir, orgs)
        
    for i, (org1,org2) in enumerate(choose2(orgs)): 
        outfile = os.path.join(reciprocal_results_dir,
                               'reciprocals_%s_%s.csv' % (org1, org2))
        print "collating reciprocals for %s,%s at %s" % (org1, org2, time.ctime())
        filename1 = os.path.join(blast_results_dir, "results_%s_%s.txt" % (org1, org2))
        filename2 = os.path.join(blast_results_dir, "results_%s_%s.txt" % (org2, org1))
        print "parsing results for %s, %s" % (org1, org2)
        print "using:", filename1, filename2
        hits1 = parse_results(filename1)
        print "parsing results for %s, %s" % (org2, org1)
        hits2 = parse_results(filename2)
        print "finding reciprocals"
        reciprocals = find_reciprocals(hits1, hits2)

        # write reciprocals
        with open(outfile, 'w') as f:
            f.write('\t'.join(['org1_gi', 'org2_gi', 'org1_locus_tag', 'org2_locus_tag',
                               'org1_annotation', 'org2_annotation']))
            f.write('\n')
            for gi1, gi2 in reciprocals:
                locus_tag_1 = locus_tag_dicts[org1][gi1]
                locus_tag_2 = locus_tag_dicts[org2][gi2]
                annotation_1 = annotation_dicts[org1][locus_tag_1]
                annotation_2 = annotation_dicts[org2][locus_tag_2]
                f.write('\t'.join((gi1, gi2, locus_tag_1, locus_tag_2,
                                   annotation_1, annotation_2)) + '\n')

                
                
def parse_results(filename):
    """Accept a file containing the results of blasting query_org against
    target_org and return a dictionary of the form {protein in query_org: first
    blast hit in target_org}"""
    hits = {}
    query_def_pattern = re.compile(r"""<Iteration_query-def>
                                       gi\|(.*)\|ref.*
                                       </Iteration_query-def>""", re.X)
    hit_def_pattern = re.compile(r"<Hit_def>gi\|(.*?)\|ref.*?</Hit_def>")
    evalue_pattern =  re.compile(r"<Hsp_evalue>(.*?)</Hsp_evalue>")
    cutoff = 1e-10
    found_best_hit = False
    with open(filename) as f:
        for line in f:
            query_def_match = query_def_pattern.search(line)
            hit_def_match = hit_def_pattern.search(line)
            evalue_match = evalue_pattern.search(line)
            if query_def_match:
                query_name = query_def_match.group(1)
                found_best_hit = False
            elif hit_def_match and not found_best_hit: 
                hit_name = hit_def_match.group(1)
            elif evalue_match and not found_best_hit:
                evalue = float(evalue_match.group(1))
                if evalue < cutoff:
                    hits[query_name] = (hit_name, cutoff)
                    found_best_hit = True 
    return hits

def find_reciprocals(d1, d2):
    """Given two dictionaries d1 and d1 collecting the matches
    for org1 against org2 and vice versa, return a list of tuples
    [(prot1, prot2)] such that prot1:prot2 is in hits1 and prot2:prot1
    is in hits2"""
    hits1 = {k:v[0] for k, v in d1.iteritems()}
    hits2 = {k:v[0] for k, v in d2.iteritems()}
    reciprocals = [(h1, h2) for h1 in hits1 for h2 in hits2
                   if h2 in hits1[h1] and h1 in hits2[h2]]
    return reciprocals

def load_reciprocals(orgs):
    """Read the reciprocal blast results and return a dictionary of form
    {results_filename:[(geneA, geneB)]}"""
    def get(org1,org2):
        filename = "reciprocals_%s_%s.csv" % (org1,org2)
        lines = open(os.path.join("reciprocal_results", filename)).readlines()
        results = [re.search(r"(.*)\t(.*)\n", line).groups() for line in lines]
        return results

    all_results = {}
    for org1,org2 in choose2(orgs):
        #for filename in dir_contents:
        print org1,org2
        try:
            results = get(org1,org2)
        except IOError:
            reversed_results = get(org2,org1)
            results = map(reverse,reversed_results)
        filename = "reciprocals_%s_%s.csv" % (org1,org2)
        all_results[filename] = results
    return all_results

def get_genome(data_dir,org):
    genome_file_name = os.path.join(data_dir, org+'.gbk')
    return SeqIO.read(genome_file_name, 'genbank')

def make_locus_tag_dict(gbk_filename):
    print gbk_filename
    genome = SeqIO.read(gbk_filename, 'genbank')
    print "finding CDSs"
    CDSs = ([feature for feature in genome.features if feature.type == 'CDS'])
    print "findings gis"
    gis = [cds.qualifiers['db_xref'][0][3:] for cds in CDSs]
    print "finding locus tags"
    locus_tags = [cds.qualifiers['locus_tag'][0] for cds in CDSs]
    return {gi:locus_tag for (gi, locus_tag) in zip(gis, locus_tags)}        
    
def make_annotation_dicts(data_dir, orgs):
    """return a dictionary of form anno_dict[org][locus_tag] == annotation"""
    anno_dict = {}
    for org in orgs:
        anno_dict[org] = {}
        genome = get_genome(data_dir, org)
        CDSs = ([feature for feature in genome.features
                 if feature.type == 'CDS'])
        for cds in CDSs:
            description = (head(cds.qualifiers['product'])
                           if 'product' in cds.qualifiers
                           else 'None')
            anno_dict[org][head(cds.qualifiers['locus_tag'])] = description
    return anno_dict

def reverse_complement(seq):
    table = string.maketrans("ACGTacgt", 'TGCAtgca')
    return seq.translate(table)[::-1]

def download_genome(accession, data_dir):
    print "Downloading the genome %s from NCBI database..." % accession,
    file_path = os.path.join(data_dir, accession+'.gbk')
    net_handle = Entrez.efetch(db='nuccore', id=accession,
                               retmode='gbwithparts', rettype='text')
    out_handle = open(file_path, 'w')
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print "done."

ORGS = ['NC_000913', 'NC_002516']
BLAST_BIN_PATH = "/home/sefa/ncbi-blast-2.2.29+/bin/"
DB_PATH = 'blastdb'
DATA_PATH = 'data'
BLAST_RESULTS_PATH = 'blast_results'
RECIPROCAL_BLAST_RESULTS_PATH = 'reciprocal_results'
epsilon = 10**-10

def main():
    populate_dbs(ORGS, DATA_PATH, BLAST_BIN_PATH, DB_PATH)
    reciprocal_blasts(ORGS)
    collate_reciprocal_blasts(ORGS, DATA_PATH,BLAST_RESULTS_PATH,
                              RECIPROCAL_BLAST_RESULTS_PATH)
