#!/usr/bin/env python

import sys
import datatable as dt
import wget
from metapub import PubMedFetcher
import os
import logging
import csv

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

## Should be taken as input arguments
doi = "6811941"
version = 47
prev_version =  version - 1

url_collated_tsv = "https://zenodo.org/record/" + str(doi) + "/files/cancermine_collated.tsv?download=1"
url_sentences_tsv = "https://zenodo.org/record/" + str(doi) + "/files/cancermine_sentences.tsv?download=1"

local_dir = basedir

local_collated_tsv = "data-raw/cancermine_collated.v" + str(version) + ".tsv"
local_sentences_tsv = "data-raw/cancermine_sentences.v" + str(version) + ".tsv"
local_citations_tsv = "output/cancermine_citations.v" + str(version) + ".tsv"
local_citations_prev_tsv = "output/cancermine_citations.v" + str(prev_version) + ".tsv"


if not os.path.exists(local_collated_tsv):
    wget.download(url_collated_tsv, local_collated_tsv)
if not os.path.exists(local_sentences_tsv):
    wget.download(url_sentences_tsv, local_sentences_tsv)

os.system('gzip -f')

fetch = PubMedFetcher()

## Load dictionary with already retrieved PMIDs
retrieved_pmids = {}
if os.path.exists(local_citations_prev_tsv):
    f = open(local_citations_prev_tsv, 'r')
    reader = csv.DictReader(f, fieldnames = ['pmid','citation','html_link'], delimiter='\t')
    for row in reader:
        retrieved_pmids[str(row['pmid'])] = str(row['pmid']) + '\t' + str(row['citation']) + '\t' + str(row['html_link'])
    f.close()

## read sentences, convert to list
list_cm = dt.fread(local_sentences_tsv, sep = "\t", fill = True).to_list()
## get unique PMIDs
pmids = sorted(set(list_cm[1]))

f = open(local_citations_tsv, "a")
i = 0
total_pmids = len(pmids)
for pmid in pmids:
    i = i + 1

    if str(pmid) in retrieved_pmids.keys():
        f.write(retrieved_pmids[str(pmid)] + '\n')
        continue
    try:
        article = fetch.article_by_pmid(str(pmid))
    except Exception as e:
        print("Oops!", e.__class__, "occurred. Could not find " + str(pmid))
        continue

    first_author_surname = "Unknown"
    if len(article.authors) > 0:
        first_author_surname = article.authors[0].split(' ')[0]
    citation = first_author_surname + " et al., " + str(article.journal) + ", " + str(article.year)
    html_link = "<a href='https://www.ncbi.nlm.nih.gov/pubmed/" + str(pmid) + "' target='_blank'>" + citation + "</a>"
    if i % 100 == 0 and i > 0:
        pct_completed = '{:3.1f}'.format((i / total_pmids) * 100)
        logging.info("Retrieved citation data for n = " + str(i) + " PMIDs ( Completed " + str(pct_completed) + "% of all PMIDs .. )")
    f.write(str(pmid) + '\t' + str(citation) + '\t' + str(html_link) + '\n')
    #i = i + 1

f.close()
