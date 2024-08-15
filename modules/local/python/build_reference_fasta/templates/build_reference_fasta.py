#!/usr/bin/env python

import os
import sys
import errno
import argparse
import platform

from seq_sim.io.fasta import Fasta

def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")

def build_fasta_from_top_blast_hits(refs, blast, output):
    # parse the blast file
    with open(blast, "r") as in_f:
        blast_lines = in_f.readlines()

    # if the blast file is empty, error out
    if len(blast_lines) == 0:
        sys.exit("Error: No top hits detected in the blast file.")

    # build list of the second column of the blast file
    top_blast_hits = []
    for line in blast_lines:
        top_blast_hits.append(line.split("\t")[1])

    # if there are duplicates, error out
    if len(top_blast_hits) != len(set(top_blast_hits)):
        sys.exit("Error: Duplicate hits in the top blast hits.")

    # parse the reference fasta file
    fasta = Fasta()
    fasta = fasta.read_file(refs)
    for id, seq in fasta.items():
        fasta[id] = ''.join([byte.decode('utf-8') for byte in seq])

    # write the top blast hits to the output fasta file
    with open(output, "w") as out_f:
        for hit in top_blast_hits:
            if hit in fasta:
                out_f.write(">" + hit + "\n")
                out_f.write(fasta[hit] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--refs", default="!{refs}")
    parser.add_argument("--blast", default="!{blast}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    build_fasta_from_top_blast_hits(args.refs, args.blast, args.output)
    dump_versions(args.process_name)