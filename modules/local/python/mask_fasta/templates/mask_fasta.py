#!/usr/bin/env python

import sys
import argparse
import platform

from seq_sim.io.fasta import Fasta

def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")

def read_fasta(path):
    fasta = Fasta()
    fasta = fasta.read_file(path)
    for id, seq in fasta.items():
        fasta[id] = ''.join([byte.decode('utf-8') for byte in seq])
    return fasta

def write_fasta(path, fasta, column_width=70):
    with open(path, "w") as out_f:
        for id, seq in fasta.items():
            out_f.write(">" + id + "\n")
            for i in range(0, len(seq), column_width):
                out_f.write(seq[i:i+column_width] + "\n")

def mask_fasta(fasta, ref, output):
    # parse the fasta files
    ref_fasta = read_fasta(ref)
    source_fasta = read_fasta(fasta)

    # Replace 'N' in the source sequence with the reference sequence's base
    output_fasta = {}
    for id, seq in source_fasta.items():
        if id in ref_fasta:
            ref_seq = ref_fasta[id]
            assert len(seq) == len(ref_seq), f"Reference Sequences must be of the same length: {id}"
            new_seq = ''.join(ref if con == 'N' else con for con, ref in zip(seq, ref_seq))
            output_fasta[id] = new_seq
        else:
            sys.exit("Error: " + id + " not found in the reference fasta file.")

    # write the masked fasta file
    write_fasta(output, output_fasta)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--fasta", default="!{fasta}")
    parser.add_argument("--ref", default="!{ref}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    mask_fasta(args.fasta, args.ref, args.output)
    dump_versions(args.process_name)
