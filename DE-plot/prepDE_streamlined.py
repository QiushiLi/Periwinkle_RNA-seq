
#!/usr/bin/env python3
import re, csv, sys, os, glob
import argparse

# Setup argparse for command line options
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', default='.')
parser.add_argument('-g', '--gene_output', default='gene_count_matrix.csv')
parser.add_argument('-t', '--transcript_output', default='transcript_count_matrix.csv')
parser.add_argument('-l', '--length', type=int, default=75)
parser.add_argument('-p', '--pattern', default=".")
parser.add_argument('-c', '--cluster', action="store_true")
parser.add_argument('-s', '--string', default="MSTRG")
parser.add_argument('-v', '--verbose', action="store_true")
parser.add_argument('--min_reads', type=int, default=10)

args = parser.parse_args()

def load_samples(input_path):
    samples = []
    if os.path.isfile(input_path):
        with open(input_path, 'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    continue
                sample_id, gtf_path = line.strip().split()
                if not os.path.exists(gtf_path):
                    samples.append((sample_id, gtf_path))
    else:
        for dirname in glob.glob(f"{input_path}/*"):
            if os.path.isdir(dirname) and re.search(args.pattern, dirname):
                gtf_files = glob.glob(f"{dirname}/*.gtf")
                samples.extend([(os.path.basename(dirname), gtf) for gtf in gtf_files])
    return samples

def process_samples(samples):
    pass

if __name__ == "__main__":
    samples = load_samples(args.input)
    process_samples(samples)
