import datetime
import time
import csv 
import argparse


def main(manifest, results_dir):
	
	sample_ids = []
	with open(manifest) as csvDataFile:
		cr = csv.DictReader(csvDataFile, delimiter='\t')
		for row in cr:
			sample = row['sample']
			if sample not in sample_ids:
				sample_ids.append(sample)

	ts = time.time()
	ts = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d_%H:%M:%S')
	with open(f"{results_dir}/{ts}.gvcf_entry_manifest.tsv", 'w+') as output_manifest:

		for sample in sample_ids:
			output_manifest.write(f"{sample}\t{results_dir}/{sample}.bam\t{results_dir}/{sample}.bam.bai\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--manifest')
	parser.add_argument('--results_dir')
	args = parser.parse_args()

	main(args.manifest, args.results_dir)
