#/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## imports
import os
import sys
import argparse
import logging as log
from .__backend import check_input
from .__backend import mkdir_force
from .__backend import grouped

## objects
from .__backend import indvSNP
from .__backend import indvAllele
from .__backend import indvSample
from .__backend import ChromosomeSNPMap

## terminal colouring
class clr:
	def __init__(self):
		pass

	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

## actual script
class SNPMatch:
	def __init__(self):
		"""
		SNPMatch docstring goes here lmao
		"""

		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='snpmatch', description='SNPMatch: Split *.PED data by chromosome while maintaining *.MAP SNP order.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		self.parser.add_argument('-p', '--ped', help='PED file. Contains your samples and their respective SNP data.', nargs=1, required=True)
		self.parser.add_argument('-r', '--report', help='FinalReport file from your GIGAMUGA run.', nargs=1, required=True)
		self.parser.add_argument('-m', '--map', help='MAP file containing SNPs and positions.', nargs=1, required=True)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()

		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'snpm__ ', clr.end, 'SNPMatch: Split PED/MAP by chromosome.'))
			log.info('{}{}{}{}'.format(clr.bold, 'snpm__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		## Input files, check exist/format
		self.infiles = [(self.args.ped[0], "PED"), (self.args.report[0], "TXT"), (self.args.map[0], "MAP")]
		check_input(self.infiles)

		## Check output path, make if doesn't exist
		self.output_target = self.args.output[0]
		if not os.path.exists(self.output_target):
			mkdir_force(self.output_target)

		##
		## Begin processing
		## First stage: form SNP order objects (per chromosome in MAP file)
		self.ordered_snpmap = self.snp_map_order()
		## Second stage: split map file into chromosomes we just created
		self.split_orderedmap()
		## Third stage: scrape report data so we can join everything together
		self.processed_alleles = self.scrape_report()
		## Fourth stage: identify all the SNP data from our PED file
		self.processed_samples = self.scrape_mutations()
		## Fifth stage: combine information so we can split our PED file into chromosome
		self.match_chromosome_snp()
		## Sixth stage: convert stored information into our split PED files
		self.split_mutation_data()


	def snp_map_order(self):

		## temporary storage before splitting into respective vectors
		preprocess_storage = []

		## create a hierarchy object for all SNPs to be stored in
		postprocess_storage = ChromosomeSNPMap()

		## open our map file, read each line and split
		## once split, make indvSNP object for that SNP
		## append to temp storage
		map_file = self.infiles[2][0]
		with open(map_file, 'r') as mapfi:
			for map_entry in mapfi.read().splitlines():
				dat_split = map_entry.split()
				current_snp = indvSNP(chromosome=dat_split[0], snp_name=dat_split[1],
									  col3=dat_split[2], col4=dat_split[3])
				preprocess_storage.append(current_snp)

		## all SNP processed, split into respective chromosome vector
		## within our vector-wide map object
		for indv_snp in preprocess_storage:
			postprocess_storage.append(indv_snp)

		return postprocess_storage

	def split_orderedmap(self):

		## Iterate over our dictionary of mappings
		processed_snpcount = 0
		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Splitting *.MAP into individual chromosomes...'))
		for chromosome, snp_list in self.ordered_snpmap.mapping.iteritems():

			## Create MAP file for this chromosome
			target_name = 'mapping_{}.map'.format(chromosome)
			target_output = os.path.join(self.output_target, target_name)

			## Populate this chromosome's MAP file with all SNPs found from input
			with open(target_output, 'w') as chrmapfi:
				for snp in snp_list:
					chrmapfi.write('{}\t{}\t{}\t{}\n'.format(snp.get_chr(), snp.get_snpname(),
															 snp.get_col3(), snp.get_col4()))
			processed_snpcount += len(snp_list)
		log.info('{}{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Split *.MAP SNP total: ', processed_snpcount))

	def scrape_report(self):

		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Combining sample SNP/allele information...'))
		## I/O data for reportfile
		report_file = self.infiles[1][0]

		## List of objects for each SNP/Allele entry in report file
		processed_alleles = []

		## Get raw input from report, only keep 'data' section
		with open(report_file, 'r') as repfi:
			untrimmed_input = repfi.read().splitlines()

		## Get to the 'data' section, only process from index of 'data'+2
		## index+2 == skip header lines in ReportFile.txt
		## only keep first four columns of trimmed input
		target_element = '[Data]'
		target_index = untrimmed_input.index(target_element)
		trimmed_input = untrimmed_input[target_index+2:]
		split_input = [x.split('\t')[0:4] for x in trimmed_input]

		## for each data entry, assign it into an allele object for easy retreival
		## add allele to processed_alleles for this instance
		for entry in split_input:
			allele_object = indvAllele(snp_name=entry[0], sample_id=entry[1],
									   allele1_fw=entry[2], allele2_fw=entry[3])
			processed_alleles.append(allele_object)

		return processed_alleles

	def scrape_mutations(self):

		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Gathering sample mutation information...'))
		## I/O for pedfile
		ped_file = self.infiles[0][0]

		## List of tuples: (sample ID/background information, mutation list)
		preprocessed_samples = []

		## List of sample objects for each entry in the PED file
		processed_samples = []

		## Get information
		with open(ped_file, 'r') as pedfi:
			data_samples = pedfi.read().splitlines()

		## Split current sample information into 'data' list and 'mutation' list
		## Data list == family, sample, mum, dad, sex, phenotype
		## Mutation list == SNP mutations present in sample
		mutation_lists = [x.split('\t')[6:] for x in data_samples]
		data_samples = [x.split('\t')[:6] for x in data_samples]

		## Fix mutation list so it is a list of tuples, rather than a single list of strings
		## as information is (allele1_fw, allele2_fw) in input... maintain structure
		grouped_mutations = []
		for sample in mutation_lists:
			if not len(sample)%2 == 0:
				log.error('{}{}{}{}'.format(clr.red,'snpm__ ',clr.end, ' The length of a specified sample SNP list in your *.PED file is not divisible by 2.'))
				sys.exit(2)
			else:
				grouped_mutations.append(list(grouped(sample, 2)))

		## Combine this data into a tuple for each sample
		for x, y in zip(data_samples, grouped_mutations):
			preprocessed_samples.append((x, y))

		## Turn tuple into an object for the sample, append to list of all samples in instance
		for current_individual in preprocessed_samples:
			sample_object = indvSample(family_id=current_individual[0][0], sample_id=current_individual[0][1],
									   mother=current_individual[0][2], father=current_individual[0][3],
									   sex=current_individual[0][4], phenotype=current_individual[0][5],
									   mutation_list=current_individual[1])
			processed_samples.append(sample_object)

		return processed_samples

	def match_chromosome_snp(self):

		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Gathering sample mutation information...'))

		## Iterate over all chromosomes we have data for...
		# for chromosome, snp_list in self.ordered_snpmap.mapping.iteritems():
		# 	print '\n\n\n\n\n\n\n\n\nWorking on Chromosome: {}'.format(chromosome)
		# 	print 'number of SNPs in this chr: {}'.format(len(snp_list))
		#
		# 	## Loop over all alleles (with associated sample_ID and mutation values)
		# 	## Loop over all snps present in the current chromosome
		# 	## If the current SNPs match, this SNP is present in this chromosome
		# 	## Identify all samples where this is the case, and append information to mapping
		# 	## i.e. create a list of SNP on <curr_chromosome> present in <sample_id>, and their value
		# 	for mutation in self.processed_alleles:
		# 		for chr_snp in snp_list:
		# 			if mutation.get_snpname() in chr_snp.get_snpname():
		# 				for individual in self.processed_samples:
		# 					if individual.get_sampleid() == mutation.get_sampleid():
		# 						target_info = (mutation.get_snpname(),
		# 									   [mutation.get_allele1_fw(), mutation.get_allele2_fw()])
		# 						individual.append(chromosome, target_info)

		## testing only on chr Y
		chromosome = 'chrY'
		snp_list = self.ordered_snpmap.mapping[chromosome]

		for mutation in self.processed_alleles:
			for chr_snp in snp_list:
				if mutation.get_snpname() in chr_snp.get_snpname():
					for individual in self.processed_samples:
						if individual.get_sampleid() == mutation.get_sampleid():
							target_info = (mutation.get_snpname(),
										   [mutation.get_allele1_fw(), mutation.get_allele2_fw()])
							individual.append(chromosome, target_info)

	def split_mutation_data(self):

		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Outputting chromosome split PED files...'))

		## Loop over all samples we have, output to the respective chromosome PED file...
		for sample in self.processed_samples:
			print 'hi'






















def main():
	try:
		SNPMatch()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'snpm__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)