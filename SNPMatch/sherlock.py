#/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## imports
import os
import sys
import argparse
import logging as log
import multiprocessing
from .__backend import clr
from .__backend import grouped
from .__backend import check_input
from .__backend import mkdir_force

## objects
from .__backend import indvSNP
from .__backend import indvAllele
from .__backend import indvSample
from .__backend import indvChromosome
from .__backend import ChromosomeSNPMap

## global
THREADS = multiprocessing.cpu_count()
PROC_ALLELE = []
PROC_SAMPLE = []
MAPPING_OUTPUTS = {}
MATCHED_OBJECTS = []

# multiprocessing worker function
def worker(input_iterable):

	"""
	Function to iterate over information within chromosome entry in our SNP mapping
	Must be defined at top-level so it can be pickled/serialised for multiprocessing
	:param input_iterable: current iteration of snpmap mapping
	:return: None
	"""

	## Get input from current iteration of snpmap
	## this instance of worker() on a discrete processor will run independently
	chrom, snp_list = input_iterable
	log.info('{}{}{}Worker (ID: {}) processing {}...'.format(clr.yellow, 'snpm__ ', clr.end, os.getpid(), chrom)),

	## current chromosome results dictionary
	current_workload = {'Chromosome': chrom}

	## Populate this chromosome's dictionary with
	## key == current sample
	## value == list, which will have all SNPs appended to it
	for individual in PROC_SAMPLE:
		current_workload[individual.get_sampleid()] = []

	# Loop over all alleles (with associated sample_ID and mutation values)
	# Loop over all snps present in the current chromosome
	# If the current SNPs match, this SNP is present in this chromosome
	# Identify all samples where this is the case, and append information to mapping
	# i.e. create a list of SNP on <curr_chromosome> present in <sample_id>, and their value
	for mutation in PROC_ALLELE:
		for chr_snp in snp_list:
			if mutation.get_snpname() in chr_snp.get_snpname():
				for individual in PROC_SAMPLE:
					if individual.get_sampleid() == mutation.get_sampleid():
						target_info = (mutation.get_snpname(),
							[mutation.get_allele1_fw(), mutation.get_allele2_fw()])
						current_workload[individual.get_sampleid()].append(target_info)

	## Append results object for this chromosome
	result_object = indvChromosome(current_workload['Chromosome'])
	result_object.set_results(current_workload)
	MATCHED_OBJECTS.append(result_object)

	## Output results to individual chromosome PED file
	chromosome_sample_strings = []

	## Loop over every sample in our results dictionary:
	for sample_key, snp_list in result_object.get_results().iteritems():
		sample_header = ''; sample_mutation = ''

		## build intro details for this sample
		for sample in PROC_SAMPLE:
			if sample.get_sampleid() == sample_key:
				sample_header = '\n{}\t{}\t{}\t{}\t{}\t{}\t'.format(
					sample.get_familyid(),
					sample.get_sampleid(),
					sample.get_mother(),
					sample.get_father(),
					sample.get_sex(),
					sample.get_phenotype()
				)

		## append sample string with mutation string
		for snp_tuple in snp_list:
			try:
				sample_mutation += '{}\t{}\t'.format(snp_tuple[1][0], snp_tuple[1][1])
			except IndexError:
				pass

		sample_string = sample_header + sample_mutation
		chromosome_sample_strings.append(sample_string)

	## Target the correct directories
	chrom = result_object.get_chromosome()
	desired_key = '{}_dir'.format(chrom)
	target_output = MAPPING_OUTPUTS[desired_key]

	## output the chromosome split ped file
	with open(target_output, 'w') as chr_outfi:
		for sample in chromosome_sample_strings:
			chr_outfi.write(sample)

	log.info('{}{}{}Worker (ID: {}) finished {}!'.format(clr.green, 'snpm__ ', clr.end, os.getpid(), chrom))

## actual program logic
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
		self.parser.add_argument('-t', '--threads', help='Number of CPU threads to use when mapping chromosome data. Default is system max.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()

		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'snpm__ ', clr.end, 'SNPMatch: Split PED/MAP by chromosome.'))
			log.info('{}{}{}{}'.format(clr.bold, 'snpm__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		## CPU threads
		self.threads = self.args.threads

		## Input files, check exist/format
		self.infiles = [(self.args.ped[0], "PED"), (self.args.report[0], "TXT"), (self.args.map[0], "MAP")]
		check_input(self.infiles)

		## Check output path, make if doesn't exist
		self.output_target = self.args.output[0]
		if not os.path.exists(self.output_target):
			mkdir_force(self.output_target)

		## Empty outputs
		self.mapping_outputs = {}
		## Create files for each chromosome
		for i in range(0,23):
			self.mapping_outputs['chr{0}_dir'.format(i)]=os.path.join(self.output_target,'chr{0}.ped'.format(i))
		self.mapping_outputs['chrX_dir'] = os.path.join(self.output_target,'chrX.ped')
		self.mapping_outputs['chrY_dir'] = os.path.join(self.output_target,'chrY.ped')
		global MAPPING_OUTPUTS; MAPPING_OUTPUTS = self.mapping_outputs

		##
		## Begin processing
		## First stage: form SNP order objects (per chromosome in MAP file)
		self.ordered_snpmap = self.snp_map_order()
		## Second stage: split map file into chromosomes we just created
		self.split_orderedmap()
		## Third stage: scrape report data so we can join everything together
		self.processed_alleles = self.scrape_alleles()
		## Fourth stage: identify all the SNP data from our PED file
		self.processed_samples = self.scrape_samples()
		## Fifth stage: combine information so we can split our PED file into chromosome
		## Each individual worker will output results to an individual chromosome PED file
		self.matched_resultobjects = []
		self.match_chromosome_snp()

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

	def scrape_alleles(self):

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

	def scrape_samples(self):

		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Gathering individual sample information...'))
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

		## inform how multi-threaded this instance will be
		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Gathering sample mutation information...'))
		log.info('{}{}{}{}'.format(clr.yellow, 'snpm__ ', clr.end, 'Launching {} process worker jobs...'.format(self.threads)))

		## set this instance's processed_alleles/samples to globals so the top-level worker function can access
		## cant pass as arguments due to requirement of passing iterator for chromosome(s)
		global PROC_ALLELE; global PROC_SAMPLE
		PROC_ALLELE = self.processed_alleles
		PROC_SAMPLE = self.processed_samples

		## wrap in KeyboardInterrupt exception
		## launch a pool of processes with the specified amount of thread (or default == max)
		## run the job through worker function, pass iterator of our chromosome/mapping data
		processor_pool = multiprocessing.Pool(self.threads)
		try:
			processor_pool.imap(worker, self.ordered_snpmap.mapping.iteritems())
		except KeyboardInterrupt:
			log.info('{}{}{}{}'.format(clr.red,'snpm__ ',clr.end,'Caught KeyboardInterrupt. Killing worker {}'.format(os.getpid())))
			processor_pool.terminate()
			processor_pool.join()
		else:
			processor_pool.close()
			processor_pool.join()

		## inform user we have closed the worker pool
		log.info('{}{}{}{}'.format(clr.green, 'snpm__ ', clr.end, 'Done, closing processor worker pool!'))

def main():
	try:
		SNPMatch()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'snpm__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)