#/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## imports
import os
import errno
import itertools

## functions
def check_input(infiles):
	for infi_tuple in infiles:
		# does current file exist?
		if os.path.isfile(infi_tuple[0]): pass
		else: return False
		# is current file correct format?
		if infi_tuple[0].endswith(infi_tuple[1]): pass
		else: return False
	return True

def mkdir_force(path):
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise

def grouped(iterable, n):
	"""
	s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
	"""
	out_list = itertools.izip(*[iter(iterable)] * n)
	return out_list

## classes
class ChromosomeSNPMap:
	def __init__(self):
		"""
		This object contains lists for every chromosome in our input samples.
		Each discrete SNP from our MAP file will be sorted into dict key/val.
		"""
		self.mapping = {}
		## Create an entry for each chromosome
		for i in range(0,23):
			self.mapping['chr{0}'.format(i)]=[]
		self.mapping['chrX'] = []
		self.mapping['chrY'] = []

	def append(self, indv_snp):
		key = 'chr{0}'.format(indv_snp.get_chr())
		self.mapping[key].append(indv_snp)

class indvSNP:
	def __init__(self, chromosome, snp_name, col3, col4):
		"""
		Individual SNP entry from our MAP file
		:param chromosome: the chromosome
		:param snp_name: what the SNP is called
		:param col3: unknown variable but requested to be included
		:param col4: unknown variable but requested to be included
		"""
		self.chr = chromosome
		self.snp_name = snp_name
		self.col3 = col3
		self.col4 = col4

	def get_chr(self):
		return self.chr
	def get_snpname(self):
		return self.snp_name
	def get_col3(self):
		return self.col3
	def get_col4(self):
		return self.col4

class indvAllele:
	def __init__(self, snp_name, sample_id, allele1_fw, allele2_fw):
		"""
		Individual allele entry from our report file
		:param snp_name: the SNP present
		:param sample_id: the sample it was found in
		:param allele1_fw: the original allele value 
		:param allele2_fw: the mutated allele value
		"""
		self.snp_name = snp_name
		self.sample_id = sample_id
		self.allele1_fw = allele1_fw
		self.allele2_fw = allele2_fw

	def get_snpname(self):
		return self.snp_name
	def get_sampleid(self):
		return self.sample_id
	def get_allele1_fw(self):
		return self.allele1_fw
	def get_allele2_fw(self):
		return self.allele2_fw

class indvSample:
	def __init__(self, family_id, sample_id, mother, father, sex, phenotype, mutation_list):
		"""
		Individual sample from our PED file
		:param family_id: 
		:param sample_id: 
		:param mother: 
		:param father: 
		:param sex: 
		:param phenotype: 
		:param mutation_list: 
		"""
		self.family_id = family_id
		self.sample_id = sample_id
		self.mother = mother
		self.father = father
		self.sex = sex
		self.phenotype = phenotype
		self.mutation_list = mutation_list

		self.mapping = {}
		for i in range(0,23):
			self.mapping['chr{0}'.format(i)]=[]
		self.mapping['chrX'] = []
		self.mapping['chrY'] = []

	def append(self, chromosome, data):
		key = '{0}'.format(chromosome)
		self.mapping[key].append(data)

	def get_familyid(self):
		return self.family_id
	def get_sampleid(self):
		return self.sample_id
	def get_mother(self):
		return self.mother
	def get_father(self):
		return self.father
	def get_sex(self):
		return self.sex
	def get_phenotype(self):
		return self.phenotype
	def get_mutationlist(self):
		return self.mutation_list
