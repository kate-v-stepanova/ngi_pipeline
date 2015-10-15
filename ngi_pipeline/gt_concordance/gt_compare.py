import os
import re
import sys
import subprocess
import vcf
import os

from config import BAM_FILE, INTERVAL_FILE, GATK_REF_FILE, GATK_VAR_FILE, GATK_HOME, SNPS_FILE, OUTPUT_DIR, CONC_FILE



def genotype_sample(sample):
	project = sample.split('_')[0]
	gatk_home = os.environ.get('GATK_HOME', GATK_HOME)
	command = "{gatk_home}/GenomeAnalysisTK.jar".format(gatk_home=gatk_home)
	bamfile = BAM_FILE.format(sample=sample, project=project)
   
	output_file = os.path.join(OUTPUT_DIR.format(sample=sample, project=project), "{sample}.vcf".format(sample=sample))
 
	if not os.path.exists(os.path.dirname(output_file)):
		os.makedirs(os.path.dirname(output_file))
	options = """-T UnifiedGenotyper  -I {bamfile} -R {gatk_ref_file} -o {sample}  -D {gatk_var_file} -L {interval_file} -out_mode EMIT_ALL_SITES """.format(bamfile=bamfile,
											sample=output_file,
						interval_file=INTERVAL_FILE,
											gatk_ref_file=GATK_REF_FILE,
											gatk_var_file=GATK_VAR_FILE)
	full_command = 'java -Xmx6g -jar ' + command + ' ' + options
	subprocess.call(full_command.split())


def generate_gt_file(sample, maf_filename):
	maf_file = open(maf_filename, 'r')
	lines = maf_file.readlines()
	header = lines[0].split(',')
	# find a substring that contains sample_id, e.g: ID0-P1234_123, where P1234_123 is sample
	sample_id = re.search(',.*?{sample}'.format(sample=sample), lines[0]).group(0).split(',')[-1]
	index = header.index(sample_id)
	maf_data = {}
	for line in lines:
		data = line.split(',')
		# rs10144418 - position number
		rs_id = data[0]
		maf_data[rs_id] = data[index]

	print maf_data

	result = []
	snps_file = open(SNPS_FILE, 'r')
	lines = snps_file.readlines()
	for line in lines:
		genotype = {}
		data = line.split()
		chromosome = data[0]
		position = data[1]
		rs_position = data[2]
		a1 = data[3]
		a2 = data[4]

		reference, alternative = maf_data[rs_position].split() # from csv file

		# convert to int - for sorting
		genotype['chromosome'] = chromosome if chromosome == "X" else int(chromosome)
		genotype['position'] = position
		genotype['rs_position'] = rs_position
		genotype['ref_sample'] = reference # from csv file
		genotype['alt_sample'] = alternative
		genotype['ref_genotype'] = a1 # from snps file
		genotype['alt_genotype'] = a2


		result.append(genotype)

	result = sorted(result, key=lambda x:x['chromosome'])

	project = sample.split('_')[0]
	# open or create file
	filename = os.path.join(OUTPUT_DIR.format(project=project, sample=sample), sample + '.gt')
	if not os.path.exists(os.path.dirname(filename)):
		os.makedirs(os.path.dirname(filename))


	output_file = open(filename, 'w+')
	output_lines = []
	for genotype in result:
		line = "{chromosome} {position} {rs_position} {ref_genotype} {alt_genotype} {ref_sample} {alt_sample}\n".format(
			chromosome=genotype['chromosome'],
			position=genotype['position'],
			rs_position=genotype['rs_position'],
			ref_sample=genotype['ref_sample'], # from csv file
			alt_sample=genotype['alt_sample'],
			ref_genotype=genotype['ref_genotype'], # from snps file
			alt_genotype=genotype['alt_genotype']
		)
		output_lines.append(line)
	output_file.writelines(output_lines)
	output_file.close()

	header_file = open(filename + '.header', 'w+')
	header_file.write('chromosome position rs_position reference alternative a1_gt a2_gt\n')
	header_file.close()


def check_concordance(sample):
	project = sample.split('_')[0]
	gt_file = os.path.join(OUTPUT_DIR.format(project=project, sample=sample), sample + '.gt')
	sample_file = os.path.join(OUTPUT_DIR.format(project=project, sample=sample), sample + '.vcf')
	vcf_reader = vcf.Reader(open(sample_file, 'r'))

	gt_file = open(gt_file, 'r')
	lines = gt_file.readlines() # rename to gt_lines

	hits = 0
	snps_number = 0
	lost = 0

	for record in vcf_reader:
		reference = str(record.REF[0])
		alternative = str(record.ALT[0])
		chromosome = str(record.CHROM)
		position = str(record.POS)

		genotype = str(record.genotype(sample)['GT'])
		a1, a2 = genotype.split('/')
		a1 = reference if a1 == '0' else alternative
		a2 = reference if a2 == '0' else alternative

		for line in lines:
			split_line = line.strip().split() # remove \n and split by ' '
			b1 = split_line[5]
			b2 = split_line[6]

			if split_line[0] == chromosome and split_line[1] == position:
				concordance = set([a1, a2]) == set([b1, b2])
				if concordance:
					hits += 1
					snps_number += 1
				else:
					if b1 != '0' and b2 != '0':
						snps_number += 1
					else:
						lost += 1
				gt_file.close()
				break


	percent_hits=(float(hits)/float(snps_number))*100
	output = [
		'Lost snps: %d\n' % lost,
		'Total number of hits: %d / %d\n' % (hits, snps_number),
		'Percent hits %d%%\n' % percent_hits
	]

	for line in output:
		print line

	# todo: do we really need these files?
	conc_filename = CONC_FILE.format(sample=sample, project=project)
	if not os.path.exists(os.path.dirname(conc_filename)):
		os.makedirs(os.path.dirname(conc_filename))
	conc_file = open(conc_filename, 'w+')
	conc_file.writelines(output)

	return  output


def run_concordance_check(sample, maf_gt_file):
	# sample='P1858_165': split -> project_id=P1858
	project = sample.split('_')[0]
	genotype_sample(sample)
	generate_gt_file(sample, maf_gt_file)
	output = check_concordance(sample)
	return output


if __name__ == '__main__':
	try:
		sample = sys.argv[1]
		maf_gt_file = sys.argv[2]
	except IndexError:
		print """Usage:
		gt_compare.py <sample> <maf-gt-filename>

		Arguments:
		<sample>
				- eg: P1234_101
		<maf-gt-filename>
				- eg: ID10_FM3_150601.csv"""
		exit()
	run_concordance_check(sample, maf_gt_file)

