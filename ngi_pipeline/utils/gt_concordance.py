import os
import shutil
import re
import subprocess

from pyexcel_xlsx import get_data

from ngi_pipeline.utils.classes import with_ngi_config
from ngi_pipeline.utils.charon import CharonSession


def parse_xlsx_files(xl_files_path):
    xl_files = [os.path.join(xl_files_path, filename) for filename in os.listdir(xl_files_path) if '.xlsx' in filename]
    genotype_data = {}
    for xl_file in xl_files:
        data = get_data(xl_file)
        sheet_name = 'HaploView_ped_0'
        data = data.get(sheet_name)
        header = data[0]

        for row in data[1:]:
            if not re.match(r"^ID\d+-P\d+_\d+", row[1]):
                continue
            # if the same sample occurs in a different 'xlsx' files, will be overwriten
            sample_id = row[1].split('-')[-1]
            project_id = sample_id.split('_')[0]
            for rs_id in header[9:]:
                rs_index = header.index(rs_id)
                if project_id not in genotype_data:
                    genotype_data[project_id] = {}
                if sample_id not in genotype_data[project_id]:
                    genotype_data[project_id][sample_id] = {}
                reference, alternative = row[rs_index].split()
                genotype_data[project_id][sample_id][rs_id] = {'gt_allele1': reference, 'gt_allele2': alternative}

        # todo: add try-except
        new_path = os.path.join(xl_files_path, 'archived')
        if not os.path.exists(new_path):
            os.makedirs(new_path)

        filename = os.path.basename(xl_file)
        # shutil.move(xl_file, os.path.join(new_path, filename))

    return genotype_data

# todo: need a config

def run_gatk(sample_id, config):
    project = sample_id.split('_')[0]
    project_analysis_path = config.get('project_analysis_path')
    bamfile_path = config.get('bamfile_path')
    bamfile = os.path.join(project_analysis_path, project, bamfile_path, "{}.clean.dedup.bam".format(sample_id))
    output_dir = config.get('sample_analysis_output_path')

    output_file = os.path.join(project_analysis_path, output_dir, "{}.vcf".format(sample_id))

    gatk_ref_file=config.get('gatk_ref_file')
    if gatk_ref_file is None:
        logging.error('Config must contain gatk_ref_file')
        raise RuntimeError('Config must contain gatk_ref_file')
    gatk_var_file=config.get('gatk_var_file')
    if gatk_var_file is None:
        logging.error('Config must contain gatk_var_file')
        raise RuntimeError('Config must contain gatk_var_file')
    interval_file = config.get('gatk_var_file')
    if interval_file is None:
        logging.error('Config must contain interval_file')
        raise RuntimeError('Config must contain interval_file')
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    options = """-T UnifiedGenotyper  -I {bamfile} -R {gatk_ref_file} -o {output_file}  -D {gatk_var_file} -L {interval_file}
        -out_mode EMIT_ALL_SITES """.format(bamfile=bamfile,
                                            output_file=output_file,
                                            gatk_ref_file=gatk_ref_file,
                                            gatk_var_file=gatk_var_file,
                                            interval_file=interval_file
                                           )
    gatk_path = config.get('gatk_path')
    if gatk_path is None:
        logging.error('Config must contain gatk_path')
        raise RuntimeError('Config must contain gatk_path')

    full_command = 'java -Xmx6g -jar {} {}'.format(gatk_path, options)
    subprocess.call(full_command.split())
    # todo: if not error: return output file, else: don't know
    return output_file


def run_concordance(vcf_filename, rs_positions):
    with open(vcf_filename) as vcf_file:
        vcf_data = vcf.Reader(vcf_file)
        sample = vcf_data.samples[0]
        for record in vcf_data:
            reference = str(record.REF[0])
            alternative = str(record.ALT[0])
            chromosome = str(record.CHROM)
            position = str(record.POS)
            chrom_position = "{}-{}".format(chromosome, position)
            rs_position = rs_positioins[chrom_position]

            genotype = str(record.genotype(sample)['GT'])
            allele1, allele2 = genotype.split('/')
            allele1 = reference if allele1 == '0' else alternative
            allele2 = reference if allele2 == '0' else alternative

            gt_allele1 = genotype_data[sample][rs_position]['gt_allele1'] # b1
            gt_allele2 = genotype_data[sample][rs_position]['gt_allele2'] # b2

            concordance = set([allele1, allele2]) == set([gt_allele1, gt_allele2])
            if not concordance:
                print chromosome, position, gt_allele1, gt_allele2, allele1, allele2

@with_ngi_config
def genotype_stockholm(config=None, config_file_path=None):
    print('analyse_all')
    print config
    print config_file_path
    # sample_id = args.__dict__.get('sample_id') or args.__dict__.get('project_id')
    # xlsx_path = "/Users/kate/work/gt_concordance"
    # 1. parse xl files

    gt_config = config.get('gt_concordance', {})
    if not gt_config:
        logging.error('config must have gt_concordance parameters')
        raise RuntimeError('Config must have gt_concordance parameters')


    rs_positions = {}
    snps_filepath = gt_config.get('snps_file')

    with open(snps_filepath) as snps_file:
        for line in snps_file.readlines():
            chromosome, position, rs_position, reference, alternative = line.split()
            chromosome_position = "{}-{}".format(chromosome, position)
            rs_positions[chromosome_position] = rs_position
    snps_file.close()
    print rs_positions


    xlsx_path = gt_config.get('xlsx_path')
    print xlsx_path
    genotype_data = parse_xlsx_files(xlsx_path)
    charon_session = CharonSession()
    for project_id in genotype_data.keys():
        print project_id

        samples = charon_session.project_get_samples(project_id)
        samples = samples.get('samples')
        for sample in samples:
            sample_id = sample.get('sampleid')
            print '\nsample: ', sample_id
            analysis_status = sample.get('analysis_status')
            genotype_status = sample.get('genotype_status')
            print 'analysis status: ', analysis_status
            print 'genotype status: ', genotype_status
            # if analysis_status == 'ANALYSED': # should be like this
            if analysis_status != 'ANALYSED':
                # if Genotype Status != 'AVAILABLE':
                if genotype_status != 'AVAILABLE':
                    print('run gatk')
                    # run gatk genotyping & change status to available
                    # save status in charon
                    vcf_file = run_gatk(sample_id, gt_config)
                    sample['genotype_status'] = 'AVAILABLE'

                else: # if genotype_status == 'AVAILABLE':
                    vcf_file = os.path.join(gt_config.get('project_analysis_path'), project_id, gt_config.get('sample_analysis_output_path'), '{}.vcf'.format(sample_id))

                # 4. run the concordance script
                run_concordance(vcf_filename=vcf_file, rs_positions=rs_positions)
            else:
                # dont move xlsx file to nosync,
                pass

    # todo: 2. check sample in charon: if analysis status != ANALISED: skip
    # todo: 3. if Genotype Status != 'AVAILABLE': run gatk genotyping
    # todo: 4. run the concordance script

# genotype_stockholm()