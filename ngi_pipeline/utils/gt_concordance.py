import os
import shutil
import re

from pyexcel_xlsx import get_data


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
    config = config.get('gt_concordance', {})
    if config:
        project = sample_id.split('_')[0]
        project_analysis_path = config.get('project_analysis_path')
        bamfile_path = config.get('bamfile_path')
        bamfile = os.path.join(project_analysis_path, project, bamfile_path, sample_id, bamfile_ext)


        output_file = os.path.join(project_path, output_dir, "{}.vcf".format(sample_id))

        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))

        options = """-T UnifiedGenotyper  -I {bamfile} -R {gatk_ref_file} -o {output_file}  -D {gatk_var_file} -L {interval_file}
            -out_mode EMIT_ALL_SITES """.format(bamfile=bamfile,
                                                output_file=output_file,
                                                gatk_ref_file=gatk_ref_file,
                                                gatk_var_file=gatk_var_file,
                                                interval_file=os.path.join(static_files_path, interval_filename)
                                              )

        full_command = 'java -Xmx6g -jar {} {}'.format(gatk_path, options)
        subprocess.call(full_command.split())
        return output_file


def genotype_stockholm():
    print('analyse_all')
    sample_id = args.__dict__.get('sample_id') or args.__dict__.get('project_id')
    xlsx_path = "/Users/kate/work/gt_concordance"
    # 1. parse xl files
    genotype_data = parse_xlsx_files(xlsx_path)
    for project_id in genotype_data.keys():
        print project_id

    # todo: 2. check sample in charon: if analysis status != ANALISED: skip
    # todo: 3. if Genotype Status != 'AVAILABLE': run gatk genotyping
    # todo: 4. run the concordance script