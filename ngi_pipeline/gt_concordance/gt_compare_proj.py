import sys
import glob
import os

from config import BAM_FILE, MAIN_CONC_FILE, OUTPUT_DIR
from gt_compare import run_concordance_check


def run_project_concordance(maf_gt_filename, project):
    # open maf file
    maf_gt_file = open(maf_gt_filename, 'r')

    # get a list of lines
    lines = maf_gt_file.readlines()
    header = lines[0].split(',')

    # get genotyped samples:
    # check if header of .csv file contains substrings with project id
    # e.g.: ID10-P1234_101, where P1234 is project id
    samples = []
    for sample in header:
        if project in sample:
	    # sample has format ID9-P1234_101
	    # need to get rid of ID9-
            samples.append(sample.split('-')[-1])


    output = []
    if samples:
        print "{number_of_samples} samples MAF-typed".format(number_of_samples=len(samples))
        seq_samples = []
        for sample in samples:
            # Check if sample has been sequenced
	    print sample
            bamfile = BAM_FILE.format(project=project, sample=sample)
	    print bamfile
            if os.path.isfile(bamfile):
                seq_samples.append(sample)
                concordance = run_concordance_check(sample, maf_gt_filename)
                output.append(concordance)
            else:
                print '{sample} not sequenced'.format(sample=sample)

	bamfiles = glob.glob(BAM_FILE.format(project=project, sample='P*'))

        print "---------------------------------"
        print "Number of Sequenced Samples %s" % len(bamfiles)
        print "Number of Samples MAF-genotyped: %s" % len(samples)
        print "Number of Samples Genotyped and Sequenced %s" % len(seq_samples)
        print "---------------------------------"
        print "Sample %Concordance"
        conc_file = open(MAIN_CONC_FILE.format(project=project), 'w+')
        for item in output:
            # take the first line and the last part of the last line
            line = item[0] + item[2].split()[-1]
            conc_file.write(line)
            print line
        print "Summarized results can be found in: %s" % MAIN_CONC_FILE.format(project=project)
        print "Results for each sample can be found in %s" % OUTPUT_DIR.format(sample="<sample>",project=project)


    else:
        print "No samples are genotyped"


if __name__ == '__main__':
    try:
        project = sys.argv[1]
        maf_gt_filename = sys.argv[2]
        run_project_concordance(maf_gt_filename, project)
    except IndexError:
        print """Usage:
        gt_compare.py <project> <maf-gt-filename>

        Arguments:
        <project>
                - eg: P1234
        <maf-gt-filename>
                - eg: ID10_FM3_150601.csv"""
        exit()
