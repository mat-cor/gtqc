import hail as hl
from hail import hadoop_open
import argparse
import sys
import subprocess
import pkg_resources

required = {'matplotlib'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing])

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

hl.init(default_reference='GRCh38')

def mahalanobis(x=None, data=None, cov=None):
    """Compute the Mahalanobis Distance between each row of x and the data
    x       : vector or matrix of data with, say, p columns.
    data    : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
    cov     : covariance matrix (p x p) of the distribution. If None, will be computed from data.
    return  : Mahalanobis distance in a numpy array
    """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
    inv_covmat = np.linalg.inv(cov)
    x_minus_mu_arr = x_minus_mu.to_numpy() # convert the x_minus_mu df to a numpy array
    x_minus_mu_bm = hl.linalg.BlockMatrix.from_numpy(x_minus_mu_arr)
    inv_covmat_bm = hl.linalg.BlockMatrix.from_numpy(inv_covmat)
    left_term = x_minus_mu_bm @ inv_covmat
    mahal = left_term @ x_minus_mu_bm.T
    return mahal.diagonal().to_numpy().ravel() # to_numpy returns ndarray and .ravel() converts the ndarray to array

def mahalanobis_distance(mt, dirname, cohortname, ancestry):
    """Compute the Mahalanobis Distance between each row of x and the data
    mt          : Hail MatrixTable
    dirname     : path to where the outputs will be saved
    cohortname  : the name of a cohort e.g. Sweden, ita etc.
    return      : filtered Hail MatrixTable
    """

    # PART 1
    print("using {} ancestry as the gnomad reference".format(ancestry))
    if ancestry == 'eur':
        AF_table_path = 'gs://dsge-covid19-data/gnomad/gnomad_nfe_AF_100p.ht'
    elif ancestry == "afr":
        AF_table_path = 'gs://dsge-covid19-data/gnomad/gnomad_afr_AF.ht'
    else:
        AF_table_path = 'gs://dsge-covid19-data/gnomad/gnomad_amr_AF.ht'

    gnomad_ht = hl.read_table(AF_table_path)
    print("The are {} SNPs in gnomAD".format(gnomad_ht.count()))

    mt_input = mt
    ht = mt_input.rows()
    print("The are {} SNPs in the data".format(ht.count()))

    # this will join the the two tables by locus, so there'll be SNPs with the same locus but different alleles (we remove them later on)
    table_joined = ht.key_by('locus').join(gnomad_ht.key_by('locus'))
    #print("The are {} SNPs common between the data and gnomAD".format(table_joined.count()))

    # remove the key so we can be able to select only a few columns
    # because both the gnomad and our data have the 'locus' and 'allele' cols, gnomad will be recoded as locus_1 and alleles_1
    table_result = table_joined.key_by()
    if ancestry == 'eur':
        table_result = table_result.select(table_result.locus, data_allele=table_result.alleles,
                                           gnomad_alleles=table_result.alleles_1,
                                           data_rsid=table_result.rsid, data_AF=table_result.variant_qc.AF[1],
                                           gnomad_AF=table_result.nfe_AF)
    elif ancestry == "afr":
        table_result = table_result.select(table_result.locus, data_allele=table_result.alleles,
                                           gnomad_alleles=table_result.alleles_1,
                                           data_rsid=table_result.rsid, data_AF=table_result.variant_qc.AF[1],
                                           gnomad_AF=table_result.afr_AF)
    else:
        table_result = table_result.select(table_result.locus, data_allele=table_result.alleles,
                                           gnomad_alleles=table_result.alleles_1,
                                           data_rsid=table_result.rsid, data_AF=table_result.variant_qc.AF[1],
                                           gnomad_AF=table_result.amr_AF)

    #flip = hl.case().when(table_result.data_allele[1] == table_result.gnomad_alleles[0], True).when(mt.allele == mt.alleles[1], False).or_missing()
    table_result = table_result.annotate(to_swap=hl.cond((table_result.gnomad_alleles[0]==table_result.data_allele[1]) & (table_result.gnomad_alleles[1]==table_result.data_allele[0]), True, False))
    #table_result = table_result.annotate(to_swap=hl.cond((str(table_result.gnomad_alleles[0])==str(table_result.data_allele[1]) and str(table_result.gnomad_alleles[1])==str(table_result.data_allele[0])), True, False))
    raw_ht = table_result.annotate(final_AF=hl.cond(table_result.to_swap == True, 1-table_result.gnomad_AF, table_result.gnomad_AF))

    # match keep only snps with either: (1) matching allele (our data and gnomad); or (2) 'swapped' alleles (to_swap is True)
    filtered_ht = raw_ht.annotate(keep=hl.cond((raw_ht.data_allele==raw_ht.gnomad_alleles) | (raw_ht.to_swap == True), True, False))
    #filtered_ht.show()
    filtered = filtered_ht.filter(filtered_ht.keep == True)
    filtered_not_common = filtered_ht.filter(filtered_ht.keep == True, keep = False)
    print("The are {} SNPs common between the data and gnomad".format(filtered.count()))
    print("The are {} SNPs found in the data but not in gnomad".format(filtered_not_common.count()))
    filtered.export(dirname + 'qc_step_4/{}/'.format(cohortname) + cohortname + '.CommonAlleles.tsv')
    filtered_not_common.export(dirname + 'qc_step_4/{}/'.format(cohortname) + cohortname + '.DataOnlyAlleles.tsv')

    # PART 2 (Mahalanobis Distance)
    file = pd.read_csv(dirname + 'qc_step_4/{}/'.format(cohortname) + cohortname + '.CommonAlleles.tsv', sep='\t')
    file = file.dropna() # remove any NAs in the df

    df_x = file[['data_AF', 'final_AF']] # select only the allele frequency columns
    df_x['mahala'] = mahalanobis(x=df_x, data=file[['data_AF', 'final_AF']]) # comput MD and create new col
    md30 = df_x[df_x['mahala'] > 30] # select variants with MD > 30

    # plot
    fig = plt.figure(figsize=(17, 11))
    ax = fig.add_subplot()
    ax.scatter(file['data_AF'], file['final_AF'], color = 'black', s = 0.3)
    ax.scatter(md30['data_AF'], md30['final_AF'], color = 'red', s = 0.3)
    ax.set_title('AF checks for {}'.format(cohortname), fontsize=20)
    ax.set_xlabel(xlabel='data_AF', fontsize=15)
    ax.set_ylabel(ylabel='gnomad_AF', fontsize=15)
    start_line = min(ax.get_xlim()[0], ax.get_ylim()[0])
    x_max = ax.get_xlim()[1]
    y_max = ax.get_ylim()[1]
    end_line = min(x_max, y_max)
    ax.plot([start_line, end_line], [start_line, end_line], ls="--")
    ax.set_xlim([-0.02, x_max])
    ax.set_ylim([-0.05, y_max])
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig('/tmp/{}.alleleFreqCheck.png'.format(cohortname), dpi=300, facecolor='w')
    hl.hadoop_copy('file:///tmp/{}.alleleFreqCheck.png'.format(cohortname), dirname + 'qc_step_4/{}/'.format(cohortname) + cohortname + '.alleleFreqCheck.png')
    plt.close()
    #plt.savefig(dirname + 'alleleFreqCheck.png', dpi=300)

    rsids = pd.merge(md30, file, how='inner', on=['data_AF','final_AF']) # merge the outliers with the original dataset wo we can get rsids for outliers (MD > 30)
    rsids_list = list(set(rsids['data_rsid'])) # convert to list and remove any duplicates using set()

    #mt = mt_input.filter_rows(hl.literal(rsids_list).contains(mt_input['rsid']), keep=False) # filter out the outliers
    #print("The are {} SNPs in the final data Aafter filtering allele checks".format(mt.count_rows()))

    #return mt
    return rsids_list


def hwe_allele_checks(mt, dirname, samplestokeep, cohorts_list, phenofile, ancestry, hwe_cas, hwe_con):
    """
        :param mt: # the full path to the mt used here is the one from the first QC step
        :param dirname: path to the bucket
        :param basename: file basename
        :param samplestokeep: list of samples from PCA to be used in HWE filtering
        :param phenoType: Case or Control
        :param cohortsFile: text file containing the name of cohorts you want to run HWE on (one cohort name each row)
        :return: writes out a hwe qced matrix table
        """
    # Annotate cohort information
    pheno = hl.import_table(phenofile, impute=True)
    pheno = pheno.key_by('ID')
    mt = mt.annotate_cols(pheno=pheno[mt.s])

    # mt = hl.read_matrix_table(mt)
    keep = hl.import_table(samplestokeep, no_header=True)
    # keep = hl.import_table(samplestokeep, key='01947TG2')
    keep = keep.key_by('f0')
    # mt = mt.filter_cols(hl.is_defined(keep.index(mt['s'])))
    mt = mt.semi_join_cols(keep)

    for i in cohorts_list:
        print("<============================== Allele frequency checks and HWE filtering for {} ==============================>".format(i))
        mtco = mt.filter_cols(mt.pheno.Cohort == i, keep=True)

        mtco = hl.variant_qc(mtco)
        # there are SNPs on chr3 around the 3:45848429 locus that we don't want to remove, so we get their rsids and remove them from
        # list of rsids already flagged to be removed from the mahalanobis_distance() function

        chr3_snps = hl.filter_intervals(mtco, [hl.parse_locus_interval('chr3:45348429-46348429')]).rsid.collect()
        allele_check_rsids = mahalanobis_distance(mtco, dirname, i, ancestry)
        rsids_to_remove = list(set(allele_check_rsids) - set(chr3_snps))
        # these are the SNPs in the chr3 locus that were selected during AF check, we flag them
        chr3_flagged_snps = list(set(chr3_snps).intersection(allele_check_rsids))

        mtco = mtco.filter_rows(hl.literal(rsids_to_remove).contains(mtco['rsid']), keep=False) # filter out the outliers
        # print("The are {} SNPs in the final data Aafter filtering allele checks".format(mt.count_rows()))

        initial_vars = mtco.count_rows()  # initial # of variants
        initial_samples = mtco.count_cols()  # initial # of samples
        print("The initial number of samples to be used in HWE filtering is {}".format(initial_samples))
        print("The initial number of variants to be used in HWE filtering is {}".format(initial_vars))

        mt_autosome = hl.variant_qc(mtco.filter_rows(mtco.locus.in_autosome_or_par()))
        n_aut = mt_autosome.count_rows()  # number of autosomal SNPs

        mt_X = mtco.filter_rows(mtco.locus.in_x_nonpar())
        mt_XF = hl.variant_qc(
            mt_X.filter_cols(mt_X.is_female == True))  # only choose females, males are haploid on chromX
        # if you do HWE considering males, the freq will be affected
        n_xf = mt_XF.count_rows()

        # SNPs_aut_HWE = mt_autosome.filter_rows(mt_autosome.variant_qc.p_value_hwe < hwe).rsid.collect()
        SNPs_aut_HWE_cases = hl.variant_qc(mt_autosome.filter_cols(mt_autosome.is_case == True))
        SNPs_aut_HWE_controls = hl.variant_qc(mt_autosome.filter_cols(mt_autosome.is_case == False))
        SNPs_aut_HWE_cas = SNPs_aut_HWE_cases.filter_rows(SNPs_aut_HWE_cases.variant_qc.p_value_hwe < hwe_cas).rsid.collect()
        SNPs_aut_HWE_con = SNPs_aut_HWE_controls.filter_rows(SNPs_aut_HWE_controls.variant_qc.p_value_hwe < hwe_con).rsid.collect()

        # SNPs_XF_hwe = mt_XF.filter_rows(mt_XF.variant_qc.p_value_hwe < hwe).rsid.collect()
        SNPs_XF_hwe_cases = hl.variant_qc(mt_XF.filter_cols(mt_XF.is_case == True))
        SNPs_XF_hwe_controls = hl.variant_qc(mt_XF.filter_cols(mt_XF.is_case == False))
        SNPs_XF_hwe_cas = SNPs_XF_hwe_cases.filter_rows(SNPs_XF_hwe_cases.variant_qc.p_value_hwe < hwe_cas).rsid.collect()
        SNPs_XF_hwe_con = SNPs_XF_hwe_controls.filter_rows(SNPs_XF_hwe_controls.variant_qc.p_value_hwe < hwe_con).rsid.collect()

        snps_to_remove_hwe = SNPs_aut_HWE_cas + SNPs_aut_HWE_con + SNPs_XF_hwe_cas + SNPs_XF_hwe_con
        print("Number of autosomal SNPs that fail HWE: " + str(len(SNPs_aut_HWE_cas + SNPs_aut_HWE_con)))
        print("Number of SNPs in females that fail HWE: " + str(len(SNPs_XF_hwe_cas + SNPs_XF_hwe_con)))
        print("Total number of SNPs that fail HWE and are to be removed: " + str(len(snps_to_remove_hwe)))

        if len(snps_to_remove_hwe) > 0:
            mtco = mtco.filter_rows(hl.literal(snps_to_remove_hwe).contains(mtco['rsid']), keep=False)

        # SNPs around the 3:45848429 locus after qc
        chr3_snps_after_qc = hl.filter_intervals(mtco, [hl.parse_locus_interval('chr3:45348429-46348429')]).rsid.collect()

        # logging
        original_stdout = sys.stdout  # Save a reference to the original standard output
        #qc_filename = '/{}'.format(i) + i + '.hwe-filter.log'
        with hadoop_open(dirname + 'qc_step_4/{}/'.format(i) + i + '.hwe-filter.log', 'w') as f:
            sys.stdout = f
            print("\n<============== Allele frequency check information ==============>")
            print("SNPs around the 3:45848429 locus flagged for allele frequency check: {}".format(len(chr3_flagged_snps)))
            print("<============== HWE filter information ==============>")
            print("{} SNPs to be included in the HWE filtering".format(initial_vars))
            print("{} samples to be included in the HWE filtering".format(initial_samples))
            print("Number of autosomal SNPs: {}".format(n_aut))
            print("Number of autosomal SNPs that fail HWE: " + str(len(SNPs_aut_HWE_cas + SNPs_aut_HWE_con)))

            print("Number of SNPs in females: {}".format(n_xf))
            print("Number of SNPs in females that fail HWE: " + str(len(SNPs_aut_HWE_cas + SNPs_aut_HWE_con)))
            print("Total number of SNPs (aut+chromX in females) removed due to HWE filtering: " + str(len(snps_to_remove_hwe)))
            print("Number of SNPs remaining after HWE filtering: {}".format(mtco.count_rows()))
            print("\nSNPs around the 3:45848429 locus after QC: {}".format(len(chr3_snps_after_qc)))
            sys.stdout = original_stdout  # Reset the standard output to its original value i.e. to printing out to the screen

        # write out the HWE filtered information into a mt
        ### NOTE: Remove this if statement because for some cohorts, you will have a length of 0 and no mt will be written
        #if len(snps_to_remove_hwe) > 0:
        name = '.' + i + '.hwe_filtered_qced'
        mt_name = name + '.mt'
        mtco.write(dirname + 'qc_step_4/{}/'.format(i) + i + '.hwe_filtered_qced.mt', overwrite=True)
        # filtered.export(dirname + '{}/'.format(cohortname) + cohortname + '.AlleleFreqCheck.tsv')
        hl.export_plink(mtco, dirname + 'qc_step_4/{}/'.format(i) + i + '.hwe_filtered_qced')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, required=True)
    # parser.add_argument('--basename', type=str)

    parser.add_argument('--mt', type=str, help="Path to the matrix table you want to use for HWE")
    parser.add_argument('--pcasamples', type=str,
                        help="File from PCA containing list of samples to be used in HWE filtering")
    parser.add_argument('--cohort', action='append', help="Name of the cohort you want to run the filters on")
    parser.add_argument('--phenofile', type=str,
                        help="Path to the phenotype file to be used for annotating genotype info with phenotype info")
    parser.add_argument('--ancestry', type=str, default="eur")
    parser.add_argument('--hwe-th-con', type=float, default=1e-6, help="HWE_controls < NUM")
    parser.add_argument('--hwe-th-cas', type=float, default=1e-10, help="HWE_cases < NUM")

    args = parser.parse_args()

    mt = hl.read_matrix_table(args.mt)

    hwe_allele_checks(mt, args.dirname, args.pcasamples, args.cohort, args.phenofile, args.ancestry, args.hwe_th_cas, args.hwe_th_con)


if __name__ == '__main__':
    main()


# EXAMPLE
# Running Allele Frequency and HWE checks
# hailctl dataproc submit hail 4allele_hwe_checks.py
# --dirname gs://dsge-covid19-data/15042021/
# --mt gs://dsge-covid19-data/15042021/cov_15042021.merged.qc.mt
# --pcasamples gs://dsge-covid19-data/15042021/pca_samples_2021_04_12
# --phenofile gs://dsge-covid19-data/phenotypes_main_2021_04_11_wo_mismatch.tsv
# --cohort italy --cohort belgium --cohort egypt --cohort iran --cohort sweden --cohort germany
