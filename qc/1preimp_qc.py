import hail as hl
import sys
import io
import argparse

from hail import hadoop_open
from typing import Tuple, Any, Dict


def summary_stats(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    results = {}

    counts = mt.aggregate_cols(hl.struct(is_case=hl.agg.counter(mt.is_case),
                                         is_female=hl.agg.counter(mt.is_female)))

    mapping = {
        True: 'case',
        False: 'control',
        None: 'unknown'
    }
    is_case_counts = {mapping[pheno]: count for pheno, count in counts['is_case'].items()}
    results['is_case_counts'] = is_case_counts

    mapping = {
        True: 'female',
        False: 'male',
        None: 'unknown'
    }
    is_female_counts = {mapping[pheno]: count for pheno, count in counts['is_female'].items()}
    results['is_female_counts'] = is_female_counts

    n_variants = mt.count_rows()
    n_samples = mt.count_cols()

    results['n_variants'] = n_variants
    results['n_samples'] = n_samples

    pheno_status = ['case', 'control', 'unknown']
    sex_status = ['female', 'male', 'unknown']

    for i in pheno_status:
        if i not in results['is_case_counts']:
            results['is_case_counts'][i] = 0

    for i in sex_status:
        if i not in results['is_female_counts']:
            results['is_female_counts'][i] = 0

    return mt, results


def read_input(dirname, basename, reference, phenofile, inputType) -> hl.MatrixTable:
    if inputType == "plink":
        hl.init(default_reference=reference)
        mt = hl.import_plink(bed=dirname + basename + '.bed',
                             bim=dirname + basename + '.bim',
                             fam=dirname + basename + '.fam')

        mt = mt.repartition(n_partitions=100, shuffle=True)

        mt.write(dirname + 'qc_step_1/' + basename + '.plink.mt', overwrite=True)

        mt = hl.read_matrix_table(dirname + 'qc_step_1/' + basename + ".plink.mt")

    elif inputType == "hail":
        mt = hl.read_matrix_table(dirname + basename + ".mt")

    # Annotate sex information
    pheno = hl.import_table(phenofile,
                            types={'Sex': hl.tint32, 'Analysis_C2': hl.tint32},
                            missing="NA")
    pheno = pheno.key_by('ID')
    mt = mt.annotate_cols(pheno=pheno[mt.s])

    # # Annotate SNP ID as chr:pos
    mt = mt.annotate_rows(SNPID=hl.str(mt.locus))

    # when sex is defined in the phenofile, annotate with sex from there as discrepancies as been already corrected,
    # otherwise keep the sex from the .fam file
    mt = mt.annotate_cols(is_female=hl.cond(hl.is_defined(mt.pheno.Sex),
                                            hl.cond(mt.pheno.Sex == 2, True, False),
                                            mt.is_female))

    # add last input to this function so that the user can specify the analysis type (provided it's in the phenofile)
    mt = mt.annotate_cols(is_case=hl.cond(mt.pheno.Analysis_C2 == 1,
                                          True,
                                          False))

    return mt


def recode_allele_n(dirname, basename, reference, phenofile, inputType) -> hl.MatrixTable:
    """
    Recode SNPs with the first allele equal to "N"
    :param dirname: directory name where the input is
    :param basename: basename of the input
    :param reference: reference (e.g. GRCh38, GRCh37)
    :param phenofile: phenotype file
    :param inputType: plink or hail
    :return: MatrixTable with alleles recoded
    """
    mt = read_input(dirname, basename, reference, phenofile, inputType)
    snps_3alleles = mt.filter_rows(mt.alleles.length() > 2).rsid.collect()
    if len(snps_3alleles) > 0:
        mt = mt.filter_rows(hl.literal(snps_3alleles).contains(mt.rsid), keep=False)
    print("Number of SNPs deleted because with 3 alleles: " + str(len(snps_3alleles)))
    return mt


def pre_qc(dirname, basename, reference, phenofile, inputType):
    """
    pre-QC steps:
    - recode SNPs with the first allele equal to N
    - remove SNPs with call rate < 0.95
    :param dirname: directory name where the input is
    :param basename: basename of the input
    :param reference: reference (e.g. GRCh38, GRCh37)
    :param phenofile: phenotype file
    :param inputType: plink or hail
    :return: MatrixTable with SNPs recoded and qc'ed by call rate
    """
    mt = recode_allele_n(dirname, basename, reference, phenofile, inputType)
    mt, pre_sum_stats = summary_stats(mt)

    # this is required to generate the log file
    old_stdout = sys.stdout  # everything printed from here will be saved to a var
    new_stdout = io.StringIO()
    sys.stdout = new_stdout

    print("<============== Sample and SNP summary (pre-qc) ==============>")
    print("{} SNPs to be included in pre-QC filtering".format(pre_sum_stats['n_variants']))
    print("{} samples to be included in pre-QC filtering".format(pre_sum_stats['n_samples']))
    print("Number of females: {}".format(pre_sum_stats['is_female_counts']['female']))
    print("Number of males: {}".format(pre_sum_stats['is_female_counts']['male']))
    print("Number of samples with unspecified sex: {}".format(pre_sum_stats['is_female_counts']['unknown']))
    print("Number of cases: {}".format(pre_sum_stats['is_case_counts']['case']))
    print("Number of controls: {}".format(pre_sum_stats['is_case_counts']['control']))
    print("Number of samples with unknown phenotype: {}".format(pre_sum_stats['is_case_counts']['unknown']))

    # SNPs on contigs chr1-X
    # Remove contigs other than chr1-22, X (some data contains contigs like 'chr1_KI270706v1_random': 4)
    contigs = ['chr%s' % i for i in range(1, 23)]
    contigs = contigs + ["chrX"]
    mt = mt.filter_rows(hl.literal(contigs).contains(mt.locus.contig))
    snps_filt_contigs = mt.count_rows()  # of snps remaining after filtering by contigs
    print("\n<============== SNP-QC part1 ==============>")
    print("Number of SNPs on contigs chr1-X: {}\t Filtered SNPs: {}"
          .format(snps_filt_contigs, (pre_sum_stats['n_variants'] - snps_filt_contigs)))

    # SNPs with position != 0
    mt = mt.filter_rows(mt.locus.position != 0)
    snps_filt_pos = mt.count_rows()  # of snps remaining after filtering by pos
    print("Number of SNPs with position != 0: {}\t Filtered SNPs: {}".
          format(snps_filt_pos, (snps_filt_contigs - snps_filt_pos)))

    # SNPs with a1 != a2
    mt = mt.filter_rows(mt.alleles[1] != mt.alleles[0])
    snps_filt_alleles = mt.count_rows()  # of snps remaining after filtering by allele
    print("Number of SNPs with a1 != a2: {}\t Filtered SNPs: {}".
          format(snps_filt_alleles, (snps_filt_pos - snps_filt_alleles)))

    # Remove SNPs with missingness > 5%
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)
    snps_filt_call = mt.count_rows()  # of snps remaining after filtering by call rate
    print("Number of SNPs with call rate >= 95%: {}\t Filtered SNPs: {}".
          format(snps_filt_call, (snps_filt_alleles - snps_filt_call)))

    print("Total SNPs removed: {}".format(pre_sum_stats['n_variants'] - snps_filt_call))
    print("Total samples removed: {}".format(pre_sum_stats['n_samples'] - mt.count_cols()))
    print("After pre-QC, there are {} SNPs and {} samples".format(snps_filt_call, mt.count_cols()))

    pre_qc_info = new_stdout.getvalue()
    sys.stdout = old_stdout  # Reset the standard output to its original value i.e. to printing out to the screen

    return mt, pre_qc_info


def get_samples_to_remove(mt, maf, miss, mfthr, dirname):
    """
    TODO: complete
    :param mt:
    :param maf:
    :param miss:
    :param mfthr: male/female threshold like m_f_thr in qc_filter()
    :return: list of samples to remove
    """
    # Consider SNPs with MAF > 0.01 for the sample QC
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    mt.filter_rows(mt.maf > maf)
    mt = hl.sample_qc(mt)

    miss_samples = mt.filter_cols(mt.sample_qc.call_rate < 1 - miss).s.collect()

    imputed_sex = hl.impute_sex(mt.GT)
    discordant_sex_controls = mt.filter_cols(
        (mt.is_case == False) & ((imputed_sex[mt.s].f_stat < mfthr[0]) & (mt.is_female == False) |
        (imputed_sex[mt.s].f_stat > mfthr[1]) & (mt.is_female == True))).s.collect()

    discordant_sex_cases = mt.filter_cols(
        (mt.is_case == True) & ((imputed_sex[mt.s].f_stat < mfthr[0]) & (mt.is_female == False) |
        (imputed_sex[mt.s].f_stat > mfthr[1]) & (mt.is_female == True))).s.collect()

    # discordant_samples = discordant_sex_controls + discordant_sex_cases

    # for tracking purposes, we only take note of ONLY CASES with sex direpancies and DON'T REMOVE THEM (only controls)
    if len(discordant_sex_cases) > 0:
        imputed_sex.filter(hl.literal(discordant_sex_cases).contains(imputed_sex['s'])).export(
            dirname + 'qc_step_1/imputed_sex_discordant.tsv')

    return miss_samples, discordant_sex_controls, discordant_sex_cases


def get_snps_to_remove(mt, miss):
    """

    :param mt: Hail MatrixTable
    :param miss: missingness threshold value
    :return: list of SNPs to remove
    """

    # there is no need to separate between autosomes and X, and females and females here (only in HWE filtering)
    mt = hl.variant_qc(mt)
    SNPs_to_remove = mt.filter_rows(mt.variant_qc.call_rate < 1 - miss).rsid.collect()

    return SNPs_to_remove


def remove_related_samples(mt, dirname):
    """
    :param mt: Hail MatrixTable
    :param dirname: path to the bucket where a file containing sample IDs for removed samples
    :return: writes out list of IDs of samples to be removed
    """

    # select only common variants (ALT AF > 10%)
    mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.1)

    # perform sampleQC so we can get call rates for comparison in the end
    mt = hl.sample_qc(mt)

    # STEP 1: Filter out regions with high LD
    high_ld_regions = hl.import_locus_intervals('gs://dsge-covid19-data/high-ld-regions_b38.tsv',
                                                reference_genome='GRCh38')
    mt = mt.filter_rows(hl.is_missing(high_ld_regions[mt.locus]))

    # STEP 2: Remove SNPs that are in LD
    # multiallelic variants must be filtered out or split before being passed to ld_prune()
    biallelic_dataset = mt.filter_rows(hl.len(mt.alleles) == 2)
    pruned_variant_table = hl.ld_prune(biallelic_dataset.GT, r2=0.2, bp_window_size=250000)
    mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

    # STEP 3: Run IBD
    # if the field with maf is not specified, Hail will compute it, resulting in more computation
    mt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    relatedness_check = hl.identity_by_descent(mt, maf=mt['maf']) # this returns a Hail Table with the sample pairs
    samples_to_remove = relatedness_check.filter(relatedness_check.ibd.PI_HAT > 0.98)

    # _localize=False means don't put this in Python, keep it as a Hail expr
    call_rate_dict = mt.aggregate_cols(hl.dict(hl.agg.collect((mt.s, mt.sample_qc.call_rate))), _localize=False)
    samples_to_remove = samples_to_remove.annotate(cr_s1=call_rate_dict[samples_to_remove.i],
                                                   cr_s2=call_rate_dict[samples_to_remove.j])
    samples_list = samples_to_remove.annotate(sample_to_remove=hl.cond(
        samples_to_remove.cr_s1 >= samples_to_remove.cr_s2, samples_to_remove.i, samples_to_remove.j))
    samples = samples_list.sample_to_remove.collect()

    if len(samples) > 0:
        with hadoop_open(dirname + 'qc_step_1/relatedness_removed_samples.tsv', 'w') as f:
            for sample in samples:
                f.write(sample + "\n")

    return samples


def qc_filter(dirname, basename, reference, maf, sample_miss, snp_miss, m_f_thr, phenofile, inputType):
    """

    :param dirname:
    :param basename:
    :param reference:
    :param maf:
    :param sample_miss:
    :param snp_miss:
    :param m_f_thr: male/female threshold like mfthr in get_samples_to_remove()
    :return: writes out a qced Hail matrix table
    """
    if m_f_thr is None:
        m_f_thr = [0.5, 0.5]
    mt, pre_qc_log = pre_qc(dirname, basename, reference, phenofile, inputType)
    initial_vars = mt.count_rows()  # initial # of variants
    initial_samples = mt.count_cols()  # initial # of samples

    sample_missingness, sample_sex_check_controls, sample_sex_check_cases = get_samples_to_remove(mt, maf, sample_miss,
                                                                                                  m_f_thr, dirname)
    related_samples = remove_related_samples(mt, dirname)
    samples_to_remove = sample_missingness + sample_sex_check_controls + related_samples

    # remove samples with high missingness
    snps_to_remove = get_snps_to_remove(mt, snp_miss)

    # Remove samples and variants that fail QC
    if len(samples_to_remove) > 0:
        mt = mt.filter_cols(hl.literal(samples_to_remove).contains(mt['s']), keep=False)
    if len(snps_to_remove) > 0:
        mt = mt.filter_rows(hl.literal(snps_to_remove).contains(mt['rsid']), keep=False)

    n_filt = mt.count_rows()  # number of SNPs after filtering out for missingness and sample

    # Remove SNPs that end up to be monomorphic (MAC = 0)
    mt = mt.annotate_rows(MAC=hl.min(mt.variant_qc.AC))
    mt = mt.filter_rows(mt.MAC > 0)
    n_monomorphic = n_filt - mt.count_rows()

    original_stdout = sys.stdout  # Save a reference to the original standard output
    with hadoop_open(dirname + 'qc_step_1/qc-filter.log', 'w') as f:
        f.write(pre_qc_log)
        sys.stdout = f
        print("\n<============== Sample-QC ==============>")
        print("{} SNPs to be included in the QC".format(initial_vars))
        print("{} samples to be included in the QC".format(initial_samples))
        print("Samples removed for missingness: " + str(len(sample_missingness)))
        print("Samples removed for duplicate sample check (IBD>0.98): " + str(len(set(related_samples))))
        print("Control samples removed for sex check: " + str(len(sample_sex_check_controls)))
        print("Number of Case samples with sex disrepancies: " + str(len(sample_sex_check_cases)))
        print("\n<============== SNP-QC part2 ==============>")
        print("Number of SNPs removed for missingness: " + str(len(snps_to_remove)))
        print("Monomorphic SNPs removed: {}".format(n_monomorphic))
        print("\n<============== Summary ==============>")
        print("Total SNPs removed: {}".format(initial_vars - mt.count_rows()))
        print("Total samples removed: {}".format(initial_samples - mt.count_cols()))
        print("After QC-filtering, there are {} SNPs and {} samples".format(mt.count_rows(), mt.count_cols()))
        sys.stdout = original_stdout

    mt.write(dirname + 'qc_step_1/' + basename + '.qc.mt', overwrite=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str, required=True)
    parser.add_argument('--inputType', type=str, required=True, choices=['hail','plink'])

    # required for QC
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--maf', default=0.01, type=float)
    parser.add_argument('--sampleMiss', default=0.02, type=float)
    parser.add_argument('--snpMiss', default=0.02, type=float)
    parser.add_argument('--mfthr', default=None)
    parser.add_argument('--phenofile', type=str,
                        help="Path to the phenotype file to be used for annotating genotype info with phenotype info")

    args = parser.parse_args()

    print("<==================== Running QC filtering without an HWE filters ====================>")
    qc_filter(args.dirname, args.basename, args.reference, args.maf, args.sampleMiss,
                args.snpMiss, args.mfthr, args.phenofile, args.inputType)


if __name__ == '__main__':
    main()

# hailctl dataproc submit CLUSTER_NAME SCRIPT [optional args to your python script...]

# EXAMPLES
# (1) Running qc without any HWE
# hailctl dataproc submit hail 1preimp_qc.py \
# --dirname gs://lindo/test/ \
# --basename plus_cov_illumina_14082020.chr0.pos0.removed \
# --phenofile gs://lindo/test/cases_controls_phenotypes_ita_be_brazil_swe_ger_27_10_2020
# --inputType plink
