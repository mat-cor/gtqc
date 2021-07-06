import hail as hl
import argparse
import pandas as pd

hl.init(default_reference="GRCh38")

def filter_snps(mt, maf):
    print("Initial number of SNPs before filtering: {}".format(mt.count_rows()))
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    mt = mt.filter_rows(mt.maf > maf)
    print("\nNumber of SNPs after MAF filtering: {}".format(mt.count_rows()))

    # MHC chr6:25-35Mb
    # chr8.inversion chr8:7-13Mb
    intervals = ['chr6:25M-35M', 'chr8:7M-13M']
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in intervals], keep=False)
    print("\nNumber of SNPs after MHC and chr8 inversions filtering: {}".format(mt.count_rows()))

    mt_ld_prune = hl.ld_prune(mt.GT, r2=0.8, bp_window_size=500000)
    mt_ld_pruned = mt.filter_rows(hl.is_defined(mt_ld_prune[mt.row_key]))
    print("\nNumber of SNPs after LD filtering: {}".format(mt_ld_pruned.count_rows()))

    return mt_ld_pruned


def load_ref(dirname, basename):
    """
    Loads a reference plink dataset, writes out a matrix table
    :param dirname: plink file directory name
    :param basename: plink base filename
    :return:
    """
    ref = hl.import_plink(bed=dirname + basename + '.bed',
                          bim=dirname + basename + '.bim',
                          fam=dirname + basename + '.fam',
                          min_partitions=100)
    print("Filtering reference")
    ref = filter_snps(ref, 0.1)
    ref.describe()

    print('sites in ref data: ' + str(ref.count()))  #
    ref.write(dirname + 'qc_step_2/' + basename + '.mt', args.overwrite)


def load_data_plink(dirname, basename):
    """
    Loads data plink dataset, return matrix table
    :param dirname: plink file directory name
    :param basename: plink base filename
    :return: MatrixTable
    """
    mt = hl.import_plink(bed=dirname + basename + '.bed',
                         bim=dirname + basename + '.bim',
                         fam=dirname + basename + '.fam',
                         min_partitions=100)
    print("Filtering data")
    mt = filter_snps(mt, 0.1)
    return mt


def load_data_mt(dirname, basename):
    """
    Loads data plink dataset, return matrix table
    :param dirname: plink file directory name
    :param basename: plink base filename
    :return: MatrixTable
    """
    mt = hl.read_matrix_table(dirname + basename + ".mt")
    print("Filtering data")
    mt = filter_snps(mt, 0.1)
    return mt


def intersect_ref(refdirname, refbasename, data, databasename, outprefix):
    """
    Intersects reference panel with the data and writes intersections as matrix tables
    :param dirname: directory name to put reference and data file intersections
    :param basename: base filename for reference data
    :param data: data
    :return:
    """
    print(refdirname + refbasename)
    this_ref = hl.read_matrix_table(refdirname + refbasename + '.mt')

    # filter data to sites in ref & array data
    data_in_ref = data.filter_rows(hl.is_defined(this_ref.rows()[data.row_key]))
    print('\nsites in ref and data, inds in data: ' + str(data_in_ref.count()))  #

    ##
    data_in_ref.write(outprefix + 'qc_step_2/' + databasename + '_intersect_1000G.mt', args.overwrite)

    # filter ref to data sites
    ref_in_data = this_ref.filter_rows(hl.is_defined(data.rows()[this_ref.row_key]))
    print('\nsites in ref and data, inds in ref: ' + str(ref_in_data.count()))  #
    ##
    ref_in_data.write(outprefix + 'qc_step_2/' + '1000G_intersect_' + databasename + '.mt', args.overwrite)


def run_pca(mt: hl.MatrixTable, out_prefix: str, overwrite: bool = True):
    """
    Run PCA on a dataset
    :param mt: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'qc_step_2/' + '1000G_scores.ht', overwrite)
    pca_scores = hl.read_table(out_prefix + 'qc_step_2/' + '1000G_scores.ht')
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    pca_scores.export(out_prefix + 'qc_step_2/' + '1000G_scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'qc_step_2/' + '1000G_loadings.ht', overwrite)  # PCA loadings


def pc_project(mt: hl.MatrixTable, loadings_ht: hl.Table, loading_location: str = "loadings", af_location: str = "pca_af"):
    """
    Projects samples in `mt` on pre-computed PCs.
    :param MatrixTable mt: MT containing the samples to project
    :param Table loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param str loading_location: Location of expression for loadings in `loadings_ht`
    :param str af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    :rtype: Table
    """
    n_variants = loadings_ht.count()
    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location]
    )
    mt = mt.filter_rows(hl.is_defined(mt.pca_loadings) & hl.is_defined(mt.pca_af) &
                        (mt.pca_af > 0) & (mt.pca_af < 1))
    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(n_variants * 2 * mt.pca_af * (1 - mt.pca_af))
    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))
    return mt.cols().select('scores')


def project_individuals(pca_loadings, project_mt):
    """
    Project samples into predefined PCA space
    :param pca_loadings: existing PCA space
    :param project_mt: matrix table of data to project
    :param project_prefix: directory and filename prefix for where to put PCA projection output
    :return:
    """
    ht_projections = pc_project(project_mt, pca_loadings)
    ht_projections = ht_projections.transmute(**{f'PC{i}': ht_projections.scores[i - 1] for i in range(1, 21)})
    return ht_projections


def main(args):
    if args.load_ref:
        load_ref(args.ref_dirname, args.ref_basename)

    if args.intersect_ref:
        intersect_ref(args.ref_dirname, args.ref_basename, load_data_mt(args.data_dirname, args.data_basename), args.data_basename, args.out_prefix)

    if args.pca_project:
        """
        Compute PCA in global reference panel, project data individuals into PCA space
        """
        ref_in_data = hl.read_matrix_table(args.out_prefix + 'qc_step_2/' + '1000G_intersect_' + args.data_basename + '.mt')
        print('\nComputing reference PCs')
        run_pca(ref_in_data, args.out_prefix)

        # project data
        pca_loadings = hl.read_table(args.out_prefix + 'qc_step_2/' + '1000G_loadings.ht')
        project_mt = hl.read_matrix_table(args.out_prefix + 'qc_step_2/' + args.data_basename + '_intersect_1000G.mt')
        ht = project_individuals(pca_loadings, project_mt)
        ht.export(args.out_prefix + 'qc_step_2/' + args.data_basename + '_scores.tsv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--load_ref', action='store_true')
    parser.add_argument('--ref_dirname', default='gs://dsge-covid19-data/1KG_HGDP/')
    parser.add_argument('--ref_basename', default='1kg_hgdp_100_partitions')
    parser.add_argument('--data_dirname', required=True)
    parser.add_argument('--data_basename', required=True)
    parser.add_argument('--out_prefix', required=True)
    parser.add_argument('--intersect_ref', action='store_true')
    parser.add_argument('--pca_project', action='store_true')

    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)


# # Run from command line
# hailctl dataproc submit hail ancestry_pca.py --intersect_ref --pca_project --overwrite
# --data_dirname gs://dsge-covid19-data/COV_ILLUMINA_15102020/
# --data_basename COV_ILLUMINA_15102020.chr0.pos0.removed.qc
# --out_prefix gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/
