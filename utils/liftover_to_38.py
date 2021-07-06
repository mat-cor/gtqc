import hail as hl
import argparse


def plink_to_mt(dirname, basename, reference):
    mt = hl.import_plink(bed=dirname + basename + '.bed',
                         bim=dirname + basename + '.bim',
                         fam=dirname + basename + '.fam',
                         reference_genome=reference)

    mt = mt.repartition(n_partitions=100, shuffle=True)

    mt.write(dirname + basename + '.mt', overwrite=True)


def main():
    parser = argparse.ArgumentParser()
    # reuired for all functions
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str, required=True)

    args = parser.parse_args()

    plink_to_mt(args.dirname, args.basename, "GRCh37")

    mt = hl.read_matrix_table(args.dirname + args.basename + ".mt")

    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

    mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, 'GRCh38', include_strand=True), old_locus=mt.locus)
    mt = mt.filter_rows(hl.is_defined(mt.new_locus) & ~mt.new_locus.is_negative_strand)

    print(mt.describe())

    mt = mt.key_rows_by(locus=mt.new_locus.result, alleles=mt.alleles)

    print(mt.describe())
    print(mt.count())

    mt.write(args.dirname + args.basename + "_GRCh38.mt", overwrite=True)

    hl.export_plink(mt, args.dirname + args.basename + "_GRCh38", varid=mt.rsid)


if __name__ == '__main__':
    main()

# hailctl dataproc submit hail liftover_to_38.py \
# --dirname gs://dsge-covid19-data/controls_DE/MIMICS_HC_autosomes_X/ \
# --basename mimics_HC_autosom_X
