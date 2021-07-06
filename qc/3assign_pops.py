import sys
import subprocess
import pkg_resources
import argparse
import hail as hl


required = {'sklearn', 'matplotlib'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing])

import pandas as pd
import numpy as np
from collections import defaultdict, namedtuple, OrderedDict
from typing import *
import random
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def load_merge_data(refscores, ref_info, datascores):

    print('\nMerging data with ref')
    ref = pd.read_table(refscores, header=0, sep='\t', compression='gzip')
    ref_info = pd.read_table(ref_info, header=0, sep='\t', low_memory=False)
    ref_info = ref_info[['project_meta.sample_id', 'hgdp_tgp_meta.Population', 'hgdp_tgp_meta.Genetic.region']]
    data = pd.read_table(datascores, header=0, sep='\t')

    # rename columns in ref_info
    ref_info.columns = ['s', 'Population', 'SuperPop']

    ref_merge = pd.merge(left=ref, right=ref_info, left_on='s', right_on='s', how='inner')

    data_ref = pd.concat([ref_merge, data], sort=False)
    print('Done merging data with ref')

    return data_ref


def assign_population_pcs(
        pop_pc_pd: pd.DataFrame,
        num_pcs: int,
        known_col: str = 'SuperPop',
        fit: RandomForestClassifier = None,
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth'
) -> Tuple[pd.DataFrame, RandomForestClassifier]:
    """
    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.
    :param Table pop_pc_pd: Pandas dataframe containing population PCs as well as a column with population labels
    :param str known_col: Column storing the known population labels
    :param str pcs_col: Columns storing the PCs
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int num_pcs: number of population PCs on which to train the model
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Dataframe containing sample IDs and imputed population labels, trained random forest model
    :rtype: DataFrame, RandomForestClassifier
    """

    # Expand PC column
    # pop_pc_pd = expand_pd_array_col(pop_pc_pd, pcs_col, num_pcs, 'PC')
    pc_cols = ['PC{}'.format(i + 1) for i in range(num_pcs)]
    print(pc_cols)
    # pop_pc_pd[pc_cols] = pd.DataFrame(pop_pc_pd[pcs_col].values.tolist())[list(range(num_pcs))]
    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

    N = len(train_data)
    print(train_data.shape)

    # Split training data into subsamples for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(list(range(0, N)), int(N * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit['s']]
        evaluate_fit = train_data.loc[~train_data['s'].isin(fit_samples)]

        # Train RF
        training_set_known_labels = train_fit[known_col].values
        training_set_pcs = train_fit[pc_cols].values
        evaluation_set_pcs = evaluate_fit[pc_cols].values

        pop_clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print('Random forest feature importances are as follows: {}'.format(pop_clf.feature_importances_))

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(len(predictions))
        print('Estimated error rate for RF model is {}'.format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    print('Classifying data')
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].values)
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].values)
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])
    print('probs shape ' + str(probs.shape))
    print('pop_pc_pd shape ' + str(pop_pc_pd.shape))
    print(probs.iloc[:3,])
    print(pop_pc_pd.iloc[:3,])
    pop_pc_pd = pd.concat([pop_pc_pd.reset_index(drop=True), probs.reset_index(drop=True)], axis=1)
    print(pop_pc_pd.shape)
    print(pop_pc_pd.iloc[:3, ])
    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label
    #pop_pc_pd = pop_pc_pd.drop(pc_cols, axis='columns')
    print(pop_pc_pd.shape)

    return pop_pc_pd, pop_clf


def pca_plot(pcs_data: str = None, ref_data: str = None, ref_info_data: str = None, include_ref: bool = None,
             cohort: str = None, out_dir: str = None, prob: float = None, phenotypes_file: str = None, pheno: str = None):

    ref = pd.read_table(ref_data, header=0, sep='\t', compression='gzip')
    ref_info = pd.read_table(ref_info_data, header=0, sep='\t', low_memory=False)
    ref_info = ref_info[['project_meta.sample_id', 'hgdp_tgp_meta.Population', 'hgdp_tgp_meta.Genetic.region']]
    ref_info.columns = ['s', 'Population', 'SuperPop']
    ref_update = pd.merge(ref, ref_info, how='left', on=['s'])

    pcs = pd.read_table(pcs_data, header=0, sep='\t')

    pheno_file = pd.read_table(phenotypes_file, header=0, sep='\t')
    pheno_file.rename(columns={'ID': 's'}, inplace=True)

    if cohort:
        pcs_df = pd.merge(pcs, pheno_file, how='left', on=['s'])
        pcs_df = pcs_df[pcs_df.Cohort == cohort]
        # Only retain samples with phenotype
        pcs_df = pcs_df[pcs_df[pheno].notna()]
    else:
        pcs_df = pcs
        pcs_df = pd.merge(pcs_df, pheno_file, how='left', on=['s'])
        # Only retain samples with phenotype
        pcs_df = pcs_df[pcs_df[pheno].notna()]

        # ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    cbPalette = {'AFR': "#999999", 'AMR': "#E69F00", 'EAS': "#56B4E9", 'EUR': "#009E73",
             'oth': "#F0E442", 'CSA': "#0072B2", 'MID': "#EEA9B8", 'OCE': "#0072B2", 'SAS': "#0072B2"}

    if include_ref is True:
        cbPalette = cbPalette
        # ref_update['col']=ref_update['SuperPop'].map(cbPalette)
        ref_update=ref_update.dropna()

        if cohort:
            plt_title = f'{cohort} projected against 1KG (p > {prob})'
            outfile = f'{cohort}.1KG.PC1-PC6.prob_{prob}.png'
        else:
            plt_title = f'All cohorts projected against 1KG (p > {prob})'
            outfile = f'all.cohorts.1KG.PC1-PC6.prob_{prob}.png'
    else:
        plt_title = f'{cohort} (p > {prob})'
        outfile = f'{cohort}.PC1-PC6.prob_{prob}.png'
        cbPalette = {k: cbPalette[k] for k in list(pcs_df['pop'].value_counts().index)}
        # pcs_df['col']=pcs_df['pop'].map(cbPalette)
        # pcs_df = pcs_df.dropna()

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(25, 25))

    for i in range(1,6):
        if i == 1 or i == 2:
            axisx = 0
        elif i == 3 or i == 4:
            axisx = 1
        else:
            axisx = 2

        if i == 1 or i == 3 or i == 5:
            axisy = 0
        else:
            axisy = 1

        if include_ref:
            axs[axisx, axisy].scatter(ref_update[f'PC{i}'], ref_update[f'PC{i+1}'], c=ref_update['SuperPop'].map(cbPalette),
                                      s=5, alpha=0.1)
        axs[axisx, axisy].scatter(pcs_df[f'PC{i}'], pcs_df[f'PC{i+1}'], c=pcs_df['pop'].map(cbPalette), s=5, alpha=1)
        axs[axisx, axisy].set_xlabel(xlabel=f'PC{i}', fontsize=15)
        axs[axisx, axisy].set_ylabel(ylabel=f'PC{i+1}', fontsize=15)
        axs[axisx, axisy].tick_params(axis='both', which='major', labelsize=12)
    handles = []

    # get population counts so we can add them to legend
    pop_counts = (pcs_df['pop'].value_counts(sort=True)).to_dict()

    for key in cbPalette:
        # if the key is not in the dict, add it
        if key not in pop_counts:
            pop_counts[key] = 0
        # manually define a new patch
        data_key = Line2D([0], [0], marker='o', color='w', label='{} (n={})'.format(key, pop_counts.get(key)),
                        markerfacecolor=cbPalette[key], markersize=10)
        handles.append(data_key)

    fig.legend(handles=handles, title='Populations', bbox_to_anchor=(0.7, 0.2), loc='lower left', frameon=False)
    fig.delaxes(axs[2][1])
    fig.suptitle(plt_title, fontsize=25)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.close()
    fig.savefig(f'/tmp/{outfile}', dpi=300, facecolor='w')
    hl.hadoop_copy(f'file:///tmp/{outfile}', f'{out_dir}plots/{outfile}')


def main(args):
    # error handling -----------------------------------------------------------
    if not args.ref_scores:
        raise Exception("Reference scores file not selected")

    if not args.data_scores:
        raise Exception("Data scores file not selected")

    if not args.dirname:
        raise Exception("Output directory where files will be saved is not specified")

    ref_scores = args.dirname + 'qc_step_2/1000G_scores.txt.bgz'
    ref_info = args.ref_info
    outdirectory = args.dirname + 'qc_step_3/'

    data_ref = load_merge_data(args.ref_scores, args.ref_info, args.data_scores)

    rf_thresh = [0.5, 0.8]
    for threshold in rf_thresh:
        pcs_df, clf = assign_population_pcs(pop_pc_pd=data_ref, num_pcs=20, min_prob=threshold)

        data_pops = pcs_df.loc[pcs_df['SuperPop'].isnull()]
        # print(data_pops['pop'].value_counts())
        superpops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
        # ["EUR", "EAS", "AMR", "AFR", "SAS"]
        cols = ['s', 'pop'] + [f'prob_{i}' for i in superpops] + [f'PC{i}' for i in range(1, 21)]
        data_pops_df = data_pops[cols]
        print(data_pops_df)

        data_pops_df.to_csv('{}pca_sup_pops_{}_probs.txt'.format(outdirectory, threshold),
                            sep='\t', index=False)

        cohorts = args.cohort

        # PCA plot for all cohorts projected against 1KG
        print("Generating PCA plots")
        pca_plot(pcs_data = '{}pca_sup_pops_{}_probs.txt'.format(outdirectory, threshold), ref_data=ref_scores,
        ref_info_data=ref_info, include_ref=True, out_dir=outdirectory, prob=threshold,
        phenotypes_file=args.phenotype_file, pheno=args.phenotype)

        for i in cohorts:
            # PCA for each cohort separately projected against 1KG
            pca_plot(pcs_data = '{}pca_sup_pops_{}_probs.txt'.format(outdirectory, threshold), ref_data=ref_scores,
            ref_info_data=ref_info, include_ref=True, cohort = i, out_dir=outdirectory, prob=threshold,
            phenotypes_file=args.phenotype_file, pheno=args.phenotype)

            # # PCA for each cohort separately
            pca_plot(pcs_data = '{}pca_sup_pops_{}_probs.txt'.format(outdirectory, threshold), ref_data=ref_scores,
            ref_info_data=ref_info, include_ref=False, cohort = i, out_dir=outdirectory, prob=threshold,
            phenotypes_file=args.phenotype_file, pheno=args.phenotype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_scores', type=str, required=True)
    parser.add_argument('--ref_info', default='gs://dsge-covid19-data/1KG_HGDP/gnomad_meta_hgdp_tgp_v1.txt')
    parser.add_argument('--data_scores', type=str, required=True)
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--phenotype_file', type=str, required=True)
    parser.add_argument('--phenotype', default='Analysis_C2')
    parser.add_argument('--cohort', action='append', help="Name of the cohort you want to generate PCA plots for")

    args = parser.parse_args()
    main(args)


# EXAMPLE
# hailctl dataproc submit hail assign_pops.py
# --ref_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/1000G_scores.txt.bgz
# --data_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/COV_ILLUMINA_15102020.chr0.pos0.removed.qc_scores.tsv
# --out_dir gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/
