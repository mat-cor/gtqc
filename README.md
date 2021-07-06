## QC
1\. Run **QC** on the whole dataset using the `1preimp_qc.py` script. Here's an example below on how to run it:
```
hailctl dataproc submit hail 1preimp_qc.py \
--dirname gs://dsge-covid19-data/COV_ILLUMINA_15102020/ \
--basename COV_ILLUMINA_15102020.chr0.pos0.removed \
--phenofile gs://dsge-covid19-data/COV_ILLUMINA_15102020/cases_controls_phenotypes_ita_be_brazil_swe_ger_29_10_2020 \
--inputType plink [or hail]
```
 * **Input(s):** PLINK .bed, .bim, and .fam files OR a Hail MatrixTable 
 * **Output(s):** (1) One or Two (if you use plink files as input) Hail MatrixTables, use the one with suffix .qc.mt for steps 2 and after; (2) log file; (3) a tsv file with all the samples removed due to relatedness, if there are any; (4) a tsv file containing of samples that failed sex check, if there are any  
 * Some of the QC checks included in step 1 above are: (1) sample/ variant missingness rate; (2) MAF checks; (3) sample relatedness (IBD); (4) sample sex check  

2\. Run **ancestry PCA** using the `2ancestry_pca.py` script
```
hailctl dataproc submit hail 2ancestry_pca.py \
--intersect_ref --pca_project --overwrite \
--data_dirname gs://dsge-covid19-data/COV_ILLUMINA_15102020/ \
--data_basename COV_ILLUMINA_15102020.chr0.pos0.removed.qc \
--out_prefix gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/
```
 * **Input(s):** qced Hail MatrixTable (from step1 above) 
 * **Output(s):** (1) Two PCA scores files, one reference (.txt.bgz) and one for the data (.tsv)
 
3\. Run **ancestry assignment** using the `3assign_pops.py` script. The script assigns population labels based on the results of PCA. Here is an example below:
 ```
hailctl dataproc submit hail 3assign_pops.py \
--ref_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/COV_ILLUMINA_15102020.chr0.pos0.removed.qc_data_scores.txt.bgz \
--data_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/data_COV_ILLUMINA_15102020.chr0.pos0.removed.qc_cases_controls_scores.tsv \
--out_dir gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/
```
* **Input(s):** Reference and data PCA scores files from step 2 above
* **Output(s):** (1) Two scores files, one at p > 0.5 and another at p > 0.8, which can then be used for PCA plots.

4\. Run **Allele frequency check against gnomAD and HWE filtering** on each cohort separately using the `4allele_hwe_checks.py` script. Here is an example below:
 ```
hailctl dataproc submit hail 4allele_hwe_checks.py \
--dirname gs://lindo/test/ \
--basename plus_cov_illumina_14082020.chr0.pos0.removed \
--mt gs://lindo/test/plus_cov_illumina_14082020.chr0.pos0.removed.qc.mt \
--pcasamples gs://lindo/test/eur.samples.to.keep \
--phenoType Control [or Case] \
--cohortsFile gs://lindo/test/cohorts.txt
```
* **Input(s):** (1) QCed Hail MatrixTable (from step 1); (2) file, without a header, containing list of samples from PCA to be used in HWE filtering (one sampleID per line); (3) text file containing the name of cohorts you want to run HWE on (one cohort name each row)
* **Output(s):** (1) HWE QCed Hail MatrixTable for each cohort specified in the cohorts file input; (2) a directory named after the cohort with log file and a png file showing the filtered SNPs during allele frequency check against gnomAD.


## Preparing data for imputation with TopMed server
https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/

The tools necessary for these steps are installed on a VM called **pre-imp**, in my home dir (/home/cordioli/). Add the following lines to your `~/.bashrc` file so that you can run the commands (e.g. plinks) without having to specify the whole path:

`alias plink='/home/cordioli/plink'`   
`alias bcftools='/home/cordioli/bcftools/bcftools'`   
`alias qctool='/home/cordioli/qctool/build/release/qctool_v2.0.7'`   

**Input:** .fam,.bed,.bim files for the Qc'ed data  
**Output:** VCFs (one per each chromosome) to submit to the imputation server

0. Rename chr names (from chrN to N, as in the ref panel):  
(if duplicates:`cut -f 2 <> | sort | uniq -d > dups`)
`plink --bfile <> --update-chr /home/cordioli/ucsc2ensembl.txt --exclude dups --make-bed --out <>`
1. Create frequency file:  
`plink --freq --bfile <input> --out <freq-file>` 
2. Execute check script:  
`perl /home/cordioli/HRC-1000G-check-bim.pl -b <bim file> -f <freq-file> -r /home/cordioli/PASS.Variantsbravo-dbsnp-all.tab.gz -n -h`   
3. The perl script ran in the previous step generates a bash file with a set of plink commands to update or remove SNPs, split the bfile and create a VCF per each chromosome.
4. Sort VCFs  
`for i in {1..23}; do bcftools sort <PREFIX>-chr${i}.vcf -Oz -o <PREFIX>-chr${i}.vcf.gz & done`
5.  _"If your input data is GRCh38/hg38 please ensure chromosomes are encoded with prefix 'chr' (e.g. chr20)"_ :  
`for i in {1..23}; do bcftools annotate -Oz --rename-chrs /home/cordioli/ensembl2ucsc.txt <PREFIX>-chr${i}.vcf.gz > <PREFIX>_chr${i}.vcf.gz & done`
6. Submit to TopMed server, selecting "rsq Filter 0.3" to have the output already filtered for imputation score.
