# This script is for downloading gnomad AF data
# examples on how to check AF index for populations
# (1) gnomad_ht.freq_index_dict[‘gnomad_afr’].show() # index 6
# (2) gnomad_ht.freq_index_dict[‘gnomad_nfe’].show() # index 2
# (3) gnomad_ht.freq_index_dict[‘gnomad_amr’].show() # index 8

# important: use hailctl dataproc start hail --region=europe-west1 --zone=europe-west1-b --requester-pays-allow-all to start the
# cluster as you should pay for accessing the gnomAD data

import hail as hl
# write HailTable with gnomAD data selecting only locus, alleles and population specific AF columns columns
gnomad_ht = hl.experimental.load_dataset(name='gnomad_genome_sites',
                                          version='2.1.1',
                                          reference_genome='GRCh38',
                                          region='us',
                                          cloud='gcp')
#print(gnomad_ht.count())
gnomad_ht.select(amr_AF=gnomad_ht.freq[8].AF).repartition(100).write('gs://dsge-covid19-data/gnomad/gnomad_amr_AF.ht', overwrite = True)
