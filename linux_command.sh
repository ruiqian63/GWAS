#path of GWAS 
cd /neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian_oldgwas/
cd /neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas
#source
source /programs/biogrids.shrc
source ~/.bashrc

#add alias
nano ~/.bashrc
alias folder='cd /neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian'

#require e3 cluster
ssh ch258782@e3-login.tch.harvard.edu
srun -A bch -p bch-interactive -n 6 --mem=56GB --qos=interactive --pty /bin/bash

#scp
scp rui.qian@navarro:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/test/rawdata/fail-het-qc.txt  /Users/Rui/downloads
scp ch258782@e3-login.tch.harvard.edu:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/gwas_prepost/sig_locus_mt_r2_no_chr6x.ld /Users/Rui/documents/GWAS/redo_no_chr6X/pre
scp ch258782@e3-login.tch.harvard.edu:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/rawdata/merged_1.vcf.gz /Users/Rui/documents/GWAS/redo_no_chr6X

scp /Users/Rui/downloads/1_qc_gwas/inversion.txt rui.qian@navarro:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/test/rawdata
scp /Users/Rui/documents/GWAS/redo_no_chr6X/pre/sig_locus.snplist ch258782@e3-login.tch.harvard.edu:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/gwas_prepost


scp -r for folder
