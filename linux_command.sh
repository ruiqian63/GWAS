#path of GWAS 
cd /neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/

#source
source /programs/biogrids.shrc
source ~/.bashrc

#add alias
nano ~/.bashrc
alias folder='cd /neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian'

#require e3 cluster
ssh ch258782@e3-login.tch.harvard.edu
srun -A bch -p bch-interactive -n 6 --mem=24GB --qos=interactive --pty /bin/bash

#scp
scp rui.qian@navarro:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/test/rawdata/fail-het-qc.txt  /Users/Rui/downloads
scp ch258782@e3-login.tch.harvard.edu:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/test_all/gwas_logistic_results.assoc.logistic /Users/Rui/documents/gwas_all
scp /Users/Rui/downloads/1_qc_gwas/inversion.txt rui.qian@navarro:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/test/rawdata
scp /Users/Rui/downloads/folate_pass_pcgc_ctd_srv.vcf ch258782@e3-login.tch.harvard.edu:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/ruiqian/data

scp -r for folder
