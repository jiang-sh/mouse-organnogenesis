#!/bin/sh
#cat organ.list | while read line; do sh bamtosnap.sh ${line}; done

ID=scATAC
batch=core40
mem=90
ppn=12
thread=12

#Sample=ForeBrain_E13.5
#Rscriptid=SnapATAC_Eye_E13_5.R
Sample=$1
resample=${Sample/./_}
Organ=${Sample%_*}
Stage=${Sample#*_}
Rscriptid=scATAC_${resample}
fileid=SnapATAC_${resample}
bash=${Sample}_ArchR_bap2_v2_execute

dsub_script=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/Scripts/dsub_09_bap2/${bash}.sh
barcode_multi=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_barcode_multiplet_ArchR_bap2_v2
bam=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/02_cellranger_atac_output/${Rscriptid}/outs/possorted_bam.bam
barcode=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/08_scATAC_ArchR/scATAC_ArchR_first_QC_tss4frag200_v2/cellsSampleBarcode/tss4_frag200_${Rscriptid}.txt
#barcode=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/02_cellranger_atac_output/${Rscriptid}/outs/raw_peak_bc_matrix/barcodes.tsv
#barcode=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/02_cellranger_atac_output/${Rscriptid}/outs/filtered_peak_bc_matrix/barcodes.tsv

cat >${dsub_script} <<EOF
#!/bin/sh
#PBS -N ${bash}
#PBS -o ${bash}.o
#PBS -e ${bash}.e
#PBS -q ${batch}
#PBS -l mem=${mem}gb,walltime=999:00:00,nodes=1:ppn=${ppn}
#HSCHED -s s+S+S
#PPN limit ${ppn}

### bap2
export PATH=/software/biosoft/software/python/python2019/bin:\$PATH
cd /software/biosoft/software/python/python2019
source venv3/bin/activate

cd ${barcode_multi}
bap2 bam -i ${bam} -o ${Rscriptid}_bap -n ${Rscriptid} -r mm10 -c ${ppn} -bt CB -w ${barcode}

EOF
dsub ${dsub_script} 
