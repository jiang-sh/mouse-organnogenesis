#!/bin/sh
#cat scATAC_organ.list | while read sample; do sh scATAC-seq_cellranger-atac_mm.sh ${sample}; done
ID=scATAC
batch=core40
mem=90
ppn=12
thread=12

#Sample=ForeBrain_E13.5
Sample=$1
Organ=${Sample%_*}
Stage=${Sample#*_}
fq_name=${ID}_${Organ}
outputid=${ID}_${Sample/./_}
bash=${ID}_${Sample}

fastqs=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_rawdata/scATAC-seq/${Sample}
ref=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_reference/refdata-cellranger-atac-mm10-1.2.0
result=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/02_cellranger_atac_output

dsub_script=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/Scripts/dsub_01_cellranger_atac/${bash}.sh

cat >${dsub_script} <<EOF
#!/bin/sh
#PBS -N ${bash}
#PBS -o ${bash}.o
#PBS -e ${bash}.e
#PBS -q ${batch}
#PBS -l mem=${mem}gb,walltime=999:00:00,nodes=1:ppn=${ppn}
#HSCHED -s Organogenesis+cellrangeratac+mm10
#PPN limit ${ppn}

# make environment
export PATH=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_software/cellranger-atac-1.2.0:\$PATH

cd ${result}
cellranger-atac count --id ${outputid} \
                      --fastqs ${fastqs} \
                      --sample ${fq_name} \
                      --reference ${ref}  \
                      --localcores ${thread} \
                      --localmem ${mem}
EOF

dsub ${dsub_script}
