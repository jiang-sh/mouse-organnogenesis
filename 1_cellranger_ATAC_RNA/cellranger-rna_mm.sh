#!/bin/bash
#cat scRNA_organ.list | while read sample; do sh 01_cellranger_rna.sh ${sample}; done

ID=scRNA
batch=core24
mem=60
ppn=12
species=mm10

#Sample=Gonad-male_E11.5
Sample=$1
Organ=${Sample%_*}
Stage=${Sample#*_}
fq_name=${ID}_${Organ}_5exp
outputid=scRNA_${Sample/./_}

fastqs=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_rawdata/scRNA-seq/${Sample}
ref=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_reference/refdata-cellranger-mm10-3.0.0
dsub_script=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/Scripts/dsub_01_cellranger_rna/scRNA_${Sample}.sh
result=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/02_cellranger_rna_output

cat >${dsub_script} <<EOF
#!/bin/bash
#PBS -N 10X_scRNA_${Sample}
#PBS -q ${batch}
#PBS -l mem=${mem}gb,walltime=999:00:00,nodes=1:ppn=${ppn}
#PBS -o 10X_scRNA_${Sample}.001.log
#PBS -e 10X_scRNA_${Sample}.001.err
#PBS -V
#HSCHED -s Organogenesis+cellranger+mm10

cd ${result}
cellranger count --id ${outputid} \
		--fastqs ${fastqs} \
		--sample ${fq_name} \
		--transcriptome ${ref} \
 		--localcores ${ppn} \
		--localmem ${mem}

EOF

dsub ${dsub_script}
