################# 
#################
#################
#################
#################
#################
#  1. Use celescope to convert to 10x format
#  2. Use cell ranger to count
#  3. RNA velocyte analysis
# change the fastq.gz file to fastq_prefix1_1.fq.gz or fastq_prefix1_2.fq.gz

# mapfile refer to https://github.com/singleron-RD/CeleScope/blob/master/docs/rna/multi_rna.md
# convert10x refer to https://github.com/singleron-RD/CeleScope/blob/master/docs/convert10X/multi_convert10X.md
# conda activate celescope
multi_convert10X \
    --mapfile  ./rna.mapfile.txt \
    --thread 30 \
    --ref_path "/media/byron/Reference/refdata-gex-GRCh38-2020-A" \
    --soft_path "/home/Bioinfo/Software/cellranger-6.1.2/cellranger" \
  	--mem 150 \
    --mod shell 
# it generate the code in shell fold

# then run this
for f in ./shell/*.sh; do
	bash "$f"
done

multi_rna	\
	--mapfile ./rna.mapfile.txt \
	--genomeDir /media/byron/Reference/hs_ensembl_99 \
	--thread 15\
	--mod shell

for i in *;
do
echo $i
sh ./shell/$i.sh
done


for i in *;
do
echo ${i}
/home/Bioinfo/Software/cellranger-6.1.2/cellranger count \
--id=${i} \
--transcriptome=/media/byron/Reference/refdata-gex-GRCh38-2020-A \
--fastqs=/media/Pierro_lab/hrr/HRA004266/${i}/02.convert \
--sample=${i} \
--localcores=20 \
--localmem=100
done

for i in HRR*;
do
echo ${i}
velocyto run10x -m /media/byron/Reference/hg38_rmsk.gtf \
	/media/Pierro_lab/hrr/HRAXXXX/cellranger/${i} \
	/media/byron/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf &
done

