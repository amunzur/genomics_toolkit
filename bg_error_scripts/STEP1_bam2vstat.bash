#!/bin/bash
#SBATCH --job-name=vstat
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --export=all
#SBATCH --chdir=/groups/wyattgrp/log
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log
#SBATCH --array=1-26

#printf "SLURM_JOB_ID=$SLURM_JOB_ID\n";

ref="/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
outputdir=/groups/wyattgrp/users/amunzur/error_rates
vstatdir=${outputdir}/vstat;
chrom=$(printf "$(seq 1 22)\nX\nY\nM\n" | sed 's/^/chr/g');
conda_profile_path="/home/amunzur/anaconda3/etc/profile.d/conda.sh"
min_depth=20

shopt -s extglob;
mkdir -p ${outputdir};
mkdir -p ${vstatdir};
tmp_dir=$(mktemp -d -p ${vstatdir});
cd ${tmp_dir};

#bam_list should be realpath of files
bam_list="/groups/wyattgrp/users/amunzur/error_rates/vip_normals_wbc_bam_list.txt"
all_bams=($(cat ${bam_list}))
bam=${all_bams[$SLURM_ARRAY_TASK_ID - 1]}

sample=$(basename ${bam} | sed 's/.bam//g');

printf "#########################################################################\n";
echo "# start running pysamstats on $sample file at $(date +"%D %H:%M")";
echo "# from $bam"; 
printf "#########################################################################\n";

source $conda_profile_path
conda activate pysamstats
echo "$chrom" | parallel --env conda_profile_path --env ref --env bam --env tmp_dir -j 0 "source $conda_profile_path; conda activate pysamstats; pysamstats --type variation -f $ref --no-dup -c {} $bam | grep -v -e '_K' -e '_G' -e '_J' -e 'chrM' > $tmp_dir/{}.vstat";

#concatenate the chroms while removing bases 
head -n1 $(ls *.vstat | head -n1) > allchr.vstat;
for i in $(echo $chrom); do
	echo "concatenating chromosome $i at $(date +"%D %H:%M")";
	tail -n +2 $i.vstat | awk -v mind=$min_depth '$4 > mind' >> allchr.vstat;
done;

echo "start adding sample name and linux formatting on allchr.vstat tumor at $(date +"%D %H:%M")";
paste <(head -n1 allchr.vstat | tr -d '\r') <(echo sample) > allchr_tmp.vstat;
cat allchr.vstat | tr -d '\r' | tail -n +2 | awk -F $'\t' -v sample=$sample '{print $0,sample}' FS="\t" OFS="\t" >> allchr_tmp.vstat;

echo "$chrom" | parallel "rm ./{}.vstat";
mv allchr_tmp.vstat ${vstatdir}/${sample}.vstat;

printf "#########################################################################\n";
echo "# finished making vstat for $sample at $(date +"%D %H:%M")";
vstat=$(readlink -ve ${vstatdir}/${sample}.vstat); 
printf "# $vstat\n"
printf "#########################################################################\n";
printf "\n"
done;

rm -r $tmp_dir;