#!/bin/bash

#for SAMPLE in HG002 HG003 HG005 HG006; do
#    echo start $SAMPLE
#    cp ~/MNT/analysis_results/ubuntu/BETA/hg38/daylily/results/day/hg38/RIH0_ANA0-${SAMPLE}-19_DBC0_0/align/strobe/snv/deep/RIH0_ANA0-${SAMPLE}-19_DBC0_0.strobe.deep.snv.sort.vcf.{gz,gz.tbi} inputs &
#    cp ~/MNT/analysis_results/ubuntu/BETA/hg38/daylily/results/day/hg38/RIH0_ANA0-${SAMPLE}-19_DBC0_0/align/strobe/RIH0_ANA0-${SAMPLE}-19_DBC0_0.strobe.mrkdup.sort.{bam,bam.bai} inputs &
#    echo end $SAMPLE
#done
#echo waiting
#wait
#echo done waiting
#sleep 100
#return 1

JOB_ID=$1


#for SAMPLE in HG001 HG004 HG007; do
#    mkdir -p "$PWD/output2/tfrecords/${SAMPLE}"

    #bin/util/submit_make_examples.sh $SAMPLE "--regions chr21 --regions chr20"
    
#    sleep 1.25
#    JOB_ID=$(sbatch --dependency=afterany:$JOB_ID  --comment RandD --partition i192 bin/util/submit_make_examples.sh $SAMPLE  | awk '{print $4}')

#done


mkdir -p output3/training_output

docker run \
  --cpus=16 \
  --memory=224g \
  -v $PWD/output2/tfrecords:/tfrecords \
  -v $PWD/output3:/output3 \
  -v $PWD/output2:/output2 \
  -v /fsx:/fsx \
  daylilyinformatics/deepvariant-avx512:1.5.0 \
  /opt/deepvariant/bin/model_train \
    --model_type=WGSSTROBE \
    --train_dir=/output3/training_output \
    --dataset_config_pbtxt=/output2/dataset_config.pbtxt \
    --batch_size=64 \
    --num_training_steps=50000

#CALL VARS
#docker run \
#  --cpus=192 \
#  --memory=324g \
#  -v /path/to/output:/output \
#  -v /home/ubuntu/MNT/analysis_results:/mnt:ro \
#  daylilyinformatics/deepvariant-avx512:1.5.0 \
#  /opt/deepvariant/bin/call_variants \
#    --outfile /output2/HG001_deepvariant.vcf.gz \
#    --examples /mnt/tfrecords/HG001/examples.tfrecord.gz \
#    --checkpoint /output2/training_output/model.ckpt-50000
