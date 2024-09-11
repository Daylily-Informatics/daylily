#!/bin/bash

# https://github.com/google/deepvariant/blob/r1.6/docs/deepvariant-training-case-study.md


mkdir -p output2/model

docker run \
  --cpus=16 \
  --memory=224g \
  -v $PWD/output2:/output2 \
  -v /fsx:/fsx \
  daylilyinformatics/deepvariant-avx512:1.5.0 \
    /opt/deepvariant/bin/model_train \
    --model_name=WGSSTROBE \
    --batch_size=4096 \
    --dataset_config_pbtxt=/output2/tfrecords/training/training_set.pbtxt \
    --eval_config_pbtxt=/output2/tfrecords/validation/validation_set.pbtxt\ 
    --number_of_steps=8000000 \
    --train_dir=/output2/model \
    --model_name=inception_v3 \
    --max_checkpoints_to_keep=10 \
    --start_from_checkpoint=model_default




