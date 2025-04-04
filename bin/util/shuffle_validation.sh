#!/bin/bash


python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
--input_pattern_list="output2/tfrecords/validation/HG00?/examples.tfrecord.gz"  \
--output_pattern_prefix="output2/tfrecords/validation/training_set.with_label.shuffled" \
--output_dataset_config="output2/tfrecords/validation/validation_set.pbtxt" \
--output_dataset_name="validation" \
--direct_num_workers=192 \
--step=-1



exit 0
