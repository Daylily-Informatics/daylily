#!/bin/bash



#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/training/HG00?/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/training/HG00?/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/training/training_set.pbtxt" \
#--output_dataset_name="training" \
#--direct_num_workers=192 \
#--step=-1


python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
--input_pattern_list="output2/tfrecords/validation/HG00?/examples.tfrecord.gz"  \
--output_pattern_prefix="output2/tfrecords/validation/HG00?/training_set.with_label.shuffled" \
--output_dataset_config="output2/tfrecords/validation/validation_set.pbtxt" \
--output_dataset_name="validation" \
--direct_num_workers=192 \
--step=-1



exit 0
#
#https://github.com/GuillaumeHolley/TFrecordShuffler

#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG002/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG002/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG002/training_set.pbtxt" \
#--output_dataset_name="HG002" \
#--direct_num_workers=192 \
#--step=-1


#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG003/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG003/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG003/training_set.pbtxt" \
#--output_dataset_name="HG003" \
#--direct_num_workers=192 \
#	--step=-1


#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG005/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG005/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG005/training_set.pbtxt" \
#--output_dataset_name="HG005" \
#--direct_num_workers=192 \
#	--step=-1


#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG006/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG006/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG006/training_set.pbtxt" \
#--output_dataset_name="HG006" \
#--direct_num_workers=192 \
#--step=-1



#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG004/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG004/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG004/training_set.pbtxt" \
#--output_dataset_name="HG004" \
#--direct_num_workers=192 \
#--step=-1

#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG001/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG001/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG001/training_set.pbtxt" \
#--output_dataset_name="HG001" \
#--direct_num_workers=192 \
#--step=-1


#python3 ../TFrecordShuffler/shuffle_tfrecords_lowmem.py \
#--input_pattern_list="output2/tfrecords/HG007/examples.tfrecord.gz"  \
#--output_pattern_prefix="output2/tfrecords/HG007/training_set.with_label.shuffled" \
#--output_dataset_config="output2/tfrecords/HG007/training_set.pbtxt" \
#--output_dataset_name="HG007" \
#--direct_num_workers=192 \
#--step=-1
