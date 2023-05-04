# Daylily Directory Structure

The `results` [directory structure](../images/assets/day_files_tree.png) is designed to cleanly organize the output of daylily.  The tree of the directories if running both the `SBS` and `B2ertD` pipelines is:
  ![](../images/assets/day_tree_dir.png)
  
The file naming convention is designed to allow straight forward integration of new aligners, deduplication, SNV/SV callers, and QC tools. Filenames are are unique w/in a batch, allowing working with files w/out namespace collision. Filenames will be unique across batches if run & experiment identifiers are manages appropriately when creating `analysis_manifese.csv` files to begin analysis batches. 
  

  
