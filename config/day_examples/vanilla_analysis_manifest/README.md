# Analysis Manifest Details

**THIS IS ALL ACTIVELY BEING DEPRECATED IN FAVOR OF MOVING TO SNAKEMAKE v8 PATTERNS**

> This info applies only to the standard analysis manifest, not the special bcl2fastq, whcih is described elsewhere.

## Analysis Manifest Details
The manifest instructs DAY as to which starting data to use. There must be at least  one data row under the header row of the file.  Eaach time the pipeline is run, the manifest is validated vs. its companion yaml schema located in workflow/shcema/ dir.  If any validateable cols fail validation, a clear error is thrown.  Note-  any number of unspecified colums may be appended to the core validated set of columns.

The manifest should be saved in the config dirwith the follwing name: 'analysis_manifest.csv'

## Format

  * 9 Required colmns which must have non-null values
  * 8 Required columns which may have null(na) values
  *  supported columns, 16 of which must be present and at leaset have an 'na' entry (no empty string fields are permittted in any colum found in the manifest.) :

sample,sample_lane,SQ,RU,EX,LANE,r1_path,r2_path,biological_sex,iddna_uid,concordance_control_path,is_positive_control,is_negative_control,sample_type,merge_single,external_sample_id,remote_bam_file_path,rclone_handle
  - _the columns should be in this order as well_

  - The last 2 columns are specifically for fetching BAM files which your user can access with a pre-configured backend handle (ie: if rclone ls someS3Acct:/abucket/and/a/path/to/a/sample.bam returns the ls info for the bam, then it can be fetched w/DAY by enterining `someS3Acct:` in the `rclone_handle` column, and in teh remote_bam_fiqle_path you would enter `/abucket/and/a/path/to/a/sample.bam`.  If these fields are populated, you will also be expected to provice on the command line the path to your crclone config file via ` --config rclone_config_file=/full/path/to/rclone/rclone.conf ` as part of the call to mod-run.

#### Additional Columns
  - may be added, but may not have spaces in their column headers and should follow the final supported column.  Null fields must be na, and values should  not have tabs or space characters.


* sample is ideally a name in the form RU#_EX#_SQ#_#, where the Run is the seq runid, the EX# is of the EXdir you're in (thouh no validation on the name happens, you could just as easily put an AT instead of an EX.  The SQ is the SQ the reads this row point to, and the last # is the lane #.... *if* the lane fastqs are not combined.  If they have been combined, this is set to 0.
* SQ is the SQuid from the sample name
* EX is the EXuid from the sample name
* RU is trhe RUuid from the sample name
  * *In our case, these are LIMS IDs, but the pipeline has no concept of lims IDS- these could just as easily be strings* Our primary use of thes IDs is to make a good attempt at making the samplke name in each report unique'ish.  Meaning we could combine several runs worth of QC data and not have sample or file name clashes.
* r1_path - /full/ path from root to r1 file
* r2_path - /full/ path from root to r2 file

## Required Columns, 'na' Values Allowed

* biological_sex=male or female or empty (an entry not mathcing make or femalke will be treated as unknown).  Mostly this information is useful as a QC check to confirm the stated sample biological sex matches what is observed.
M,male,F,Female,female
* iddna_uid -ie: TA12343 or if there are 2+ expected, connect them with : with no spaces.  ie: TA1234:TA9876
* concordance_control_path - empty or a  path to a directory with a vcf file named 'sample.vcf' and a bed file named 'sample.bed'. These will mostly be coriel controls. The directory these files live in will be used to name the subdirs the concoirdance results are stored in . SNVs only.   *NOT IMPLEMENTED YET*
* is_positive_control - empty or yes or no (no and empty treated the same. To display on reports.) *NOT USED YET*
* is_negative_control - empty or yes or no (no and empty treated the same. To display on reports.) *NOT USED YET*
* sample_type - empty, or: blood, saliva, buccal swab, etc.  May be displayed in QC reports or used to guide pass/fail decisions. *NOT USED YET*
* ... all *NOT USED YET* fields must still be present in the sample sheet

#### Optional, Informational Columns
  - As many additional columns as you wish may be present in the sample sheet as long as there is no collision with any known column names.  Names should also be limited to the simple alphanumeric character and not contain spaces. The informational colmns should follow the last of the supported column.  Filed values, if null, must be 'na', and values should  not contain tabs or space characters.

# Creating sample sheets
You may make them manually, but there are helper scripts in bin to generate them for common use cases (like from an AT bcl2fq dir, or from just generic fastq files of any provenance (which need to have expected extensions however). The core README has a detail page on the Manifest for more info.

## From a BCL2FASTQ AT Directory.
* For any run through our prod pipeline, you can generate a compliant sample sheet with the script `bin/at2mgsamp_sheet.py`

  > This scripts takes 6 positional arguments:
        1) /full/path/to/bcl2fq/at/output/
        2) RU#
        3) EX#
        4) (split is not fully tested)either "split" or "merge".  This informs mod if lane specific data should be merged or split
        5) Path to create output links w/in.
        6) manifest_file_name.csv
this will leave you a sample sheet named as entered in #6.   Copy that file to the config dir names  as such:
`cp newly_creates_manifest.csv config/analysis_manifest.csv`

* with this saved to the expected location, it will not be used when executinf any commands.

## A script to generate a manifest from Metro Park Reads
> main]$bin/mp_reads_to_manifest.py

This scripts takes 4 positional arguments:
        1) /full/path/to/fqdir/
        2) EX#
        3) Path to create output links w/in.
        4) sample_sheet_file_name.csv
* The path in #1 should contin files with the extension R1.fastq.gz or R2.fastq.gz

## Script to take arbitrary reads and build a minimal manifest.
##### Well, the script was written, but appears un-checkedin somewhere out there...
* It retuired the end of the reads follow wht convention of bcl2fq reads(so R1/R2 could be found,and pairs identified bt the string preceding R1R2, Uniqe integer were assigned as sample IDs and a mapping table was produced to key the originall sample/fike names to the new ones.)   R0 was used for the run, and whichever working experiment directory I had created with Vera was used as the EX (Runs begin associated with one EX, but EXs are not designed only for runs.  EXs are just containers to uniquely Identify and help locate analysis results,  which commonly is just for one run, but an ex may hold just a few samples from a run, data from several runs, collaborator data, etc.  So, this activity requires a EX# always. EX#s also help (but do not fully protect) avoid reanalyses on the same experiment from beig confused, as each re-analysis, is in ita own EX, can be unambiguously described.)   And, this is also not enforced.  You may choose to what extent you wish to leverage EXs (if you have a system you like which keeps everything within one EX, this is not prohibited.,)

# Example Manifest Files

## Comprehensive

  > clia_wgs_validation_giab.csv
  This manifest contains all GIAB samples from a WGS validation run.   The fields are annotated to test IDDNAs, sex, and concordance, and will be the most comprehensive test case.

## Three low coverage examples you may choose to test run with, with just the basic info.

* 0.01xwgs_HG002.samplesheet.csv
* 0.05xwgs_HG002.samplesheet.csv


### CSVs
* I notoriously loathe CSVs, but it's the standard here and I went with it (until I switch to using PEPs!).  So, no commas please in the optional columns. And please do not quote the characters betwen the commas.


### More on rclone config

Rclone is a very feature rich tool.  You can find the docs here.  In addition to connecting and working with all manner of data backends and transfer protocols, it also supports the non-glamorous 'local' storage, and in the same pattern as any remote storage. So this example should translate to any valid remote you define.

```
cd $DAY_ROOT
source dyinit  --project <PROJECT>
mg-a local
rclone -L --config config/day_example_analysis_manifests/rclone_xxx_local_example.conf cat localxxx:/fsx/data/external_data/research_experiments/ManualRuns/EX3636/b37/wgs_data/SQ42092_LS1585266_RE2698_G3/Lib_for_SequenceIndex_SQ42092.aligned.deduped.sort.bam | samtools view -H
```

  - This should return the header from the BAM being referenced.  If this works for you, and it should if your DAY env built, then the BAM manifest features should work.
  - You can initialize new rclone remote handles by running `rclone config`.  The entry for the local is simply:

  ```
  [localxxx]
  type = local

  ```
  - You can test run with the included rcloone.conf from XXX with: (it runs a few hours):
  ```
  cp .tiny_data/data/0.01xwgs_HG002.samplesheet.csv config/analysis_manifest.csv

   mg-r seqqc -p -k --config rclone_conf_file=config/day_example_analysis_manifests/rclone_xxx_local_example.conf -n

   mg-r seqqc -p -k --config rclone_conf_file=config/day_example_analysis_manifests/rclone_xxx_local_example.conf
  ```


# Optional, Functional Manifest Columns

This is a paradigm not really disucsssed in Snakemake-land (correction:  this is advised against... though I'm not convinced yet it needs to be avoided).  Snakemake argues this kind of workflow relates config should all reside in teh profile yaml files.  I have found that for this kind of routine param changing, those files are a small inconvenience to edit. And the YAML files are largely controlling sample agnostic behavior or sample irrelevant behavior. All which leave it possible but in my expereince, awkward to build w/in workflow sample spoecific behavior.

The analysis manifest feels quite natural for setting sample specific behaviour, though it is arguably expanding the intention of the manifest. So, please consider this an exploratory approach and not necessarily a strong vote in favor of this approach (yet at least).

## For The Most Part, The Following Must Be Manually Annotated At This Time
  * Note- the handful of manifest build scripts available presnetly only fill in the fields wich can not contain na.  All of these fields will be blank in manifests generated with those scripts.
  * This said, the information to populate much of the below is knowable, but reqiores a tigheter coupling to LIMS and whatnot, which was not a priority.
  * It also feels like there is an opportunity to build a experiment google sheet to a manifest, as those tables contain a lot of this info, but in a format not yet reproducible to allow programatic use if them as data sources.

#### Run, Add More Samples & Run Again, Define a New Tool, Run Again....

* Like all workflow managers, a pipeline can be re-run anytime.  If the pipe has complete successfully, and nothing new data or code wise has been added-  you'll be informed there is nothing to do.

##### however

*  If you add new samples to the manifest, it should note that fact and run those samples, leaving the already complete samples as-is.
*  If a new tool is defined,  this should be noted and if re-run, just this tool will be executed for all samples.
*  So-  you may always run with a barebones manifest, and re-run in the future
*  ^^^ Many caveats apply.... always verify a re-run appears to be doing the right thing by using the `-n ` flag.

---

required columns ([examples](../config/day_example_analysis_manifests))

---

### col namne == `concordance_control` and `external_id`
> This colmns cells may contain 'na' to skip running a concordance check, or a contain  a path to a directory which is expected to contain 1-many sub directories(of any name, but ideally brief and desriptive of the sample specific concordance data w/in, ie: may_varset_v1 ).   Each subdir is scanned for 2 files, with names `external_id.vcf.gz(and .tbi)` and `external_id.bed` where `external_id` is dran from the external ID column.  Each directory with these files will have the sample snv.vcf compared against the concordance directory `external_id.vcf`, with the usual concordance stats calculated w/in the regions defined in the `external_id.bed`.  These resporrs are saves in the sample results path per- subdirectory dataset.   For the GIAB samples, the is defined the entire dataset for each sample w/int the genome wide HC regions, then subset to AY20(for comparison with the validation) and AY84, the new exome.  NA12878 has a few additional datasets included.  The genome wide 'giab' results are compiles across all samples and is presented in the multiqc report.  But there will be many more granular reports saved w.the sample results data.

*Concordance With Existing Sample Results*
_From Any Source_
May be internal or publoc datasets~
This approach was kept minimal to be flexible. You can run concordance checks on any sample, not only the well known controls.   The sample vcf produced by mod will be compared to all matching bed+vcf pairs in the search directory (and agin the names of the subdirs they are in will be used to tag the results, so try to be concise and descriptive) on any vcf and bed file that is detected in a subdirectory of the given seach dir.
If you would like to check the concordnce of a clinical sample you are processing with results from a different pipeline- or different version of the current pipeline, you simply define a subdir  holding links (or copies) to the VCF and BED to use for the analysis.
!! It is worth noting, concordance against random VCFs is going to look much less clean.  This is b/c the truth sets are hihgly processed to eliminate regions too hard to be useful to include.   You may wish to use a more constrined bed region to caclulate within for this kind of analysis.

### col namne == `neg_control`

* Nothing implemented for this field in the pipeline yet. It is presently serving as metadata stored with he analysis results in a consistent format, which will facilitate analysis of these specific data but there being  standard and very clear place to these data, as well as being of value when analysis across runs happens more commonly.
* One use in the pipeline coule be to mark samples in multiqc of this class when the report builds.

### col namne == `pos_control`


* Nothing implemented for this field in the pipeline yet. It is presently serving as metadata stored with he analysis results in a consistent format, which will facilitate analysis of these specific data but there being  standard and very clear place to these data, as well as being of value when analysis across runs happens more commonly.
* One use in the pipeline coule be to mark samples in multiqc of this class when the report builds.


### col namne == `sample_type`

* Nothing implemented for this field in the pipeline yet. It is presently serving as metadata stored with he analysis results in a consistent format, which will facilitate analysis of these specific data but there being  standard and very clear place to these data, as well as being of value when analysis across runs happens more commonly.
* One use in the pipeline could be to mark samples in multiqc of this class when the report builds.
* Another might be to modulate pipeline params for samples of various types, even in manifests with mixed types.

### col namne == `biological_sex`
* Values may be : na,female,none,unknown,male
* This is usesul as metadata to have saved with the analysis, and is also used in a QC check ( presented in the final multiqc, which flags samples with a sex derived from alignment that does not match the sex entered in this column)

### col namne == `iddna_uid`
* The TA# >>(look up syntax to define multiple)<< which will be used to generate iddna metrics.  If null, metrics are still scanned, no expexted numbers will be generated, but you will still get a report of any iddnas identified (all classified as unexpected)

### col namne == `merge_single`

* Is being deprecated. To merge all sample lane fastqs or process them as dinsrinct pairs is inferreed from the names set in the sample and sample_lane cols.  And for the moment, singl lane processing has not been worked on in some time, and should be considered experimental.






## Optional Columns May Be Specified (and should be carefully checked as these cols are not validated as of yet.)

## Pre-Alignment Downsampling Options

  Note-  the downsampling options described below will involve a nested transformation of the raw fast which importantly for bwa, sends it a single interleaved stream of fastq. In the case bwa receives interleaved fastq's, you must set the -p flag, otherwise, alignment proceed w/not error, and you will find yourself with subly very strange results.   I mention this in case you hack anything up-- the 43`-p` flag should be set automatically for the code supporting it presently.

### col namne == `subsample_to_pct`

> Randomly (with preset seed) Align this % of the input fastq reads. Values must be between 0-1.0, which will downsample by PERCENTAGE from the raw fastq, during alignment.  This is a different method of downsampling, and will not be appropriate in all cases- but may also be superior to time saving downsampling of already processed artifacts, like from an already aligned and deduped BAM.  Again, it will depend on your desired outcome, but this is a very unbiased way to downsample.

### col name == `subsample_to_cov`

> This also applies subsampling reads during alignment, so the coverage of the read set has not been precisely determined by an aligner.  Coverge is estimates from the input fasstq reads, and that coverage estimate is used to deteremine the subsample rate.Thid approach will almost cetainly not give you on the nose post laignment coverage, as specified.  It should reflect the pre=duplicate screeening coveage more closely, but i the library has a substantial amount og non-humsan reads, theose will br considered to count towards coverage and then now align, leving the final coverage below what is reuqested.   It is feasible to plut in a few new tools which can sample a set of randome reads from the fastqs, use ultra rapid alignment to estimate the % of raw reasd whcih will align.  And the dedup % could be estimated directly from the fastqs as well- numerous tools can accomplish this.  With that added here, whe resulting sampling rate can be much better turned and yilerd a final coverage closer to the requested #.  But, for the momenht, it's a rough implementation.
#### Allowed Col values
  - `na` to no downsample, or an `integer from 0-n`, if you request a coverage downsample value that appears greater than the likley coverage, an error will the trown that this request appears unable to be satisfied.
  - This feature operates per sample, so in the case you are combining multiple fastqs to a single bam, you must specify for each lane-row to be combines the same value, and it should be  in the context of the merged set of fastqs, not per contributing fastq.
  - Also, b/c this is per sample, you may enable this for some samples and skip `na` others w/in the samne manifest.

---
---

## Raw FastQ Ingestion Methods

### col name == `fq_ingestion_approach`
> There str mulyiple ways for program to ingest their input data.  Aligners and a a handful of other programs all largely follow the same pattern, so this allows conrol over which method to use to feed fastq data to, say, bwa-mem2 (proc sub is advised).


#### Allowed Col values
*raw reads*
  * Commonly, and most approachly, the program is given files to read.  This has the benefit of being clear (if the file being read from it not messed with) - but has overhead in the opening and consumin gthe file. Limitstion str thst msny programs with flags to accept file inputs can not take multiple files. In this mode, the files are presented in order raw:
        -`R1.fq.gz R2.fqgz`?


*process_substitution*
Process substitution allows running an arbitrart command and presenting it's output as a file.  This is very powerful and can get a little gnarly, but the primary use I make of it is pretty straight forward.  ie:
  - 4 lane sequencer run completes, generating 4R1 and 4R2 fastq files.  Often the next step of a pipeline is to concatenate each all of the R1 and all R2 files into a sigle file each, which then bwa can read in and align.  Using process substitution, you can `cat, or rcat, or unpigz` the list of 4R1 inside this construct `<(cat R1L1.fq R1L2.fq R1L3.fq R1L4.fq )` and specify it as the R1 input file, and bwa will begin reading from it immediately.  There is no copy step needed.  Also, using `unpigz` is advised, it outperforms the others by a small but meaningful margin.
  - So, the bwa command, in toy form, would loook like:
        `bwa mem huref.fast <(unpigx -c -q -- R1.fast.gz )   <(unpigz -c -q -- R2.fastq.gz )`
        **There seems to even be a gain when dealing with local compressed files as unpigz is faster at reading these compressed files than most everythiing else.


## Exotic Source Data Utilization

### Local, atypical fatq's

* These are well handled. A few rules should be followed for best results, but you can fill out the required w/required value columns as follows:
1) sample = can be any string of non-controversial characters and not including any '.''s.  Try to keep the names reasonably short.    The name should have 4 '_', with the information between them up to you to decide on--- though whatever you chooose to appears between the '_' should match the column for that bit ie: RU_EX_SQ_0
2) sample_lane = should be sample_1 (of if multiple lanes are merging, sample_[1,2,3,4] as appropriate.).  The string in this column must be unique and the the string 'sample' begins with.
3,4,5) RU/EX/SQ --- these do not need to actually be these IDs, but need to match the corresponding field in the sampe and sample_lane strings.
6)Lane - 0 implies 'artificial' lane and most commonly the result of merging lanes. the *_0 lane in the 'sample' will become the name used in pipeline files.  so, when merging fastqs, the sample filed string should end in a _0, and the sample_lane string should have the lane number at the end _1.  When not combinging, the sample_# and sample_lane_# may be the same.
7,8) path to a R1 fastq and R2 fastq for a lane, or subset or combined fastq file pair.
....then the remaining columns, the documentation above all applies.
#### Experimental External Data Manifest Builder scripts
    * A few `use at own risk, and absolutely sanity check results` scripts live in bin/ to process atypical reads and produce manifests. Those may be useful out of the box.



---
---

## Not Built, Woth Investigating

*store in unaligned bam, then align bam* +++ NOT BUILT
  - An option to explore.  There are a pretty compelling list of reasons to consider choosing SAM/BAM as the primary read storage format.


---
---
---
# BCL2FASTQ Sample Sheet

  * The format of these is dictated by bcl2fastq. And there are a variety of styles in active use depending on the library prep and sequencer type.  (read: they can not be used interchngibly per se.  Confirm you have the correct format with a lab guru, and if re-processing a LIMS managed run, then you're in luck-- the top of the BCL2FQ AT directory will have a SampleSheet.csv which can be used and certain fields modified as-is if you copy iy locslly.  There is an example HowTo Run BCL2Fastq.md with specifics.
  > It is worth mentioning that running bcl2fastq, unlike most of the tools in DAY, *requires* you edit the profile rule_config.yaml in the `config/day_profiles/PROFILENAME/` dir.  <small>this is the set of core config files whch defines the profiled command arguments and whatnot.  Typically you won't need to edit anything here... and please do not edit the files in `template` as tgese are the template fils ysed to build new config files for every newe deployment of mod.
