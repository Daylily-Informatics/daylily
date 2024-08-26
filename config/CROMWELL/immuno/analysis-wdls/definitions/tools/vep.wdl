version 1.0

import "../types.wdl"  # !UnusedImport

# Do not use this task directly, use the vep workflow below.
task vepTask {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    Array[String] plugins
    Boolean coding_only = false
    Array[String] custom_args
    # Required files are necessary to force localization. The call itself uses them
    # via the custom_args field, which is a string and won't localize its parts, but
    # does need to be pointed to the right inputs dir after localization
    Array[File] required_files # !UnusedDeclaration
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
    File? synonyms_file
  }

  Float cache_size = 2*size(cache_dir_zip, "GB")  # doubled to unzip
  Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 10 + round(reference_size + vcf_size + cache_size + size(synonyms_file, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "64GB"
    bootDiskSizeGb: 30
    cpu: 4
    docker: "mgibio/vep_helper-cwl:vep_105.0_v1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String annotated_path = basename(basename(vcf, ".gz"), ".vcf") + "_annotated.vcf"
  String cache_dir = basename(cache_dir_zip, ".zip")

  command <<<
    mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
    #custom vep inputs (required files) get localized and we have to define this variable
    #pointing to their current path so that the custom string works as expected
    ~{if length(required_files) > 0 then "custom_inputs_dir=$(dirname ~{required_files[0]})" else ""}

    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf \
    --vcf \
    --fork 4 \
    --term SO \
    --transcript_version \
    --offline \
    --cache \
    --symbol \
    -o ~{annotated_path} \
    -i ~{vcf} \
    ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
    ~{if coding_only then "--coding_only" else ""} \
    --~{pick} \
    --dir ~{cache_dir} \
    --fasta ~{reference} \
    ~{sep=" " custom_args} \
    ~{sep=" " prefix("--plugin ", plugins)}  \
    ~{if everything then "--everything" else ""} \
    --assembly ~{ensembl_assembly} \
    --cache_version ~{ensembl_version} \
    --species ~{ensembl_species}
  >>>

  output {
    File annotated_vcf = annotated_path
    File vep_summary = annotated_path + "_summary.html"
  }
}


# I do not like that this has to be a VM but WDL spec 1.0 has forced my hand
# In the future if WDL gets better support for Array[struct], make this call
# happen in WDL instead of a task
task parseVepCustomAnnotationIntoArg {
  input { VepCustomAnnotation obj }
  runtime {
    preemptible: 1
    maxRetries: 2 
    docker: "python:3.10" 
  }
  command <<<
    python <<CODE
    check_existing = "~{true="--check_existing" false="" obj.annotation.check_existing}"
    custom = ",".join([
        "$" + "custom_inputs_dir/~{basename(obj.annotation.file)}",
        "~{obj.annotation.name}",
        "~{obj.annotation.data_format}",
        "~{obj.method}",
        "~{true=1 false=0 obj.force_report_coordinates}",
        "~{sep="," obj.annotation.vcf_fields}"
    ])
    print(f"{check_existing} --custom {custom}")
    CODE
  >>>
  output { String custom_arg = read_lines(stdout())[0] }
}


workflow vep {
  input {
    File vcf
    File cache_dir_zip
    File reference
    File reference_fai
    File reference_dict
    Array[String] plugins
    String ensembl_assembly
    String ensembl_version
    String ensembl_species
    File? synonyms_file
    Array[VepCustomAnnotation] custom_annotations = []  # !UnverifiedStruct
    Boolean coding_only = false
    Boolean everything = true
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String pick = "flag_pick"
  }

  scatter(vca in custom_annotations) {
    # Why is it so hard to just do `Array[File]? append File`
    Array[File] required_files = flatten([select_first([vca.annotation.secondary_files, []]), [vca.annotation.file]])

    call parseVepCustomAnnotationIntoArg { input: obj=vca }
  }

  call vepTask {
    input:
    vcf=vcf,
    cache_dir_zip=cache_dir_zip,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    plugins=plugins,
    ensembl_assembly=ensembl_assembly,
    ensembl_version=ensembl_version,
    ensembl_species=ensembl_species,
    synonyms_file=synonyms_file,
    required_files=flatten(required_files),
    custom_args=parseVepCustomAnnotationIntoArg.custom_arg,
    coding_only=coding_only,
    everything=everything,
    pick=pick
  }

  output {
    File annotated_vcf = vepTask.annotated_vcf
    File vep_summary = vepTask.vep_summary
  }
}
