# Common Errors

A few common errors that we've encountered and potential solutions

## Out of Space

These errors indicate that the disk storage space has been
filled. That part is pretty straightforward. The part that's a bit
more of a pitfall is that _depending on what output_ this happened on,
you need to change different disk sizes.

If the failed write happened at `/cromwell_root` path, then `disks:
"local-disk ..."` needs to be increased. However, if the failed write
happens during to `stdin` or `stdout`, or any of the other standard
Linux-y places, then you'll need to increase the value of
`bootDiskSizeGb`. Cromwell in GCP mounts two disks, at minimum: the
boot disk, and a local-disk. Boot disk handles all the operating
system files, but local-disk is where almost all of your "work" is
going to happen, besides piping between commands.

## File missing

This applies more to newly converted files than hardened ones _but_
many runs failed because a file wasn't included in the
instance. Generally, this happens because the CWL did not specify a
secondaryFile that it assumed would exist next to the passed in
file. This works on the cluster, because the tools just look for the
file and it already sits where it's expected. This does not work on
the cloud because that file is never sent to the instance. The
solution is to add this parameter explicitly to the WDL and pass it
through, top down.

## CommandException: No URLs matched

This is one of two things. Either (A) the input is malformed or
otherwise incorrect, or (B) the specified file was not uploaded to the
bucket. These are both instances of the general version of the error,
"No file has been uploaded to the specified URL".

# Differences from CWL

Last confirmed mirror with the analysis-workflows CWL repo was commit
788bdc99c1d5b6ee7c431c3c011eb30d385c1370, PR#1063, Apr6 2022. Commits
from that point on may deviate unless compared. Update these values if
that is done.

## Directory types must be a zip file, or Array[File]

There is not yet a supported Directory type in WDL. Instances of this
like `Directory vep_cache_dir` which involve nested directory structure are
replaced with `File vep_cache_dir_zip`. Instances of this like
`Directory hla_call_files` which are just a flat collection of files are
replaced with `Array[File] hla_call_files`.


## Input files must prefix arguments with the name of the workflow

Input files must prefix each argument with the name of the workflow
they're going to run, because a WDL file can contain multiple
workflows or pass inputs over a layer if they aren't propagated
through in the definition. e.g. to call workflow `somaticExome` with
input `foo`, yaml key must be `somaticExome.foo`

If WDLs are being used leveraging the
[`cloud-workflows/scripts/cloudize-workflow.py` helper
script](https://github.com/griffithlab/cloud-workflows/tree/main/scripts),
the generated input file will have this handled already.
