# Running The Bundled Tests
Daylily comes with a test dataset of 0.01x coverage data.  This will automatically be used when running a pipeline if not supercede with an `config/analysis_manifest.csv` file.  To run, simply execute:
```bash
# from the cloned root directory

source dyinit source dyinit  --project <PROJECT>
dy-a local
dy-r produce_deduplicated_bams

```
