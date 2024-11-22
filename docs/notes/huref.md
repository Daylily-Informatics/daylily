# Overview
my notes re: human references


# Builds

## (GRCh38/hg38)
_from chatgpt 4.o_
### Chromosome Sizes (GRCh38/hg38)

| Chromosome | Size (bp)       | Notes on Mapping Difficulty                        |
|------------|-----------------|---------------------------------------------------|
| 1          | 248,956,422     | Largest chromosome, may have challenges in repetitive regions. |
| 2          | 242,193,529     | Large, with some complex regions.                 |
| 3          | 198,295,559     | Relatively easier compared to larger chromosomes. |
| 4          | 190,214,555     | Contains complex centromeric and telomeric regions. |
| 5          | 181,538,259     | Moderate complexity for mapping.                  |
| 6          | 170,805,979     | Some regions are MHC-heavy, leading to high variability. |
| 7          | 159,345,973     | Relatively easier due to moderate size.           |
| 8          | 145,138,636     | Moderate size and complexity.                     |
| 9          | 138,394,717     | Contains a large heterochromatic block, challenging to map. |
| 10         | 133,797,422     | Smaller, relatively straightforward to map.       |
| 11         | 135,086,622     | Some repetitive and GC-rich regions.             |
| 12         | 133,275,309     | Moderate size, often mapped well.                 |
| 13         | 114,364,328     | Smaller, but contains repetitive regions.         |
| 14         | 107,043,718     | Smaller, relatively easier to map.                |
| 15         | 101,991,189     | Moderate difficulty due to repetitive regions.    |
| 16         | 90,338,345      | Smaller and gene-rich, often easier to map.       |
| 17         | 83,257,441      | Smaller, with gene-dense and repetitive regions.  |
| 18         | 80,373,285      | Smaller, moderate mapping difficulty.             |
| 19         | 58,617,616      | Gene-rich and GC-rich, potentially challenging.   |
| 20         | 64,444,167      | Smaller and generally easier to map.              |
| 21         | 46,709,983      | Smallest autosome, generally easier to map.       |
| 22         | 50,818,468      | Small, with some repetitive regions.              |
| X          | 156,040,895     | Large, may have challenges in pseudoautosomal regions. |
| Y          | 57,227,415      | Smaller, male-specific, many repetitive regions.  |
| MT         | 16,569          | Mitochondrial DNA, very small, mapped easily but with care for contamination. |

#### Easiest Chromosomes to Map with BWA
- MT	16,569	Mitochondrial DNA, very small, mapped easily but with care for contamination.
- Small Autosomes: Chromosomes 21, 22, 20, and 19 are relatively small and easier to map due to fewer repetitive and complex regions (with exceptions for gene-rich areas like chromosome 19).
- X Chromosome: Easier for female samples due to lack of Y contamination; pseudoautosomal regions on the X-Y boundary may pose challenges for males.
- MT (Mitochondrial DNA): Extremely small and typically straightforward to map, provided high-quality sequencing and contamination checks are performed.

#### Considerations for Mapping with BWA
##### Repetitive Regions:

- Large chromosomes (e.g., 1, 2) and regions near centromeres or telomeres are more repetitive, making mapping difficult.
- Use the bwa mem algorithm with appropriate parameters for handling repeats (e.g., -T for mismatch penalties).
- GC-Rich and Gene-Dense Regions:

- Chromosomes 19 and 22 are gene-rich and GC-rich, which can lead to mapping biases.
- Use sequencing libraries with uniform coverage across GC content to minimize issues.
- Highly Variable Regions:

- Chromosome 6 (MHC region) and chromosome 9 (heterochromatic block) are challenging to map due to high variability.
- Mapping Parameters for Easier Mapping:

- Use higher mismatch penalties (-B) and lower gap extension penalties (-O) to improve specificity.
- Consider -t for multithreading to utilize more CPUs efficiently.
- Let me know if you need further assistance with specific mapping strategies or parameter tuning for BWA!


# For use in chunking variant callers
```
chunk_chrms: 1~1-50000000,1~50000001-100000000,1~100000001-150000000,1~150000001-200000000,1~200000001-248956422,2~1-50000000,2~50000001-100000000,2~100000001-150000000,2~150000001-200000000,2~200000001-242193529,3~1-50000000,3~50000001-100000000,3~100000001-150000000,3~150000001-198295559,4~1-50000000,4~50000001-100000000,4~100000001-150000000,4~150000001-190214555,5~1-50000000,5~50000001-100000000,5~100000001-150000000,5~150000001-181538259,6~1-50000000,6~50000001-100000000,6~100000001-150000000,6~150000001-170805979,7~1-50000000,7~50000001-100000000,7~100000001-150000000,7~150000001-159345973,8~1-50000000,8~50000001-100000000,8~100000001-145138636,9~1-50000000,9~50000001-100000000,9~100000001-138394717,10~1-50000000,10~50000001-100000000,10~100000001-133797422,11~1-50000000,11~50000001-100000000,11~100000001-135086622,12~1-50000000,12~50000001-100000000,12~100000001-133275309,13~1-50000000,13~50000001-100000000,13~100000001-114364328,14~1-50000000,14~50000001-100000000,14~100000001-107043718,15~1-50000000,15~50000001-100000000,15~100000001-101991189,16~1-50000000,16~50000001-90338345,17~1-50000000,17~50000001-83257441,18~1-50000000,18~50000001-80373285,19~1-50000000,19~50000001-58617616,20~1-50000000,20~50000001-64444167,21~1-46709983,22~1-50000000,22~50000001-50818468,23~1-50000000,23~50000001-100000000,23~100000001-150000000,23~150000001-156040895,24~1-50000000,24~50000001-57227415,25~1-16569
```



## b37
### Chromosome Sizes (b37)


| Chromosome | Size (bp)       | Notes on Mapping Difficulty                        |
|------------|-----------------|---------------------------------------------------|
| 1          | 249,250,621     | Largest chromosome, may have challenges in repetitive regions. |
| 2          | 243,199,373     | Large, with some complex regions.                 |
| 3          | 198,022,430     | Relatively easier compared to larger chromosomes. |
| 4          | 191,154,276     | Contains complex centromeric and telomeric regions. |
| 5          | 180,915,260     | Moderate complexity for mapping.                  |
| 6          | 171,115,067     | Some regions are MHC-heavy, leading to high variability. |
| 7          | 159,138,663     | Relatively easier due to moderate size.           |
| 8          | 146,364,022     | Moderate size and complexity.                     |
| 9          | 141,213,431     | Contains a large heterochromatic block, challenging to map. |
| 10         | 135,534,747     | Smaller, relatively straightforward to map.       |
| 11         | 135,006,516     | Some repetitive and GC-rich regions.             |
| 12         | 133,851,895     | Moderate size, often mapped well.                 |
| 13         | 115,169,878     | Smaller, but contains repetitive regions.         |
| 14         | 107,349,540     | Smaller, relatively easier to map.                |
| 15         | 102,531,392     | Moderate difficulty due to repetitive regions.    |
| 16         | 90,354,753      | Smaller and gene-rich, often easier to map.       |
| 17         | 81,195,210      | Smaller, with gene-dense and repetitive regions.  |
| 18         | 78,077,248      | Smaller, moderate mapping difficulty.             |
| 19         | 59,128,983      | Gene-rich and GC-rich, potentially challenging.   |
| 20         | 63,025,520      | Smaller and generally easier to map.              |
| 21         | 48,129,895      | Smallest autosome, generally easier to map.       |
| 22         | 51,304,566      | Small, with some repetitive regions.              |
| X          | 155,270,560     | Large, may have challenges in pseudoautosomal regions. |
| Y          | 59,373,566      | Smaller, male-specific, many repetitive regions.  |
| MT         | 16,569          | Mitochondrial DNA, very small, mapped easily but with care for contamination. |


#### Easiest Chromosomes to Map with BWA
- same as hg38 basically.


# For use in chunking variant callers
```
chunk_chrms: 1~1-50000000,1~50000001-100000000,1~100000001-150000000,1~150000001-200000000,1~200000001-249250621,2~1-50000000,2~50000001-100000000,2~100000001-150000000,2~150000001-200000000,2~200000001-243199373,3~1-50000000,3~50000001-100000000,3~100000001-150000000,3~150000001-198022430,4~1-50000000,4~50000001-100000000,4~100000001-150000000,4~150000001-191154276,5~1-50000000,5~50000001-100000000,5~100000001-150000000,5~150000001-180915260,6~1-50000000,6~50000001-100000000,6~100000001-150000000,6~150000001-171115067,7~1-50000000,7~50000001-100000000,7~100000001-150000000,7~150000001-159138663,8~1-50000000,8~50000001-100000000,8~100000001-146364022,9~1-50000000,9~50000001-100000000,9~100000001-141213431,10~1-50000000,10~50000001-100000000,10~100000001-135534747,11~1-50000000,11~50000001-100000000,11~100000001-135006516,12~1-50000000,12~50000001-100000000,12~100000001-133851895,13~1-50000000,13~50000001-100000000,13~100000001-115169878,14~1-50000000,14~50000001-100000000,14~100000001-107349540,15~1-50000000,15~50000001-100000000,15~100000001-102531392,16~1-50000000,16~50000001-90354753,17~1-50000000,17~50000001-81195210,18~1-50000000,18~50000001-78077248,19~1-50000000,19~50000001-59128983,20~1-50000000,20~50000001-63025520,21~1-48129895,22~1-50000000,22~50000001-51304566,23~1-50000000,23~50000001-100000000,23~100000001-150000000,23~150000001-155270560,24~1-50000000,24~50000001-59373566,25~1-16569
```