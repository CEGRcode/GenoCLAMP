# GenoCLAMP: consolidation of sequence motifs from genomic-scale data

### Justin S Cha<sup>1</sup>, B Franklin Pugh<sup>1*</sup>, William KM Lai<sup>1,2*</sup>

<sup>1</sup>Department of Molecular Biology and Genetics, Cornell University, USA
<br>
<sup>2</sup>Department of Computational Biology, Cornell University, USA
<br>
<sup>*</sup>Co-corresponding authors

### Correspondence: wkl29@cornell.edu

## Abstract
Genome-wide protein–DNA mapping studies have produced thousands of enriched sequence motifs, many of which are redundant or highly similar because transcription factors (TFs) often bind cooperatively as distinct complexes to shared sites. Consolidating these motifs into representative profiles is essential for reducing redundancy and improving the characterization of TF–DNA interactions. Here, we present Genomic CLustered Alignment of Motif Profiles (GenoCLAMP), a new algorithm that unifies motif-clustering and alignment to collapse redundant motifs into high-quality merged matrix representations. Using motifs generated from a large-scale Saccharomyces cerevisiae ChIP-exo project, we demonstrate that GenoCLAMP reliably recovers known TF complexes that bind common motifs and uncovers evidence for previously uncharacterized motifs. Moreover, GenoCLAMP-consolidated motifs exhibit distinct positional and combinatorial patterns across genomic features, revealing a broader DNA-sequence “grammar” that structures TF binding across the genome. Together, these results establish GenoCLAMP as a robust and efficient framework for motif consolidation, enhancing both computational analyses and biological interpretation of large-scale TF–DNA binding datasets.

## General software notes
GenoCLAMP includes a Python version and a Julia version of the software. Python is usually already installed on most operating systems, so the Python version may be easier to set up. The Julia version may be better for your needs if speed is important. Singularity/Apptainer definition files are provided for both versions.

## Software details (Python)
### Dependencies
 - Python >= 3.6
 - numpy
 - scipy
 - numba
 - openpyxl (for the output summary)
 - cairosvg (optional, for the output summary)

### Getting started
To download GenoCLAMP:
```
git clone https://github.com/CEGRcode/GenoCLAMP.git
```
Support via PyPI and conda is coming soon.

### Running GenoCLAMP
Example command on simulated motifs:
```
python genoclamp-python/run_genoclamp.py --meme simulated_motifs/SIM-1.meme --output-dest simulated_motifs/SIM-1_results
```
GenoCLAMP takes in a list of motifs using the MEME file format. This list can be provided using `--meme`:
```
python genoclamp-python/run_genoclamp.py --meme motif1.meme motif2.meme ... {other_args}
```
They can also be provided as a newline-delimited text file using `--meme-list`:
```
python genoclamp-python/run_genoclamp.py --meme-list meme_files.txt {other_args}
```
The default output location is `./genoclamp_out`. This can be changed with `--output-dest` or `-o`:
```
python genoclamp-python/run_genoclamp.py --output-dest new_genoclamp_out {other_args}
```

### Flags
 - `--meme`: Input MEME files
 - `--meme-list`: Plain text file with one MEME file path per line
   - One of `--meme` or `--meme-list` is required
 - `--nsites-thresh`: Motifs with nsites less than the threshold will be filtered out
 - `--evalue-thresh`: Motifs with E-value greater than the threshold will be filtered out
 - `--info-score-thresh`: Motifs with information score less than the threshold will be filtered out
 - `--periodicity1-thresh`: Motifs with periodicity score for period 1 greater than the threshold will be filtered out
 - `--periodicity2-thresh`: Motifs with periodicity score for period 2 greater than the threshold will be filtered out
 - `--periodicity3-thresh`: Motifs with periodicity score for period 3 greater than the threshold will be filtered out
 - `--pc`: Alpha parameters for the Dirichlet prior (one value per nucleotide)
 - `--min-base-overlap`: Minimum number of overlapping bases required for merging clusters
 - `--min-information-overlap`: Minimum bit overlap dot product required for merging clusters
 - `--max-information-overhang`: Maximum sum of absolute bit difference allowed for merging clusters
 - `--concentration`: Concentration parameter (clustering score = BLLR * cluster_size ^ concentration)
 - `--n-workers`: Number of worker threads (default is number of CPUs)
 - `--trim-thresh`: Bases on the periphery of the consensus motif with information below the threshold will be trimmed
 - `--get-sites`: When specified, outputs BED files of the union of binding sites for each cluster (input MEME files must include site information)
 - `--output-dest`: Output directory

### Output
GenoCLAMP will create folders inside the folder specified by `--output-dest` with names formatted as `clusterX`. Each folder will contain the following files:

 - Aligned Position Frequency Matrices (PFMs) in transfac format named `clusterX_aligned-motifs.transfac`
 - Consensus PFM in TRANSFAC format named `clusterX_consensus-motif.transfac`
 - Consensus PFM in MEME format named `clusterX_consensus-motif.meme`
 - Aligned stack of logos as an SVG named `clusterX_aligned-motifs.svg`
 - Consensus logo as an SVG named `clusterX_consensus-motif.svg`
 - Optional BED file of binding sites named `clusterX_binding-sites.bed`

Example output for one of the clusters in the simulated data is in `simulated_motifs/cluster404`. A summary spreadsheet of all clusters, `summary.xlsx`, is also written to the output root directory (Python version only).

## Software details (Julia)
The Julia version of the code is about 30-50% faster and also has less overhead than the Python version.
### Dependencies
- Julia >= 1.0
- [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
- [XML.jl](https://github.com/JuliaComputing/XML.jl)
- [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl)

### Running GenoCLAMP
Example command on simulated motifs:
```
julia genoclamp-julia/run_genoclamp.jl --meme simulated_motifs/SIM-1.meme --output-dest simulated_motifs/SIM-1_results
```
Arguments and output are the same as for the Python version, with one exception: the Julia version does not support `--n-workers`. To control parallelism, pass `--threads` directly to Julia: `julia --threads <n_workers> genoclamp-julia/run_genoclamp.jl ...`