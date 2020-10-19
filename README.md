# What's New in Version 2.1
- Built in UMI collapsing algorithm
- Parallel processing enabled on Mac, Linux, and Slurm HPC systems
- Automatic file detection

## Dependencies
- [R](https://www.r-project.org/)
- [GNU Parallel](https://www.gnu.org/software/parallel/)

## Setup
- All options from IronThrone v 1.0 still function as described below.
- In a single directory, place amplicon R1/R2 files, barcode whitelist file, config file, IronThrone-GoT script, and the appropriate Parallelized UMI Collapse script
- fastq files will be automatically detected if filenames contain `R1` or `R2`
- Barcode whitelist file will be automatically detected if filename contains `B/barcode`
- Config file will be automatically detected if file extension is `.config`

## New Option Parameters
| Option | Description |
| ------ | ----------- |
| `-tl/--target_lines` | desired file length for split fastq files (must be multiple of 4) |
| `-pcr/--pcr_read_threshold` | ratio above which majority of PCR reads must be in order for a UMI to be called definitively  (default: 0.5) |

<a href="https://github.com/landau-lab/IronThrone-GoT"><img src="https://github.com/landau-lab/IronThrone-GoT/blob/master/GoT_logo.png" border="0"></a>

# <a name="started"></a>IronThrone-GoT
Processes GoT amplicon data and generates a table of metrics

(created by the [Landau lab](https://www.landaulab.org/) at the New York Genome Center)


## Getting Ready
#### Amplicon R1.fastq & R2.fastq files
- R1.fastq (includes barcodes and UMIs sequences)
- R2.fastq (includes primer, shared and targeted sequences)
> `IronThrone-GoT` supports both uncompressed and GZipped formats.

#### Config file
To scan the expected sequences at specific positions in amplicon reads,\
.config file should be prepared as follows:

(linear_example.config)
```yaml
GCAGCAGAGAAACAAATGAAGGA 1 23   #(1st row)  primer.sequence start.pos end.pos
AAACAGGACGAGGAGCAGAGG 25 45    #(2nd row)  shared.sequence start.pos end.pos
AAGAAGACAAG 59 69              #(3rd row)      WT.sequence start.pos end.pos
TGAGGACAAAG 59 69              #(4th row)     MUT.sequence start.pos end.pos
```

(circ_example.config)
```yaml
GAAGAAGACAAGAAACGCAAAGAGG 1 25        #(1st row)          primer.sequence start.pos end.pos
AGGAGGAGGCAGAGGACAA 26 44             #(2nd row)          shared.sequence start.pos end.pos
GGAGGATGATG 45 55                     #(3rd row)              WT.sequence start.pos end.pos
TTGTCGGAGGA 45 55                     #(4th row)             MUT.sequence start.pos end.pos
TGAGGATGAGGAGGATGAGG 66 85            #(5th row)   PCR#2Fw_in_WT.sequence start.pos end.pos
TGAGGATGAGGAGGATGAGG 71 90            #(6th row)  PCR#2Fw_in_MUT.sequence start.pos end.pos
CACTGAGAATGTAAGAACTACAAACAA 96 122    #(7th row)   PCR#2Rv_in_WT.sequence start.pos end.pos
CACTGAGAATGTAAGAACTACAAACAA 101 127   #(8th row)  PCR#2Fv_in_MUT.sequence start.pos end.pos
```
> A config file should be prepared in tab-separated format.
> As an example, see `testdata/linear_example.config` or `testdata/circ_example.config`.


#### Whitelist barcode file
This file provides a list of barcodes that are put in preparation for both GoT and 10X libraries, so that false barcodes are distinguished.

Whitelisted barcode files provided by 10X are collected at `barcodes10X/` directory.

Choose the whitelist barcode file considering which 10X chemistry version used:

- 10X v2.chemistry: `737K-august-2016.txt` (default)
- 10X v3.chemistry: `3M-february-2018.txt`
> Download the `3M-february-2018.txt` file [here](https://drive.google.com/open?id=1-kMeT_asRhYu9dlCq6CkN49wDMrgcLXn), and save at `barcodes10X/` directory.

Also, the length of UMI depends on 10X chemistry version so set `--umilen` if needed.
- 10X v2.chemistry: `10` (default)
- 10X v3.chemistry: `12`

## Option parameters
#### Usage
```bash
IronThrone-GoT [options] --run linear --fastqR1 <in.R1.fastq> --fastqR2 <in.R2.fastq>
                --config <in.config> --sample <out.prefix> --outdir <out.path>
```

#### Required inputs
| Option | Description |
| ------ | ----------- |
| `-r/--run` | run module for processing `linear` or `circ` GoT (default: `linear`) |
| `-f1/--fastqR1` | input R1.FASTQ FILE (input file can be in GZip format with .gz extension) |
| `-f2/--fastqR2` | input R2.FASTQ FILE (input file can be in GZip format with .gz extension) |
| `-c/--config` | input CONFIG FILE   (input file should be in tab-separated) |

#### Additional options
| Option | Description |
| ------ | ----------- |
| `-m/--mmtch` | allowed mismatch ratio to grep the expected sequences (default: 0.2) |
| `-p/--postP` | cutoff for the posterior probability in the barcode replacement (default: 0.99) |
| `-d/--dupcut` | cutoff for the total number of duplication (default: 1) |
| `-b/--bclen` | length of barcode (default: 16) |
| `-u/--umilen` | length of UMI (default: 10) |
| `-w/--whitelist` | file for whitelisted barcodes (737K-august-2016.txt) |
| `-s/--sample` | prefix for outputs (myGoT) |
| `-o/--outdir` | path for outputs (./out) |
| `-l/--log` | logfile name (myGoT.log) |
| `-t/--thread` | number of threads to run in parallel (default: 4) |
| `-k/--keepouts` | if set to 1, keeping intermediate files (default: 0) |
| `-v/--verbose` | if set to 1, returning more logs (default: 0) |
| `-h/--help` | show option parameters |


## Output table
(prefix).summTable.txt is a table of aggregated (semi-colon separated) and averaged metrics sorted by unique cell barcode and UMI.

| Column | Description |
| ------ | ----------- |
| BC | barcode |
| whitelist | Y:barcode was perfectly matched with the whitelisted barcodes<br/>N:barcode was not perfectly matched with the whitelisted barcodes,<br/>&nbsp;&nbsp;&nbsp;&nbsp;and replaced with statistical significance |
| UMI | UMIs sharing the same barcode |
| num.WT.in.dups | counts of WT calls for each duplicate read (barcode + each UMI) |
| num.MUT.in.dups | counts of MUT calls for each duplicate read (barcode + each UMI) |
| num.amb.in.dups | counts of amb calls for each duplicate read (barcode + each UMI) |
| call.in.dups | calls for each duplicate read (barcode + each UMI) |
| avg.base_error.R2 | average of base errors for entire R2.fastq reads sharing the same barcode |
| avg.base_error.primer | average of base errors for the primer region in R2.fastq reads sharing the same barcode |
| avg.base_error.shared | average of base errors for the shared region in R2.fastq reads sharing the same barcode |
| avg.base_error.WT | average of base errors for the WT region in R2.fastq reads sharing the same barcode |
| avg.base_error.MUT | average of base errors for the MUT region in R2.fastq reads sharing the same barcode |
| mismatch.primer | counts of mismatches in the primer region for each duplicate read (barcode + each UMI) |
| mismatch.shared | counts of mismatches in the shared region for each duplicate read (barcode + each UMI) |
| mismatch.WT | counts of mismatches in the WT region for each duplicate read (barcode + each UMI) |
| mismatch.MUT | counts of mismatches in the MUT region for each duplicate read (barcode + each UMI) |
| WT.calls | counts of total WT call |
| MUT.calls | counts of total MUT call |
| amb.calls | counts of total amb call |


## (Run an example)
Test with example .fastq and .config files at `testdata/` directory.

(linear-GoT)
```bash
IronThrone-GoT --run linear --config linear_example.config \
               --fastqR1 linear_example.R1.fastq --fastqR2 linear_example.R2.fastq \
               --outdir game --log of --sample throne
```
(circ-GoT)
```bash
IronThrone-GoT --run circ --config circ_example.config \
               --fastqR1 circ_example.R1.fastq --fastqR2 circ_example.R2.fastq \
               --outdir game --log of --sample throne
```
