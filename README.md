# Bulk Verkko assembly pipeline

_THIS REPO IS NOT PUBLICLY SUPPORTED_

Collection of scripts to download/process/assemble large number of samples via verkko, optimized to run on Biowulf / Helix.

See [hprc.txt](hprc.txt) for info on the actual commands available and example usage.

## Dependencies

* Perl 5+
* Python
* CutAdapt
* SeqRequester
* Hifiasm
* Verkko
* MashMap
* Yak
* Compleasm
* Quaak
* Meryl
* Seqtk

## CHM13v2

For chromosome assignment using MashMap alignments, download the following files:
* [chm13v2.0.fa](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa)
* [chm13v2.0.hpc.fa](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.hpc.fa)

## Configuration

Multiple instances can be run using the same data by setting up different versions.

* `config.ini` (example: [here](config.ini)).
```bash
[1]	# configuration version
samples=['$root/cache/input.tsv']  # sample metadata including path
rasm='$root/assemblies'            # assembly output dir
rsoft='$root/software-v1'          # versioned softwares, linked under this path. Not necessarily tied with config version
data='$root/seq'                   # sequences will be soft-linked under this path or written here for pre-processing
busco='$root/hprc-cache/busco'     # path to BUSCO db file for compleasm
refn="/path/to/chm13v2.0.fa"       # chm13v2 for chromosome labeling
refc="/path/to/chm13v2.0.hpc.fa"   # homopolymer compressed chm13v2 for chromosome labeling, for the homopolymer compressed assembly graph annotation
odb="primates_odb10"               # compleasm BUSCO db
hybridCorrection="false"           # set for ont-only correction
```

## Fetch data

* Prepare `input.tsv` with local or s3 paths
* `local://` if data is local
* Reads can be cram / bam / fq / fa formats, bgzipped (preferred for fq / fa files) or gzipped
  ```sh
  # Example to fetch read paths
  ls $path.fq.gz | awk '{print "local:/"$1}' | tr '\n' ','
  ```
* See https://github.com/skoren/R3_HPRC for updated examples and use [convertToHPRCFormat.py](https://github.com/NHGRI/marbl_utils/blob/master/sequence_tools/convertToHPRCFormat.py) for a script to generate a compatible tsv from R3 inputs.

* input.tsv format (example: [here](hprc-cache/b1.tsv))

  | sample_id | isMaleSample | hifi  | ont | hic | child_ilmn | mat_id | mat_ilmn | mat_id | mat_ilmn  |
  | --------- | ------------ | ----- | --- | --- | ---------- | ------ | -------- | ------ | --------- |
  | SAMPLE_A      | true\|false  | ['path1','path2',...] |  ['path1','path2',...] |  ['path1','path2',...] | ['path1','path2',...] | SAMPLE_A_MAT | ...
  
  All fields are optional

* Fetch
  ```sh
  module load perl/5.36
  perl hprc.pl list --v1
  perl hprc.pl fetch --v1
  ```

## Read filtering

* HiFi cutadapt
  ```sh
  perl hprc.pl read-filter --v2
  # Confirm "READY-TO-COMPUTE" before submission
  perl hprc.pl read-filter --v2 --submit
  ```

## Read stats

* This is the only step done on Helix for heavy IO
* Takes usually long unless manually parallized
* Note that `read-stats` does not have a `--submit` option. It runs directly

  ```sh
  # On tmux, or make a read-stats.sh to run with nohup
  perl hprc.pl read-stats --v1
  # For nohup:
  nohup ./read-stats.sh ?>> read-stats.log &
  ```

* Confirm readlen:ok
  ```sh
  perl hprc.pl read-stats --v1
  ```

* Yakmers (optional, in case we have Ilmn reads)
  ```sh
  perl hprc.pl yakmers --v1
  perl hprc.pl yakmers --v1 --submit
  ```

## Read correction

* Hifiasm read-correction
  ```sh
  perl hprc.pl read-correct --v1
  # Confirm the .sh scripts under data

  # This is to submit all samples under input.tsv
  perl hprc.pl read-correct --v1 --submit

  # To submit specific samples
  perl hprc.pl read-correct --v1 --submit --sample A B C ...
  ```

## Verkko assembly

* Script auto-detects what's available, what's not and will submit in 4 modes unless specified

  ```sh
  perl hprc.pl assemble --verkko-full --v3
  ```

* Here are Verkko modes that will be submitted with --verkko-full

  The script auto-detects what's available (READY-TO-COMPUTE) and shows the status. Scripts will be generated to submit on biowulf.

  Modes can be submitted separately (e.g. `--verkko-base`) if needed.

  | Hi-C | Parental-Ilmn | Verkko Mode |
  | :---: | :------------: | ----------- |
  |   O   |        O       | verkko-base, verkko-hi-c, verkko-trio,  verkko-thic |
  |   X   |        O       | verkko-base, verkko-trio |
  |   O   |        X       | verkko-base, verkko-hi-c |

## Analyze

* Analyze submits chr-assign and analysis for each verkko-modes available
  ```sh
  perl hprc.pl analyze --v1
  # Check what's ready
  perl hprc.pl analyze --v1 --submit
  ```

* KNOWN ISSUE: Chromosome reorienting and renaming is currently using a custom script. This will be updated soon.


