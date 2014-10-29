---
title: MethyC-Seq Analysis Notes
layout: post
---

## Preliminary steps

### Bcl-conversion

1. Enter a screen

```bash
screen -S bcl-conversion
```
2. Navigate to directory with run (eg. 130909_SNL119_0105_AC2GYKACXX)
3. Check sample sheet configured correctly. If only one adapter in lane, remove adapter sequence from sample sheet.
4. Run the following (modified with correct run name, sample sheet etc). Can change final value to change number of reads in files:

```bash
/usr/local/packages/CASAVA_v1.8.2/bcl2fastq/build/bin/configureBclToFastq.pl --input-dir /dd_rundata/hiseq/Runs/130909_SNL119_0105_AC2GYKACXX/Data/Intensities/BaseCalls/ --sample-sheet /dd_rundata/hiseq/Runs/130909_SNL119_0105_AC2GYKACXX/SampleSheet.csv --fastq-cluster-count 50000000
```

5. Navigate to newly created Unaligned directory (under top run directory) and enter:

```bash
nohup make -j 12
```

### Moving and renaming files

1. Copy run files from run directory to working directory:

```bash
cp Project_E_grandis ~/working_data
```

2. Rename fastq files to s_1_sequence.txt, s_2_sequence.txt etc.
3. Store sequence files a separate directory, eg. sequences. If you have data from the same library but multiple runs, store in separate directories.

## Mapping

Mapping and postmap needs to be done on the high memory group:

```
echo "alias limit='cgexec -g *:sys_limits/high'” >> /home/<username>/.bashrc
```

Can do multiple samples at a time

Use `map.php` (0 mismatches, or use `map_1mm.php` or `map_2mm.php` for 1 or 2 mismatches) to map all the reads to the genome (follow instructions):

```bash
limit php /home/lister/working_data/php/methpipe_se/map.php | tee -a log.txt
```

For PE data

```bash
limit php /home/lister/working_data/php/methpipe_pe/map.php | tee -a log.txt
```

## Post-map

Do one sample at a time. This step will generate all the tables in mySQL (used for AnnoJ and DMR script).

This script can take mapped read runs, merge sets, convert to .slam format, sort reads,  collapse, trim, split, stack, hammer, import reads, stacks and mC to MYSQL.

Start with a mapped dir containing the subdir that contain the `\*_final` mapped files.

Navigate to directory above mapped run data and start a screen:

```bash
screen -S postmap
```

Start the postmap script as follows:

```bash
limit php /home/lister/working_data/php/methpipe_se/post_map.php | tee -a log.txt
```

You will see the following prompts:

```
  -Do you want to perform stage 1 (merge mapped runs, convert to .slam, sort, collapse, trim reads) (y/n): y
  -Do you want to perform stage 2 (import mapped reads into MYSQL) (y/n): y
  -Do you want to perform stage 3 (stack and hammer) (y/n): y
  -Do you want to perform stage 4 (import stacks into MYSQL) (y/n): y
  -Do you want to perform stage 5 (import mC's into MYSQL) (y/n): y
  -Do you want to perform stage 6 (correct mammalian mCH for genotype) (y/n): n
  -Do you want to perform stage 7 (make and import allC tables) (y/n): y
  -Do you want to perform stage 8 (identify partially methylated domains) (y/n): n
```

Enter the path to mapped folder when prompted.

Shows a summary of all options, then:

```
- Based on these settings, do you want to proceed (y/n): y
- Number of libraries that make up the sample: 1
- Enter run folder names in library 1 (space delim): sequences  <-- name of sample folder
```

# Methylpy DMR finder

Add the following to your `.bashrc`:

```bash
alias methylenv='source /usr/local/virtualenv/methylenv/bin/activate; export PYTHONPATH=/usr/local/packages/methylpy:/usr/local/packages/methylpy/methylpy'
```

All Methylpy steps must be done in using the methylenv. To exit the methylenv, type `deactivate`.

To clone methylpy from repository:

```bash
git clone https://ryanlister@bitbucket.org/schultzmattd/methylpy.git
```

To update, from the directory created by the clone:

```bash
git pull
```

To test methylpy:

```bash
python methylpy_test.py
```

There are 3 steps to the DMR finding algorithm:

1. Perform a root mean square test (you can think of it like a chisquare test) on each site across all samples. P-values are simulated (i.e., randomize the data a bunch of times and see if you get a significant result), which adjusts for multiple testing.

2. Calculate threshold p-value for a desired FDR.

3. Aggregate any significant sites within X bp and showing changes in the same direction (e.g, sample A is methylated and sample B is unmethylated) into a window.


## Generation of allC files

Edit the allC generating script: `create_allc_file_template_hs.py`

Sample base names:

```python
["sample_1_name", "sample_2_name"]
```

Database names:

```python
["database_name","database_name"]
```

MySQL server:

```python
["localhost","localhost"]
```

### Run allC generating script
Move to folder called “allC”:

```bash
$ methylenv
(methylenv) $ python create_allc_file_egrandis.py
```

## DMRfind

Run all samples at the same time within an experiment.

Edit the python script named `DMR_find.py` with your sample names and parameters.

Run script, in a folder named “DMR”:

```bash
python dmr_find.py > dmr_find_ouput.txt
```

## Histogram correction

Modify `histogram_correction.py` script with name of `_rms_results.tsv` file from allC step.

Run script:

```bash
python histogram_correction.py >> histogram_correction_output.txt
```

Use the p-value determined by histogram correction for the collapse step.

## Collapse

Edit the collapse.py script with your sample names and parameters. This may need to be changed and run several times to find the right parameters.
Run the script on the `_rms_results.tsv` file.

```bash
python collapse.py
```
 
