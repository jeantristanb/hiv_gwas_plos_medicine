# Leave-One-Study-Out Meta-analysis

A leave-one-study-out sensitivity analysis was performed from the original Africa METAL meta-analysis.

The procedure automatically creates a separate METAL analysis for **each `PROCESS` file** in the original meta-analysis. In each analysis, one summary-statistics file is removed while all other files and meta-analysis settings remain unchanged.

This provides a sensitivity analysis to determine whether an association is strongly influenced by any individual contributing GWAS dataset.

## Input METAL file

The original METAL analysis file is:

```bash
metalfile=$PWD/metal_all.meta
```

The script first checks that every summary-statistics file referenced by a `PROCESS` statement exists.

## Leave-one-study-out script

```bash
DirI=$PWD

metalfile=$DirI/metal_all.meta

# Check that all summary-statistics files referenced by
# PROCESS statements exist before starting the analyses.
for File in `grep PROCESS $metalfile | awk '{print $NF}'`
do
    if [ ! -f $File ]
    then
        echo "$File not exist"
        exit
    fi
done


# Create one meta-analysis for each PROCESS file,
# excluding that file from the analysis.
for File in `grep PROCESS $metalfile | awk -F"[ /]" '{print $NF}'`
do

    # Remove the extension to create a short analysis name.
    NewFile=${File%.*}

    # Create a dedicated directory for the leave-one-out analysis.
    newdir="del_"$NewFile
    mkdir -p $newdir

    cd $newdir

    # Create a new METAL input file with the selected
    # summary-statistics file removed.
    #
    # The OUTFILE is also changed so that each analysis
    # produces a separate output file.
    cat $metalfile | \
        grep -v $File | \
        sed "s/OUTFILE meta_hiv_africa_x .tbl/OUTFILE $newdir .tbl/" \
        > newmeta.in


    # Create an SLURM submission script for this
    # leave-one-study-out meta-analysis.
    echo """#!/bin/bash -l
#SBATCH -J $NewFile
#SBATCH -o $NewFile"_%j.txt"
#SBATCH -e $NewFile"_%j.txt"
#SBATCH -t 60:20:00
#SBATCH -c 4
#SBATCH -p batch
#SBATCH --mem=20GB

~/bin/METAL/build/metal/metal newmeta.in
""" > run.bash


    # Submit the meta-analysis to SLURM.
    sbatch run.bash

    cd $DirI

done
```

## How the script works

The script uses the original:

```text
metal_all.meta
```

as the template.

It identifies every file specified by a `PROCESS` statement:

```bash
grep PROCESS $metalfile
```

For each file, it:

1. Checks that the input file exists.
2. Creates a directory named `del_<summary_statistics_file>`.
3. Removes that one `PROCESS` line from the original METAL file.
4. Keeps all other `PROCESS` statements unchanged.
5. Changes the METAL output name to the leave-one-out directory name.
6. Creates an SLURM submission script.
7. Submits the analysis using `sbatch`.

## Example

If the original METAL file contains:

```text
PROCESS meta_x_fem_bc.txt
PROCESS meta_x_fem_cc.txt
PROCESS meta_x_fem_oc.txt
```

the script creates separate analyses such as:

```text
del_meta_x_fem_bc/
del_meta_x_fem_cc/
del_meta_x_fem_oc/
```

For example, the analysis in:

```text
del_meta_x_fem_bc/
```

contains all the original datasets **except**:

```text
meta_x_fem_bc.txt
```

The analysis in:

```text
del_meta_x_fem_cc/
```

contains all the original datasets **except**:

```text
meta_x_fem_cc.txt
```

and so on.

## Important distinction

This is a **leave-one-summary-statistics-file-out** analysis.

It does **not** remove an entire cohort at a time.

For example, because JCS contributes several separate GWAS summary-statistics files, the script will independently perform analyses leaving out:

```text
meta_x_fem_bc.txt
meta_x_fem_cc.txt
meta_x_fem_oc.txt
meta_x_male_oc.txt
meta_x_fem_nc.txt
meta_x_male_nc.txt
```

one at a time.

Likewise, ARK, AWI-Gen and Uganda files are removed **individually**, according to their individual `PROCESS` entries.

Therefore, the number of leave-one-out analyses equals the number of `PROCESS` statements in `metal_all.meta`.

## SLURM resources

Each leave-one-out meta-analysis is submitted independently with:

```text
CPU:        4 cores
Memory:     20 GB
Time:       60:20:00
Partition:  batch
```

METAL is executed using:

```bash
~/bin/METAL/build/metal/metal newmeta.in
```

## Output structure

Each analysis is stored in its own directory:

```text
del_<summary_statistics_file>/
```

with:

```text
newmeta.in
run.bash
<summary_statistics_file>_<SLURM_JOB_ID>.txt
<leave-one-out METAL output>
```

The resulting meta-analysis outputs can be compared with the original full Africa meta-analysis to assess the influence of each individual summary-statistics dataset.


## Forestplot
* `forestplot_sensitivity.r` : extract information from each study and plot by rsid

