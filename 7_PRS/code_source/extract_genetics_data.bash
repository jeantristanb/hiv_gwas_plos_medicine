### Genotype data preparation

The following commands prepare the PLINK genotype datasets used for the PRS and PCA analyses.

First, genomic positions are extracted from the reference PLINK `.bim` file:

```bash
awk '{print $1"\t"$4"\t"$4}' \
    "$dirdata/all.bim" > all_pos.bim
```

The corresponding variants are then extracted from the reference dataset:

```bash
plink \
    -bfile "$dirdata/all" \
    --extract range all_pos.bim \
    --make-bed \
    -out all_plk
```

The same genomic positions are extracted from the cleaned genotype dataset used for the analysis:

```bash
plink \
    -bfile "$dirdatsa/all_clean" \
    --extract range all_pos.bim \
    --make-bed \
    -out array_plk
```

This produces the PLINK dataset:

```text
array_plk
```

which is subsequently used for principal component calculation and downstream PRS analyses.

