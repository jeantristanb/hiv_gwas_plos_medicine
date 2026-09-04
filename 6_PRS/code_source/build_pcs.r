```r
# Calculate principal components (PCs) by study site and sex
# ------------------------------------------------------------

# Create output directory
dirout <- "pcs_res/"
dir.create(dirout, showWarnings = FALSE)

# Read PRS and phenotype data
data1 <- read.table(
    "all/hiv_eur/PRScs/eur/eur.profile",
    header = TRUE
)

data3 <- read.table(
    "all/hiv_eur/pheno/allsite_5_format.pheno",
    header = TRUE
)

# Merge PRS and phenotype data
all <- merge(
    data1,
    data3,
    by = c(1, 2)
)

# Calculate PRS
all$sumscore <- all$SCORE * all$CNT


# Function to calculate PCs
# ------------------------

calculate_pcs <- function(data, suffix) {

    for (site in unique(data$Study.Site)) {

        site_data <- data[data$Study.Site == site, ]

        fileout <- paste0(
            dirout,
            site,
            suffix
        )

        # Create PLINK individual list
        write.table(
            site_data[, c(1, 2)],
            file = paste0(fileout, ".ind"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )

        # LD pruning
        system(
            paste(
                "~/bin/plink2",
                "-bfile array_plk",
                "--keep", paste0(fileout, ".ind"),
                "--indep-pairwise 100 10 0.1",
                "-maf 0.01",
                "--out", fileout
            )
        )

        # Calculate 20 principal components
        system(
            paste(
                "~/bin/plink2",
                "-bfile array_plk",
                "--keep", paste0(fileout, ".ind"),
                "--extract", paste0(fileout, ".prune.in"),
                "--pca 20",
                "--out", fileout
            )
        )
    }
}


# PCs for all participants
# ------------------------

calculate_pcs(
    all,
    "_all"
)


# PCs for men
# -----------

all_men <- all[all$sex == 0, ]

calculate_pcs(
    all_men,
    "_men"
)


# PCs for women
# -------------

all_women <- all[all$sex == 1, ]

calculate_pcs(
    all_women,
    "_women"
)
```

