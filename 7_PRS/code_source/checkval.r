```r
# PRS statistical analysis
# ========================

library(pROC)


# -------------------------------------------------------------------------
# 1. Read and merge data
# -------------------------------------------------------------------------

prs_data <- read.table(
    "random1/hiv_eur/PRScs/eur/eur.profile",
    header = TRUE
)

phenotype <- read.table(
    "../data/pheno/allsite.tsv",
    header = TRUE
)

phenotype_format <- read.table(
    "random1/hiv_eur/pheno/allsite_5_format.pheno",
    header = TRUE
)

# Merge PRS and phenotype data
all_data <- merge(
    prs_data,
    phenotype,
    by = c(1, 2)
)

all_data <- merge(
    all_data,
    phenotype_format,
    by = c(1, 2)
)

# Calculate PRS
all_data$sumscore <- all_data$SCORE * all_data$CNT


# -------------------------------------------------------------------------
# 2. Function for PRS analysis
# -------------------------------------------------------------------------

run_prs_analysis <- function(data, residual_covariates,
                             regression_covariates, output_file) {

    results <- list()

    for (site in unique(data$Study.Site)) {

        # Select individuals from the current study site
        site_data <- data[data$Study.Site == site, ]

        # ---------------------------------------------------------------
        # Residualize PRS
        # ---------------------------------------------------------------

        residual_formula <- as.formula(
            paste(
                "sumscore ~",
                paste(residual_covariates, collapse = " + ")
            )
        )

        site_data$sumscore_residuals <- residuals(
            glm(
                residual_formula,
                data = site_data
            )
        )

        # ---------------------------------------------------------------
        # ROC / AUC analysis
        # ---------------------------------------------------------------

        roc_auc <- roc(
            site_data$hiv,
            site_data$sumscore_residuals,
            direction = "<"
        )

        auc_ci <- ci(roc_auc)

        auc_results <- data.frame(
            auc = as.numeric(roc_auc$auc),
            auc_ci5 = as.numeric(auc_ci[1]),
            auc_ci95 = as.numeric(auc_ci[3])
        )

        # ---------------------------------------------------------------
        # Logistic regression
        # ---------------------------------------------------------------

        regression_formula <- as.formula(
            paste(
                "hiv ~",
                paste(regression_covariates, collapse = " + ")
            )
        )

        model <- glm(
            regression_formula,
            data = site_data,
            family = binomial()
        )

        model_summary <- summary(model)

        # Extract PRS coefficient
        prs_result <- model_summary$coefficients["sumscore", ]

        # ---------------------------------------------------------------
        # Store results
        # ---------------------------------------------------------------

        results[[site]] <- c(
            study.site = site,
            as.list(prs_result),
            auc_results,
            n = nrow(site_data)
        )
    }

    # Combine results
    results <- do.call(
        rbind,
        lapply(results, as.data.frame)
    )

    rownames(results) <- NULL

    # Write results
    write.csv(
        results,
        file = output_file,
        row.names = FALSE
    )

    return(results)
}


# -------------------------------------------------------------------------
# 3. Analysis of all participants
# -------------------------------------------------------------------------

all_results <- run_prs_analysis(
    data = all_data,
    residual_covariates = c(
        "Pcs_1",
        "Pcs_2",
        "Pcs_3",
        "Pcs_4",
        "sex",
        "age"
    ),
    regression_covariates = c(
        "sumscore",
        "sex",
        "age",
        "Pcs_1",
        "Pcs_2",
        "Pcs_3",
        "Pcs_4"
    ),
    output_file = "all_prsrs.csv"
)


# -------------------------------------------------------------------------
# 4. Analysis of men
# -------------------------------------------------------------------------

all_men <- all_data[all_data$sex == 0, ]

men_results <- run_prs_analysis(
    data = all_men,
    residual_covariates = c(
        "Pcs_1",
        "Pcs_2",
        "Pcs_3",
        "Pcs_4",
        "age"
    ),
    regression_covariates = c(
        "sumscore",
        "sex",
        "age",
        "Pcs_1",
        "Pcs_2",
        "Pcs_3",
        "Pcs_4"
    ),
    output_file = "all_men_prsrs.csv"
)


# -------------------------------------------------------------------------
# 5. Analysis of women
# -------------------------------------------------------------------------

all_women <- all_data[all_data$sex == 1, ]

women_results <- run_prs_analysis(
    data = all_women,
    residual_covariates = c(
        "Pcs_1",
        "Pcs_2",
        "Pcs_3",
        "Pcs_4",
        "age"
    ),
    regression_covariates = c(
        "sumscore",
        "sex",
        "age",
        "Pcs_1",
        "Pcs_2",
        "Pcs_3",
        "Pcs_4"
    ),
    output_file = "all_women_prsrs.csv"
)
```

