# HIV-1 GWAS Meta-analysis Without Self-Identified Information

This analysis evaluates the robustness of the HIV-1 GWAS meta-analysis when **self-identified information is not used for AWI-Gen and ARK participants**.

The analysis was performed in a separate working directory:

```text
/spaces/jeantristan/isabel/hiv/hiv_noselfid/awigen
```

## Rationale

Self-identified information was available for some participants in the AWI-Gen and ARK datasets.

For this sensitivity analysis, self-identified information was removed from the AWI-Gen and ARK analysis where applicable, and the GWAS and meta-analysis were repeated.

The other contributing datasets already contained the relevant information required for their analyses and were therefore not modified in the same way.

The purpose of this analysis is to determine whether removing self-identified information changes the GWAS or meta-analysis results.

---

# Workflow

The analysis consists of five main steps:

1. **Check and extract individuals according to the self-ID information**
2. **Run the sex/phenotype-specific SAIGE GWAS**
3. **Prepare the METAL configuration**
4. **Run the new meta-analysis**
5. **Compare the new meta-analysis with the original analysis**

```text
Original phenotype/genotype data
             |
             v
Check self-ID information
             |
             v
Extract eligible individuals
             |
             v
       run_saige_xx
             |
             v
      SAIGE GWAS results
             |
             v
       metal.config
             |
             v
        run_metal
             |
             v
   No-self-ID meta-analysis
             |
             v
Compare with original
meta-analysis results
```

---

# 1. Check and extract individuals using self-ID

The first step identifies the individuals that should be included in the no-self-ID analysis.

The purpose of this step is to verify the available self-identified information and generate the corresponding individual lists for the subsequent GWAS analyses.

The relevant script is:

```bash
check_extract_individual_with_selfid
```

This step is performed before running SAIGE.

The resulting participant lists are then

