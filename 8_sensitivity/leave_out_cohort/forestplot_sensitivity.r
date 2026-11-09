cmt<-1
headfile<-read.table('../del_meta_x_ark_fem/del_meta_x_ark_fem1.tbl', header=T,comment.char="", nrows=1)
rowname_file<-names(headfile)

for(File in system('ls *.tbl' ,intern=T)){
Data<-read.table(File, header=F,comment.char="")
Data$File<-File
if(cmt==1) DataF<-Data
else DataF<-rbind(DataF,Data)
cmt<-cmt+1
}
names(DataF)<-c(rowname_file,'cohort_Del')
allbp<-read.table('bp_snp.bp')
allbp$Chromosome<-paste('chr',sapply(strsplit(allbp$V1,split=':'),function(x)x[1]),sep='')
allbp$Position<-as.integer(sapply(strsplit(allbp$V1,split=':'),function(x)x[2]))
all_filter<-merge(allbp,DataF,all.y=F, all.x=T, by=c('Chromosome','Position'))

library(dplyr)
library(ggplot2)

plot_data <- all_filter %>%
  mutate(
    OR = exp(Effect),
    CI_low = exp(Effect - 1.96 * StdErr),
    CI_high = exp(Effect + 1.96 * StdErr),
    cohort_Del = sub("\\.tbl$", "", cohort_Del)
  )

plot_data$old_cohort_Del<-plot_data$cohort_Del

plot_data <- plot_data %>%
  mutate(
    cohort_label = cohort_Del %>%
      sub("^del_meta_x_", "", .) %>%
      sub("^del_mungestats_uganda_", "uganda_", .) %>%
      sub("^meta_hiv_africa_x", "hiv_africa", .) %>%
      sub("1$", "", .) %>%
      sub("hiv_africa_jcs1.tbl", "jcs_all", .) 
  )
cohort_labels_cord <- c(
  "fem_bc"              = "BC (F)",
  "ark_fem"             = "ARK (F)",
  "ark_male"            = "ARK (M)",
  "fem_oc"              = "OC (F)",
  "male_oc"             = "OC (M)",
  "fem_cc"              = "CC (F)",
  "fem_nc"              = "NC (F)",
  "male_nc"             = "NC (M)",
  "jcs_all"             = "JCs",
  "awi_east_fem"        = "AWI-East (F)",
  "awi_east_male"       = "AWI-East (M)",
  "awi_south_fem"       = "AWI-South (F)",
  "awi_south_male"      = "AWI-South (M)",
  "uganda_fem_chrprefix" = "Uganda (F)",
  "uganda_male_chrprefix" = "Uganda (M)",
  "hiv_africa"        = "HIV (Full)"
)

cohort_order <- c(
  "HIV (Full)",
  "JCs",
  # JCS
  "BC (F)",
  "CC (F)",
  "NC (F)",
  "NC (M)",
  "OC (F)",
  "OC (M)",

  # ARK
  "ARK (F)",
  "ARK (M)",

  # AWI-Gen
  "AWI-East (F)",
  "AWI-East (M)",
  "AWI-South (F)",
  "AWI-South (M)",

  # Uganda
  "Uganda (F)",
  "Uganda (M)"
)

plot_data <- plot_data %>%
  mutate(
    cohort_discarded = cohort_labels_cord[cohort_label]
  )

rsid_mapping <- c(
  "chr1:216797017" = "rs112358207",
  "chr3:46367240"  = "rs2040388",
  "chr3:46370051"  = "rs2227010",
  "chr8:28875139"  = "rs113298445"
)

plot_data <- plot_data %>%
  mutate(
    rsID = rsid_mapping[MarkerName]
  )

plot_data <- plot_data %>%
  mutate(
    cohort_discarded = factor(
      cohort_discarded,
      levels = rev(cohort_order)
    )
  )
x <- plot_data

ggplot(x, aes(y = cohort_discarded, x = OR)) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +

  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      colour = cohort_discarded == "HIV (Full)"
    ),
    orientation = "y",
    height = 0.2
  ) +

  geom_point(
    aes(
      colour = cohort_discarded == "HIV (Full)"
    ),
    size = 2
  ) +

  scale_colour_manual(
    values = c(
      "FALSE" = "black",
      "TRUE"  = "red"
    ),
    guide = "none"
  ) +

  scale_x_log10() +

  facet_wrap(~ rsID,scales = "free_x") + labs(x = "Odds ratio (95% CI)",y = "Cohort excluded") +theme_classic()
ggsave('comparison_model_discarded.tif')
