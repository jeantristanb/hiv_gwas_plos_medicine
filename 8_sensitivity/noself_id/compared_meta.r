library(data.table)
library(dplyr)
DataI<-fread('/dataG/cancer/isabel_work/pGWAS/metal/HIV/meta_hiv_africa_x/meta_hiv_africa_x1.tbl')
Data_nohivself<-fread('../meta/meta_hiv_africa_x1.tbl')
allmerge<-merge(DataI,Data_nohivself, by='MarkerName',all.x=T,suffixes=c('.inital','.nohivself'))
allmerge_below<-allmerge[ allmerge$"P-value.inital"<5*10**-6 | allmerge$"P-value.nohivself"<5*10**-6, ]

write.csv(allmerge_below, file='comparison_pvalue.csv')
allmerge_below<-read.csv('comparison_pvalue.csv')
plot(-log10(allmerge_below$P.value.inital), -log10(allmerge_below$P.value.nohivself))
Data<-allmerge_below[ allmerge_below$P.value.inital<=5*10**-8,c('MarkerName','Chromosome.inital','Position.inital',"Allele1.inital","Allele2.inital",'Freq1.inital','Effect.inital','StdErr.inital','P.value.inital','Freq1.nohivself','Effect.nohivself','StdErr.nohivself','P.value.nohivself')]


Data_out <- Data %>%
  mutate(
    across(
      c(Allele1.inital, Allele2.inital),
      toupper
    ),
    across(
      c(Freq1.inital, Effect.inital, StdErr.inital,
        Freq1.nohivself, Effect.nohivself, StdErr.nohivself),
      ~ formatC(., format = "f", digits = 2)
    ),
    across(
      c(P.value.inital, P.value.nohivself),
      ~ formatC(., format = "e", digits = 1)
    )
  )
write.csv(Data_out, "comparison_pvalue_sig.csv", row.names = FALSE, quote = FALSE)
