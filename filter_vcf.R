library(dplyr)

vep_workir="/Documents/2024/"

## SNPs
snps <- read.table("snps.pass.info.genotypes.txt")
snps_df <- as.data.frame(snps)
#snps_df$V3 == "1/1"
FN <- snps_df %>%
  filter(if_all(V3:V7, `==`, c("1/1", "1|1")))

# NN <- snps_df %>% 
#   filter(if_all(V8:V12, `==`, "0/1")) 

FN_NN <- snps_df %>% 
  filter(if_all(V3:V7, `==`, "1/1") & if_any(V8:V12, `==`, "0/1") & if_any(V8:V12, `==`, "0/0"))

FN_NN1 <- snps_df %>% 
  filter(if_all(V3:V7, `==`, "1/1") & if_any(V8:V12, `==`, "0/1") & if_any(V8:V12, `==`, "0/0") & !if_any(V8:V12, `==`, "1/1"))
  
FN_NN2 <- snps_df %>% 
  filter(if_all(V3:V7, `==`, "1/1")  & !if_any(V8:V12, `==`, "1/1"))

FN_NN3 <- snps_df %>% 
  filter(if_all(V3:V7, `==`, "1|1")  & !if_any(V8:V12, `==`, "1|1"))

FN_NN4 <- snps_df %>% 
  filter(if_all(V3:V7, `==`, c("1/1", "1|1"))  & !if_any(V8:V12, `==`, c("1/1", "1|1")))

write.table(FN_NN2, "/Documents/2024/filtered.snp2.txt", sep = "\t")
write.table(FN_NN3, "/Documents/2024/filtered.snp3.txt", sep = "\t")
write.table(FN_NN4, "/Documents/2024/filtered.snp4.txt", sep = "\t")
