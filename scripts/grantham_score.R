library(tidyverse)
# Calculate grantham score

grantham <- readr::read_tsv("https://gist.githubusercontent.com/danielecook/501f03650bca6a3db31ff3af2d413d2a/raw/5583a134b36b60762be6cd54002a0f4044338cd0/grantham.tsv") %>%
  tidyr:: gather(SECOND,SCORE, -FIRST) %>% dplyr::filter(SCORE > 0)  

calculate_grantham <- function(a1, a2) {
  (grantham %>% dplyr::filter(FIRST == a1, SECOND == a2))$SCORE
}

grantham

sample_id <- c("ERR10513574")

df <- read.csv(paste("../data/samples/PRJEB56841/dRNA/", sample_id, "_APOBEC-like.vcf", sep=""), sep = "\t", header = TRUE)
df_CDS <- df %>% filter(mutation_category != "intergenic") %>% filter(motif == 'TC' | motif == 'GA')
OUTFILE <- paste("../data/samples/PRJEB56841/dRNA/", sample_id, "_APOBEC-like_CDS.vcf", sep="")

score <- mapply(calculate_grantham, df_CDS$parent_aa, df_CDS$mutated_aa, SIMPLIFY = FALSE)
score <- as.numeric(score)
df_CDS$grantham_score <- score
df_CDS[is.na(df_CDS)] <- 0

#0–50 conservative, yellow; 51–100 moderately conservative, pink; 101–150 moderately radical, blue; >150 radical, brown; +synonymous, green1
df_CDS <- df_CDS %>% mutate(grantham_rank = case_when(
  grantham_score > 50 & grantham_score < 101 ~ "moderately conservative",
  grantham_score > 100 & grantham_score < 151 ~ "moderately radical",
  grantham_score > 149 ~ "radical",
  grantham_score < 51 & mutation_category != "synonymous" ~ "conservative",
  mutation_category == "synonymous" ~ "synonymous",
))
df_CDS <- df_CDS %>% mutate(grantham_rank_color = case_when(
  grantham_score > 50 & grantham_score < 101 ~ pink,
  grantham_score > 100 & grantham_score < 151 ~ blue,
  grantham_score > 149 ~ brown,
  grantham_score < 51 & mutation_category != "synonymous" ~ yellow,
  mutation_category == "synonymous" ~ green1,
))

df_CDS$parent_aa[df_CDS$mutation_category == 'synonymous'] <- '='
df_CDS$mutated_aa[df_CDS$mutation_category == 'synonymous'] <- '='

write.csv(df_CDS, OUTFILE, row.names = FALSE)
