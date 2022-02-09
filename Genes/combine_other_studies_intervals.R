# Make a list of intervals found in other studies and note article name for use in phenogram
library(tidyverse)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")

# load gene list of other studies
turpin<- read_xlsx("../other_studies/Other_studies.xlsx", sheet="Turpin_2016", skip=1)
turpin <- turpin %>% mutate(Start=BP-5e5, End=BP+5e5) %>% dplyr::select(CHR, Start, End)
colnames(turpin) <- c("Chr", "Start", "End")
turpin$article <- "Turpin et al. 2016"

wang<- read_xlsx("../other_studies/Other_studies.xlsx", sheet="Wang_2015")
wang <- wang %>%  dplyr::select(chromsome, start, end)
colnames(wang) <- c("Chr", "Start", "End")
wang$article <- "Wang et al. 2015"
wang$Start<- wang$Start*1e6
wang$End <- wang$End*1e6

leamy<- read_xlsx("../other_studies/Other_studies.xlsx", sheet="Leamy_2014")
leamy <- leamy %>% dplyr::select(Ch, `CI (Mb)`) %>% separate(`CI (Mb)`, sep="-", into=c("Start", "End")) %>% drop_na(Ch)
colnames(leamy)<- c("Chr", "Start", "End")
leamy$article <- "Leamy et al. 2014"
leamy$Start <- as.numeric(as.character(leamy$Start))*1e6
leamy$End <- as.numeric(as.character(leamy$End))*1e6

mcknite<- read_xlsx("../other_studies/Other_studies.xlsx", sheet="McKnite_2012")
mcknite <- mcknite %>% dplyr::select(chr,CI, ...6)
colnames(mcknite) <- c("Chr", "Start", "End")
mcknite$article <- "McKnite et al. 2012"
mcknite$Start <- mcknite$Start*1e6
mcknite$End <- mcknite$End*1e6

benson<- read_xlsx("../other_studies/Other_studies.xlsx", sheet="Benson_2010")
benson <- benson %>% dplyr::select(chromosome,`95% CI, Mb`) %>% drop_na(chromosome) %>% separate(`95% CI, Mb`, sep="-", into=c("Start", "End"))
colnames(benson)<- c("Chr", "Start", "End")
benson$article <- "Benson et al. 2010"
benson$Start <- as.numeric(as.character(benson$Start))*1e6
benson$End <- as.numeric(as.character(benson$End))*1e6

snijders<- read_xlsx("../other_studies/Other_studies.xlsx", sheet="Snijders_2017")
snijders <- snijders %>% dplyr::select(chr, interval_start.mm10, interval_stop.mm10)
snijders$chr<- gsub("chr", "", snijders$chr)
colnames(snijders)<- c("Chr", "Start", "End")
snijders$article <- "Snijders et al. 2017"

kemis <- read_xlsx("../other_studies/Other_studies.xlsx", sheet="Kemis_2019", skip=2)
kemis <- kemis %>% dplyr::select(Chr, ci_lo, ci_hi)
colnames(kemis)<- c("Chr", "Start", "End")
kemis$article <- "Kemis et al. 2019"
kemis$Start <- as.numeric(as.character(kemis$Start))*1e6
kemis$End <- as.numeric(as.character(kemis$End))*1e6
all_articles <- rbind(wang,leamy)

all_articles <- rbind(all_articles,mcknite)
all_articles <- rbind(all_articles,benson)
all_articles <- rbind(all_articles,snijders)
all_articles <- rbind(all_articles,kemis)
all_articles <- rbind(all_articles,turpin)
all_articles$Chr <- paste0("chr",trimws(all_articles$Chr))
# all_articles$Start<- as.numeric(as.character(trimws(all_articles$Start)))
# all_articles$End<- as.numeric(as.character(trimws(all_articles$End)))

# overlaps$startA<- as.numeric(as.character(overlaps$startA))


overlaps_with_article <- left_join(overlaps, all_articles, by=c("chrA"="Chr","startA"="Start")) %>% distinct(chrA, startA, endA, article)
overlaps_with_article_end <- left_join(overlaps, all_articles, by=c("chrA"="Chr", "endA"="End" )) %>% distinct(chrA, startA, endA, article)
overlaps_with_article$end_article <- overlaps_with_article_end$article
overlaps_with_article_total <- overlaps_with_article %>% mutate(all_articles=ifelse(is.na(article), end_article, article)) 

write_delim(overlaps_with_article_total, "overlaps_with_articles.csv", delim=";")
write_delim(all_articles, "all_articles.csv", delim=";")
# for phenogram

# in overlaps: a is other studies, b is our study
sum_freq_our_intervals <- overlaps %>% 
  group_by(chrB, startB, endB) %>% 
  count()