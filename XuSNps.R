


setwd("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/Xu_GRS_Mods")
SNPs <- vector()

for (mod in list.files("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/Xu_GRS_Mods")){

  df <- read.csv(mod, sep = '\t')  
  df <- df$rsid
  SNPs <- c(SNPs, df)
      
}

SNPs <- unique(SNPs)

write.csv(SNPs, file = 'XuSNPs.txt', row.names = F)
