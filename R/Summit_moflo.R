#loading and bouncing BioRad-S3e
#----------------------------------
library(dplyr)
library(ggplot2)
library(sp)

inputdir <- '~/DATA/BGU_Israel_2022/FC/MoFlo_/220828'
outdir <-paste0(inputdir,'/out')
dir.create(outdir)

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

myfiles <- list.files(path= inputdir, pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles[5], path= inputdir, alter.names=TRUE)
fs <- read.flowSet(myfiles[3], path= inputdir, alter.names=TRUE)


for(i in 1:length(fs)){
  df <- fs[[i]]@exprs %>% data.frame()

  #duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  #timeticks=fs[[i]]@description$TIMETICKS
  #patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id= gsub('\\.fcs','',fs[[i]]@description$'GUID')
  n_events <- nrow(df)
  long_string <- paste0(sample_id," (n=",n_events,")")

  write.table(df, file=paste0(outdir,'/data.00',i,'_',sample_id,'.txt'),
              quote = FALSE,
              sep='\t',
              row.names = FALSE)

  fl4_fl2 <- df %>% filter(FL1.Height > 5) %>% ggplot(aes(x = `FL4.Height`, y =  `FL2.Height`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 350)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(10001,10000000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(10001,10000000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)

  ggplot2::ggsave(paste0(outdir,'/',i,'_',sample_id, '_density-plots_fl4_fl2','.pdf'), fl4_fl2)
  medium_string <- paste0(patient_id," ", sample_id)

  ssc_fl2 <- df %>% filter(FL1.Height > 5) %>% ggplot(aes(x = `SSC.Height`, y =  `FL2.Height`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 350)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(100001,10000000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(10001,10000000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)

  ggplot2::ggsave(paste0(outdir,'/',i,'_',sample_id, '_density-plots_ssc_fl2','.pdf'), ssc_fl2)
  medium_string <- paste0(patient_id," ", sample_id)

  p.comb <- gridExtra::grid.arrange(ssc_fl2,fl4_fl2, nrow=1)

  ggplot2::ggsave(paste0(outdir,'/combined_',i,'_',sample_id, '_','.pdf'), p.comb,height = 15, width = 15)
}
