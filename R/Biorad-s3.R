#loading and bouncing BioRad-S3e
#----------------------------------
library(dplyr)
library(ggplot2)
library(sp)

inputdir <- '~/DATA/BioRad_S3e/10082022'
outdir <-paste0(inputdir,'/out')
dir.create(outdir)

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

myfiles <- list.files(path= inputdir, pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path= inputdir, alter.names=TRUE)


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

  gp <- df %>% filter(FL1.HEIGHT > 5) %>% ggplot(aes(x = `FSC.HEIGHT`, y =  `FL1.HEIGHT`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(101,1000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,1000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)

  ggplot2::ggsave(paste0(outdir,'/',i,'_',sample_id, '_density-plots_','.pdf'), gp)
  medium_string <- paste0(patient_id," ", sample_id)

  fsc_fl1 <- df %>% filter(FL1.HEIGHT > 5) %>% ggplot(aes(x = `FSC.HEIGHT`, y =  `FL1.HEIGHT`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(101,1000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,1000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  fsc_fl2 <- df %>% filter(FL1.HEIGHT > 5) %>% ggplot(aes(x = `FSC.HEIGHT`, y =  `FL2.HEIGHT`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(101,1000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(101,1000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  fsc_fl3 <- df %>% filter(FL1.HEIGHT > 5) %>% ggplot(aes(x = `FSC.HEIGHT`, y =  `FL3.HEIGHT`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(101,1000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,1000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  fsc_fl4 <- df %>% filter(FL1.HEIGHT > 5) %>% ggplot(aes(x = `FSC.HEIGHT`, y =  `FL4.HEIGHT`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(101,1000001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,1000001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  p.comb <- gridExtra::grid.arrange(fsc_fl1, fsc_fl2, fsc_fl3, fsc_fl4, nrow=2)

  ggplot2::ggsave(paste0(outdir,'/all_',i,'_',sample_id, '_','.pdf'), p.comb,height = 10, width = 10)
}
