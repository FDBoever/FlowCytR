#loading and bouncing FACsSortData

# run this in terminal
#----------------------------------

#cd ~/DATA/FACSort/2207/2707
#for f in Data.*; do mv "$f" "$f.fcs"; done


library(dplyr)
library(ggplot2)
library(sp)

inputdir <- '/Users/sa01fd/DATA/FACSort/2208/080822'
inputdir <- '~/DATA/FACSort/2211/1129/'
inputdir <- '~/DATA/FACSort/2212/1206/'


outdir <-paste0(inputdir,'/out')
dir.create(outdir)

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

myfiles <- list.files(path= inputdir, pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path= inputdir, alter.names=TRUE)

for(i in 1:length(fs)){
  df <- fs[[i]]@exprs %>% data.frame()

  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")

  write.table(df, file=paste0(outdir,'/data.00',i,'_',sample_id,'.txt'),
              quote = FALSE,
              sep='\t',
              row.names = FALSE)

  gp <- df %>% dplyr::filter(FL1.H > 5) %>%
    ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL1.H`)) +
    ggplot2::scale_fill_gradientn(colours = myColor) +
    ggplot2::stat_binhex(bins = 200)+
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio=1)+
    ggplot2::annotation_logticks() +
    ggplot2::ggtitle(long_string)

  ggplot2::ggsave(paste0(outdir,'/',i,'_',sample_id, '_density-plots_','.pdf'), gp)
  medium_string <- paste0(patient_id," ", sample_id)

  ssc_fl1 <- df %>% dplyr::filter(FL1.H > 5) %>%
    ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL1.H`)) +
    ggplot2::scale_fill_gradientn(colours = myColor) +
    ggplot2::stat_binhex(bins = 200)+
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio=1)+
    ggplot2::annotation_logticks() +
    ggplot2::ggtitle(medium_string)

  ssc_fl2 <- df %>% dplyr::filter(FL1.H > 5) %>%
    ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL2.H`)) +
    ggplot2::scale_fill_gradientn(colours = myColor) +
    ggplot2::stat_binhex(bins = 200)+
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio=1)+
    ggplot2::annotation_logticks() +
    ggplot2::ggtitle(medium_string)

  ssc_fl3 <- df %>% dplyr::filter(FL1.H > 5) %>%
    ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL3.H`)) +
    ggplot2::scale_fill_gradientn(colours = myColor) +
    ggplot2::stat_binhex(bins = 200)+
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio=1)+
    ggplot2::annotation_logticks() +
    ggplot2::ggtitle(medium_string)

  fl3_fl1 <- df %>%
    dplyr::filter(FL1.H > 5) %>%
    ggplot2::ggplot(ggplot2::aes(x = `FL3.H`, y =  `FL1.H`)) +
    ggplot2::scale_fill_gradientn(colours = myColor) +
    ggplot2::stat_binhex(bins = 200)+
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio=1)+
    ggplot2::annotation_logticks() +
    ggplot2::ggtitle(medium_string)

  p.comb <- gridExtra::grid.arrange(ssc_fl1, ssc_fl2, ssc_fl3, fl3_fl1, nrow=2)

  ggplot2::ggsave(paste0(outdir,'/all_',i,'_',sample_id, '_','.pdf'), p.comb,height = 10, width = 10)
}




df.stat <- read.delim('~/DATA/FACSort/2212/1206/stat_final_1206') %>% data.frame()

outdf <- c()
for(i in 1:length(fs)){
  fname <- myfiles[i]
  df <- fs[[i]]@exprs %>% data.frame()
  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")

  outdf <- rbind(outdf,
                 c('FileName' = fname,
                   'duration' = duration,
                   'timetick' = timeticks,
                   'patient_id' = patient_id,
                   'sample_id' = sample_id,
                   'n_events' =  n_events))
}

outdf <- outdf %>% data.frame() %>%
  dplyr::mutate(File.Name = sub('\\.fcs','',FileName))


df.stat <- df.stat %>% left_join(outdf, by='File.Name')


write.table(df.stat, file=paste0(outdir,'/stat_out1206.txt'),
            quote = FALSE,
            sep='\t',
            row.names = FALSE)










df %>% dplyr::filter(FL1.H > 5) %>%
  dplyr::slice(1:10000) %>%
  ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL1.H`)) +
  ggplot2::scale_fill_gradientn(colours = myColor) +
  ggplot2::geom_point(size=0.1)+
  ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         limits = c(1,10001)) +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         limits = c(4,10001)) +
  ggplot2::theme_classic() +
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::annotation_logticks() +
  ggplot2::ggtitle(long_string)


df %>% dplyr::filter(FL1.H > 5) %>%
  dplyr::slice(1:10000) %>%
  ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL1.H`)) +
  ggplot2::scale_fill_gradientn(colours = myColor) +
  ggplot2::geom_point(size=0.1)+
  ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         limits = c(1,10001)) +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         limits = c(4,10001)) +
  ggplot2::theme_classic() +
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::annotation_logticks() +
  ggplot2::ggtitle(long_string)



greys <- grep("gr[ea]y", colours(), value = TRUE)
df %>% dplyr::filter(FL1.H > 5) %>%
  dplyr::slice(1:50000) %>%
  ggplot2::ggplot(ggplot2::aes(x = `SSC.H`, y =  `FL1.H`)) +
  ggplot2::stat_density2d(ggplot2::aes(fill=..level..), geom="polygon",bins=75)+
  ggplot2::geom_point(size=0.05,color='red',alpha=0.7)+
  ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         limits = c(1,10001)) +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         limits = c(4,10001)) +
  ggplot2::theme_classic() +
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::annotation_logticks() +
  ggplot2::ggtitle(long_string) +
  ggplot2::scale_fill_gradientn(colours = greys)









for(i in 1:length(fs)){

  df <- fs[[i]]@exprs %>% data.frame()


  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")

  write.table(df, file=paste0(outdir,'/data.00',i,'_',sample_id,'.txt'),
              quote = FALSE,
              sep='\t',
              row.names = FALSE)

  gp <- df %>% filter(FL1.H > 5) %>% slice(1:10000) %>% ggplot(aes(x = `SSC.H`, y =  `FL1.H`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)

  ggplot2::ggsave(paste0(outdir,'/10k_',i,'_',sample_id, '_density-plots_','.pdf'), gp)


  medium_string <- paste0(patient_id," ", sample_id)

}






df %>% filter(FL1.H > 5) %>% slice(1:10000) %>% ggplot(aes(x =  `FL1.H`)) +
  geom_histogram(bins = 200) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  theme_classic() + theme(aspect.ratio=0.25)+ annotation_logticks(sides='b') + ggtitle(long_string)





df %>% filter(FL1.H > 5) %>% ggplot(aes(x =  `FL1.H`)) +
  geom_histogram(bins = 200) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  theme_classic() + theme(aspect.ratio=0.25)+ annotation_logticks(sides='b') + ggtitle(long_string)




df %>% filter(FL1.H > 5) %>% ggplot(aes(x =  `FL1.H`)) +
  geom_histogram(bins = 200) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  theme_classic() + theme(aspect.ratio=0.25)+ annotation_logticks(sides='b') + ggtitle(long_string)





i=169
df <- fs[[i]]@exprs %>% data.frame()


duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
timeticks=fs[[i]]@description$TIMETICKS
patient_id=fs[[i]]@description$'PATIENT ID'
sample_id=fs[[i]]@description$'SAMPLE ID'
n_events <- nrow(df)

medium_string <- paste0(patient_id," ", sample_id)

ssc_fl1 <- df %>% filter(FL1.H > 2) %>% ggplot(aes(x = `SSC.H`, y =  `FL1.H`)) +
  scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(4,10001)) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(medium_string)


ssc_fl3 <- df %>% filter(FL1.H > 2) %>% ggplot(aes(x = `SSC.H`, y =  `FL3.H`)) +
  scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(4,10001)) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)



fl3_fl1 <- df %>% filter(FL1.H > 2) %>% ggplot(aes(x = `FL3.H`, y =  `FL1.H`)) +
  scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(4,10001)) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)


grid.arrange(ssc_fl1, ssc_fl3, fl3_fl1, nrow=1)



polygon.x <- c(600,80,80,600)
polygon.y <- c(80,80,600,600)




df %>% filter(FL1.H > 2) %>% ggplot(aes(x = `SSC.H`, y =  `FL1.H`)) +
  geom_point(size=0.5, aes(color=gate)) +
  #scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(4,10001)) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string) +
  geom_polygon(data=data.frame(cbind(polygon.x,polygon.y)),aes(x=polygon.x, y=polygon.y),alpha=0,color='black')

df %>% filter(FL1.H > 2) %>% ggplot(aes(x = `SSC.H`, y =  `FL1.H`)) +
  scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1,10001)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(4,10001)) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)+
  geom_polygon(data=data.frame(cbind(polygon.x,polygon.y)),aes(x=polygon.x, y=polygon.y),alpha=0,color='black')



polygon.x <- c(200,150, 75, 60, 60, 75, 150, 200)
polygon.y <- c(75, 60, 60, 75 , 150, 200, 200, 150)


for(i in 1:length(fs)){

  df <- fs[[i]]@exprs %>% data.frame()
  df$gate <- point.in.polygon(df$SSC.H, df$FL1.H, polygon.x,polygon.y)


  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")

  write.table(df, file=paste0(outdir,'/data.00',i,'_',sample_id,'.txt'),
              quote = FALSE,
              sep='\t',
              row.names = FALSE)

  gp <- df %>% filter(FL1.H > 2) %>% ggplot(aes(x = `SSC.H`, y =  `FL1.H`)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(1,10001)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  limits = c(4,10001)) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + ggtitle(long_string)+
    geom_polygon(data=data.frame(cbind(polygon.x,polygon.y)),aes(x=polygon.x, y=polygon.y),alpha=0,color='black')


  ggplot2::ggsave(paste0(outdir,'/gated_',i,'_',sample_id, '_density-plots_','.pdf'), gp)

}
