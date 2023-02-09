

myfiles <- list.files(path="~/DATA/FACSort/data/test", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="~/DATA/FACSort/data/test")

fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)


#
fs_comp <- fs

#Quality control and cleanup if needed
# anomalous cells detected via flow rate checks, aqcuisition checks and dynamic range checks
# THIS IS GREAT IF YOU WANT TO STANDARDISE THE CLEANUP PROCESS, BUT IF IT IS CLEAN OR CLEANED, OR NO PROBLEMS WITH DATA< THIS STEP IS NOT NECCESARY
# best served when one has cloggy samples that may have pulsing of the fluidics and will cut out the bad clogs, and stitch the data back
# -- it does not replace traditional gating algorithms, or manual gating
# -- if you use flowAI or flowClean, still do traditional cleanup heirarchy
fs_comp_clean <- flowAI::flow_auto_qc(fs_comp)


trans <- flowCore::estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,3:7]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)






library("Phenoflow")

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)

library(dplyr)


#Load data
myfiles <- list.files(path="C:/Users/hallc/Downloads/FlowRepository_FR-FCM-ZZZV_files", pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles[1:10], path="C:/Users/hallc/Downloads/FlowRepository_FR-FCM-ZZZV_files", alter.names=TRUE)

#Load data
myfiles <- list.files(path="/Users/sa01fd/DATA/FACSort/2110/21/fcs", pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="/Users/sa01fd/DATA/FACSort/2110/21/fcs", alter.names=TRUE)

myfiles <- list.files(path="/Users/sa01fd/DATA/FACSort/2111/25/fcs", pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="/Users/sa01fd/DATA/FACSort/2111/25/fcs", alter.names=TRUE)




myfiles <- list.files(path="/Users/sa01fd/DATA/FACSort/2111/25/fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="/Users/sa01fd/DATA/FACSort/2111/25/fcs", alter.names=TRUE)

myfiles <- list.files(path="/Users/sa01fd/DATA/FACSort/2205/0405/Rock", pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="/Users/sa01fd/DATA/FACSort/2205/0405/Rock", alter.names=TRUE)





#----------------------- QUICKBOUNCE MIKE MOFLO
myfiles <- '3.1.16_170m_5L_0.4um_conc_4%PFA_Hst_TEM.fcs'
fs <- read.flowSet(myfiles, path="/Users/sa01fd/Downloads/", alter.names=TRUE)



df <- fs[[1]]@exprs %>% data.frame()

df %>% head(n=10000) %>% ggplot(aes(FSC.Height, FL1.Height)) + geom_point(size=0.1,color='darkgrey') + geom_density_2d(aes(colour = ..level..), contour_var = "count") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))

df %>% dplyr::filter(FSC.Height > 1000000) %>%
  dplyr::filter(FL1.Height > 1000000) %>%
  ggplot(aes(FSC.Height, FL1.Height)) + geom_point(size=0.1,color='darkgrey') + # geom_density_2d(aes(colour = ..level..), contour_var = "count") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+theme(aspect.ratio=1)



for(ct in c(1000000,2500000,3000000,4000000,5000000,10000000)){
  df.cutoff <- df %>% dplyr::filter(FSC.Height > ct) %>%
    dplyr::filter(FL1.Height > ct)
  message(dim(df.cutoff))
  write.table(df.cutoff, file=paste0('~/DATA/3.1.16_170m_5L_0.4um_conc_4%PFA_Hst_TEM_FSC_FL1_cutoff_at_',ct/1000000,'mil.txt'),
              quote = FALSE,
              sep='\t',
              row.names = FALSE)
}



df %>% dplyr::filter(FSC.Height > 1000000) %>%
  dplyr::filter(FL1.Height > 1000000) %>%
  ggplot(aes(FSC.Height, FL1.Height)) + geom_point(size=0.1,color='darkgrey') + # geom_density_2d(aes(colour = ..level..), contour_var = "count") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+theme(aspect.ratio=1)



write.table(df, file='~/DATA/3.1.16_170m_5L_0.4um_conc_4%PFA_Hst_TEM.txt',
            quote = FALSE,
            sep='\t',
            row.names = FALSE)

write.table(df %>% head(n=100000), file='~/DATA/3.1.16_170m_5L_0.4um_conc_4%PFA_Hst_TEM_100000.txt',
            quote = FALSE,
            sep='\t',
            row.names = FALSE)
#-----------------------

#BASH
#for f in fltr_bac.*; do mv "$f" "$f.fcs"; done

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}






myfiles <- list.files(path="/Users/sa01fd/DATA/FACSort/2112/13",pattern='.fcs', ignore.case = TRUE)
#myfiles <- myfiles[grepl('fltr_bac',myfiles)]
myfiles <- myfiles[1:length(myfiles)]
fs <- read.flowSet(myfiles, path="/Users/sa01fd/DATA/FACSort/2112/13", alter.names=TRUE)


i=2
df <- fs[[i]]@exprs %>% data.frame()
df$density <- get_density(df$SSC.H, df$FL1.H, n = 100)

duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
timeticks=fs[[i]]@description$TIMETICKS
patient_id=fs[[i]]@description$'PATIENT ID'
sample_id=fs[[i]]@description$'SAMPLE ID'
n_events <- nrow(df)
long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")

#
# df %>% head(n=100000) %>%
#   ggplot(aes(SSC.H,  FL1.H)) +
#   geom_point(size=0.1,color='darkgrey') +theme_classic() +
#   #geom_density_2d(aes(colour = ..level..), contour_var = "count") +
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   #coord_cartesian(xlim=c(1,10000),ylim=c(1,10000)) +
#                 ggtitle(long_string)


# non-binned dot plot
fr <- fs[[i]]
df <- fortify(fr)

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))



x_thr <- 140
y_thr <- 140

ggplot(df, aes(x = `SSC.H`, y =  `FL1.H`)) +
  scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()+
  geom_hline(yintercept = y_thr)+
  geom_vline(xintercept = x_thr)+
  geom_hline(yintercept = 10)




df %>% dplyr::mutate(G1 = ifelse(SSC.H>x_thr & FL1.H>y_thr,'euk','bac')) %>%
  dplyr::mutate(G1 = ifelse(FL1.H< 10 ,'debris',G1)) %>%
  group_by(G1) %>% tally()








df %>% #head(n=10000) %>%
  ggplot(aes(SSC.H,  FL1.H)) +
  theme_classic() +
  #geom_density_2d(aes(colour = ..level..), contour_var = "count") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
 # coord_cartesian(xlim=c(1,10000),ylim=c(1,10000)) +
  ggtitle(long_string)





for(i in st:end){

}


#=====
# data is stored in fs[[1]]@exprs
i=5
n_experiments = 14

st<-1
end<-4

df.all <- c()
for(i in st:end){
  df.a <-  fs[[i]]@exprs %>% data.frame()
  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df.a)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")
  df.a$duration = duration
  df.a$patient_id = patient_id
  df.a$sample_id = sample_id
  df.a$n_events = n_events
  df.a$i = i
  df.all <- rbind(df.all, df.a)

}





df.all %>% dplyr::mutate(G1 = ifelse(SSC.H>x_thr & FL1.H>y_thr,'euk','bac')) %>%
  dplyr::mutate(G1 = ifelse(FL1.H< 10 ,'debris',G1)) %>%
  group_by(G1,patient_id,sample_id,i) %>% tally() %>%
  ggplot(aes(x=i,y=n))+geom_bar(aes(fill=G1),position="stack", stat="identity")


df.all %>% dplyr::mutate(G1 = ifelse(SSC.H>x_thr & FL1.H>y_thr,'euk','bac')) %>%
  dplyr::mutate(G1 = ifelse(FL1.H< 10 ,'debris',G1)) %>%
  group_by(G1,patient_id,sample_id,i) %>% tally() %>%
  ggplot(aes(x=i,y=n))+geom_bar(aes(fill=G1),position="fill", stat="identity")



plist<-list()
for(grpl in c('fdb33','fdb36','ij30','ij29','ij31','ij159')){
  df.sel <- df.all %>%
    tibble() %>%
    mutate(sample_id=as.character(sample_id)) %>%
    dplyr::filter(grepl(grpl,sample_id))

  plist<-list()
  for(i in unique(df.sel$i)){
    message(i)
    plist[[i]] <- df.sel %>% ggplot(aes(SSC.H, FL1.H))  +
      scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + facet_wrap(~i,ncol=2) + ggtitle(grpl)

  }
  p1 <- arrangeGrob(grobs = plist, ncol = 2)

  ggplot2::ggsave(paste0(grpl, 'density-plots_','.pdf'), p1)

  }

multiplot(plotlist=braylist,cols=10)




for(i in st:end){
  df.1 <- fs[[i]]@exprs %>% data.frame()
  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df.1)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")


  tryCatch({
  p1 <- df.1 %>% ggplot(aes(SSC.H, FL3.H)) + geom_point(size=0.1,color='darkgrey') + geom_density_2d(aes(colour = ..level..), contour_var = "count") +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_color_continuous(type = "viridis")

  p2 <- df.1 %>% ggplot(aes(FL1.H, FL3.H)) + geom_point(size=0.1,color='darkgrey') + geom_density_2d(aes(colour = ..level..), contour_var = "count") +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_color_continuous(type = "viridis")

  p3 <- df.1 %>% ggplot(aes(SSC.H, FL1.H)) + geom_point(size=0.1,color='darkgrey') + geom_density_2d(aes(colour = ..level..), contour_var = "count") +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_color_continuous(type = "viridis")


  #p4 <- df.1 %>% ggplot(aes(SSC.H, FL2.H)) + geom_point(size=0.1,color='darkgrey') + geom_density_2d(aes(colour = ..level..), contour_var = "count") +
  #  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_color_continuous(type = "viridis")


 #p.all <- cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
  p.all <- cowplot::plot_grid(p1, p2,p3, labels = "AUTO")

  # now add the title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      long_string,
      x = 0,
      hjust = 0,
      size=14
    ) +
    ggplot2::theme(
      plot.margin = margin(0, 0, 0, 7)
    )

  p.final <- cowplot::plot_grid(
    title, p.all,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )


  ggplot2::ggsave(paste0(sample_id," ", patient_id,'_','density-plots_',i,'.pdf'), p.final)

  },
 error = function(e)
   print("You can't calculate the log of a character"))


  tryCatch({


  p1 <-  df.1 %>% ggplot(aes(SSC.H, FL1.H))  +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  p2 <- df.1 %>% ggplot(aes(FL1.H, FL3.H)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  p3 <- df.1 %>% ggplot(aes(SSC.H, FL1.H)) +
    scale_fill_gradientn(colours = myColor) + stat_binhex(bins = 200)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()

  #p4 <- df.1 %>% ggplot(aes(SSC.H, FL2.H)) + geom_point(size=0.1,color='darkgrey') + geom_density_2d(aes(colour = ..level..), contour_var = "count") +
  #  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_color_continuous(type = "viridis")


  #p.all <- cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
  p.all <- cowplot::plot_grid(p1, p2,p3, labels = "AUTO")

  # now add the title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      long_string,
      x = 0,
      hjust = 0,
      size=14
    ) +
    ggplot2::theme(
      plot.margin = margin(0, 0, 0, 7)
    )

  p.final <- cowplot::plot_grid(
    title, p.all,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )


  ggplot2::ggsave(paste0(sample_id," ", patient_id,'_','density-plots_',i,'.pdf'), p.final)

  #end of trycatch
  },
error = function(e)
  print("You can't calculate the log of a character"))

}








for(i in 1:n_experiments){
  df.1 <- fs[[i]]@exprs %>% data.frame()
  duration=gsub("\\.\\)","",gsub("Time \\(","",fs[[i]]@description$'$P8S'))
  timeticks=fs[[i]]@description$TIMETICKS
  patient_id=fs[[i]]@description$'PATIENT ID'
  sample_id=fs[[i]]@description$'SAMPLE ID'
  n_events <- nrow(df.1)
  long_string <- paste0(patient_id," ", sample_id," (n=",n_events," ,",duration,")")

  p1 <- df.1 %>% ggplot(aes(SSC.H, FL3.H)) + geom_hex(bins=100) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_fill_continuous(type = "viridis")

  p2 <- df.1 %>% ggplot(aes(FL1.H, FL3.H))  + geom_hex(bins=100) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_fill_continuous(type = "viridis")

  p3 <- df.1 %>% ggplot(aes(SSC.H, FL1.H)) + geom_hex(bins=100) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_fill_continuous(type = "viridis")


  p4 <- df.1 %>% ggplot(aes(SSC.H, FL2.H))  + geom_hex(bins=100) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_fill_continuous(type = "viridis")


  p.all <- cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
 # p.all <- cowplot::plot_grid(p1, p2,p3, labels = "AUTO")

  # now add the title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      long_string,
      x = 0,
      hjust = 0,
      size=14
    ) +
    ggplot2::theme(
      plot.margin = margin(0, 0, 0, 7)
    )

  p.final <- cowplot::plot_grid(
    title, p.all,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )


  ggplot2::ggsave(paste0('density_plots_',patient_id," ", sample_id,'.pdf'), p.final)

}












gridExtra::grid.arrange(p1,p2,ncol=2)


df.1 %>% select(SSC.H, FL3.H, FL1.H) %>%
  tidyr::pivot_longer(cols = c(FL3.H, FL1.H),names_to='fl')


#Events per time
df.1 %>%
  group_by(Time) %>%
  tally() %>%
  ggplot(aes(Time,n)) + geom_point(size=0.1,color='darkgrey') +
  theme_classic() + theme(aspect.ratio=0.2)+
  scale_color_continuous(type = "viridis") +
  ylab('n events')











df.1 %>% ggplot(aes(SSC.H, FL3.H)) + geom_point(size=0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()


df.1 %>% ggplot(aes(SSC.H, FL1.H)) + geom_point(size=0.3) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()


+




df.1 %>% ggplot(aes(SSC.H, FL3.H)) + geom_hex(bins = 200)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks() + scale_fill_continuous(type = "viridis")

geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis")




df.1 %>% ggplot(aes(SSC.H, FL3.H)) + geom_density_2d() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_classic() + theme(aspect.ratio=1)+ annotation_logticks()




library(scales) # to access break formatting functions
# x and y axis are transformed and formatted
p2 <- ggplot(Animals, aes(x = body, y = brain)) + geom_point() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()
# log-log plot without log tick marks
p2
# Show log tick marks
p2 + annotation_logticks()














fs_comp_clean <- flow_auto_qc(fs)


#matrix<-spillover(fs[[1]])$SPILL
#colnames(matrix)<-c("X.FITC.A.", "X.Pacific.Blue.A.", "X.Alexa.680.A.", "X.APC.A.", "X.PE.Cy7.A.", "X.PE.Cy55.A.", "X.PE.Tx.RD.A.", "X.PE.Green.laser.A.")
#fs_comp <-compensate(fs,matrix)

fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,c(4,6:12)]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#create the empty gating set
auto_gs<-GatingSet(fs_comp_clean_trans)

#cell gate
fs_data<- gs_pop_get_data(auto_gs)
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC.A","SSC.A")))
gs_pop_add(auto_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(auto_gs)
autoplot(auto_gs, x="FSC.A", y="SSC.A", "noneDebris_gate", bins=256)

#Singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC.A", "FSC.H")))
gs_pop_add(auto_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(auto_gs)
autoplot(auto_gs, x = 'FSC.A', y = 'FSC.H', "singlets", bins = 256)

#===============================









library(ggcyto)
dataDir <- system.file("extdata",package="flowWorkspaceData")
gs_orig <- load_gs(list.files(dataDir, pattern = "gs_bcell_auto",full = TRUE))


dataDir<-'/Users/sa01fd/DATA/FACSort/2110/21/fcs/'
gs_orig <- load_gs(list.files(dataDir, pattern = "fcs",full = TRUE))


gs <- gs_clone(gs_orig)

autoplot(gs, "CD3")
p <- autoplot(gs, "CD3", bins = 261)
p

#first sample
gs[[1]]@data

plot(gs[[1]])




library(dplyr)

df <- gs[[2]]@data[[1]]@exprs %>% data.frame()

df %>% ggplot2::ggplot(aes(FSC.A,SSC.A)) +
  geom_hex(bins=216)

df %>% ggplot2::ggplot(aes(SSC.A,CD38)) +
  geom_hex(bins=216)


CD38
Live



df <- df[,c('FSC.A','SSC.A','Live','CD38')]

df <- df.1[,c('SSC.H','FL1.H','FL2.H' )]

df <- cbind(metadata, data[metadata$name, ])


umap_rec <- recipes::recipe(~.,data=df) %>%
  #update_role(name, new_role = 'id' ) %>% #descide what it is...
  recipes::step_normalize(recipes::all_predictors()) %>% #  center and scale
  embed::step_umap(recipes::all_predictors())

umap_prep <- recipes::prep(umap_rec)
umap.out <-  recipes::juice(umap_prep)

umap.kmer.data <- umap.out %>%
  dplyr::select(umap_1, umap_2) %>%
  data.frame()


umap.kmer.data %>%
  ggplot(aes(umap_1,umap_2))+geom_hex(bins=60)





rownames(umap.kmer.data) <- umap.out$name

message(':: running dbscan')
res <- dbscan::dbscan(umap.kmer.data, eps = 1, minPts = 20)
n_clusters <- length(unique(res$cluster))
message(paste0('-- ', n_clusters, ' clusters detected'))

umap.kmer.data$cluster <- res$cluster
umap.kmer.plot <- umap.kmer.data %>% tibble::rownames_to_column(var='name')

p.umap.dbscan <- umap.kmer.data %>%
  ggplot(aes(umap_1,umap_2,color=as.factor(cluster)))+
  ggplot2::geom_hline(yintercept = 0, size = 0.25, colour = "#bdbdbd") +
  ggplot2::geom_vline(xintercept = 0, size = 0.25, colour = "#bdbdbd") +
  ggplot2::xlab(paste("UMAP 1",sep="")) +
  ggplot2::ylab(paste("UMAP 2",sep="")) +
  geom_point() +
  ggplot2::ggtitle(paste0('UMAP - CLR transformed 5-mer Binned dbscan')) +
  ggplot2::scale_color_manual(name='DBSCAN',values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_clusters+1))+
  #fdb_style()+
  theme(plot.title = element_text(size = 9))

p.umap.dbscan





