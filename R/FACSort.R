df <- read.delim('~/DATA/FACSort/FACSort_beads_210921.txt')

df.flowrate <- cbind('NR'=c(2,4,8,10,12,16,20),'flowrate'=c(9.5,19,38,47.5,57,76.11,95.14)) %>% data.frame

df.long <- df %>% tidyr::pivot_longer(cols = c('r1','r2','r3'))
df.long

my.formula <- y ~ x
p = df.long %>%
  ggplot2::ggplot(ggplot2::aes(time,value,color=as.factor(NR)))+
  ggplot2::geom_point(ggplot2::aes()) +
  ggplot2::geom_smooth(method='lm', formula = my.formula,fullrange=TRUE)+
  ggplot2::ylab('number of events') +
  ggplot2::facet_wrap(~uLbeads,scales='free') +
  ggpmisc::stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  expand_limits(x = 0, y = 0) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999"))#+
  #fdb_style()



ggplot2::ggsave(filename=paste0('events_time','.pdf'),
                plot=p,
                width = 8,
                height = 8,
                unit='in')


p = df.long %>% ggplot(aes(NR,value,color=as.factor(time)))+
  geom_point(aes()) +
  geom_smooth(method='lm', formula = my.formula,fullrange=TRUE)+
  facet_wrap(~uLbeads,scales='free') +
  ylab('number of events') +
  ggpmisc::stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  expand_limits(x = 0, y = 0) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999")) #+ fdb_style()


ggplot2::ggsave(filename=paste0('events_NR','.pdf'),
                plot=p,
                width = 8,
                height = 8,
                unit='in')





RColorBrewer::brewer.pal(10,name='Set1')

my.formula <- y ~ x
df.long %>% filter(uLbeads == 200) %>%
  ggplot(aes(time,value,color=as.factor(NR)))+
  geom_point() +
  geom_smooth(method='lm', formula = my.formula,fullrange=TRUE)+
  facet_wrap(~NR,scales='free') +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  fdb_style() +expand_limits(x = 0, y = 0) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999"))



p = df.long %>% filter(uLbeads == 200) %>%
  ggplot(aes(time,value,color=as.factor(NR)))+
  geom_point() +
  geom_smooth(method='lm', formula = my.formula,fullrange=TRUE)+
  facet_wrap(~NR,scales='free') +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  fdb_style() +expand_limits(x = 0, y = 0) +
  xlab('number of events') +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999"))

ggplot2::ggsave(filename=paste0('facetted','.pdf'),
                plot=p,
                width = 12,
                height = 12,
                unit='in')


#-------------
p <- df.long %>% select(-name) %>%
  dplyr::group_by(mLMQ, uLbeads, NR, time) %>%
  dplyr::summarise(mean=mean(value), sd=sd(value)) %>%
  ggplot(aes(mean,sd, color= as.factor(NR))) +
  geom_point() +
  fdb_style() +   scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999")) +
  facet_wrap(~uLbeads, scales='free')

ggplot2::ggsave(filename=paste0('heteroskedasticity','.pdf'),
                plot=p,
                width = 11,
                height = 5,
                unit='in')


df.long  %>% filter(uLbeads == 200) %>% select(-name) %>%
  dplyr::group_by(mLMQ, uLbeads, NR, time) %>%
  dplyr::summarise(mean=mean(value), se=sd(value)/mean(value)) %>%
  ggplot(aes(mean,se, color= as.factor(NR))) +
  geom_point() +
  fdb_style()




overall.mean <- df.calculated   %>%
  filter(!is.na(events_per_mL_original)) %>%
  summarise(mean=mean(events_per_mL_original)) %>% pull()

overall.mean_plus_100 <- df.calculated  %>%
  filter(!is.na(events_per_mL_original)) %>%
  filter(value>1000) %>%
  summarise(mean=mean(events_per_mL_original)) %>% pull()


class_mean <- df.calculated %>%
  group_by(uLbeads) %>% filter(!is.na(events_per_mL_original)) %>%
  summarise(mean = mean(events_per_mL_original))



# CUTS


df.calculated <- df.long %>% left_join(df.flowrate, by ='NR') %>%
  mutate(volume = (flowrate/60)*time,
         events_per_uL_sample = value/volume,
         sample_volume = mLMQ*1000+uLbeads,
         dilution_factor = sample_volume/uLbeads,
         events_per_uL_original = events_per_uL_sample * dilution_factor,
         events_per_mL_original = events_per_uL_original * 1000)

df.calculated %>%
  ggplot(aes(NR,events_per_mL_original,group=NR)) + geom_boxplot(outlier.size=NULL) +geom_jitter() + facet_wrap(~uLbeads)


df.calculated %>%
  ggplot(aes(NR,events_per_mL_original,group=NR)) + geom_boxplot(outlier.size=NULL) +geom_jitter() + facet_wrap(~uLbeads)

df.calculated %>% mutate(bins = cut(value, c(0, 125,250,500,750,1000,1500,2000,2500,3000,3500,4000,4500, 5000,7500,10000,12500 , Inf))) %>%
  ggplot(aes(bins,events_per_mL_original,group=bins)) + geom_boxplot(outlier.size=NULL) +geom_jitter() + facet_wrap(~uLbeads) +
  geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df.calculated %>% mutate(bins = cut(value, c(0, 100,250,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000, Inf))) %>%
  ggplot(aes(bins,events_per_mL_original,group=bins)) + geom_boxplot(outlier.size=NULL) +geom_jitter() +
  geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed") +
  facet_wrap(~uLbeads)

 p <- df.calculated %>% mutate(bins = cut(value, c(0, 100,250,500,750,1000,2000,2500,5000,7500,10000, Inf))) %>%
  ggplot(aes(bins,events_per_mL_original,group=bins)) + geom_boxplot(outlier.size=NULL) +geom_jitter() + facet_wrap(~uLbeads) +
  geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed") +
   xlab('nr of measured events') +
   ylab('estimated events in original sample') +
  fdb_style(aspect.ratio=0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0('Binned_estimates_2','.pdf'),
                plot=p,
                width = 10,
                height = 5,
                unit='in')

 p <-  df.calculated %>% mutate(bins = cut(value, c(0, 1000,2500,5000,10000, Inf))) %>%
    ggplot(aes(bins,events_per_mL_original,group=bins)) + geom_boxplot(outlier.size=NULL) +
    geom_jitter() +
    xlab('nr of measured events') +
    ylab('estimated events in original sample') +
    facet_wrap(~uLbeads) +
    geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed") +
    fdb_style(aspect.ratio=0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  ggplot2::ggsave(filename=paste0('Binned_estimates_3','.pdf'),
                  plot=p,
                  width = 10,
                  height = 5,
                  unit='in')



#!!!!!!!!!!!!!!!
  p <- df.calculated %>%
    ggplot(aes(value,events_per_mL_original))  +
    geom_point(aes(color=NR)) +
    geom_smooth()+
    facet_wrap(~uLbeads) +
    xlab('nr of measured events') +
    ylab('estimated events in original sample') +
    geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed") + fdb_style(aspect.ratio=0.5)

  ggplot2::ggsave(filename=paste0('_events_estimate_value_','.pdf'),
                  plot=p,
                  width = 10,
                  height = 5,
                  unit='in')



  df.calculated %>% ggplot(aes(value,events_per_mL_original))  + geom_point(aes(color=NR)) + geom_smooth()+
    facet_wrap(~uLbeads,scales='free_x') +
    geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed")



  df.calculated %>% ggplot(aes(value,events_per_mL_original))  + geom_point() + geom_smooth()+ facet_grid(NR~uLbeads) +
    geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed")

  df.calculated %>% filter(uLbeads==200) %>% ggplot(aes(value,events_per_mL_original))  + geom_point() + geom_smooth()+
    facet_grid(NR~time,scales='free') +
    geom_hline(data = class_mean, aes(yintercept = mean), linetype = "dashed")



df.calculated %>%
  ggplot(aes(as.character(uLbeads),events_per_mL_original)) + geom_boxplot(outlier.size=NULL) +geom_jitter(width = .3)




# with mapped means

df.calculated %>%
  ggplot(aes(events_per_mL_original, color=as.factor(uLbeads),fill=as.factor(uLbeads))) +
 geom_histogram() +
  geom_vline(xintercept=overall.mean,color='blue') +
  geom_vline(xintercept=overall.mean_plus_100,color='red') +
  geom_vline(data = class_mean, aes(xintercept = mean), linetype = "dashed") +
  facet_wrap(~uLbeads,ncol=1)

p = df.calculated %>%
  ggplot(aes(events_per_mL_original, color=as.factor(uLbeads),fill=as.factor(uLbeads))) +
  geom_density() +
  geom_vline(xintercept=overall.mean,color='blue') +
  geom_vline(xintercept=overall.mean_plus_100,color='red') +
  geom_vline(data = class_mean, aes(xintercept = mean), linetype = "dashed") +
  facet_wrap(~uLbeads,ncol=1)

ggplot2::ggsave(filename=paste0('density_','.pdf'),
                plot=p,
                width = 6,
                height = 6,
                unit='in')


#2 levels

class_mean <- df.calculated %>%
  group_by(uLbeads,NR) %>%
  filter(!is.na(events_per_mL_original)) %>%
  summarise(mean = mean(events_per_mL_original))

p <- df.calculated %>%
  ggplot(aes(events_per_mL_original, color=as.factor(uLbeads),fill=as.factor(uLbeads))) +
  geom_histogram() +
  geom_vline(data = class_mean, aes(xintercept = mean), linetype = "dashed") +
  geom_vline(xintercept=overall.mean,color='blue') +
  geom_vline(xintercept=overall.mean_plus_100,color='red') +
  facet_grid(NR~uLbeads)

ggplot2::ggsave(filename=paste0('density__per_NR','.pdf'),
                plot=p,
                width = 6,
                height = 6,
                unit='in')


class_mean <- df.calculated %>%
  group_by(Class, Year) %>%
  summarise(Mean = mean(events_per_mL_original))



ggplot(class_comparison, aes(x = Income, fill = as.factor(Year))) +
  geom_histogram() +
  facet_grid(Year ~ Class, scales = "free")





df.calculated %>% mutate(bins = cut(value, c(0, 100,250,500,750,1000,2000,2500,5000,7500,10000, Inf))) %>%
  ggplot(aes(bins,events_per_mL_original,group=bins)) + geom_boxplot(outlier.size=NULL) +geom_jitter() + facet_wrap(~uLbeads) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

