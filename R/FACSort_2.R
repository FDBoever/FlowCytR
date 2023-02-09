
df <- read.delim('~/DATA/FACSort/FACSort_beads_210921.txt')
df <- df %>% mutate(syringe = ifelse(uLbeads == 25, 'mike','fdb'))
df <-  df %>% mutate(time_min = time /60)

analysis.name <- '210921_bds'

df <- read.delim('~/DATA/FACSort/FACSort_beads_220921.txt')
df <-  df %>% mutate(time_min = time /60)
analysis.name <- '220921_bds'

#===================================================================
# ASSESSING FLOW RATE

# flowrate calculated by Mike by measuring the inside of the syringe
df.flowrate <- cbind('NR'=c(2,4,8,10,12,16,20),'flowrate'=c(9.5,19,38,47.5,57,76.11,95.14)) %>% data.frame

# Extracting the slope to predict unknown values
# estimate the slope of the regression and use that to transform NR to more interpretable flow rates
lm.f <- lm(flowrate ~ NR,df.flowrate)
conversion_factor <- unname(lm.f$coefficients[2])

#Checkng linearity
# this is a linear relationship
p <- df.flowrate %>%
  ggplot(aes(NR,flowrate/1000)) +
  geom_point()+
  fdb_style() + geom_smooth(method='lm',fullrange=TRUE,se=FALSE) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab(bquote('NR')) +
  ylab(bquote('Flow Rate (mL'~min^-1*')'))+
  ggtitle(paste0('slope=',conversion_factor))
p

ggplot2::ggsave(filename=paste0('220921_flowrate_vs_NR','.pdf'),
                plot=p,
                width = 4,
                height = 4,
                unit='in')


#Lets assume this is linear, with a slope stored in conversion_factor
df.long <- df %>%
  dplyr::mutate(flowrate = NR * conversion_factor) %>%
  tidyr::pivot_longer(cols = c('r1','r2','r3'))


#==============================================================================
# particles per minute based on raw counts (howmany counted in a minute)
#---------------------------------------------------------
# in the third syringe, I played with flow rates
# I dont claculate the slope to get an estimate of beads per minute, but litteraly counted in per minute (triplicates)
# allowing me to do more measurments up the range into higher flow rates


df.model <- df.long %>% filter(syringe=='mike3')
#dilution facor
dil.fac <- 4150/150

#df.model <- df.long %>% filter(syringe=='fdb') %>% filter(time==61)
#dil.fac <- 4200/200

# since the machine is counting 61 seconds instead of 60 we need to convert this back
df.model <- df.model %>%
  dplyr::mutate(value = (value/61)*60)


model = lm(value ~ flowrate, df.model)
slope <- unname(model$coefficients[2])


#--- initial concentration (per 1000 beads)
conc.init <- slope * dil.fac * 1000


df.model$predicted = predict(
  object = model,
  newdata = df.model
)


# predict 95% confidence interval
df.model$CI = predict(
  object = model,
  newdata = df.model,
  se.fit = TRUE
)$se.fit * qnorm(1 - (1-0.95)/2)



p.2<-df.model %>% #filter(syringe=='mike3') %>%
  ggplot(aes(flowrate/1000,value/1000)) +
  fdb_style() + geom_smooth(method='lm',fullrange=TRUE,se=FALSE,size=.5,color='black',fill='lightgrey') +
  geom_point()+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab(bquote('Flow Rate (mL'~min^-1*')')) +
  ylab(bquote(''~x10^-3*' beads'~min^-1*'')) +
  geom_line(aes(y = predicted/1000 + CI/1000), color = "grey",linetype='dashed',size=.5)+ #+ # upper
  geom_line(aes(y = predicted/1000 - CI/1000), color = "grey",linetype='dashed',size=.5) +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  ggtitle(paste0('raw per minute measurments (',conc.init * 10^-6 ,')'))+
  theme(plot.title = element_text(size = 6))
p.2

ggplot2::ggsave(filename=paste0(analysis.name,'220921_raw_per_minute','.pdf'),
                plot=p.2,
                width = 5,
                height = 5,
                unit='in')


df.model %>% #filter(syringe=='mike3') %>%
  ggplot(aes(flowrate/1000,value/1000)) +
  fdb_style() + geom_smooth(method='lm',fullrange=TRUE,se=FALSE,size=.5,color='black',fill='lightgrey') +
  geom_point()+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab(bquote('Flow Rate (mL'~min^-1*')')) +
  ylab(bquote(''~x10^-3*' beads'~min^-1*'')) +
  geom_line(aes(y = predicted/1000 + CI/1000), color = "grey",linetype='dashed',size=.5)+ #+ # upper
  geom_line(aes(y = predicted/1000 - CI/1000), color = "grey",linetype='dashed',size=.5) +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  ggtitle(paste0('raw per minute measurments (',conc.init * 10^-6 ,')'))+
  theme(plot.title = element_text(size = 6))






#==============================================================================
# Caculate based on slopes of counts ~ time
#≠======================================

df.long.original <- df.long
df.long <- df.long.original %>% filter(syringe %in% c('mike1','mike2'))



#overall quick assesment of the relationships
my.formula <- y ~ x
p.3 = df.long %>% ggplot(aes(time_min,value,color=as.factor(NR)))+
  geom_point(aes()) +
  geom_smooth(method='lm', formula = my.formula,fullrange=TRUE,se=FALSE)+
  ylab('number of events') +
  facet_wrap(~syringe,scales='free') +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  fdb_style() +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999")) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab('elapsed time (min)')

p.3

ggplot2::ggsave(filename=paste0(analysis.name,'_events_per_time','.pdf'),
                plot=p.3,
                width = 10,
                height = 10,
                unit='in')




p.4 <- df.long %>%
  #filter(syringe == 'mike2') %>%
  ggplot(aes(time_min,value,color=as.factor(NR)))+
  geom_point(aes()) +
  geom_smooth(method='lm', formula = my.formula,fullrange=TRUE,se=FALSE)+
  ylab('number of events') +
  facet_wrap(~syringe+NR,scales='free') +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  fdb_style() +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999")) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab('elapsed time (min)')

ggplot2::ggsave(filename=paste0(analysis.name,'_events_per_time_split','.pdf'),
                plot=p.4,
                width = 11,
                height = 11,
                unit='in')



#=======================================
# linear models per 'replicate'
# !!!! warning, not super convinced about this, as the replicates are not related - basically randomly combined, pitch this to mike
# !!!! Alternativelly I need to randomly pull, and iterate a lot, similar to tidymodels bootrstrap sampling??


library(broom)
library(tidyr)

lm.out <- df.long %>%
  filter(!is.na(value)) %>%
  group_by(name,NR,uLbeads,mLMQ,syringe,flowrate) %>%
  do(lm = tidy(lm(value ~ time_min, data = .))) %>%
  unnest(lm)

#similarly, we dont do this right????
p.5 <- lm.out %>%
  filter(term=='time_min') %>%
  ggplot(aes(flowrate/1000,estimate/1000,color=name)) +
  geom_point() +
  fdb_style() +
  geom_smooth(method='lm',fullrange=TRUE,se=FALSE) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  facet_wrap(~syringe,scales='free') +
  xlab(bquote('Flow Rate (mL'~min^-1*')')) +
  ylab(bquote(''~x10^-3*' (beads'~min^-1*')'))

ggplot2::ggsave(filename=paste0(analysis.name,'_beads_mL_per_flowrate_split','.pdf'),
                plot=p.5,
                width = 8,
                height = 8,
                unit='in')


p.6 <- lm.out %>%
  filter(term=='time_min') %>%
  ggplot(aes(flowrate/1000,estimate/1000)) +
  geom_point() +
  fdb_style() +
  geom_smooth(method='lm',fullrange=TRUE,se=FALSE,size=.5,color='black',fill='lightgrey') +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  facet_wrap(~syringe,scales='free') +
  xlab(bquote('Flow Rate (mL'~min^-1*')')) +
  ylab(bquote(''~x10^-3*' (beads'~min^-1*')'))+
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE)

ggplot2::ggsave(filename=paste0(analysis.name,'_beads_mL_per_flowrate','.pdf'),
                plot=p.6,
                width = 8,
                height = 8,
                unit='in')












#--------------------------------------

library(dplyr)
library(raster)

#MIND YOU ALL OVERLAYED NOW
p.7 <- df.long %>%
  group_by(NR,uLbeads,mLMQ,syringe,flowrate, time, time_min) %>%
  summarise(mean = mean(value) ,
            sd = sd(value),
            cv= cv(value),
            n=n()) %>%
#  filter(cv<40) %>%
  dplyr::filter(!is.na(mean)) %>%
  ggplot(aes(mean,cv,fill=syringe)) +
    geom_point(shape=21) +
    fdb_style() +
    xlab('mean counted particles') +
    ylab('Coefficient of variance (%)')+
  scale_fill_manual(values=c('black','grey'))


ggplot2::ggsave(filename=paste0(analysis.name,'_CV_vs_count','.pdf'),
                plot=p.7,
                width = 6,
                height = 6,
                unit='in')

















p

ggplot2::ggsave(filename=paste0('220921_pipette_events_time','.pdf'),
                plot=p,
                width = 8,
                height = 8,
                unit='in')


p = df.long %>% ggplot(aes(NR,value,color=as.factor(time)))+
  geom_point(aes()) +
  geom_smooth(method='lm', formula = my.formula,fullrange=TRUE)+
  facet_wrap(~syringe,scales='free') +
  ylab('number of events') +
  stat_poly_eq(formula = my.formula ,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  fdb_style() +expand_limits(x = 0, y = 0) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999"))


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

