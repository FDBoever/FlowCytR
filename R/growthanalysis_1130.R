library(ggplot2)
library(dplyr)
library(nlraa)

df <- read.delim('~/DATA/FACSort/stat_out_exp1130.txt')
head(df)

conversionfactor = 87.719

df <- df %>%
  mutate(cellsml = Events * conversionfactor) %>%
  mutate(culture = gsub("(.+?)(\\_.*)", "\\1", Sample.ID))


df %>% dplyr::filter(grepl('KAT',sample_id)) %>%
  dplyr::filter(Label=='R3') %>%
  ggplot2::ggplot(ggplot2::aes(time,cellsml/1000,
                               group=culture,
                               color=culture)) +
  ggplot2::geom_point()+
  theme_classic()+
  theme(aspect.ratio=1)


df_av <- df %>% group_by(culture, time,Label) %>% summarize(mean=mean(cellsml), sd=sd(cellsml))

df_av %>%
  dplyr::filter(Label=='R3') %>%
  ggplot2::ggplot(ggplot2::aes(time,mean/1000, group=culture,color=culture)) +
  ggplot2::geom_point()+
  ggplot2::geom_line()+
  ggplot2::theme_classic()+
  ggplot2::theme(aspect.ratio=1)+
  facet_wrap(~culture)

predict<-data.frame(x,y)

df_av %>%
  dplyr::filter(grepl('SL48',culture)) %>%
  dplyr::filter(Label=='R1') %>%
  ggplot2::ggplot(ggplot2::aes(time,mean/1000, group=culture)) +
  ggplot2::geom_point()+
  #ggplot2::geom_line()+
  ggplot2::theme_classic()+
  ggplot2::geom_point(data=df %>%
                        dplyr::filter(grepl('SL48',sample_id)) %>%
                        dplyr::filter(Label=='R1'),
                      aes(time,cellsml/1000, group=culture),shape=21)+
  ggplot2::theme(aspect.ratio=1)


df %>%
  dplyr::filter(Label=='R5') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml/1000, group=culture,color=culture)) +
  ggplot2::geom_point()+
  ggplot2::theme_classic()+
  ggplot2::theme(aspect.ratio=1)+
  facet_wrap(~culture)

  ggplot2::scale_color_manual(values=c('black','grey'))

df %>% dplyr::filter(grepl('BA4C3',sample_id)) %>%
  dplyr::filter(Label=='R1') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml)) +
  ggplot2::geom_point()


df %>% dplyr::filter(grepl('SL48',sample_id)) %>%
  dplyr::filter(Label=='R1') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml/1000)) +
  ggplot2::geom_point()+
  ggplot2::theme_classic()+
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::scale_color_manual(values=c('black','grey'))


df %>% dplyr::filter(grepl('OSN25E8',sample_id)) %>%
  dplyr::filter(Label=='R3') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml/1000, group=culture,color=culture)) +
  ggplot2::geom_point()+
  ggplot2::theme_classic()+
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::scale_color_manual(values=c('black','grey'))

df %>% dplyr::filter(grepl('SWF3D3',Sample.ID)) %>%
  dplyr::filter(Label=='R3') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml/1000, group=culture,color=culture)) +
  ggplot2::geom_point()+
  ggplot2::theme_classic()+
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::scale_color_manual(values=c('black','grey'))


df %>% dplyr::filter(grepl('SWF3B3',sample_id)) %>%
  dplyr::filter(Label=='R3') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml/1000, group=culture,color=culture)) +
  ggplot2::geom_point()+
  ggplot2::theme_classic()+
  ggplot2::theme(aspect.ratio=1)+
  ggplot2::scale_color_manual(values=c('black','grey'))

df %>% dplyr::filter(grepl('BA4C3',sample_id)) %>%
  dplyr::filter(Label=='R1') %>% ggplot2::ggplot(ggplot2::aes(time,cellsml/1000, group=culture,color=culture)) +
  ggplot2::geom_point()+theme_classic()+theme(aspect.ratio=1)+
  ggplot2::scale_color_manual(values=c('black','grey'))




analysis.grid <- rbind(
  c("SL48","R1","R2"),
  c("BA4C3","R1","R2"),
  c("SWF2B3nephr","R3","R2"),
  c("OSN25E8","R3","R2"),
  c("SWF3B3micr","R5","R6"),
  c("SWF3D3mam","R3","R2"),
  c("SWF2KAT","R10","R3"),
  c("SWF2KATcoll","R10","R3")) %>% data.frame()

colnames(analysis.grid) = c('culture','predator','prey')

for(i in 1:nrow(analysis.grid)){
  errchecker='FALSE'

  message(i)
  sel.culture <- analysis.grid$culture[i]
  sel.gate <- analysis.grid$predator[i]

  df.growth <- df %>%
    dplyr::filter(culture==sel.culture) %>%
    dplyr::filter(Label==sel.gate) %>%
    dplyr::select(time, cellsml)

  fit <- tryCatch(
    expr = {
      nls(cellsml ~ SSlogis(time, Asym, xmid, scal), df.growth)
    },
    error = function(e){
      message('forced exponental fit')
      fit <- nls(cellsml ~ SSexpf(time, a, c), data = df.growth)
      return(fit)
    }
  )

  p.pred <- ggplot(df.growth, aes(time, cellsml)) +
    geom_point(shape=21) +
    geom_function(
      fun = function(x) {
        predict(fit, newdata = data.frame(time = x))
      },
      colour = "grey",
    ) +
    theme_classic() + theme(aspect.ratio=1)+
    ggplot2::ggtitle(paste(sel.culture, sel.gate))+
    ggplot2::xlab('Time (h)') +
    ggplot2::ylab('Cells mL-1')

  #-------------------- PREY
  sel.gate <- analysis.grid$prey[i]

  df.growth <- df %>%
    dplyr::filter(culture==sel.culture) %>%
    dplyr::filter(Label==sel.gate) %>%
    dplyr::select(time, cellsml)

  fit.prey <- tryCatch(
    expr = {
      nls(cellsml ~ SSlogis(time, Asym, xmid, scal), df.growth)
    },
    error = function(e){
      message('forced exponental fit')
      fit.prey <- nls(cellsml ~ SSexpf(time, a, c), data = df.growth)
      return(fit.prey)
    }
  )

  p.prey <- ggplot2::ggplot(df.growth, ggplot2::aes(time, cellsml)) +
    geom_point(shape=21) +
    geom_function(
      fun = function(x) {
        predict(fit.prey, newdata = data.frame(time = x))
      },
      colour = "grey",
    ) +
    theme_classic() + theme(aspect.ratio=1)+
    ggplot2::ggtitle(paste(sel.culture, sel.gate))+
    ggplot2::xlab('Time (h)') +
    ggplot2::ylab('Cells mL-1')

  p.comb <- gridExtra::grid.arrange(p.pred, p.prey, nrow=2)

  ggplot2::ggsave(paste0('~/DATA/FACSort','/',i,'_testing','.pdf'), p.comb,height = 10, width = 10)

  }



