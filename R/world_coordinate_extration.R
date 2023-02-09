library('ggmap')

library(ggmap)
library(ggplot2)


df_world <- map_data("world")  %>% data.frame()

#df_world %>% data.frame() %>% ggplot(aes(x,y)) + geom_point(size=.2)

#df_world %>% ggplot(aes(long,lat)) + geom_point(size=.2)


dfcoast %>%
  mutate(long_c= ifelse(long<0,180+(180-abs(long)),long)) %>% ggplot(aes(long_c,lat)) + geom_point(size=.05)

c_df_coast <-
  dfcoast %>%
  mutate(long_c= ifelse(long<0,180+(180-abs(long)),long))




coast <- rnaturalearth::ne_coastline(scale='medium')
dfcoast <- ggplot(data=coast)

#data from
dfcoast$data



BC.df%>% ggplot(aes(long,lat)) + geom_point(size=.2)



library(ggmap)
pngMAP_df = get_map(location = c(142.5, -8.10, 143.1, -7.9), source = "osm", zoom = 12)
ggmap(pngMAP_df)


BC.df

library(raster)
library(maps)
library(mapdata)
library(maptools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(ggsn)
library(tidyverse)
library(here)


BC.shp <- readOGR('~/Downloads/hakai_guide_to_r-master/data/2_Shapefile/COAST_TEST2.shp')
BC.shp <- readOGR('/Users/sa01fd/Downloads/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp')

BC.df <- fortify(BC.shp)

write.table(BC.df %>% dplyr::select(lat,long,order),'~/DATA/coastline_gshhg2.3.7_longlat_hires.txt',quote = FALSE,sep='\t',row.names = FALSE)





### if you want to make a larger Central coast map, just change the extent selected
CCoast <- extent(-128.48, -127.9, 51.5, 52.1)
# crop the map to the new extent
CC.shp2 <- crop(BC.shp,CCoast)
# fortify
CC.df <- fortify(CC.shp2)

# Jenn graph
fig1 <- ggplot()+ theme_bw()+
  geom_polygon(data= CC.df, aes(x=long,y=lat,group= group),
               colour= "black", size=0.1, fill='grey85')+
  coord_cartesian(xlim = c(-128.5, -127.95), ylim=c(51.63, 52.05)) +
  #geom_point(data=EXPTsites, aes(x=long, y=lat, shape=otter), size=3.3, colour="blue", stroke=1.3)+  #add this to plot site locations
  #scale_shape_manual(values=c(21,24))+       #this makes different shapes for otter "yes" and otter "no" sites
  scale_x_continuous(breaks=c(-128.4, -128.2, -128.0))+
  scalebar(CC.df, dist = 5, st.size=3.5, height=0.014, dd2km = TRUE, model = 'WGS84', anchor = c(x = -128.33, y = 51.64))+
  north(data = CC.df, scale = 0.07, symbol = 3, anchor= c(x = -128.465, y = 52.056)) +
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10),
        legend.position = "none"); fig1

write.table(df_world %>% select(long,lat,region),'~/DATA/borders_longlat.txt',quote = FALSE,sep='\t',row.names = FALSE)
write.table(dfcoast$data %>% select(long,lat),'~/DATA/coastline_longlat.txt',quote = FALSE,sep='\t',row.names = FALSE)


write.table(c_df_coast,'~/DATA/coastline_longlat_glued.txt',quote = FALSE,sep='\t',row.names = FALSE)


nrow(dfcoast$data)
nrow(df_world)





coast <- rnaturalearth::ne_coastline(scale=scale)
dfcoast <- ggplot(data=coast)

#data from
dfcoast$data



coast <- rnaturalearth::ne_coastline(scale='medium')
dfcoast <- ggplot(data=coast)
dfcoast <- dfcoast$data %>% select(long,lat)

rivers <- rnaturalearth::ne_download(scale='medium', type='rivers',category='physical')
dfrivers <- ggplot(data=rivers)

lakes <- rnaturalearth::ne_download(scale='medium', type='lakes',category='physical')
dfrivers <- ggplot(data=rivers)


write.table(dfcoast,'~/DATA/coastline_longlat.txt',quote = FALSE,sep='\t',row.names = FALSE)


coast <- rnaturalearth::ne_coastline(scale='large')
dfcoast <- ggplot(data=coast)
dfcoast <- dfcoast$data %>% select(long,lat)
write.table(dfcoast,'~/DATA/coastline_longlat_hires.txt',quote = FALSE,sep='\t',row.names = FALSE)



df_uk <- map_data("lakes")  %>% data.frame()

