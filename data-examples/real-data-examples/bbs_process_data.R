################################################################################
### Prepare BBS data 
### Note that the bbsAssistant::bbs_obs function (line 26) has been depreciated at the time
### when the scripts were uploaded to GitHub. The BBS data used to fit models 
### are available in the Data folder. 
### The R script here provides reference about data manipulation.
################################################################################
rm(list = ls())

library(dplyr)
library(readr)
library(bbsAssistant)

outpath <- "/path/to/output/"


#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
routes <- read_csv("/path/to/load/bbs_routes.csv")
routes <- routes %>% arrange(StateNum, Route) %>% select(StateNum, Route, Latitude, Longitude)

#-------------------------------------------------------------------------------
# Process data
#-------------------------------------------------------------------------------
bbs <- bbsAssistant::bbs_obs %>% 
  filter(CountryNum == 840, Year == 2019) %>% 
  select(StateNum, Route, AOU, SpeciesTotal)

bbs_bird <- bbs %>% filter(AOU == "05930")  %>% arrange(StateNum, Route)

bbs_bird$StateNum <- as.character(bbs_bird$StateNum)
routes$StateNum <- as.character(as.integer(routes$StateNum))

bird <- left_join(bbs_bird, routes, by = c("StateNum", "Route"))
dat <- data.frame(lon = bird$Longitude, 
                  lat = bird$Latitude, 
                  obs = bird$SpeciesTotal)
bound <- -102
dat <- dat %>% arrange(lon, lat) %>% filter(lon > bound) %>% 
  group_by(lon, lat) %>% mutate(obs = median(obs))
dat <- dat[!duplicated(dat[, 1:2]),]

#-------------------------------------------------------------------------------
# Write CSV
#-------------------------------------------------------------------------------
write_csv(dat, file = "/path/to/save/bbs_dat.csv")

#-------------------------------------------------------------------------------
# Plot data
#-------------------------------------------------------------------------------
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

usa <- rnaturalearth::ne_states(returnclass = "sf", country = "United States of America")
bbox <- c(minlon = bound, maxlon = -66, minlat = 24, maxlat = 50)
base <- ggplot(data = usa) + geom_sf() +
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]),
           ylim = c(bbox["minlat"], bbox["maxlat"]),
           expand = FALSE)

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

plot_dat_save <- base + 
  geom_point(aes(lon, lat, color = obs, size = obs), data = dat) + theme_bw() + 
  scale_size_continuous(guide = "none") + 
  labs(color = "", x = "Longitude", y = "Latitude") + 
  scale_colour_gradientn(colours = myPalette(100), limits = c(0, max(dat[,3]))) + 
  theme(plot.margin = unit(c(0.5,0,0,0), "cm")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.height= unit(2.2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent")) 

setEPS()
postscript(paste0(outpath, "05930_counts", ".eps"), width = 7, height = 6, bg = "transparent")
print(plot_dat_save)
dev.off()