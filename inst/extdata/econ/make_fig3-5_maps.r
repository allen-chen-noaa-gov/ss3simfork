# library(raster)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(dplyr)
# library(haven)
library(sf)
library(geosphere)
library(ggforce)
library(ggpubr)

ymin <- 0
ymax <- 3
xmin <- -3
xmax <- 0
yin <- c(seq(ymin, ymax, by = 0.75))
horizlines <- data.frame(cbind(x = rep(c(xmin, xmax),
        length(yin)), 
    y = rep(yin, each = 2), 
    layer = rep(yin, each = 2),
    newx = rep(c(xmin, xmax), 
        length(yin))))
xin <- seq(xmin, xmax, by = 0.75) 
vertlines <- data.frame(cbind(x = rep(xin, each = 2), 
    y = rep(c(ymin, ymax), 
        length(xin)), 
    layer = seq(xmin, xmax, by = 0.75)+10, newx = seq(xmin, xmax, by = 0.75)))
bindbdnew <- rbind(horizlines, vertlines)
bindbdnew$x <- as.numeric(bindbdnew$x)
bindbdnew$y <- as.numeric(bindbdnew$y)
bindbdnew$newx <- as.numeric(bindbdnew$newx)
pts_sf <- bindbdnew %>% 
    st_as_sf(coords = c("newx","y")) %>% 
    st_set_crs(4326)
pts_sf_line <- pts_sf %>% group_by(layer) %>% 
    summarize(m = mean(x), do_union=FALSE) %>% 
    st_cast("LINESTRING")
sf_poly <- st_collection_extract(st_polygonize(st_union(pts_sf_line)))
csvout <- data.frame(st_coordinates(sf_poly))
csvout$L2 <- as.factor(csvout$L2)
# csvout <- csvout[!(csvout$L2 %in% c(2,9,10,13,14,15,16)),]
csvout <- csvout[!(csvout$L2 %in% c(2,5,9,13)),]
csvout$L2 <- droplevels(csvout$L2)
# levels(csvout$L2) <- c(9,3,6,5,8,4,7,2,1)
levels(csvout$L2) <- c(3,1,2,5,6,4,9,8,7)
csvout1 <- by(csvout[,c("X", "Y")], csvout$L2, centroid)
csvout[,c("Cent_lon","Cent_lat")] <- NA
for(i in names(csvout1)) {
    csvout[csvout$L2 == i,c("Cent_lon","Cent_lat")] <- t(matrix(
       rep(csvout1[[i]], dim(csvout[csvout$L2 == i,])[1]), 2))
}

otherdatfin <- readRDS(paste0(getwd(),
  "/inst/extdata/econ/mortup_9loc/otherdatfin_mortup_9loc.RData"))

temp <- data.frame(cbind(L2=otherdatfin$choicefin$V1,
    catch=otherdatfin$catchfin$V1))
catches <- ddply(temp, .(L2), summarize, Mean_catch = mean(catch))
catches$Mean_catch <- seq(0.75, 1.50, (1.50-0.75)/8)

test <- merge(csvout, catches, by = c("L2"), all.x=TRUE, all.y=TRUE)
choices <- data.frame(table(otherdatfin$choicefin))
names(choices)[names(choices) == 'V1'] <- 'L2'
names(choices)[names(choices) == 'Freq'] <- 'End'
test <- merge(test, choices, by = c("L2"), all.x=TRUE, all.y=TRUE)
starts <- data.frame(table(otherdatfin$startloc))
names(starts)[names(starts) == 'Var1'] <- 'L2'
names(starts)[names(starts) == 'Freq'] <- 'Start'
test <- merge(test, starts, by = c("L2"), all.x=TRUE, all.y=TRUE)

#equal
mapped <- ggplot() +
    geom_path(data = test, aes(x = X, y = Y, group = L2), size = 0.5) +
    # geom_point(aes(x = x, y = y, color = layer), size = 0.5) +
    # geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),
        # fill="grey", color="black", size=0.75) +
    geom_polygon(data = test,
      aes(x = X, y = Y, group = L2, fill = Mean_catch)) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    geom_text(data = test, aes(x = Cent_lon - .25, y = Cent_lat + .25,
        group = L2, label = Start,
        color = "C"), size = 4) +
    # geom_text(data=test, aes(x = Cent_lon-.25, y = Cent_lat+.15, group = L2,
    #     label = End,
    #     color="D"), size = 5) +
    geom_text(data = test, aes(x = Cent_lon + .25, y = Cent_lat + .25,
        group = L2,
        label = paste("Area \n", L2),
        color = "B"), size = 4, lineheight = .75) +
    geom_segment(aes(x = -1.875, y = 1.875, xend = -1.875, yend = 1.125),
        arrow = arrow(length = unit(0.75, "cm")), size = 2, color = "red",
        alpha = 0.95) +
    #arrow labels
    geom_text(aes(x = -1.825, y = 1.875),
        label = "1", size = 4, color = "red") +
    labs(fill = "Abundance") +
    scale_color_manual(values = c("Orange", "Green", "Black"),
    # labels = c("Location",
    #     "Number of observations \nthat began at this location",
    #     "Number of observations \nthat chose this location")) +
    labels = c("Area",
        "Number of observations")) +
    guides(color = guide_legend(title = NULL, order = 1)) +
    theme(text = element_text(size = 12), legend.position = "top")

base <- mapped

ggsave(mapped, file = paste0(getwd(),
  "/inst/extdata/econ/mortup_9loc.png"),
  width = 7.5, height = 4.21875, units = "in", dpi = 800)

otherdatfin <- readRDS(paste0(getwd(),
  "/inst/extdata/econ/mortup_9loc_singlestart/", 
  "otherdatfin_mortup_9loc_singlestart.RData"))

temp <- data.frame(cbind(L2=otherdatfin$choicefin$V1,
    catch=otherdatfin$catchfin$V1))
catches <- ddply(temp, .(L2), summarize, Mean_catch = mean(catch))
catches$Mean_catch <- seq(0.75, 1.50, (1.50-0.75)/8)

test <- merge(csvout, catches, by = c("L2"), all.x=TRUE, all.y=TRUE)
choices <- data.frame(table(otherdatfin$choicefin))
names(choices)[names(choices) == 'V1'] <- 'L2'
names(choices)[names(choices) == 'Freq'] <- 'End'
test <- merge(test, choices, by = c("L2"), all.x=TRUE, all.y=TRUE)
starts <- data.frame(table(otherdatfin$startloc))
names(starts)[names(starts) == 'Var1'] <- 'L2'
names(starts)[names(starts) == 'Freq'] <- 'Start'
test <- merge(test, starts, by = c("L2"), all.x=TRUE, all.y=TRUE)

#start1
mapped <- ggplot() +
    geom_path(data = test, aes(x = X, y = Y, group = L2), size = 0.5) +
    # geom_point(aes(x = x, y = y, color = layer), size = 0.5) +
    # geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),
        # fill="grey", color="black", size=0.75) +
    geom_polygon(data = test,
      aes(x = X, y = Y, group = L2, fill = Mean_catch)) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    geom_text(data = test, aes(x = Cent_lon - .25, y = Cent_lat + .25,
        group = L2, label = Start,
        color = "C"), size = 4) +
    # geom_text(data=test, aes(x = Cent_lon-.25, y = Cent_lat+.15, group = L2,
    #     label = End,
    #     color="D"), size = 5) +
    geom_text(data = test, aes(x = Cent_lon + .25, y = Cent_lat + .25,
        group = L2,
        label = paste("Area \n", L2),
        color = "B"), size = 4, lineheight = .75) +
    geom_segment(aes(x = -0.375, y = 0.375, xend = -0.375*.9, yend = 1.125*.9),
        arrow = arrow(length = unit(0.75, "cm")), size = 2, color = "red",
        alpha = 0.95) +
    geom_segment(aes(x = -0.375, y = 1.125, xend = -1.125*.9, yend = 1.875*.9),
        arrow = arrow(length = unit(0.75, "cm")), size = 2, color = "red",
        alpha = 0.95) +
    geom_arc(aes(x0 = -1.125, y0 = 1.875, r = 0.25,
               start = 0, end = -1.875 * pi),
           arrow = arrow(), size = 2, color = "red", alpha = 0.95) +
    geom_segment(aes(x = -1.125, y = 1.875,
        xend = -1.875*0.95, yend = 1.125*0.9),
        arrow = arrow(length = unit(0.75, "cm")), size = 2, color = "red",
        alpha = 0.95) +
    geom_segment(aes(x = -1.875, y = 1.125, xend = -1.875, yend = 1.875),
        arrow = arrow(length = unit(0.75, "cm")), size = 2, color = "red",
        alpha = 0.95) +
    #arrow labels
    geom_text(aes(x = -0.375*.9, y = 0.375),
        label = "1", size = 4, color = "red") +
    geom_text(aes(x = -0.375*.9, y = 1.125),
        label = "2", size = 4, color = "red") +
    geom_text(aes(x = -1.125*0.95, y = 1.875*1.1),
           label = "3", size = 4, color = "red") +
    geom_text(aes(x = -1.125*0.95, y = 1.875),
        label = "4", size = 4, color = "red") +
    geom_text(aes(x = -1.875, y = 1.125*0.9),
        label = "5", size = 4, color = "red") +
    labs(fill = "Abundance") +
    scale_color_manual(values = c("Orange", "Green", "Black"),
    # labels = c("Location",
    #     "Number of observations \nthat began at this location",
    #     "Number of observations \nthat chose this location")) +
    labels = c("Area",
        "Number of observations")) +
    guides(color = guide_legend(title = NULL, order = 1)) +
    theme(text = element_text(size = 12), legend.position = "top")

singlestart <- mapped

ggsave(mapped, file = paste0(getwd(),
  "/inst/extdata/econ/mortup_9loc_singlestart.png"),
  width = 7.5, height = 4.21875, units = "in", dpi = 800)

test <- ggarrange(base, singlestart, ncol = 2, nrow = 1, common.legend = TRUE)

ggsave(test, file = paste0(getwd(),
    "/inst/extdata/econ/fig3", ".png"),
    width = 7.5, height = 4.21875, units = "in", dpi = 800)

# #start1_few
# mapped <- ggplot() +
#     geom_path(data=test, aes(x = X, y = Y, group = L2), size = 0.5) +
#     # geom_point(aes(x = x, y = y, color = layer), size = 0.5) +
#     # geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),
#         # fill="grey", color="black", size=0.75) +
#     geom_polygon(data=test, aes(x = X, y = Y, group = L2, fill=Mean_catch)) +
#     scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
#     geom_text(data=test, aes(x = Cent_lon-.35, y = Cent_lat+.35, group = L2, 
#     label=Start, 
#         color="C"), size = 3) +    
#     geom_text(data=test, aes(x = Cent_lon-.35, y = Cent_lat+.25, group = L2, 
#     label=End, 
#         color="D"), size = 3) +
#     geom_text(data=test, aes(x = Cent_lon+.35, y = Cent_lat+.35, group = L2, 
#         label=paste("Area \n", L2), 
#         color="B"), size = 3) +        
#     geom_segment(aes(x = -0.5, y = 0.5, xend = -0.5, yend = 1.5*.9),
#         arrow = arrow(length = unit(0.75, "cm")), size=2, color="red", 
#         alpha=0.5) +
#     geom_segment(aes(x = -0.5, y = 1.5, xend = -0.5, yend = 2.5),
#         arrow = arrow(length = unit(0.75, "cm")), size=2, color="red", 
#         alpha=0.5) +  
#     #arrow labels    
#     geom_text(aes(x = -0.4, y = 0.5),
#         label="1", size=2.5, color="red") +
#     geom_text(aes(x = -0.4, y = 1.5),
#         label="2", size=2.5, color="red") +   
#     labs(fill='Abundance') +    
#     scale_color_manual(values=c("Orange", "Green", "Black"), 
#     labels = c("Location", 
#         "Number of observations \nthat began at this location",
#         "Number of observations \nthat chose this location")) +
#     guides(color=guide_legend(title=NULL, order=1)) +
# 	theme(text = element_text(size=10))
#     # , legend.spacing.x = unit(1.0, 'cm')
# ggsave(mapped, file=paste0(nicedir, "\\start1_fewobs", ".png"), 
#     width = 7.5, height = 4.21875, units = "in", dpi=800)

ymin <- 0
ymax <- 4.5
xmin <- -4.5
xmax <- 0
yin <- c(seq(ymin, ymax, by = 0.75))
horizlines <- data.frame(cbind(x = rep(c(xmin, xmax),
        length(yin)), 
    y = rep(yin, each = 2), 
    layer = rep(yin, each = 2),
    newx = rep(c(xmin, xmax), 
        length(yin))))
xin <- seq(xmin, xmax, by = 0.75) 
vertlines <- data.frame(cbind(x = rep(xin, each = 2), 
    y = rep(c(ymin, ymax), 
        length(xin)), 
    layer = seq(xmin, xmax, by = 0.75)-10, newx = seq(xmin, xmax, by = 0.75)))
bindbdnew <- rbind(horizlines, vertlines)
bindbdnew$x <- as.numeric(bindbdnew$x)
bindbdnew$y <- as.numeric(bindbdnew$y)
bindbdnew$newx <- as.numeric(bindbdnew$newx)
pts_sf <- bindbdnew %>% 
    st_as_sf(coords = c("newx","y")) %>% 
    st_set_crs(4326)
pts_sf_line <- pts_sf %>% group_by(layer) %>% 
    summarize(m = mean(x), do_union=FALSE) %>% 
    st_cast("LINESTRING")
sf_poly <- st_collection_extract(st_polygonize(st_union(pts_sf_line)))
csvout <- data.frame(st_coordinates(sf_poly))
csvout$L2 <- as.factor(csvout$L2)
csvout <- csvout[(csvout$L2 %in% c(27,23,19,18,24,29,3,4)),]
# csvout <- csvout[(csvout$L2 %in% c(21,20,19,15,16,11,17,7)),]
# csvout <- csvout[(csvout$L2 %in% c(21,20,19,15,13,11,5,6)),]
# csvout <- csvout[!(csvout$L2 %in% c(2,13,14,15,16)),]
csvout$L2 <- droplevels(csvout$L2)
# levels(csvout$L2) <- c(8,6,5,4,7,3,2,1)
levels(csvout$L2) <- c(7,8,6,3,2,5,1,4)
csvout1 <- by(csvout[,c("X", "Y")], csvout$L2, centroid)
csvout[,c("Cent_lon","Cent_lat")] <- NA
for(i in names(csvout1)) {
    csvout[csvout$L2 == i,c("Cent_lon","Cent_lat")] <- t(matrix(
       rep(csvout1[[i]], dim(csvout[csvout$L2 == i,])[1]), 2))
}

otherdatfin <- readRDS(paste0(getwd(),
  "/inst/extdata/econ/mortup_8loc_random_20/",
  "otherdatfin_mortup_8loc_random_20_1.RData"))

csvout$id <- seq_len(nrow(csvout))
temp <- data.frame(cbind(L2=otherdatfin$choicefin$V1,
    catch=otherdatfin$catchfin$V1))
catches <- ddply(temp, .(L2), summarize, Mean_catch = mean(catch))
catches <- rbind(catches, c(7,0))
catches$Mean_catch <- otherdatfin$betavar
test <- merge(csvout, catches, by = c("L2"), all.x=TRUE, all.y=TRUE)
choices <- data.frame(table(otherdatfin$choicefin))
names(choices)[names(choices) == 'V1'] <- 'L2'
names(choices)[names(choices) == 'Freq'] <- 'End'
test <- merge(test, choices, by = c("L2"), all.x=TRUE, all.y=TRUE)
starts <- data.frame(table(otherdatfin$startloc))
names(starts)[names(starts) == 'Var1'] <- 'L2'
names(starts)[names(starts) == 'Freq'] <- 'Start'
test <- merge(test, starts, by = c("L2"), all.x=TRUE, all.y=TRUE)
test[is.na(test)] <- 0
test <- test[order(test$id), ] 

#random 1
mapped <- ggplot() +
    geom_path(data = test, aes(x = X, y = Y, group = L2), size = 0.5) +
    # geom_point(aes(x = x, y = y, color = layer), size = 0.5) +
    # geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),
        # fill="grey", color="black", size=0.75) +
    geom_polygon(data = test,
      aes(x = X, y = Y, group = L2, fill = Mean_catch)) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    geom_text(data = test, aes(x = Cent_lon - .25, y = Cent_lat + .25,
        group = L2, label = Start,
        color = "C"), size = 2.3) +
    # geom_text(data=test, aes(x = Cent_lon-.25, y = Cent_lat+.15, group = L2,
    #     label = End,
    #     color="D"), size = 5) +
    geom_text(data = test, aes(x = Cent_lon + .25, y = Cent_lat + .25,
        group = L2,
        label = paste("Area \n", L2),
        color = "B"), size = 2.3, lineheight = .75) +
    #arrows        
    geom_segment(aes(x = -0.375, y = 0.375, 
        xend = -0.375*0.9, yend = 1.125*0.9),
        arrow = arrow(length = unit(0.75, "cm")), size=1, color="red", 
        alpha=0.95) +
    geom_arc(aes(x0 = -0.375, y0 = 1.125, r = 0.25,
               start = 0, end = -1.125*1.5 * pi),
           arrow = arrow(), size=1, color="red", alpha=0.95) + 
    #arrow labels    
    geom_text(aes(x = -0.375*0.85, y = 0.375),
        label="1", size=2.3, color="red") +
    geom_text(aes(x = -0.375*0.9, y = 1.125*1.1),
        label="2", size=2.3, color="red") +  
    geom_text(aes(x = -3, y = 0.75),
        label="Habitat", size=5, color="Orange") +
    geom_segment(aes(x = -3.5625, y = 0.9375, xend = -2.4375, yend = 0.9375),
        size=3, color="orange") +
    geom_segment(aes(x = -3.5625, y = 1.125, xend = -2.4375, yend = 1.125),
        size=3, color="orange") +
    geom_segment(aes(x = -3.5625, y = 0.5625, xend = -2.4375, yend = 0.5625),
        size=3, color="orange") +
    geom_segment(aes(x = -3.5625, y = 0.375, xend = -2.4375, yend = 0.375),
        size=3, color="orange") +
    labs(fill = "Abundance") +
    scale_color_manual(values = c("Orange", "Green", "Black"),
    # labels = c("Location",
    #     "Number of observations \nthat began at this location",
    #     "Number of observations \nthat chose this location")) +
    labels = c("Area",
        "Number of observations")) +
    guides(color = guide_legend(title = NULL, order = 1), keyheight = 0) +
    theme(text = element_text(size = 10), legend.position = "top",
    legend.box="vertical", legend.key.size = unit(0.4, 'cm'),
    legend.text = element_text(size=8), legend.title = element_text(size=8),
    legend.spacing.y = unit(0.0, 'cm'))

random1 <- mapped

ggsave(mapped, file = paste0(getwd(),
    "/inst/extdata/econ/mortup_8loc_random_20_1.png"),
    width = 3.2, height = 3.6, units = "in", dpi=800)

otherdatfin <- readRDS(paste0(getwd(),
  "/inst/extdata/econ/mortup_8loc_random_20/",
  "otherdatfin_mortup_8loc_random_20_2.RData"))

csvout$id <- seq_len(nrow(csvout))
temp <- data.frame(cbind(L2=otherdatfin$choicefin$V1,
    catch=otherdatfin$catchfin$V1))
catches <- ddply(temp, .(L2), summarize, Mean_catch = mean(catch))
catches$Mean_catch <- otherdatfin$betavar
test <- merge(csvout, catches, by = c("L2"), all.x=TRUE, all.y=TRUE)
choices <- data.frame(table(otherdatfin$choicefin))
names(choices)[names(choices) == 'V1'] <- 'L2'
names(choices)[names(choices) == 'Freq'] <- 'End'
test <- merge(test, choices, by = c("L2"), all.x=TRUE, all.y=TRUE)
starts <- data.frame(table(otherdatfin$startloc))
names(starts)[names(starts) == 'Var1'] <- 'L2'
names(starts)[names(starts) == 'Freq'] <- 'Start'
test <- merge(test, starts, by = c("L2"), all.x=TRUE, all.y=TRUE)
test[is.na(test)] <- 0
test <- test[order(test$id), ] 

#random 2
mapped <- ggplot() +
    geom_path(data = test, aes(x = X, y = Y, group = L2), size = 0.5) +
    # geom_point(aes(x = x, y = y, color = layer), size = 0.5) +
    # geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),
        # fill="grey", color="black", size=0.75) +
    geom_polygon(data = test,
      aes(x = X, y = Y, group = L2, fill = Mean_catch)) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    geom_text(data = test, aes(x = Cent_lon - .25, y = Cent_lat + .25,
        group = L2, label = Start,
        color = "C"), size = 2.3) +
    # geom_text(data=test, aes(x = Cent_lon-.25, y = Cent_lat+.15, group = L2,
    #     label = End,
    #     color="D"), size = 5) +
    geom_text(data = test, aes(x = Cent_lon + .20, y = Cent_lat + .25,
        group = L2,
        label = paste("Area \n", L2),
        color = "B"), size = 2.3, lineheight = .75) +
    #arrows        
    geom_arc(aes(x0 = -0.375, y0 = 0.375, r = 0.10,
          start = 0, end = -1.75 * pi),
          arrow = arrow(length = unit(0.10, "inches")), size=1.0,
          color="red", alpha=0.95) + 
    geom_arc(aes(x0 = -0.375, y0 = 0.375, r = 0.20,
          start = 0, end = -1.75 * pi),
          arrow = arrow(length = unit(0.10, "inches")), size=1.0,
          color="red", alpha=0.95) + 
    #arrow labels    
    geom_text(aes(x = -0.375*0.90, y = 0.375*1.35),
        label="1", size=2.3, color="red") +
    geom_text(aes(x = -0.375*0.85, y = 0.375*1.55),
        label="2", size=2.3, color="red") +  
    geom_text(aes(x = -3, y = 0.75),
        label="Habitat", size=5, color="Orange") +
    geom_segment(aes(x = -3.5625, y = 0.9375, xend = -2.4375, yend = 0.9375),
        size=3, color="orange") +
    geom_segment(aes(x = -3.5625, y = 1.125, xend = -2.4375, yend = 1.125),
        size=3, color="orange") +
    geom_segment(aes(x = -3.5625, y = 0.5625, xend = -2.4375, yend = 0.5625),
        size=3, color="orange") +
    geom_segment(aes(x = -3.5625, y = 0.375, xend = -2.4375, yend = 0.375),
        size=3, color="orange") +
    labs(fill = "Abundance") +
    scale_color_manual(values = c("Orange", "Green", "Black"),
    # labels = c("Location",
    #     "Number of observations \nthat began at this location",
    #     "Number of observations \nthat chose this location")) +
    labels = c("Area",
        "Number of observations")) +
    guides(color = guide_legend(title = NULL, order = 1), keyheight = 0) +
    theme(text = element_text(size = 10), legend.position = "top",
    legend.box="vertical", legend.key.size = unit(0.4, 'cm'),
    legend.text = element_text(size=8), legend.title = element_text(size=8),
    legend.spacing.y = unit(0.0, 'cm'))

random2 <- mapped

ggsave(mapped, file = paste0(getwd(),
    "/inst/extdata/econ/mortup_8loc_random_20_2.png"),
    width = 3.2, height = 3.6, units = "in", dpi = 800)

test <- ggarrange(random1, random2, ncol = 2, nrow = 1)

ggsave(test, file = paste0(getwd(),
    "/inst/extdata/econ/fig4", ".png"),
    width = 6.4, height = 3.6, units = "in", dpi = 800)

################################################################################

trendsave <- list()
for (i in c(1, 5, 9, 13)) {

otherdatfin <- readRDS(paste0(getwd(),
  "/inst/extdata/econ/mortup_8loc_trend_20/",
  "otherdatfin_mortup_8loc_trend_20_", i, ".RData"))

csvout$id <- seq_len(nrow(csvout))
temp <- data.frame(cbind(L2=otherdatfin$choicefin$V1,
    catch=otherdatfin$catchfin$V1))
catches <- ddply(temp, .(L2), summarize, Mean_catch = mean(catch))
if (i == 13) {
  catches <- rbind(catches, c(7,0))
}
catches$Mean_catch <- otherdatfin$betavar
test <- merge(csvout, catches, by = c("L2"), all.x=TRUE, all.y=TRUE)
choices <- data.frame(table(otherdatfin$choicefin))
names(choices)[names(choices) == 'V1'] <- 'L2'
names(choices)[names(choices) == 'Freq'] <- 'End'
test <- merge(test, choices, by = c("L2"), all.x=TRUE, all.y=TRUE)
starts <- data.frame(table(otherdatfin$startloc))
names(starts)[names(starts) == 'Var1'] <- 'L2'
names(starts)[names(starts) == 'Freq'] <- 'Start'
test <- merge(test, starts, by = c("L2"), all.x=TRUE, all.y=TRUE)
test[is.na(test)] <- 0
test <- test[order(test$id), ] 

#trend
mapped <- ggplot() +
    geom_path(data = test, aes(x = X, y = Y, group = L2), size = 0.5) +
    # geom_point(aes(x = x, y = y, color = layer), size = 0.5) +
    # geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),
        # fill="grey", color="black", size=0.75) +
    geom_polygon(data = test,
      aes(x = X, y = Y, group = L2, fill = Mean_catch)) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43",
      limits=c(0.65,1.90)) +
    geom_text(data = test, aes(x = Cent_lon - .2, y = Cent_lat + .25,
        group = L2, label = Start,
        color = "C"), size = 2.8) +
    # geom_text(data=test, aes(x = Cent_lon-.25, y = Cent_lat+.15, group = L2,
    #     label = End,
    #     color="D"), size = 5) +
    geom_text(data = test, aes(x = Cent_lon + .25, y = Cent_lat + .25,
        group = L2,
        label = paste("Area \n", L2),
        color = "B"), size = 2.8, lineheight = .75) +
    geom_text(aes(x = -3, y = 0.75),
        label = "Habitat", size = 5, color = "Orange") +
    geom_segment(aes(x = -3.5625, y = 0.9375, xend = -2.4375, yend = 0.9375),
        size = 3, color = "orange") +
    geom_segment(aes(x = -3.5625, y = 1.125, xend = -2.4375, yend = 1.125),
        size = 3, color = "orange") +
    geom_segment(aes(x = -3.5625, y = 0.5625, xend = -2.4375, yend = 0.5625),
        size = 3, color = "orange") +
    geom_segment(aes(x = -3.5625, y = 0.375, xend = -2.4375, yend = 0.375),
        size = 3, color = "orange") +
    labs(fill = "Abundance") +
    scale_color_manual(values = c("Orange", "Green", "Black"),
    # labels = c("Location",
    #     "Number of observations \nthat began at this location",
    #     "Number of observations \nthat chose this location")) +
    labels = c("Area",
        "Number of observations")) +
    guides(color = guide_legend(title = NULL, order = 1)) +
    theme(text = element_text(size = 12), legend.position = "top")

trendsave[[i]] <- mapped

ggsave(mapped, file = paste0(getwd(),
    "/inst/extdata/econ/mortup_8loc_trend_20_", i, ".png"),
    width = 7.5, height = 4.21875, units = "in", dpi = 800)

}

test <- ggarrange(trendsave[[1]], trendsave[[5]], trendsave[[9]],
  trendsave[[13]], ncol = 2, nrow = 2, common.legend = TRUE)

ggsave(test, file = paste0(getwd(),
    "/inst/extdata/econ/fig5", ".png"),
    width = 7.5, height = 4.21875, units = "in", dpi = 800)

    # scale_x_continuous(breaks = round(seq(minlon, maxlon, by = 1),0)) +
    # scale_y_continuous(breaks = round(seq(minlat,maxlat, by = 1),0)) +
    # coord_cartesian(ylim=c(minlat, maxlat), xlim=c(minlon,maxlon)) +
    # scale_colour_manual(values = my.cols, name = "Depth") +
    # geom_polygon(aes(x = X, y = Y, group = L2, fill=Relative)) +
    # facet_grid(Fleet~Year) +
    # xlab("Longitude") +
    # ylab("Latitude") +
    # scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    # labs(fill='Abundance (Log biomass)') +
    # theme(text = element_text(size=25))

# testplot <- data.frame(cbind(otherdatfin$choicefin$V1, probmovesave))
# testplot$Location <- testplot$V1

# histprobs = ggplot(data=testplot) + 
#     geom_histogram(aes(x=probmovesave), bins=10) +
#     facet_wrap(~Location, scales = "free_y", ncol = 3) +
#     theme(text = element_text(size=8), axis.title.y=element_text(vjust=1.5)) + 
#     xlab("Probability of location choice") + 
#     ylab("Number of observations")+
#     coord_cartesian(xlim =c(0,1)) +
#     theme(text = element_text(size=10))
# ggsave(histprobs, file=paste0(nicedir, "\\start1_fewobs_support", ".png"), 
#     width = 7.5, height = 4.21875, units = "in", dpi=800)