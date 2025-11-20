library(plyr)
library(ggplot2)
library(RColorBrewer)
library(ss3simfork)

dirsin <- c("mortup_9loc",
  "mortup_9loc_singlestart",
  "mortup_8loc_random_20",
  "mortup_8loc_trend_20")

itmax <- 200

linegraphsave <- NULL
for (i in 1) {

# rawoutputdir for files (iterations and scenarios)
rawfolder <- "C:\\Users\\Allen.Chen\\Desktop\\size_select\\"
rawoutputdir <- paste0(rawfolder,
  dirsin[i])

abundout <- list()
for (it in 1:itmax) {

# The default naming convention for each abundance index .csv is below, if
# the file name is left NULL, but can be specified if desired.
casenamecor <- paste0(rawoutputdir, "\\D4-E1-F2-M0-cod\\", it,
  "\\cor-abund-D4-E1-F2-M0-cod-")
casenamebias <- paste0(rawoutputdir, "\\D3-E1-F2-M0-cod\\", it,
  "\\bias-abund-D3-E1-F2-M0-cod-")
casenamebias <- paste0(rawoutputdir, "\\D3-E1-F2-M0-cod\\", it,
  "\\bias-abund-D3-E1-F2-M0-cod-")
casenamesamp <- paste0(rawoutputdir, "\\D2-E1-F2-M0-cod\\", it,
  "\\bias-abund-D2-E1-F2-M0-cod-")

if (file.exists(paste(casenamecor, it, ".csv", sep =  "")) &&
  file.exists(paste(casenamebias, it, ".csv", sep = "")) &&
  file.exists(paste(casenamesamp, it, ".csv", sep = ""))) {
abundcor <- readabund(filename = paste(casenamecor, it, ".csv", sep = ""),
  type = "Corrected")

abundtrue <- readabund(filename = paste(casenamebias, it, ".csv", sep = ""),
  type = "True", istruecpue = TRUE)

abundbias <- readabund(filename = paste(casenamebias, it, ".csv", sep = ""),
  type = "Uncorrected")

abundsamp <- readabund(filename = paste(casenamesamp, it, ".csv", sep = ""),
  type = "Randomly sampled")

tempabundout <- rbind(abundtrue, abundbias, abundcor, abundsamp)

tempabundout$iter <- it

abundout[[it]] <- tempabundout
}

}

totabundout <- do.call(rbind, abundout)

totabundout <- totabundout[totabundout$CPUE > 0, ]

linegraph <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), mean,
  na.rm = TRUE)
names(linegraph)[names(linegraph) == "x"] <- "RelDiff"

linegraphse <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), sd)
names(linegraphse)[names(linegraphse) == "x"] <- "sd"

linegraphlen <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), length)
names(linegraphlen)[names(linegraphlen) == "x"] <- "len"

linegraph <- merge(linegraph, linegraphse, by = c("Year", "Index"))
linegraph <- merge(linegraph, linegraphlen, by = c("Year", "Index"))

linegraph$se <- linegraph$sd/sqrt(linegraph$len)

linegraph$Index <- revalue(linegraph$Index, c("Corrected" = "Corrected     ",
  "Randomly sampled" = "Randomly sampled     ",
  "True" = "True     ",
  "Uncorrected" = "Uncorrected     "))

linegraph <- linegraph[order(linegraph$Year), ]

pd <- position_dodge(0.05)

lineout <- ggplot(data = linegraph,
  aes(x = Year, y = RelDiff)) +
  geom_line(aes(colour = Index, group = Index), size = 2, position = pd) +
  geom_point(aes(shape = Index), size = 4, fill = "white", stroke = 1.5,
    position = pd) +
  xlab("Year") + ylab("Relative Differences in Abundance") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "orange",
    size = 2) +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(legend.direction = "horizontal", legend.position = "top",
    legend.key = element_rect(size = 4),
    legend.text = element_text(size = 22),
    legend.key.size = unit(2.5, "lines"), legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 30),
    axis.title.y = element_text(vjust = 1.5)) +
  scale_shape_manual(values = c(21, 22, 25, 24)) +
  scale_colour_manual(values = brewer.pal(4, "Set1"))

###############################################################################
totabundout <- do.call(rbind, abundout)

# totabundout <- totabundout[totabundout$CPUE > 0, ]
# totabundout <- totabundout[is.na(totabundout$CPUE) == FALSE, ]

totabundtemp <- reshape(
  data = totabundout[c("year", "iter", "Index", "RelDiff")],
    idvar = c("year", "iter"),
    v.names = c("RelDiff"),
    timevar = "Index",
    direction = "wide")
totabundtemp <- totabundtemp[is.na(totabundtemp$RelDiff.Corrected) == FALSE, ]

totabundtemp$RelDiff.Uncorrected <- totabundtemp$RelDiff.Uncorrected -
  totabundtemp$RelDiff.True

totabundtemp$RelDiff.Corrected <- totabundtemp$RelDiff.Corrected -
  totabundtemp$RelDiff.True

totabundtemp[c("RelDiff.Randomly sampled")] <-
  totabundtemp[c("RelDiff.Randomly sampled")] -
  totabundtemp$RelDiff.True

totabundout <- reshape(data = totabundtemp,
  idvar = c("year", "iter"),
  v.names = c("RelDiff"),
  timevar = "Index",
  direction = "long")

linegraph <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), mean,
  na.rm = TRUE)
names(linegraph)[names(linegraph) == "x"] <- "RelDiff"

linegraphse <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), sd, na.rm = TRUE)
names(linegraphse)[names(linegraphse) == "x"] <- "sd"

linegraphlen <- aggregate(totabundout$RelDiff,
  list(Year = totabundout$year, Index = totabundout$Index), length)
names(linegraphlen)[names(linegraphlen) == "x"] <- "len"

linegraph <- merge(linegraph, linegraphse, by = c("Year", "Index"))
linegraph <- merge(linegraph, linegraphlen, by = c("Year", "Index"))

linegraph$se <- linegraph$sd/sqrt(linegraph$len)

linegraph$Index <- revalue(linegraph$Index, c("Corrected" = "Corrected     ",
  "Randomly sampled" = "Randomly sampled     ",
  "True" = "True     ",
  "Uncorrected" = "Uncorrected     "))

linegraph <- linegraph[!(linegraph$Index == "True     "), ]

linegraph$Scenario <- i

linegraphsave <- rbind(linegraphsave, linegraph)

pd <- position_dodge(0.05)

lineout <- ggplot(data = linegraph,
  aes(x = (Year), y = RelDiff)) +
  geom_line(aes(colour = Index, group = Index), size = 1, position = pd) +
  geom_point(aes(shape = Index), size = 2, fill = "white", stroke = 0.75,
    position = pd) +
  geom_ribbon(aes(x = Year, ymin = RelDiff - (2 * se),
    ymax = RelDiff + (2 * se),
    fill = Index, alpha = Index)) +
  xlab("Year") + ylab("Relative error in abundance") +
  scale_x_continuous(breaks = c(88:100)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.direction = "horizontal", legend.position = "top",
    legend.key = element_rect(size = 1.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.75, "lines"), legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 10),
    axis.title.y = element_text(vjust = 1.5)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_alpha_discrete(range = c(0.45, 0.25)) +
  scale_colour_manual(values = brewer.pal(4, "Set1")[c(1, 2, 4)],
  aesthetics = c("colour", "fill"))
print(lineout)

}

lineout <- ggplot(data = linegraphsave,
  aes(x = (Year), y = RelDiff)) +
  geom_line(aes(colour = Index, group = Index), size = 2, position = pd) +
  geom_point(aes(shape = Index), size = 3, fill = "white", stroke = 1.2,
    position = pd) +
  geom_ribbon(aes(x = Year, ymin = RelDiff - (2 * se),
    ymax = RelDiff + (2 * se),
    fill = Index, alpha = Index)) +
  facet_wrap(~Scenario) +
  xlab("Year") + ylab("Relative error in abundance") +
  scale_x_continuous(breaks = c(88:100)) +
  coord_cartesian(ylim = c(-0.019, 0.172)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.direction = "horizontal", legend.position = "top",
    strip.text = element_text(size = 10, margin = margin()),
    legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.key = element_rect(size = 2.5),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.75, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 20),
    axis.title.y = element_text(vjust = 1.5)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_alpha_discrete(range = c(0.45, 0.25)) +
  scale_colour_manual(values = brewer.pal(4, "Set1")[c(1, 2, 4)],
  aesthetics = c("colour", "fill"))
print(lineout)

ggsave(lineout, filename = paste0(getwd(), "/inst/extdata/econ/fig6.png"),
  width = 8, height = 8, dpi = 800)