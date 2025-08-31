library(ggplot2)

samples <- read.csv("AllCoord.csv", stringsAsFactors = FALSE)
cn <- tolower(names(samples))
lon_idx <- which(cn %in% c("longitude","lon","long","x"))
lat_idx <- which(cn %in% c("latitude","lat","y"))
samples$lon <- as.numeric(samples[[lon_idx[1]]])
samples$lat <- as.numeric(samples[[lat_idx[1]]])

ggplot(samples, aes(x = lon, y = lat, color = Samples)) +
  borders("world", colour = "gray85", fill = "gray90") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_minimal() +
  labs(title = "Sample Locations on World Map",
       x = "Longitude", y = "Latitude", color = "Region")
