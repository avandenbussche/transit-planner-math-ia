library(dbscan)
library(ggplot2)
library(mclust)
library(geosphere)
library(proxy)
library(igraph)
library(ggmap)

source("funcs.R")



# Load data into dataframe
POIs <- read.csv("pois.csv", header = FALSE)
colnames(POIs) <- c("id", "google_places_id", "name", "latitude", "longitude", "weight_determining_type", "weight")
POIs$id <- NULL # Delete the column for google_places_id
POIs$google_places_id <- NULL # Delete the column for google_places_id



# Only selected weighted POIs and change the domain of their weight
weighted.POIs <- POIs[POIs$weight > 0,]
weighted.POIs$weight <- fit.domain(weighted.POIs$weight, 0, max(weighted.POIs$weight), 0, 10)



# 0.0006 80
# 0.0004 50 A bit too dense imo
# 0.0004 60 perf with stdev 500
# 0.0003 65
# 0.0004 40 actually very interesting as an extended network but downtown is underserved
# 0.0005 40 cool but still weird
# 0.00053 50 getting there
# 0.0003 65 GOOD BUT SMALL
# Perform DBSCAN and save cluster to POIs dataframe
x <- as.matrix(weighted.POIs[, 2:3])
dbscan.results <- dbscan(x, eps = 0.0003, minPts = 65, weights = as.numeric(weighted.POIs$weight))
weighted.POIs$cluster <- as.factor( dbscan.results$cluster ) # So ggplot does not interpolate colors



# Cool plot of the city (looks kinda old!)
plot(weighted.POIs$longitude, weighted.POIs$latitude, main = "POIs in Montreal", xlab = "Longitude", ylab = "Latitude", col = adjustcolor("black", alpha=0.15), cex = 0.1)



# Plot and save map of first round of clusters
ggplot(weighted.POIs, aes(x = longitude, y = latitude, color = cluster)) +
  geom_point(shape = 16, size = 0.1, alpha = 0.15) +
  theme(panel.background = element_blank(), legend.position = "none") +
  coord_fixed() + labs(x = "Longitude", y = "Latitude")

ggplot() +
  geom_point(data = weighted.POIs[weighted.POIs$cluster %in% c(0),], aes(x = longitude, y = latitude), shape = 16, size = 0.1, alpha = 0.05, color = "black") +
  geom_point(data = weighted.POIs[!(weighted.POIs$cluster %in% c(0)),], aes(x = longitude, y = latitude, color = cluster), shape = 16, size = 0.1, alpha = 1) +
  theme(panel.background = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  coord_fixed() + labs(x = "Longitude", y = "Latitude") + ggtitle("Clusters of POIs Found by DBSCAN Algorithm")
#ggsave("yul.pdf")



# Loop through all clusters to find their centroids and weights
cluster.centroids <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(cluster.centroids) <- c("cluster", "centroid_latitude", "centroid_longitude")
num.clusters <- max(dbscan.results$cluster) # Find the total number of clusters
acceptable.cluster.values <- seq(1, num.clusters)
for (this_cluster in acceptable.cluster.values) {
  centroid_latitude <- weighted.mean( weighted.POIs[weighted.POIs$cluster %in% this_cluster,]$latitude, weighted.POIs[weighted.POIs$cluster %in% this_cluster,]$weight )
  centroid_longitude <- weighted.mean( weighted.POIs[weighted.POIs$cluster %in% this_cluster,]$longitude, weighted.POIs[weighted.POIs$cluster %in% this_cluster,]$weight )
  row_to_add <- data.frame(cluster = this_cluster, centroid_latitude = centroid_latitude, centroid_longitude = centroid_longitude)
  cluster.centroids <- rbind(cluster.centroids, row_to_add)
}
rm(row_to_add, centroid_latitude, centroid_longitude, this_cluster)




# For each centroid:
#  Find all the POIs in a given radius and their weights
search.radius <- 1600 # in metres
cluster.centroids["weight"] <- 0
cluster.centroids["line"] <- 0
for (this_cluster in acceptable.cluster.values) {
  this_centroid <- cluster.centroids[cluster.centroids$cluster %in% this_cluster,]
  
  lat_north <- distance.to.latitude(this_centroid$centroid_latitude, search.radius, "N")
  lat_south <- distance.to.latitude(this_centroid$centroid_latitude, search.radius, "S")
  lon_east <- distance.to.longitude(this_centroid$centroid_longitude, this_centroid$centroid_latitude, search.radius, "E")
  lon_west <- distance.to.longitude(this_centroid$centroid_longitude, this_centroid$centroid_latitude, search.radius, "W")
  
  relevant.POIs.weights <- subset(weighted.POIs, latitude <= lat_north & latitude >= lat_south & longitude >= lon_west & longitude <= lon_east, select = weight)
  relevant.POIs.coords <- subset(weighted.POIs, latitude <= lat_north & latitude >= lat_south & longitude >= lon_west & longitude <= lon_east, select = c(latitude, longitude))
  relevant_POIs_distance_from_centroid <- distHaversine(as.vector(this_centroid[,2:3]), relevant.POIs.coords)
  
  # STDEV was 500
  sum_of_weights <- sum( weight.given.normal.distribution(relevant.POIs.weights, 300, relevant_POIs_distance_from_centroid) )
  cluster.centroids[cluster.centroids$cluster %in% this_cluster,]$weight <- sum_of_weights
}; rm(sum_of_weights, this_cluster, this_centroid, lat_north, lat_south, lon_east, lon_west, relevant_POIs_distance_from_centroid)

region.westernmost.lon <- min(cluster.centroids$centroid_longitude)
region.easternmost.lon <- max(cluster.centroids$centroid_longitude)

ggplot(cluster.centroids, aes(x = centroid_longitude, y = centroid_latitude, color = weight)) + geom_point(shape = 16, size = 0.8, alpha = 1) + theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed() + labs(x = "Longitude", y = "Latitude") + ggtitle("Centroids of Clusters of POIs in Greater Montreal \n as Found by DBSCAN Algorithm (Noise Removed)")
ggsave("cool.pdf")




# Populate distance and weight matrices
sum.of.weights <- function(x, y) { x + y }
distance.matrix <- distm(cluster.centroids[,2:3], fun = distHaversine)
combined.weight.matrix <- as.matrix( dist(cluster.centroids[,4], sum.of.weights) )

cost.matrix <- matrix(Inf, nrow = num.clusters ^ 2, ncol = 3)
i <- 1
for (x in acceptable.cluster.values) {
  for (y in acceptable.cluster.values) {
    # 300 was 500
    cost.matrix[i,] <- c(x, y, path.cost.3(distance.matrix[x, y], combined.weight.matrix[x, y], 300, 10^(-5), 100))
    i <- i + 1
  }
}; rm(i, x, y)




# Find the central station (the closest possible station to the centroid of the entire network)
central.station.index <- cluster.centroids[ cluster.centroids$weight == max(cluster.centroids$weight), ]$cluster




# Find the transit lines in the network
gmm <- Mclust( rev(cluster.centroids[,2:3]) )
cluster.centroids$gmm_cluster <- gmm$classification
pdf("axes.pdf", width=11, height=8.5) 
plot(gmm, what = "classification")
dev.off()

ggplot(cluster.centroids, aes(x = centroid_longitude, y = centroid_latitude, color = as.factor(gmm_cluster))) + geom_point(shape = 16, size = 0.2, alpha = 1) + theme(panel.background = element_blank()) + coord_fixed() + labs(x = "Longitude", y = "Latitude")

# Prepare GMM cluster data
gmm.clusters <- list(clusters = 1:gmm$G, cluster.means = gmm$parameters$mean,
                     covariance.matrices = gmm$parameters$variance$sigma, eigenvalues = NULL,eigenvectors = NULL,
                     termini = NULL, type.determining.magnitude = NULL, average.magnitude = 0, types = NULL)

# Find eigenvalues and eigenvectors of each GMM cluster
for (G in gmm.clusters$clusters) {
  eigenstuff <- eigen(gmm.clusters$covariance.matrices[,,G])
  gmm.clusters$eigenvalues[G] <- list(eigenstuff$values)
  gmm.clusters$eigenvectors[G] <- list(eigenstuff$vectors)
  
  # Calculate the vectors
  scaled.eigen.vec.m <- eigenstuff$vectors[,1] * sqrt( eigenstuff$values[1] )
  scaled.eigen.vec.n <- eigenstuff$vectors[,2] * sqrt( eigenstuff$values[2] )
  vector.a <- as.vector(gmm.clusters$cluster.means[,G]) + scaled.eigen.vec.m
  vector.b <- as.vector(gmm.clusters$cluster.means[,G]) + scaled.eigen.vec.n
  vector.c <- as.vector(gmm.clusters$cluster.means[,G]) - scaled.eigen.vec.m # More for plotting purposes
  vector.d <- as.vector(gmm.clusters$cluster.means[,G]) - scaled.eigen.vec.n # More for plotting purposes
 
  # Plot eigenvectors
  plot( c(gmm.clusters$cluster.means[,G][1], vector.a[1], vector.b[1], vector.c[1], vector.d[1]), c(gmm.clusters$cluster.means[,G][2], vector.a[2], vector.b[2], vector.c[2], vector.d[2]) )
  arrows( gmm.clusters$cluster.means[,G][1], gmm.clusters$cluster.means[,G][2], vector.a[1], vector.a[2] )
  arrows( gmm.clusters$cluster.means[,G][1], gmm.clusters$cluster.means[,G][2], vector.b[1], vector.b[2] )
  arrows( gmm.clusters$cluster.means[,G][1], gmm.clusters$cluster.means[,G][2], vector.c[1], vector.c[2] )
  arrows( gmm.clusters$cluster.means[,G][1], gmm.clusters$cluster.means[,G][2], vector.d[1], vector.d[2] )
  
  # Find relevant slope
  mag.a <- magnitude( gmm.clusters$cluster.means[,G] - vector.a )
  mag.b <- magnitude( gmm.clusters$cluster.means[,G] - vector.b )
  if (mag.a >= mag.b) {
    slope <- vector.a[2] / vector.a[1]
    gmm.clusters$type.determining.magnitude[G] <- list(mag.a)
  } else {
    slope <- vector.b[2] / vector.b[1]
    gmm.clusters$type.determining.magnitude[G] <- list(mag.b)
  }
  
  # Find points where line leaves region of interest (point-slope form)
  lat.westernmost.endpoint <- slope * (region.westernmost.lon - gmm.clusters$cluster.means[,G][1]) + gmm.clusters$cluster.means[,G][2]
  lat.easternmost.endpoint <- slope * (region.easternmost.lon - gmm.clusters$cluster.means[,G][1]) + gmm.clusters$cluster.means[,G][2]
  
  # Find index of closest stations to endpoints (the termini)
  current.closest.westernmost.index <- -1
  current.closest.westernmost.distance <- Inf
  current.closest.easternmost.index <- -1
  current.closest.easternmost.distance <- Inf
  
  current.searchable.clusters <- cluster.centroids[cluster.centroids$gmm_cluster == G,]
  for (this_cluster in current.searchable.clusters$cluster) {
    this_centroid <- current.searchable.clusters[current.searchable.clusters$cluster %in% this_cluster,]
    western.distance <- distHaversine(this_centroid[,3:2], c(region.westernmost.lon, lat.westernmost.endpoint))
    eastern.distance <- distHaversine(this_centroid[,3:2], c(region.easternmost.lon, lat.easternmost.endpoint))
    if ( western.distance < current.closest.westernmost.distance ) {
      current.closest.westernmost.index <- this_centroid$cluster
      current.closest.westernmost.distance <- western.distance
    }
    if ( eastern.distance < current.closest.easternmost.distance ) {
      current.closest.easternmost.index <- this_centroid$cluster
      current.closest.easternmost.distance <- eastern.distance
    }
  }
  
  gmm.clusters$termini[G] <- list( c(current.closest.westernmost.index, current.closest.easternmost.index) )
  
}; rm(G, scaled.eigen.vec.m, scaled.eigen.vec.n, vector.a, vector.b, vector.c, vector.d, mag.a, mag.b, slope,
      lat.westernmost.endpoint, lat.easternmost.endpoint, current.closest.westernmost.index,
      current.closest.westernmost.distance, current.closest.easternmost.index, current.closest.easternmost.distance,
      current.searchable.clusters, this_centroid, this_cluster, western.distance, eastern.distance)



# Prepare igraph data
cost.matrix.filtered <- cost.matrix[ cost.matrix[,3] > 0, ]
entire.possible.network <- graph.data.frame(cost.matrix.filtered[, 1:2], directed = FALSE) %>%
  set_vertex_attr("gmm_cluster", value = cluster.centroids$gmm_cluster) %>%
  set_vertex_attr("cluster_id", value = cluster.centroids$cluster)
E(entire.possible.network)$weight <- cost.matrix.filtered[,3]

# Prepare line graphing data table
line.data.for.drawing <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(line.data.for.drawing) <- c("line", "x", "y", "xend", "yend")


# Get average magnitude of each line and determine whether line is metro or train
gmm.clusters$average.magnitude <- mean(unlist(gmm.clusters$type.determining.magnitude))
for (G in gmm.clusters$clusters) {
  print( (as.numeric(gmm.clusters$type.determining.magnitude[G]) - gmm.clusters$average.magnitude)^2  )
  print(G)
  if (gmm.clusters$type.determining.magnitude[G] <= gmm.clusters$average.magnitude) {
    gmm.clusters$types[G] <- list("metro")
    metro <- TRUE 
  } else {
    gmm.clusters$types[G] <- list("train")
    metro <- FALSE
  }
  
  # If the line is a metro, its two termini must connect
  if (metro) {
    western.terminus <- gmm.clusters$termini[[G]][1]
    eastern.terminus <- gmm.clusters$termini[[G]][2]
    
    this.metro.line <- induced.subgraph(entire.possible.network,
                                        which( V(entire.possible.network)$gmm_cluster %in% c(G)) )

    this.shortest.path <- shortest_paths(this.metro.line, from = match(as.character(western.terminus),
                                                                       V(this.metro.line)$cluster_id),
                                         to = match(as.character(eastern.terminus),
                                                    V(this.metro.line)$cluster_id), mode = "in")
    
    
    # Minimum Spanning Tree
    this.metro.line.for.tree <- induced.subgraph(this.metro.line, 
                                                  which(V(this.metro.line) %in% this.shortest.path$vpath[[1]]) )
   
    min.spanning.tree <- minimum.spanning.tree(this.metro.line.for.tree)
    V(min.spanning.tree)$label.cex <- 0.2
    plot.igraph(min.spanning.tree,
      layout = as.matrix(cluster.centroids[cluster.centroids$cluster %in% this.shortest.path$vpath[[1]]$name, 3:2]),
      vertex.label.size = "",
      vertex.size = .2,
      rescale = FALSE,
      xlim = c(min(cluster.centroids[cluster.centroids$cluster %in% this.shortest.path$vpath[[1]]$name, 3]), max(cluster.centroids[cluster.centroids$cluster %in% this.shortest.path$vpath[[1]]$name, 3])),
      ylim = c(min(cluster.centroids[cluster.centroids$cluster %in% this.shortest.path$vpath[[1]]$name, 2]), max(cluster.centroids[cluster.centroids$cluster %in% this.shortest.path$vpath[[1]]$name, 2])))
    title(as.character(G))
    
    
    # Update clusters to line if more than 5 stations
    if (length(this.shortest.path$vpath[[1]]) > 5) {
      cluster.centroids[cluster.centroids$cluster %in% this.shortest.path$vpath[[1]]$name, ]$line <- G
      
      # Compile data for min spanning tree so it can be drawn in final map
      for (i in 1:length(E(min.spanning.tree))) {
        x <- cluster.centroids[cluster.centroids$cluster == as.integer(ends(min.spanning.tree, E(min.spanning.tree))[i,1]), 3]
        xend <- cluster.centroids[cluster.centroids$cluster == as.integer(ends(min.spanning.tree, E(min.spanning.tree))[i,2]), 3]
        y <- cluster.centroids[cluster.centroids$cluster == as.integer(ends(min.spanning.tree, E(min.spanning.tree))[i,1]), 2]
        yend <- cluster.centroids[cluster.centroids$cluster == as.integer(ends(min.spanning.tree, E(min.spanning.tree))[i,2]), 2]
        row.to.add <- data.frame(G, x, y, xend, yend)
        colnames(row.to.add) <- c("line", "x", "y", "xend", "yend")
        line.data.for.drawing <- rbind(line.data.for.drawing, row.to.add)
      }; rm(x, xend, y, yend, row.to.add, i)
    }

  } else { # Otherwise, the train line must connect to either a metro station that reaches the central station or directly to the central station
    
  }
}; rm(G, metro, this.metro.line, this.shortest.path, western.terminus, eastern.terminus)





# Plot all the lines on one map
geographical.map.data <- cluster.centroids[cluster.centroids$line > 0,]
geographical.map.average.longitude <- mean(geographical.map.data$centroid_longitude)
geographical.map.average.latitude <- mean(geographical.map.data$centroid_latitude)
geographical.map.center <- c(lon = geographical.map.average.longitude, lat = geographical.map.average.latitude)
geographical.map.constructor <- get_map(location = geographical.map.center, zoom = 12, source="google", maptype = "roadmap", color="bw")
# Add existing stations (from STM)
all.existing.stations <- read.csv("stm_stops.csv")
relevant.existing.stations <- all.existing.stations[all.existing.stations$location_type == 1,]
# Plot
overlay.map <- ggmap(geographical.map.constructor) +
  geom_point(data = relevant.existing.stations,
             aes(x = stop_lon, y = stop_lat), color = "black", shape = 17) + 
  geom_point(data = geographical.map.data,
             aes(x = centroid_longitude, y = centroid_latitude, color = as.factor(line))) +
  geom_segment(data = line.data.for.drawing, aes(x = x, y = y, xend = xend, yend = yend, color = as.factor(line))) +
  theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(x = "Longitude", y = "Latitude") +
  ggtitle("Final Suggested Metro System")
overlay.map
