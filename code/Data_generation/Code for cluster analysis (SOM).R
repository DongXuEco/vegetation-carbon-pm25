# Cluster AGB, BGB, SOC and PM2.5 adsorption using 16 SOM units and four groups.
# Input rows must contain finite numeric values; variables are z-score standardized.
# Cluster numbers are labels; interpret them from their variable profiles.
library(dplyr)
library(readxl)
library(kohonen)

input_file <- ""
output_dir <- ""
variable_columns <- c("AGB", "BGB", "SOC", "Q")  # Corresponding Excel column names.
cluster_method <- "ward.D2"  # Euclidean distances between SOM codebook vectors.

if (!nzchar(input_file) || !nzchar(output_dir)) {
  stop("Set input_file and output_dir.")
}
read_excel(input_file) -> esb
if (length(variable_columns) != 4L || anyDuplicated(variable_columns) ||
    !all(variable_columns %in% names(esb))) {
  stop("Specify four distinct columns for AGB, BGB, SOC and Q.")
}
esv <- esb %>% select(all_of(variable_columns))
if (!all(vapply(esv, is.numeric, logical(1)))) {
  stop("All four variables must be numeric.")
}
if (nrow(esv) < 16L || !all(is.finite(as.matrix(esv)))) {
  stop("Provide at least 16 complete rows with finite values.")
}
if (any(vapply(esv, sd, numeric(1)) == 0)) {
  stop("Each variable must have nonzero variance.")
}
scale(esv) -> esv

set.seed(2022)
g <- somgrid(xdim=4, ydim=4, topo="rectangular")
map <- som(esv, grid=g, radius=1)

plot(map, type="counts")
plot(map, type="codes")

# Group the 16 units, then assign each sample to its unit's cluster.
unit_tree <- hclust(dist(map$codes[[1]]), method=cluster_method)
unit_cluster <- cutree(unit_tree, k=4)
sample_cluster <- unname(unit_cluster[map$unit.classif])

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
prototypes <- data.frame(
  SOM_unit=seq_len(nrow(map$codes[[1]])),
  Cluster=unname(unit_cluster),
  map$codes[[1]], check.names=FALSE
)
write.csv(prototypes, file=file.path(output_dir, "ESB.csv"), row.names=FALSE)

output <- esb
output$SOM_unit <- map$unit.classif
output$Cluster <- sample_cluster
write.csv(output, file=file.path(output_dir, "BundleResultsSoM.csv"), row.names=FALSE)

