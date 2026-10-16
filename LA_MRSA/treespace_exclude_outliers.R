
#SKIP # copy from slack
#SKIP 
#SKIP # Load required libraries
#SKIP library(ape)
#SKIP library(treespace)
#SKIP library(Cairo)
#SKIP options(bitmapType='cairo')
#SKIP getOption("device") 
#SKIP 
#SKIP capabilities()
#SKIP 
#SKIP # save(tree_list, file = paste("_", cnt, "_tree_list.RData", sep=""))
#SKIP load("_1361_tree_list.RData")
#SKIP length(tree_list)
#SKIP 
#SKIP class(tree_list) = "multiPhylo"
#SKIP 
#SKIP ts_res <- treespace(tree_list, nf=2) # default method= treeVec for Kendall Colijn metric vector (for rooted tree)
#SKIP save(ts_res, file = paste("ts_res.RData", sep=""))
#SKIP 
#SKIP # Load or generate example tree data
#SKIP # Here we're generating a set of 10 random trees with 15 tips each
#SKIP # trees <- rmtree(10, 15)  # Replace with your own tree data
#SKIP # Run treespace analysis to get distances and MDS coordinates
#SKIP # ts_res <- treespace(trees, nf = 2)  # nf = 2 for a 2D MDS plot
#SKIP # Plot MDS coordinates of the trees
#SKIP #plot(ts_res$points, main = "Treespace MDS Plot", xlab = "MDS1", ylab = "MDS2")
#SKIP #text(ts_res$points, labels = 1:length(trees), pos = 3, cex = 0.7)  # Label each point
#SKIP 
#SKIP ts_res_df <- ts_res$pco$li
#SKIP x = ts_res_df[,1]
#SKIP y = ts_res_df[,2]
#SKIP plot(ts_res_df[,1], ts_res_df[,2])
#SKIP 
#SKIP identify(x, y, paste(1:len(x)))
#SKIP 
#SKIP 
#SKIP # Calculate the centroid of the MDS points (mean coordinates)
#SKIP centroid <- colMeans(ts_res$points)
#SKIP 
#SKIP # Calculate Euclidean distances of each tree from the centroid
#SKIP distances <- apply(ts_res$points, 1, function(x) sqrt(sum((x - centroid)^2)))
#SKIP 
#SKIP # Define a threshold for outliers (e.g., 2 standard deviations from the mean distance)
#SKIP threshold <- mean(distances) + 2 * sd(distances)
#SKIP 
#SKIP # Identify the indices of trees that are outliers
#SKIP outliers <- which(distances > threshold)
#SKIP 
#SKIP # Print outliers
#SKIP print(outliers)
#SKIP 
#SKIP plotGrovesD3(grove_results$pco)
#SKIP xx> grove_results <-findGroves(res)


# [11/13/2024_1330/Wed]
# cd ~/Output_Nov5ccir_2rpt/pt09_treespace_contour_plot
load("1361_tree_list.RData")
length(tree_list)
names(tree_list)

# run save_ts_res.R on JC-MAC

# GPT in R how to identify points outside of an 2-dim slanted ellipsoid area in x-y plot

pdf("included_genes.pdf")
load("ts_res.RData")
names(ts_res)
ts_res_df <- ts_res$pco$li
x = ts_res_df[,1] + 2000
y = ts_res_df[,2] - 1000
# plot(x, y)
range(x)
range(y)
#-------------------------------------------------------------------------------------
# Define the ellipsoid parameters
# a <- 3  # semi-major axis (along x-axis)
# b <- 2  # semi-minor axis (along y-axis)
a = 22000
b = 3400
theta <- pi / 4 / 2  # rotation angle (45/2 degrees)
# Rotation matrix for the given angle
R_old <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow=2, byrow=TRUE)
R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, byrow=TRUE)
# Define some points to check (x, y coordinates)
# points <- matrix(c(4, 1, 2, 2, 0, 0, 1, 1), ncol=2, byrow=TRUE)
points = matrix(c(x, y), ncol=2, byrow=FALSE)
# Function to check if a point is inside the rotated ellipsoid
is_inside_ellipsoid <- function(point) {
  # Apply the rotation matrix to the point
  point_rot <- R %*% point  # Rotated coordinates
  # Ellipsoid equation check
  (point_rot[1]^2 / a^2 + point_rot[2]^2 / b^2) <= 1
}
# Apply the check to each point
inside_check <- apply(points, 1, is_inside_ellipsoid)
# Print whether points are inside or outside
inside_check
# Plot the ellipsoid and the points
# plot(NA, xlim=c(-4, 4), ylim=c(-4, 4), xlab="X", ylab="Y", main="2D Rotated Ellipsoid and Points")
plot(NA, xlim=range(x), ylim=range(y), xlab="X", ylab="Y", main="2D Rotated Ellipsoid and Points")
# Draw the ellipse (axis-aligned)
theta_vals <- seq(0, 2*pi, length.out=100)
ellipse_x <- a * cos(theta_vals)
ellipse_y <- b * sin(theta_vals)
# Apply rotation to the ellipse
ellipse_rot <- t(R_old %*% rbind(ellipse_x, ellipse_y))
# ellipse_rot <- (R %*% rbind(ellipse_x, ellipse_y))
lines(ellipse_rot[,1], ellipse_rot[,2], col="blue")  # Draw the rotated ellipsoid
# Plot the points
points(points, col=ifelse(inside_check, "green", "red"), pch=16)
# Add a legend
legend("topright", legend=c("Inside", "Outside"), col=c("green", "red"), pch=16)
#-------------------------------------------------------------------------------------
dev.off()
row.names(ts_res_df)[ inside_check ]
write.table(row.names(ts_res_df)[ inside_check ], file="included_genes.txt", row.names = FALSE, quote=FALSE)

# cd ~/Output_Nov5ccir_2rpt/pt09_treespace_contour_plot; R --no-save < ~/Summer2024/ccirwork/treespace_exclude_outliers.R


