# ==========================================================================
# Simulation Effect of Noise on Categorization
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, patchwork, ggplot2, ggtext, ellipse, plyr, ggnewscale, stringr)
theme_set(theme_classic())

# ==========================================================================
nf <- 2 # number of features
s <- 1/10 # modification factor
levels <- c("low", "medium", "high") # category overlap levels
breaks <- seq(6, 16, length.out = 5)
# ==========================================================================

# ==========================================================================
# Makes the information-integration category structures
# ==========================================================================
a_low <- c(92.5, 127.5) # centroid of category A (medium-low)
b_low <- c(127.5, 92.5) # centroid of category B (medium-low)

a_medium <- c(101.5, 118.5) # centroid of category A (medium)
b_medium <- c(118.5, 101.5) # centroid of category B (medium)

a_high <- c(106.1, 113.9) # centroid of category A (medium-high)
b_high <- c(113.9, 106.1) # centroid of category B (medium-high)

means <- mget(c("a_low", "b_low", "a_medium", "b_medium", "a_high", "b_high"))
means <- lapply(means, "*", s)

var <- 162.5; covar <- 112 # variances and covariances of each category
var <- var * s^2; covar <- covar * s^2
c_cov <- matrix(c(var, covar, covar, var), nrow = nf)

# ==========================================================================
# Makes the transfer grid (same for all category structures)
# ==========================================================================
grid_mean <- colMeans(do.call("rbind", means))
range <- lapply(grid_mean, function(i) {i + c(-1, 1) * 4 * sqrt(var)})
values <- sapply(1:length(range), function(i) {seq(min(range[[i]]), max(range[[i]]), by = 2 * sqrt(var))})
trans <- as.data.table(expand.grid(x = values[, 1], y = values[, 2]))
trans[, alpha := (x %in% c(min(x), mean(x), max(x)) & y %in% c(min(y), mean(y), max(y)))]

# ==========================================================================
# Makes the information-integration plot
# ==========================================================================
contours <- as.data.table(ldply(means, ellipse, x = cov2cor(c_cov), scale = sqrt(diag(c_cov)), npoints = 1000))
contours[, c("category", "overlap") := data.table(str_split_fixed(.id, "_", 2))]

ii <- ggplot(contours, aes(x, y, group = category, color = category)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  # geom_point(data = trans, aes(group = NULL, fill = alpha), pch = 21, color = grey(.25), size = 4) +
  # scale_fill_manual(values = c(grey(.25, 0), grey(.25, 1))) +
  # new_scale("fill") +
  geom_point(data = trans, aes(group = NULL), color = "white", size = 4) +
  geom_point(data = trans, aes(group = NULL, alpha = alpha), color = grey(.25), size = 4) +
  # scale_alpha_manual(values = c(.2, .8)) +
  geom_polygon(size = 1, aes(fill = category), alpha = .2) + 
  geom_point(data = contours[, .(x = mean(x), y = mean(y)), by = .(category, overlap)], pch = 4, size = 3, stroke = 1) +
  scale_color_manual(values = c(grey(.6), grey(.1))) +
  scale_fill_manual(values = c(grey(.6), grey(.1))) +
  coord_equal(xlim = c(5, 17), ylim = c(5, 17)) +
  scale_x_continuous(breaks = breaks) +
  scale_y_continuous(breaks = breaks) +
  ggtitle("(b) Information-Integration Category Structures") +
  labs(x = "Feature 1", y = "Feature 2") +
  facet_wrap(~factor(overlap, levels = levels, labels = paste0(str_to_title(levels), " Overlap"))) +
  theme(legend.position = "none")

# ==========================================================================
# Makes the rule-based category structures (rotations of ii) and plot
# ==========================================================================
theta <- 45 * pi/180
t_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
means <- lapply(means, function(i) {(i - grid_mean) %*% t_matrix + grid_mean})
diag(c_cov) <- rev(eigen(c_cov)$values)
c_cov[upper.tri(c_cov) | lower.tri(c_cov)] <- 0

contours <- as.data.table(ldply(means, ellipse, x = cov2cor(c_cov), scale = sqrt(diag(c_cov)), npoints = 1000))
contours[, c("category", "overlap") := data.table(str_split_fixed(.id, "_", 2))]

rb <- ggplot(contours, aes(x, y, group = category)) + 
  geom_vline(xintercept = 11, lty = 2) +
  # geom_point(data = trans, aes(group = NULL, fill = alpha), pch = 21, color = grey(.25), size = 4) +
  # scale_fill_manual(values = c(grey(.25, 0), grey(.25, 1))) +
  # new_scale("fill") +
  geom_point(data = trans, aes(group = NULL), color = "white", size = 4) +
  geom_point(data = trans, aes(group = NULL, alpha = alpha), color = grey(.25), size = 4) +
  # scale_alpha_manual(values = c(.2, .8)) +
  geom_polygon(size = 1, aes(fill = category, color = category), alpha = .2) + 
  geom_point(data = contours[, .(x = mean(x), y = mean(y)), by = .(category, overlap)], aes(color = category), pch = 4, size = 3, stroke = 1) +
  scale_color_manual(values = c(grey(.6), grey(.1))) +
  scale_fill_manual(values = c(grey(.6), grey(.1))) +
  coord_equal(xlim = c(5, 17), ylim = c(5, 17)) +
  scale_x_continuous(breaks = breaks) +
  scale_y_continuous(breaks = breaks) +
  ggtitle("(a) Rule-Based Category Structures") +
  labs(x = "Feature 1", y = "Feature 2") +
  facet_wrap(~factor(overlap, levels = levels, labels = paste0(str_to_title(levels), " Overlap"))) +
  theme(legend.position = "none")

rb / ii 

ggsave("~/Projects/P1/analyses/figures/title case/fig-3.png", height = 6, width = 8)
