# ==========================================================================
# Analysis Simulation Effect of Noise on Categorization
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, purrr, Rsolnp, doRNG, MASS, stringr, patchwork, ggplot2, ggtext, ggh4x, ggrepel)
theme_set(theme_bw())

dropLeadingZero <- function(l){
  str_replace(l, '0(?=.)', '')
}

files <- list.files("~/Projects/P1/data/final/", full.names = T)
files <- grep(files, pattern = "binary", value = T)
dt <- rbindlist(lapply(files, fread))
dt[, noise := factor(noise, levels = c("none", "x", "w", "l", "s"))]
dt[, task := factor(task, levels = c("rb", "ii"))]
dt[, overlap := factor(overlap, levels = c("low", "medium", "high"))]
dt <- dt[order(task, overlap, iter, noise)]

# ==========================================================================
# Parameter recovery analysis
# ==========================================================================
par_dt <- unique(dt[, -c("id", "resps", "true.preds", "fit.preds", "gof")])
par_dt <- par_dt[!(fit.w1 == .5 & fit.l == 1 & fit.tau == 1)] 
par_dt <- melt(par_dt, id.vars = c("task", "overlap", "noise", "iter"), variable.name = "par")
par_dt[, c("true","par") := data.table(str_split_fixed(par, "\\.", 2))]
par_dt <- dcast(par_dt, task + overlap + noise + iter + par ~ true)
par_dt[, par := factor(par, levels = c("w1", "l", "tau", "sigma"))]
par_dt[, range := ifelse(par == "w1", 1, ifelse(par == "sigma", .5, 2))]

gofs <- par_dt[, .(cor = cor(fit, true),
                   MAE = mean(abs(true - fit)/range),
                   MSD = mean((fit - true)/range)), by = .(task, par, noise)]
table_res <- dcast(melt(gofs, id.vars = c("task", "par", "noise")), variable ~ task + noise + par)
xtable::xtable(table_res, type = "latex")

plot_gofs <- melt(gofs, id.vars = c("task", "par", "noise"))
plot_gofs[, label := format(round(value, 2), nsmall = 2)]
plot_gofs[, label := paste0(substr(label, 0, 1), substr(label, 3, 5))]
plot_gofs[, Condition := noise]
levels(plot_gofs$Condition) <- str_to_title(c("no noise", "noisy perception", "noisy attention", "noisy sensitivity", "noisy similarity"))
levels(plot_gofs$variable) <- c("'Correlation'", "'Mean Absolute Error (MAE)'", "'Mean Signed Difference (MSD)'")
levels(plot_gofs$task) <- c("(a) Rule-Based Category Structures", "(b) Information-Integration Category Structures")

size <- 2.7
pos_width <- .9; width <- .75

p_list <- list()
for (i in 1:length(levels(plot_gofs$task))) {
  p_list[[i]] <- ggplot(plot_gofs[task == levels(task)[i]], aes(par, value, fill = Condition, color = Condition)) +
    geom_hline(yintercept = 0) +
    geom_col(position = position_dodge(width = pos_width), linewidth = 1, width = width) +
    geom_text(aes(label=label, vjust = -sign(value) + .5), position = position_dodge(width = pos_width), hjust = .55, size = size, show.legend = F) +
    scale_fill_manual(values = c(alpha("black", .6), viridisLite::viridis(4, end = .8, alpha = .6, option = "E"))) +
    scale_color_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
    scale_x_discrete(name = "Parameter", labels = c(expression(italic("w"["1"]), italic("c"), italic(tau), italic(sigma)))) +
    scale_y_continuous(name = "Value", limits = ~ c(ifelse(min(.x) < 0, -.3, 0), max(ceiling(10*max(.x))/10, .3)),
                       breaks = ~ c(.x[1:2], mean(.x)), expand = c(0, 0)) + 
    labs(title = levels(plot_gofs$task)[i]) +
    coord_cartesian(clip = "off") + theme_classic() + 
    ggh4x::facet_nested_wrap(vars(variable), scales = "free", labeller = label_parsed, nest_line = element_line()) +
    theme(panel.grid = element_blank(), panel.spacing.y = unit(1, "cm"), 
          strip.background = element_blank(), strip.text = element_text(size = 16),
          axis.text = element_text(size = 16), axis.title = element_text(size = 16), title = element_text(size = 18),
          legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.key.size = unit(.4, 'cm'))
  
}
# (p_list[[1]] / plot_spacer() / p_list[[2]]) + plot_layout(heights = c(4, 1 ,4), guides = "collect") & theme(legend.position = "bottom")
wrap_plots(p_list, nrow = 2, guides = "collect") & theme(legend.position = "bottom")
ggsave("~/Projects/P1/analyses/figures/title case/fig-c2.png", width = 12, height = 6.5)

levels(par_dt$noise) <- str_to_title(c("'no noise'", "'noisy perception'", "'noisy attention'", "'noisy sensitivity'", "'noisy similarity'"))
levels(par_dt$par) <- c(expression(paste("i. Attention Weight ", italic("w"["1"]))), 
                        expression(paste("ii. Distance Sensitivity ", italic("c"))),
                        expression(paste("iii. Response Determinism ", italic(tau))),
                        expression(paste("iv. Standard Deviation ", italic(sigma)))) 
par_dt[, c("noise_bin", "where") := tstrsplit(noise, " ")]
par_dt[, noise_bin := ifelse(noise_bin == "'No", "'No Noise'", "'Noise in'")]
par_dt[, where := ifelse(where == "Noise'", "", paste0("'", where))]
par_dt[, where := factor(where, levels = str_to_title(c("'perception'", "'attention'", "'sensitivity'", "'similarity'")))]

p_list <- list()
for (i in 1:length(levels(par_dt$task))) {
  p_list[[i]] <- ggplot(par_dt[task == levels(task)[i]], aes(round(true, 2), round(fit, 2), color = noise)) +
    geom_point(alpha = .1) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    scale_x_continuous(expand = c(.0, .0), n.breaks = 2, labels = dropLeadingZero) +
    scale_y_continuous(expand = c(.0, .0), n.breaks = 2, labels = dropLeadingZero) +
    scale_color_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
    labs(x = "True Value", y = "Estimated Value (No Noise Assumed)", title = levels(plot_gofs$task)[i]) +
    theme_classic() +
    ggh4x::facet_nested_wrap(vars(par, noise_bin, where), ncol = 5, labeller = label_parsed, scales = "free", nest_line = element_line()) +
    theme(panel.grid = element_blank(), legend.position = "none", aspect.ratio = 1,
          strip.background = element_blank(), strip.text = element_text(size = 16),
          panel.spacing.y = unit(1, "cm"))
}
wrap_plots(p_list, ncol = 2) & theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), title = element_text(size = 16))
ggsave("~/Projects/P1/analyses/figures/title case/fig-c1.png", width = 16, height = 12) # 9 or 12 depending on whether sigma is contained

# ==========================================================================
# Probabilistic predictions analysis
# ==========================================================================
preds <- dt[, .(mean = mean(true.preds), sd = sd(true.preds)), by = .(task, overlap, id, noise, iter)]
preds_dt <- preds[, .(1 - mean(mean), median(sd)), by = .(task, id, noise)]
levels(preds_dt$noise) <- str_to_title(c("no noise", "noisy perception", "noisy attention", "noisy sensitivity", "noisy similarity"))

ids <- c(1, 3, 5, 11, 13, 15, 21, 23, 25)
preds_dt <- preds_dt[id %in% ids]

# ==========================================================================
# Rule-based category structures

levels <- c(21, 11, 1, 23, 13, 3, 25, 15, 5)
preds_dt[, id := factor(id, levels = levels[levels %in% ids])]
preds_dt[, lvl := as.numeric(apply(.SD, 1, function(i) {which(i["id"] == levels[levels %in% ids])}))]
labels <- paste0("<img src='~/Projects/P1/analyses/figures/axis-icons/9-", levels[levels %in% ids],  ".png' width='10' /><br>")

p1 <- ggplot(preds_dt[task == "rb" & id %in% ids], aes(id, V1, col = noise, fill = noise, group = noise)) +
  geom_hline(yintercept = .5, lty = 2, color = "darkgrey") +
  geom_point(size = 2) +
  geom_rect(alpha = .2, aes(xmin = lvl - .3, xmax = lvl + .3, ymin = V1 - V2, ymax = V1 + V2)) +
  facet_grid(.~noise) +
  # facet_grid(overlap~noise, labeller = labeller(overlap = c(low = "Low Overlap", medium = "Medium Overlap", high = "High Overlap"))) +
  scale_fill_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  scale_color_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  ggtitle("(a) Rule-Based Category Structures") +
  scale_y_continuous(name = expression(paste("Prediction Pr(", italic("Right"), "|", italic("i"), ")")), limits = c(0, 1), expand = c(0, 0), breaks = c(0, .5, 1)) +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_markdown(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(.5, "lines"))

p1 <- p1 + scale_x_discrete(name = expression(paste("Transfer Stimulus ", italic("i"))), labels = labels)
# ggsave("Projects/P1/analyses/figures/final/fig-predictions-rb-9.png", width = 7, height = 5.5)

# ==========================================================================
# Information-integration category structures

levels <- c(21, 11, 23, 1, 13, 25, 3, 15, 5)
preds_dt[, id := factor(id, levels = levels[levels %in% ids])]
preds_dt[, lvl := as.numeric(apply(.SD, 1, function(i) {which(i["id"] == levels[levels %in% ids])}))]
labels <- paste0("<img src='~/Projects/P1/analyses/figures/axis-icons/9-", levels[levels %in% ids],  ".png' width='10' /><br>")

p2 <- ggplot(preds_dt[task == "ii" & id %in% ids], aes(id, V1, col = noise, fill = noise, group = noise)) +
  geom_hline(yintercept = .5, lty = 2, color = "darkgrey") +
  geom_point(size = 2) +
  geom_rect(alpha = .2, aes(xmin = lvl - .3, xmax = lvl + .3, ymin = V1 - V2, ymax = V1 + V2)) +
  facet_grid(.~noise) +
  # facet_grid(overlap~noise, labeller = labeller(overlap = c(low = "Low Overlap", medium = "Medium Overlap", high = "High Overlap"))) +
  scale_fill_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  scale_color_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  ggtitle("(b) Information-Integration Category Structures") +
  scale_y_continuous(name = expression(paste("Prediction Pr(", italic("Right"), "|", italic("i"), ")")), limits = c(0, 1), expand = c(0, 0), breaks = c(0, .5, 1)) +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_markdown(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(.5, "lines"))

p2 <- p2 + scale_x_discrete(name = expression(paste("Transfer Stimulus ", italic("i"))), labels = labels)
# ggsave("Projects/P1/analyses/figures/final/fig-predictions-ii-9.png", width = 7, height = 5.5)

p1 / p2
# ggsave("~/Projects/P1/analyses/figures/title case/fig-6.png", width = 9, height = 5)
# ggsave("~/Projects/P1/analyses/figures/title case/fig-a2.png", width = 9, height = 8.5)

# ==========================================================================
# Binary responses analysis
# ==========================================================================
resps <- dt[, .(mean = mean(resps)), by = .(task, overlap, noise, iter, id)]
resps_dt <- resps[, .(1 - mean(mean), var(mean)), by = .(task, overlap, noise, id)][, .(mean(V1), mean(V2)), by = .(task, noise, id)]
levels(resps_dt$noise) <- str_to_title(c("no noise", "noisy perception", "noisy attention", "noisy sensitivity", "noisy similarity"))

ids <- c(1, 3, 5, 11, 13, 15, 21, 23, 25)
resps_dt <- resps_dt[id %in% ids]

# ==========================================================================
# Rule-based category structures

levels <- c(21, 11, 1, 23, 13, 3, 25, 15, 5)
resps_dt[, id := factor(id, levels = levels[levels %in% ids])]
resps_dt[, lvl := as.numeric(apply(.SD, 1, function(i) {which(i["id"] == levels[levels %in% ids])}))]
labels <- paste0("<img src='~/Projects/P1/analyses/figures/axis-icons/9-", levels[levels %in% ids],  ".png' width='10' /><br>")

p3 <- ggplot(resps_dt[task == "rb" & id %in% ids], aes(id, V1, col = noise, fill = noise, group = noise)) +
  geom_hline(yintercept = .5, lty = 2, color = "darkgrey") +
  geom_point(size = 2) +
  geom_rect(alpha = .2, aes(xmin = lvl - .3, xmax = lvl + .3, ymin = V1 - V2, ymax = V1 + V2)) +
  facet_grid(.~noise) +
  # facet_grid(overlap~noise, labeller = labeller(overlap = c(low = "Low Overlap", medium = "Medium Overlap", high = "High Overlap"))) +
  scale_fill_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  scale_color_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  ggtitle("(a) Rule-Based Category Structures") +
  scale_y_continuous(name = expression(paste("Category ", italic("Right"), " Responses")), limits = c(0, 1), expand = c(0, 0), breaks = c(0, .5, 1)) +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_markdown(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(.5, "lines"))

p3 <- p3 + scale_x_discrete(name = expression(paste("Transfer Stimulus ", italic("i"))), labels = labels)
# ggsave("Projects/P1/analyses/figures/final/fig-responses-rb-9.png", width = 7, height = 5.5)

# ==========================================================================
# Information-integration category structures

levels <- c(21, 11, 23, 1, 13, 25, 3, 15, 5)
resps_dt[, id := factor(id, levels = levels[levels %in% ids])]
resps_dt[, lvl := as.numeric(apply(.SD, 1, function(i) {which(i["id"] == levels[levels %in% ids])}))]
labels <- paste0("<img src='~/Projects/P1/analyses/figures/axis-icons/9-", levels[levels %in% ids],  ".png' width='10' /><br>")

p4 <- ggplot(resps_dt[task == "ii" & id %in% ids], aes(id, V1, col = noise, fill = noise, group = noise)) +
  geom_hline(yintercept = .5, lty = 2, color = "darkgrey") +
  geom_point(size = 2) +
  geom_rect(alpha = .2, aes(xmin = lvl - .3, xmax = lvl + .3, ymin = V1 - V2, ymax = V1 + V2)) +
  facet_grid(.~noise) +
  # facet_grid(overlap~noise, labeller = labeller(overlap = c(low = "Low Overlap", medium = "Medium Overlap", high = "High Overlap"))) +
  scale_fill_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  scale_color_manual(values = c("black", viridisLite::viridis(4, end = .8, option = "E"))) +
  ggtitle("(b) Information-Integration Category Structures") +
  scale_y_continuous(name = expression(paste("Category ", italic("Right"), " Responses")), limits = c(0, 1), expand = c(0, 0), breaks = c(0, .5, 1)) +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_markdown(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(.5, "lines"))

p4 <- p4 + scale_x_discrete(name = expression(paste("Transfer Stimulus ", italic("i"))), labels = labels)
# ggsave("Projects/P1/analyses/figures/final/fig-responses-ii-9.png", width = 7, height = 5.5)

p3 / p4
# ggsave("~/Projects/P1/analyses/figures/title case/fig-5.png", width = 9, height = 5)
# ggsave("~/Projects/P1/analyses/figures/title case/fig-a1.png", width = 9, height = 8.5)
