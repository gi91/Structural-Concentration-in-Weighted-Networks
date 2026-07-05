# =============================================================================
#  Network Concentration Index (NCI) - Simulation Study
#  Riso & Zoia (2026)
#
#  Implemented indices - exact order of Table 1 in the paper:
#
#   1. Baseline NCI         psi(A)    = w'Aw  / (1 - HHI)
#   2. Density-adjusted     psi_d(A)  = w'Aw  / ((1 - HHI) * delta)
#   3. Null-model NCI       psi_0(A)  = (w'Aw) / (w'Ew)
#   4. Degree-constrained   psi_d(A)  = (w'Aw) / (w'Bw)
#   5. Weighted NCI         psi_w(W)  = w'Ww  / (1 - HHI)
#   6. Transformed-data NCI psi_f(W)  = w'f(W)w / (1 - HHI)   [f = sqrt]
#   7. Multi-layer NCI      psi_ML    = psi(Sum alpha_l A_l)
#
#  Comparisons with HHI:
#   - Dedicated figures: NCI vs HHI as network structure varies
#   - Sensitivity of the indices as concentration varies (alpha -> HHI)
#   - Scatter NCI ~ HHI with random weights (random simplex)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(reshape2)
})

# Resolve namespace conflicts: dplyr takes precedence over MASS, stats, reshape2
select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
rename  <- dplyr::rename

set.seed(42)

# -----------------------------------------------------------------------------
#  GRAPHICAL THEME
# -----------------------------------------------------------------------------

theme_nci <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      strip.background  = element_rect(fill = "grey95", colour = "grey75"),
      strip.text        = element_text(face = "bold", size = base - 1),
      plot.title        = element_text(face = "bold", size = base + 1, hjust = 0),
      plot.subtitle     = element_text(size = base - 2, colour = "grey45"),
      axis.title        = element_text(size = base - 1),
      legend.position   = "bottom",
      legend.key.size   = unit(0.45, "cm"),
      legend.text       = element_text(size = base - 2),
      plot.margin       = margin(8, 12, 8, 8)
    )
}

pal3 <- c(
  "Core-periphery" = "#2166ac",
  "Random"         = "#f4a582",
  "Peripheral"     = "#d6604d"
)

# -----------------------------------------------------------------------------
#  INDEX FUNCTIONS  (Table 1 of the paper only)
# -----------------------------------------------------------------------------

# HHI  Sum w_i^2
hhi <- function(w) sum(w^2)

# Helper: common denominator  1 - HHI
denom <- function(w) {
  d <- 1 - hhi(w)
  if (d <= 0) stop("1 - HHI <= 0: degenerate weights")
  d
}

# 1. Baseline NCI  psi(A) = w'Aw / (1 - HHI)
nci_base <- function(w, A) {
  as.numeric(t(w) %*% A %*% w) / denom(w)
}

# 2. Density-adjusted NCI  psi_delta(A) = w'Aw / ((1 - HHI) * delta)
#    delta = observed density = |E| / (N(N-1))
nci_dens <- function(w, A) {
  N     <- nrow(A)
  delta <- sum(A) / (N * (N - 1))
  if (delta <= 0) return(NA_real_)
  as.numeric(t(w) %*% A %*% w) / (denom(w) * delta)
}

# 3. Weighted NCI  psi_w(W) = w'Ww / (1 - HHI)
#    W: symmetric matrix of continuous weights (diag = 0)
#    DIFFERS from NCI_base because it uses link intensities gamma_ij in (0,1]
#    instead of simple presence/absence (0/1).
nci_weighted <- function(w, W) {
  as.numeric(t(w) %*% W %*% w) / denom(w)
}

# 4. Transformed-data NCI  psi_f(W) = w'f(W)w / (1 - HHI)
#    f applied entry-wise to the continuous matrix W.
#    Default f = sqrt: compresses high intensities  -> psi_f < psi_w  when gamma_ij < 1
#    WARNING: always pass a continuous W, NOT a binary A.
#    If a binary A were passed: f(0)=0, f(1)=1 -> psi_f = psi_base  (collapse!)
nci_trans <- function(w, W, f = sqrt) {
  M <- f(W)                          # entry-wise transformation
  as.numeric(t(w) %*% M %*% w) / denom(w)
}

# 5. Multi-layer NCI  psi_ML = psi_base applied to A_comb = Sum alpha_l A_l
nci_multi <- function(w, layers, alpha) {
  stopifnot(length(layers) == length(alpha),
            abs(sum(alpha) - 1) < 1e-9)          # alpha_l sum to 1
  A_comb <- Reduce("+", Map("*", alpha, layers))
  nci_base(w, A_comb)
}

# 6. Null-model NCI  psi_0(A) = (w'Aw) / (w'Ew)
#    E: expected Erdos-Renyi matrix with parameter p
#    Property: E[psi_0] = 1 for any p
nci_null <- function(w, A, p) {
  N   <- nrow(A)
  E   <- matrix(p, N, N); diag(E) <- 0
  den <- as.numeric(t(w) %*% E %*% w)
  if (den <= 0) return(NA_real_)
  as.numeric(t(w) %*% A %*% w) / den
}

# 7. Degree-constrained NCI  psi_d(A) = (w'Aw) / (w'Bw)
#    B: configuration model  B_ij = d_i d_j / (2m)
nci_deg <- function(w, A) {
  d_vec <- rowSums(A)
  m     <- sum(d_vec)
  if (m == 0) return(NA_real_)
  B <- outer(d_vec, d_vec) / m; diag(B) <- 0
  den <- as.numeric(t(w) %*% B %*% w)
  if (den <= 0) return(NA_real_)
  as.numeric(t(w) %*% A %*% w) / den
}

# --- Labels (Table 1 order) --------------------------------------------------
# EXACT order of Table 1 in the paper:
#   Baseline -> Density-adjusted -> Null-model -> Degree-constrained ->
#   Weighted -> Transformed-data -> Multi-layer
INDEX_NAMES <- c(
  "NCI_base"  = "1. Baseline NCI",
  "NCI_dens"  = "2. Density-adjusted",
  "NCI_null"  = "3. Null-model",
  "NCI_deg"   = "4. Degree-constrained NCI",
  "NCI_wgt"   = "5. Weighted NCI",
  "NCI_trans" = "6. Transformed-data NCI",
  "NCI_multi" = "7. Multi-layer NCI"
)
IDX <- names(INDEX_NAMES)   # order automatically propagated to ALL figures

# -----------------------------------------------------------------------------
#  GLOBAL PARAMETERS
# -----------------------------------------------------------------------------

N        <- 10
w_fixed  <- c(0.30, 0.20, 0.15, 0.10, 0.08, 0.06, 0.04, 0.03, 0.02, 0.02)
stopifnot(abs(sum(w_fixed) - 1) < 1e-9)

cat(sprintf("Fixed weights: HHI = %.4f  |  1-HHI = %.4f\n\n",
            hhi(w_fixed), 1 - hhi(w_fixed)))

# -----------------------------------------------------------------------------
#  NETWORK GENERATORS
# -----------------------------------------------------------------------------

# Erdos-Renyi network with probability p
gen_er <- function(N, p) {
  A <- matrix(0L, N, N)
  for (i in seq_len(N - 1))
    for (j in (i + 1):N)
      if (runif(1) < p) { A[i, j] <- 1L; A[j, i] <- 1L }
  A
}

# Core-periphery network: dense core + sparse periphery
gen_core <- function(N, p_core = 0.9, p_periph = 0.1, core_size = 4) {
  A <- matrix(0L, N, N)
  for (i in seq_len(core_size - 1))
    for (j in (i + 1):core_size)
      if (runif(1) < p_core) { A[i, j] <- 1L; A[j, i] <- 1L }
  for (i in (core_size + 1):(N - 1))
    for (j in (i + 1):N)
      if (runif(1) < p_periph) { A[i, j] <- 1L; A[j, i] <- 1L }
  A
}

# Peripheral network (core with reversed probabilities)
gen_periph <- function(N) gen_core(N, p_core = 0.1, p_periph = 0.7, core_size = 4)

# Continuous matrix W built on A (weights on existing edges)
make_W <- function(A, lo = 0.2, hi = 1.0) {
  N <- nrow(A)
  raw <- A * matrix(runif(N * N, lo, hi), N, N)
  W   <- (raw + t(raw)) / 2
  diag(W) <- 0
  W
}

# Two layers for the multi-layer index: layer 2 = attenuated version of A
make_layers <- function(A) {
  N  <- nrow(A)
  L2 <- A * matrix(runif(N * N, 0, 0.5), N, N)
  list(A, (L2 + t(L2)) / 2)
}

# -----------------------------------------------------------------------------
#  HELPER FUNCTION: computes all 7 indices on a single network
# -----------------------------------------------------------------------------

calc_all <- function(w, A, p_ref) {
  W      <- make_W(A)
  layers <- make_layers(A)
  c(
    NCI_base  = nci_base(w, A),
    NCI_dens  = nci_dens(w, A),
    NCI_wgt   = nci_weighted(w, W),
    NCI_trans = nci_trans(w, W),          # continuous W -> different from NCI_wgt
    NCI_multi = nci_multi(w, layers, c(0.6, 0.4)),
    NCI_null  = nci_null(w, A, p_ref),
    NCI_deg   = nci_deg(w, A)
  )
}

# -----------------------------------------------------------------------------
#  PART A - DETERMINISTIC SCENARIOS (fixed weights, different structures)
# -----------------------------------------------------------------------------

# Scenario 1: Core-periphery (nodes 1-4 dense, 5-10 sparse)
A_core <- matrix(c(
  0,1,1,1,0,0,0,0,0,0,
  1,0,1,1,0,0,0,0,0,0,
  1,1,0,1,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,
  0,0,0,0,0,1,1,0,0,0,
  0,0,0,0,1,0,1,0,0,0,
  0,0,0,0,1,1,0,0,0,0,
  0,0,0,0,0,0,0,0,1,1,
  0,0,0,0,0,0,0,1,0,1,
  0,0,0,0,0,0,0,1,1,0
), 10, 10, byrow = TRUE)

# Scenario 2: Peripheral (connectivity concentrated on low-weight nodes)
A_periph <- matrix(c(
  0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,1,1,0,0,0,
  0,0,0,1,1,0,1,1,0,0,
  0,0,1,0,1,1,0,1,0,0,
  0,1,0,0,0,1,1,0,1,1,
  1,0,0,0,0,0,0,1,0,1,
  0,0,0,0,0,0,0,1,1,0
), 10, 10, byrow = TRUE)

# Scenario 3: Random (comparable density)
A_rand <- matrix(c(
  0,1,0,0,0,0,1,0,0,0,
  1,0,0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,1,0,0,0,
  0,0,1,0,0,0,0,1,0,0,
  0,0,0,0,0,1,0,1,0,0,
  0,1,0,0,1,0,0,0,1,0,
  1,0,1,0,0,0,0,0,1,0,
  0,0,0,1,1,0,0,0,0,1,
  0,0,0,0,0,1,1,0,0,1,
  0,0,0,0,0,0,0,1,1,0
), 10, 10, byrow = TRUE)

p_ref <- sum(A_core) / (N * (N - 1))   # reference density

# Detail table
det_tab <- data.frame(
  Scenario = c("Core-periphery", "Peripheral", "Random"),
  HHI      = hhi(w_fixed),
  Density  = c(sum(A_core), sum(A_periph), sum(A_rand)) / (N * (N - 1))
)

set.seed(42)
vals <- rbind(
  calc_all(w_fixed, A_core,   p_ref),
  calc_all(w_fixed, A_periph, p_ref),
  calc_all(w_fixed, A_rand,   p_ref)
)
det_tab <- cbind(det_tab, round(vals, 4))
cat("=== Table A - Deterministic scenarios ===\n")
print(det_tab, row.names = FALSE)

# Long format for ggplot
det_long <- det_tab |>
  dplyr::select(Scenario, all_of(IDX)) |>
  pivot_longer(all_of(IDX), names_to = "Index", values_to = "Value") |>
  mutate(
    Label    = INDEX_NAMES[Index],
    Label    = factor(Label, levels = rev(unname(INDEX_NAMES))),
    Scenario = factor(Scenario, levels = c("Core-periphery", "Peripheral", "Random"))
  )

# -- Figure A: barchart of all indices by scenario ----------------------------

fig_A <- ggplot(det_long, aes(x = Value, y = Label, fill = Scenario)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65,
           colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            position = position_dodge(width = 0.72),
            hjust = -0.08, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = pal3) +
  scale_x_continuous(
    limits = c(0, 2),
    breaks = seq(0, 2, 0.5)
  ) +
  labs(
    title = "",
    x = sprintf("Fixed weights, HHI = %.4f | network structure varies",
                hhi(w_fixed)),
    y = NULL,
    fill = ""
  ) +
  theme_nci() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),   # <-- force bold legend
    legend.title = element_text(face = "bold")
  )

print(fig_A)


ggsave("fig_A.jpg", fig_A, width = 7, height = 5, dpi = 600, bg = "white")


# =============================================================================
#  NETWORK PLOTS - final version
#  - plot_net_disconnected : for networks with multiple components (A_core)
#  - plot_net_connected    : for connected networks (A_periph, A_rand)
# =============================================================================

library(igraph); library(scales); library(ggplot2)

# -----------------------------------------------------------------------------
#  COMMON HELPER: builds df_nodes and df_edges given a layout
# -----------------------------------------------------------------------------

make_dfs <- function(g, lay, w, size_range = c(5, 14)) {
  N        <- igraph::vcount(g)
  df_nodes <- data.frame(
    x      = lay[, 1],
    y      = lay[, 2],
    label  = paste0("N", seq_len(N)),
    weight = w,
    size   = scales::rescale(w, to = size_range)
  )
  el <- igraph::as_edgelist(g, names = FALSE)
  df_edges <- if (nrow(el) == 0) {
    data.frame(x=numeric(0), y=numeric(0), xend=numeric(0), yend=numeric(0))
  } else {
    data.frame(
      x    = lay[el[,1], 1], y    = lay[el[,1], 2],
      xend = lay[el[,2], 1], yend = lay[el[,2], 2]
    )
  }
  list(nodes = df_nodes, edges = df_edges)
}

# -----------------------------------------------------------------------------
#  COMMON HELPER: adds label-offset columns + draws the ggplot
# -----------------------------------------------------------------------------

draw_net <- function(df_nodes, df_edges, col_vec,
                     thresh = 9, off = 0.06,
                     cx_ref = NULL, cy_ref = NULL) {
  
  # Radial direction for outside labels
  if (is.null(cx_ref)) cx_ref <- rep(mean(df_nodes$x), nrow(df_nodes))
  if (is.null(cy_ref)) cy_ref <- rep(mean(df_nodes$y), nrow(df_nodes))
  
  df_nodes$outside <- df_nodes$size < thresh
  dx <- df_nodes$x - cx_ref
  dy <- df_nodes$y - cy_ref
  nm <- pmax(sqrt(dx^2 + dy^2), 0.01)
  
  df_nodes$lx <- ifelse(df_nodes$outside, df_nodes$x + dx/nm * off, df_nodes$x)
  df_nodes$ly <- ifelse(df_nodes$outside, df_nodes$y + dy/nm * off, df_nodes$y)
  
  ggplot() +
    
    geom_segment(
      data = df_edges,
      aes(x=x, y=y, xend=xend, yend=yend),
      colour = "grey50", linewidth = 0.75, lineend = "round"
    ) +
    
    geom_point(
      data = df_nodes, aes(x=x, y=y, size = size + 2.5),
      colour = "grey12", show.legend = FALSE
    ) +
    
    geom_point(
      data = df_nodes, aes(x=x, y=y, size = size),
      colour = col_vec, alpha = 0.93, show.legend = FALSE
    ) +
    
    geom_text(
      data = df_nodes[!df_nodes$outside, ],
      aes(x=x, y=y, label=label),
      colour = "white", fontface = "bold", size = 3.3,
      vjust = 0.5, hjust = 0.5
    ) +
    
    geom_text(
      data = df_nodes[df_nodes$outside, ],
      aes(x=lx, y=ly, label=label),
      colour = "grey10", fontface = "bold", size = 3.0,
      vjust = 0.5, hjust = 0.5
    ) +
    
    scale_size_identity() +
    scale_x_continuous(limits = c(0, 1), expand = c(0.08, 0.08)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.08, 0.08)) +
    coord_fixed(ratio = 1, clip = "off") +
    theme_void() +
    theme(
      plot.margin     = margin(4, 4, 4, 4),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# -----------------------------------------------------------------------------
#  1. DISCONNECTED NETWORK  (A_core: 3 components)
# -----------------------------------------------------------------------------

plot_net_disconnected <- function(A, w, node_col, seed = 42) {
  
  N    <- nrow(A)
  g    <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)
  set.seed(seed)
  
  zone_centers <- list(
    c(0.28, 0.50, 0.38),   # largest component (4-node core)
    c(0.78, 0.74, 0.14),   # periphery cluster 1
    c(0.78, 0.26, 0.14)    # periphery cluster 2
  )
  
  comp_order <- order(table(comp$membership), decreasing = TRUE)
  lay        <- matrix(0, N, 2)
  cx_ref     <- numeric(N)
  cy_ref     <- numeric(N)
  
  for (rank in seq_along(comp_order)) {
    k   <- comp_order[rank]
    idx <- which(comp$membership == k)
    sg  <- igraph::induced_subgraph(g, idx)
    zc  <- zone_centers[[rank]]
    
    if (length(idx) == 1) {
      sub_lay <- matrix(c(0.5, 0.5), 1, 2)
    } else {
      sub_lay <- igraph::layout_in_circle(sg)
      for (d in 1:2) {
        rng <- range(sub_lay[, d])
        if (diff(rng) > 1e-9)
          sub_lay[, d] <- (sub_lay[, d] - mean(rng)) / (diff(rng) / 2)
      }
      sub_lay[, 1] <- sub_lay[, 1] * zc[3] + zc[1]
      sub_lay[, 2] <- sub_lay[, 2] * zc[3] + zc[2]
    }
    
    lay[idx, ]    <- sub_lay
    cx_ref[idx]   <- zc[1]
    cy_ref[idx]   <- zc[2]
  }
  
  dfs <- make_dfs(g, lay, w)
  draw_net(dfs$nodes, dfs$edges,
           col_vec  = rep(node_col, N),
           cx_ref   = cx_ref,
           cy_ref   = cy_ref)
}

# -----------------------------------------------------------------------------
#  2. CONNECTED NETWORK  (A_periph, A_rand)
#     KK layout normalized to [0,1]x[0,1] with generous padding
# -----------------------------------------------------------------------------
plot_net_disconnected <- function(A, w, node_col, seed = 42,
                                  push_away = NULL) {
  N <- nrow(A)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  set.seed(seed)
  
  # Manual layout: core radius increased to 0.18 -> N3 always visible
  r_core   <- 0.18
  r_periph <- 0.12
  
  lay    <- matrix(0, N, 2)
  cx_ref <- numeric(N)
  cy_ref <- numeric(N)
  
  # Core (N1-N4): circle with generous radius
  angles_core <- seq(pi/4, 2*pi + pi/4, length.out = 5)[1:4]
  lay[1:4, 1] <- 0.28 + r_core * cos(angles_core)
  lay[1:4, 2] <- 0.50 + r_core * sin(angles_core)
  cx_ref[1:4] <- 0.28
  cy_ref[1:4] <- 0.50
  
  # Periphery cluster 1 (N5-N7)
  angles_p1 <- seq(pi/6, 2*pi + pi/6, length.out = 4)[1:3]
  lay[5:7, 1] <- 0.78 + r_periph * cos(angles_p1)
  lay[5:7, 2] <- 0.74 + r_periph * sin(angles_p1)
  cx_ref[5:7] <- 0.78
  cy_ref[5:7] <- 0.74
  
  # Periphery cluster 2 (N8-N10)
  angles_p2 <- seq(-pi/6, 2*pi - pi/6, length.out = 4)[1:3]
  lay[8:10, 1] <- 0.78 + r_periph * cos(angles_p2)
  lay[8:10, 2] <- 0.26 + r_periph * sin(angles_p2)
  cx_ref[8:10] <- 0.78
  cy_ref[8:10] <- 0.26
  
  # Optional push_away (compatibility with existing calls)
  if (!is.null(push_away)) {
    for (node_name in names(push_away)) {
      ni  <- as.integer(node_name)
      fac <- push_away[[node_name]]$factor
      dx  <- lay[ni, 1] - cx_ref[ni]
      dy  <- lay[ni, 2] - cy_ref[ni]
      nm  <- sqrt(dx^2 + dy^2)
      if (nm > 1e-9) {
        lay[ni, 1] <- cx_ref[ni] + dx * fac
        lay[ni, 2] <- cy_ref[ni] + dy * fac
      }
    }
  }
  
  dfs <- make_dfs(g, lay, w)
  draw_net(dfs$nodes, dfs$edges,
           col_vec = rep(node_col, N),
           cx_ref  = cx_ref,
           cy_ref  = cy_ref)
}

plot_net_connected <- function(A, w, node_col, seed = 42) {
  
  N <- nrow(A)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  set.seed(seed)
  
  lay <- igraph::layout_with_kk(g)
  
  # Normalize to [0.08, 0.92] - generous padding avoids label clipping
  for (d in 1:2) {
    rng <- range(lay[, d])
    if (diff(rng) > 1e-9)
      lay[, d] <- (lay[, d] - rng[1]) / diff(rng) * 0.84 + 0.08
  }
  
  dfs <- make_dfs(g, lay, w)
  draw_net(dfs$nodes, dfs$edges,
           col_vec = rep(node_col, N))
}
# =============================================================================
#  RENDERING
# =============================================================================


fig_core <- plot_net_disconnected(
  A         = A_core,
  w         = w_fixed,
  node_col  = pal3["Core-periphery"],
  seed      = 42,
  push_away = list(
    "7" = list(factor = 0.3),
    "9" = list(factor = 0.3)
  )
)

fig_periph <- plot_net_connected   (A_periph, w_fixed, pal3["Peripheral"])
fig_rand   <- plot_net_connected   (A_rand,   w_fixed, pal3["Random"])

print(fig_core)
print(fig_periph)
print(fig_rand)

ggsave("fig_core.pdf",   fig_core,   width=5, height=5, device=cairo_pdf)
ggsave("fig_periph.pdf", fig_periph, width=5, height=5, device=cairo_pdf)
ggsave("fig_rand.pdf",   fig_rand,   width=5, height=5, device=cairo_pdf)

ggsave("fig_core.png",   fig_core,   width=5, height=5, dpi=600, bg="white")
ggsave("fig_periph.png", fig_periph, width=5, height=5, dpi=600, bg="white")
ggsave("fig_rand.png",   fig_rand,   width=5, height=5, dpi=600, bg="white")


# -----------------------------------------------------------------------------
#  PART B - MONTE CARLO  (R = 5000 per mechanism)
# -----------------------------------------------------------------------------

R_sim <- 5000
cat(sprintf("\nMonte Carlo: R = %d per mechanism...\n", R_sim))

sim_one_structure <- function(gen_fun, R, w, label, p_ref) {
  rows <- vector("list", R)
  for (k in seq_len(R)) {
    A <- gen_fun()
    v <- tryCatch(calc_all(w, A, p_ref), error = function(e) rep(NA_real_, 7))
    rows[[k]] <- as.list(v)
  }
  bind_rows(rows) |>
    mutate(Structure = label, HHI_w = hhi(w))
}

set.seed(42)
mc_core   <- sim_one_structure(function() gen_core(N),         R_sim, w_fixed, "Core-periphery", p_ref)
mc_rand   <- sim_one_structure(function() gen_er(N, p_ref),    R_sim, w_fixed, "Random",          p_ref)
mc_periph <- sim_one_structure(function() gen_periph(N),       R_sim, w_fixed, "Peripheral",      p_ref)

mc_all <- bind_rows(mc_core, mc_rand, mc_periph) |>
  mutate(Structure = factor(Structure,
                            levels = c("Core-periphery", "Random", "Peripheral")))

# Check difference NCI_wgt vs NCI_trans (if = 1.000 -> still collapsed)
cat(sprintf("Correlation NCI_wgt ~ NCI_trans: %.4f  [must be < 1]\n",
            cor(mc_all$NCI_wgt, mc_all$NCI_trans, use = "pairwise")))

# Long format
mc_long <- mc_all |>
  pivot_longer(all_of(IDX), names_to = "Index", values_to = "Value") |>
  mutate(
    Label = INDEX_NAMES[Index],
    Label = factor(Label, levels = unname(INDEX_NAMES))
  )

# -- Figure 1: MC distributions (kernel densities) ----------------------------
fig1 <- ggplot(mc_long, aes(x = Value, fill = Structure, colour = Structure)) +
  geom_density(alpha = 0.28, linewidth = 0.55) +
  facet_wrap(~ Label, scales = "free", ncol = 4) +
  scale_fill_manual(values = pal3) +
  scale_colour_manual(values = pal3) +
  labs(
    title    = "",
    x = "Index value", y = "Density",
    fill = "Structure", colour = "Structure"
  ) +
  theme_nci(base = 9) +
  theme(aspect.ratio = 0.75)

print(fig1)

ggsave("Density.jpg", fig1, width = 7, height = 5, dpi = 600, bg = "white")


# -- Figure 2: boxplots --------------------------------------------------------
fig2 <- ggplot(mc_long, aes(x = Structure, y = Value, fill = Structure)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25, linewidth = 0.35) +
  facet_wrap(~ Label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pal3) +
  labs(
    title = "Figure 2. Boxplots of NCI variants by network structure",
    x = NULL, y = "Index value", fill = "Structure"
  ) +
  theme_nci(base = 9) +
  theme(
    aspect.ratio = 0.75,
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

print(fig2)

# -----------------------------------------------------------------------------
#  PART C - VERIFICATION OF THE PAPER'S PROPOSITIONS
# -----------------------------------------------------------------------------

p_grid  <- seq(0.05, 0.95, by = 0.05)
R_prop  <- 5000

# -- Figure 3: Proposition P4 - E[NCI_base] = p under Erdos-Renyi -------------
cat("\nVerification P4: E[NCI_base] = p ...\n")
prop4_df <- bind_rows(lapply(p_grid, function(p) {
  vals <- replicate(R_prop, {
    A <- gen_er(N, p)
    tryCatch(nci_base(w_fixed, A), error = function(e) NA_real_)
  })
  data.frame(p    = p,
             mu   = mean(vals, na.rm = TRUE),
             se   = sd(vals,   na.rm = TRUE) / sqrt(sum(!is.na(vals))))
}))
cat(sprintf("  Max |E[NCI] - p| = %.5f\n",
            max(abs(prop4_df$mu - prop4_df$p), na.rm = TRUE)))


fig3 <- ggplot(prop4_df, aes(x = p, y = mu)) +
  geom_ribbon(aes(ymin = mu - 1.96 * se, ymax = mu + 1.96 * se),
              fill = "#2166ac", alpha = 0.18) +
  geom_line(colour = "#2166ac", linewidth = 0.9) +
  geom_point(colour = "#2166ac", size = 2.0) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  labs(
    x = expression(bold(p)),
    y = expression(bold(E*"["*psi(omega, A)*"]")),
    title = ""
  ) +
  coord_fixed(ratio = 1) +
  theme_nci() +
  theme(
    axis.title.x = element_text(face = "bold",size=20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x  = element_text(face = "bold",size=18),
    axis.text.y  = element_text(face = "bold",size=18)
  )

print(fig3)


ggsave("Pro4.jpg", fig3, width = 7, height = 5, dpi = 600, bg = "white")


# -- Figure 4: E[NCI_null] = 1 for every p ------------------------------------
cat("Verification: E[NCI_null] = 1 ...\n")
prop_null_df <- bind_rows(lapply(p_grid, function(p) {
  vals <- replicate(R_prop, {
    A <- gen_er(N, p)
    tryCatch(nci_null(w_fixed, A, p), error = function(e) NA_real_)
  })
  data.frame(p  = p,
             mu = mean(vals, na.rm = TRUE),
             se = sd(vals,   na.rm = TRUE) / sqrt(sum(!is.na(vals))))
}))

fig4 <- ggplot(prop_null_df, aes(x = p, y = mu)) +
  geom_ribbon(aes(ymin = mu - 1.96 * se, ymax = mu + 1.96 * se),
              fill = "#d6604d", alpha = 0.18) +
  geom_line(colour = "#d6604d", linewidth = 0.9) +
  geom_point(colour = "#d6604d", size = 2.0) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  annotate("text", x = 0.82, y = 1.04,
           label = "E[NCI_null] = 1  (theoretical)",
           colour = "grey40", size = 3.2) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0.7, 1.3), breaks = seq(0.7, 1.3, 0.1)) +
  labs(
    title    = "Figure 4. E[NCI_null] = 1 under Erdos-Renyi (any p)",
    subtitle = sprintf("%d replications per p  |  band: +/-1.96 SE", R_prop),
    x = "Link probability p", y = "E [ NCI_null ]"
  ) +
  theme_nci()

print(fig4)

# -- Figure 5: NCI_dens = 1 with uniform weights ------------------------------
# (Proposition P2: with w = 1/N, NCI_dens equals 1 by construction
#  because psi_delta = [w'Aw / (1-HHI)] / delta and, with uniform w,
#  w'Aw = delta(1-1/N))
cat("Verification P2: NCI_dens = 1 with uniform weights ...\n")
w_unif <- rep(1 / N, N)

prop2_df <- bind_rows(lapply(p_grid, function(p) {
  vals <- replicate(R_prop, {
    A <- gen_er(N, p)
    tryCatch(nci_dens(w_unif, A), error = function(e) NA_real_)
  })
  data.frame(p  = p,
             mu = mean(vals, na.rm = TRUE),
             se = sd(vals,   na.rm = TRUE) / sqrt(sum(!is.na(vals))))
}))
cat(sprintf("  Max |E[NCI_dens, w=unif] - 1| = %.5f\n",
            max(abs(prop2_df$mu - 1), na.rm = TRUE)))

fig5 <- ggplot(prop2_df, aes(x = p, y = mu)) +
  geom_ribbon(aes(ymin = mu - 1.96 * se, ymax = mu + 1.96 * se),
              fill = "#4dac26", alpha = 0.18) +
  geom_line(colour = "#4dac26", linewidth = 0.9) +
  geom_point(colour = "#4dac26", size = 2.0) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  annotate("text", x = 0.75, y = 1.06,
           label = "E[NCI_dens] = 1  (theoretical)",
           colour = "grey40", size = 3.2) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0.7, 1.3), breaks = seq(0.7, 1.3, 0.1)) +
  labs(
    title    = "Figure 5. Proposition P2: NCI_dens = 1 with uniform weights",
    subtitle = sprintf("w = 1/N  |  %d replications per p  |  band: +/-1.96 SE", R_prop),
    x = "Link probability p", y = "E [ NCI_dens ]"
  ) +
  theme_nci()

print(fig5)

# -----------------------------------------------------------------------------
#  PART D - NCI vs HHI COMPARISONS
# -----------------------------------------------------------------------------

# -- Figure 6: Sensitivity to the leading node weight -> HHI varies -----------
# w(alpha) = (alpha, (1-alpha)/(N-1), ...)  with alpha in [0.10, 0.90]
# All NCI curves + HHI on the same network A_core

alpha_grid <- seq(0.10, 0.90, by = 0.02)
set.seed(42)

sens_df <- bind_rows(lapply(alpha_grid, function(a) {
  w_a    <- c(a, rep((1 - a) / (N - 1), N - 1))
  W_a    <- make_W(A_core)          # W fixed for a clean comparison
  layers <- make_layers(A_core)
  data.frame(
    alpha     = a,
    HHI       = hhi(w_a),
    NCI_base  = nci_base(w_a, A_core),
    NCI_dens  = nci_dens(w_a, A_core),
    NCI_wgt   = nci_weighted(w_a, W_a),
    NCI_trans = nci_trans(w_a, W_a),
    NCI_multi = nci_multi(w_a, layers, c(0.6, 0.4)),
    NCI_null  = nci_null(w_a, A_core, p_ref),
    NCI_deg   = nci_deg(w_a, A_core)
  )
}))

# Distinct colours for 8 series
# Colours in the exact order of Table 1 (used in fig6 sensitivity)
pal8 <- c(
  "HHI"                      = "grey30",
  "1. Baseline NCI"          = "#2166ac",
  "2. Density-adjusted"      = "#4dac26",
  "3. Null-model NCI"        = "#d6604d",
  "4. Degree-constrained NCI"= "#b35806",
  "5. Weighted NCI"          = "#7b3294",
  "6. Transformed-data NCI"  = "#e66101",
  "7. Multi-layer NCI"       = "#1a9850"
)

sens_long <- sens_df |>
  pivot_longer(-alpha, names_to = "key", values_to = "Value") |>
  mutate(
    Label    = case_when(
      key == "HHI"      ~ "HHI",
      TRUE              ~ INDEX_NAMES[key]
    ),
    Label    = factor(Label, levels = names(pal8)),
    LineType = ifelse(key == "HHI", "dashed", "solid")
  )

fig6 <- ggplot(sens_long,
               aes(x = alpha, y = Value,
                   colour = Label, linetype = LineType)) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = pal8) +
  scale_linetype_identity() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  labs(
    title    = "Figure 6. Sensitivity of NCI variants to leading weight alpha (and HHI)",
    subtitle = sprintf(
      "Core-periphery network A_core  |  w1 = alpha,  w2..%d = (1-alpha)/%d",
      N, N - 1),
    x = "Leading node weight alpha", y = "Index value",
    colour = NULL
  ) +
  guides(linetype = "none") +
  theme_nci() +
  theme(legend.position = "right",
        legend.text      = element_text(size = 8))

print(fig6)

# -- Figure 7: scatter NCI ~ HHI with random weights (all three mechanisms) ---
# For each draw: random weights from the simplex + random network

cat("\nScatter NCI vs HHI (random weights)...\n")
set.seed(42)
R_scatter <- 800

hhi_scatter <- bind_rows(lapply(
  c("Core-periphery", "Random", "Peripheral"),
  function(lab) {
    gen <- switch(lab,
                  "Core-periphery" = function() gen_core(N),
                  "Random"         = function() gen_er(N, p_ref),
                  "Peripheral"     = function() gen_periph(N)
    )
    replicate(R_scatter, {
      A   <- gen()
      w_r <- diff(c(0, sort(runif(N - 1)), 1))   # random weights on the simplex
      W_r <- make_W(A)
      lyr <- make_layers(A)
      data.frame(
        Structure = lab,
        HHI       = hhi(w_r),
        NCI_base  = tryCatch(nci_base(w_r, A),             error = function(e) NA_real_),
        NCI_dens  = tryCatch(nci_dens(w_r, A),             error = function(e) NA_real_),
        NCI_wgt   = tryCatch(nci_weighted(w_r, W_r),       error = function(e) NA_real_),
        NCI_trans = tryCatch(nci_trans(w_r, W_r),          error = function(e) NA_real_),
        NCI_multi = tryCatch(nci_multi(w_r, lyr, c(.6,.4)),error = function(e) NA_real_),
        NCI_null  = tryCatch(nci_null(w_r, A, p_ref),      error = function(e) NA_real_),
        NCI_deg   = tryCatch(nci_deg(w_r, A),              error = function(e) NA_real_)
      )
    }, simplify = FALSE) |> bind_rows()
  }
))

hhi_long <- hhi_scatter |>
  pivot_longer(all_of(IDX), names_to = "Index", values_to = "NCI") |>
  mutate(
    Label     = INDEX_NAMES[Index],
    Label     = factor(Label, levels = unname(INDEX_NAMES)),
    Structure = factor(Structure,
                       levels = c("Core-periphery", "Random", "Peripheral"))
  )

fig7 <- ggplot(hhi_long, aes(x = HHI, y = NCI, colour = Structure)) +
  geom_point(alpha = 0.08, size = 0.5, shape = 16) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.5) +
  facet_wrap(~ Label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = pal3) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    title    = "",
    x = "HHI", y = "NCI value", colour = "Structure"
  ) +
  theme_nci(base = 9) +
  theme(aspect.ratio = 0.78)

print(fig7)

ggsave("HHIvsNCI.jpg", fig7, width = 7, height = 6, dpi = 600, bg = "white")


# -- Figure 8: NCI ~ HHI correlation by index (Pearson coefficient) -----------
cor_hhi <- hhi_scatter |>
  group_by(Structure) |>
  summarise(
    across(all_of(IDX),
           ~ cor(.x, HHI, use = "pairwise.complete.obs"),
           .names = "{.col}"),
    .groups = "drop"
  ) |>
  pivot_longer(all_of(IDX), names_to = "Index", values_to = "r") |>
  mutate(
    Label     = INDEX_NAMES[Index],
    Label     = factor(Label, levels = unname(INDEX_NAMES)),
    Structure = factor(Structure,
                       levels = c("Core-periphery", "Random", "Peripheral"))
  )

fig8 <- ggplot(cor_hhi, aes(x = r, y = Label, fill = Structure)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65,
           colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", r)),
            position = position_dodge(width = 0.72),
            hjust = -0.1, size = 2.7, fontface = "bold") +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
  scale_fill_manual(values = pal3) +
  scale_x_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, 0.1)) +
  labs(
    title    = "Figure 8. Pearson correlation between NCI variants and HHI",
    subtitle = "By network structure - random weights on simplex",
    x = "Pearson r (NCI, HHI)", y = NULL, fill = "Structure"
  ) +
  theme_nci() +
  theme(panel.grid.major.y = element_blank())

print(fig8)

# -- Figure 9: inter-index correlation heatmap --------------------------------
cor_mat <- cor(mc_all[, IDX], use = "pairwise.complete.obs")
colnames(cor_mat) <- rownames(cor_mat) <- unname(INDEX_NAMES)

cor_long_mat <- reshape2::melt(cor_mat, varnames = c("X", "Y"), value.name = "r")

fig9 <- ggplot(cor_long_mat, aes(X, Y, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), 
            size = 2.7, colour = "grey10", fontface = "bold") +
  scale_fill_gradient2(
    low = "#d6604d", mid = "white", high = "#2166ac",
    midpoint = 0, limits = c(-1, 1), name = "r"
  ) +
  labs(
    title    = "",
    subtitle = "",
    x = NULL, y = NULL
  ) +
  theme_nci(base = 9) +
  theme(
    aspect.ratio = 1,
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    legend.position = "none"
  )

print(fig9)
ggsave("Corr.jpg", fig9, width = 6, height = 6, dpi = 600, bg = "white")


# -----------------------------------------------------------------------------
#  SAVING
# -----------------------------------------------------------------------------

out <- "/mnt/user-data/outputs"
dir.create(out, showWarnings = FALSE)

ggsave(file.path(out, "fig0_networks.pdf"),        fig0,  width = 13,  height = 5,    dpi = 300)
ggsave(file.path(out, "figA_deterministic.pdf"),   fig_A, width = 11,  height = 6.5,  dpi = 300)
ggsave(file.path(out, "fig1_mc_density.pdf"),      fig1,  width = 14,  height = 8.5,  dpi = 300)
ggsave(file.path(out, "fig2_mc_boxplot.pdf"),      fig2,  width = 14,  height = 8.5,  dpi = 300)
ggsave(file.path(out, "fig3_prop4_enci_eq_p.pdf"), fig3,  width = 6,   height = 6,    dpi = 300)
ggsave(file.path(out, "fig4_enull_eq_1.pdf"),      fig4,  width = 7,   height = 4.5,  dpi = 300)
ggsave(file.path(out, "fig5_prop2_dens_unif.pdf"), fig5,  width = 7,   height = 4.5,  dpi = 300)
ggsave(file.path(out, "fig6_sensitivity_alpha.pdf"),fig6, width = 11,  height = 5.5,  dpi = 300)
ggsave(file.path(out, "fig7_nci_vs_hhi_scatter.pdf"),fig7,width = 13,  height = 10,   dpi = 300)
ggsave(file.path(out, "fig8_corr_nci_hhi.pdf"),    fig8,  width = 10,  height = 5.5,  dpi = 300)
ggsave(file.path(out, "fig9_corr_matrix.pdf"),     fig9,  width = 8.5, height = 8.5,  dpi = 300)

cat("\nFigures saved in:", out, "\n")







# =============================================================================
#  APPENDIX B.5 - SCALABILITY TABLE
#  N in {10, 50, 100} x 3 mechanisms x R = 1000 replications
#  Erdos-Renyi link probability p = 0.267 across N; core-periphery and
#  peripheral mechanisms retain the group-specific probabilities of
#  Appendix B.2-B.3, so their observed density varies with N
#  Random weights from the simplex (order-statistic method)
# =============================================================================

R_scale  <- 1000
N_grid   <- c(10, 50, 100)
delta_tgt <- sum(A_core) / (10 * 9)   # approx. 0.267, computed on A_core N=10

cat(sprintf("\nScalability: delta_target = %.4f\n", delta_tgt))

# -----------------------------------------------------------------------------
#  Generators adapted for arbitrary N with fixed density
# -----------------------------------------------------------------------------

# Erdos-Renyi with p = delta_tgt
gen_er_scale <- function(N, p = delta_tgt) gen_er(N, p)

# Core-periphery: core = first floor(N*0.4) nodes
# p_core and p_periph fixed at 0.9 and 0.1 as in the paper
# The expected density is:
#   delta_exp = [k(k-1)*p_core + (N-k)(N-k-1)*p_periph] / [N(N-1)]
# with k = floor(0.4*N)
gen_core_scale <- function(N, p_core = 0.9, p_periph = 0.1) {
  core_size <- max(2, floor(0.4 * N))
  gen_core(N, p_core = p_core, p_periph = p_periph, core_size = core_size)
}

# Peripheral: reversed probabilities
gen_periph_scale <- function(N) {
  core_size <- max(2, floor(0.4 * N))
  gen_core(N, p_core = 0.1, p_periph = 0.7, core_size = core_size)
}

# Random weights from the simplex via the order-statistic method (Devroye, 2006)
draw_simplex <- function(N) {
  u <- sort(runif(N - 1))
  diff(c(0, u, 1))
}

# -----------------------------------------------------------------------------
#  Main loop
# -----------------------------------------------------------------------------

set.seed(42)

scale_results <- bind_rows(lapply(N_grid, function(n) {
  
  cat(sprintf("  N = %d ...\n", n))
  
  p_n <- delta_tgt   # target density = p for ER
  
  mechs <- list(
    list(label = "Core-periphery", gen = function() gen_core_scale(n)),
    list(label = "Random",         gen = function() gen_er_scale(n, p_n)),
    list(label = "Peripheral",     gen = function() gen_periph_scale(n))
  )
  
  bind_rows(lapply(mechs, function(m) {
    
    nci_vals <- numeric(R_scale)
    hhi_vals <- numeric(R_scale)
    delta_vals <- numeric(R_scale)
    
    for (k in seq_len(R_scale)) {
      A   <- m$gen()
      w_r <- draw_simplex(n)
      
      nci_vals[k]   <- tryCatch(
        nci_base(w_r, A),
        error = function(e) NA_real_
      )
      hhi_vals[k]   <- hhi(w_r)
      delta_vals[k] <- sum(A) / (n * (n - 1))
    }
    
    data.frame(
      N         = n,
      Mechanism = m$label,
      mean_NCI  = mean(nci_vals,   na.rm = TRUE),
      sd_NCI    = sd(nci_vals,     na.rm = TRUE),
      mean_HHI  = mean(hhi_vals,   na.rm = TRUE),
      mean_delta= mean(delta_vals, na.rm = TRUE),
      n_valid   = sum(!is.na(nci_vals))
    )
  }))
}))

# -----------------------------------------------------------------------------
#  Print R table
# -----------------------------------------------------------------------------

cat("\n=== Scalability Table ===\n")
print(scale_results, digits = 4, row.names = FALSE)

# -----------------------------------------------------------------------------
#  Automatic LaTeX output
# -----------------------------------------------------------------------------

cat("\n% -- LaTeX: Scalability Table ------------------------------------------\n")
cat("\\begin{table}[htbp]\n")
cat("\\centering\\small\n")
cat("\\renewcommand{\\arraystretch}{1.3}\n")
cat("\\caption{Mean baseline NCI (standard deviation in parentheses) and mean\n")
cat("HHI across $R = 1{,}000$ replications, by network size and generating\n")
cat("mechanism. Weights drawn uniformly from $\\Delta_{N-1}$ via the\n")
cat("order-statistic method; seed~$= 42$.}\n")
cat("\\label{tab:scalability}\n")
cat("\\begin{tabular}{clccc}\n")
cat("\\hline\n")
cat("$N$ & Mechanism & Mean NCI (Std.\\ Dev.) & Mean HHI & Mean $\\delta$ \\\\\n")
cat("\\hline\n")

prev_N <- -1
for (i in seq_len(nrow(scale_results))) {
  r <- scale_results[i, ]
  
  # Separator row between blocks with different N
  if (r$N != prev_N && prev_N != -1) cat("\\hline\n")
  prev_N <- r$N
  
  cat(sprintf(
    "%d & %s & %.3f (%.3f) & %.3f & %.3f \\\\\n",
    r$N,
    r$Mechanism,
    r$mean_NCI,
    r$sd_NCI,
    r$mean_HHI,
    r$mean_delta
  ))
}

cat("\\hline\n")
cat("\\multicolumn{5}{l}{\\footnotesize \\textit{Notes}: NCI~$=$~Baseline\n")
cat("Network Concentration Index $\\psi(\\boldsymbol{w}, A)$.\n")
cat("HHI~$=$~Herfindahl--Hirschman Index. $\\delta$~$=$~observed network\n")
cat("density. All results seeded at~42.}\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
# -----------------------------------------------------------------------------
#  LaTeX TABLE
# -----------------------------------------------------------------------------

cat("\n% -- LaTeX Table 1: 7 NCI variants ---------------------------------------\n")
cat("\\begin{table}[htbp]\n\\centering\\small\n")
cat("\\caption{Variants of the NCI (Table 1, Riso \\& Zoia 2026).}\n")
cat("\\label{tab:nci_variants}\n\\renewcommand{\\arraystretch}{1.3}\n")
cat("\\begin{tabular}{clllcccc}\\hline\n")
cat("& \\textbf{Variant} & \\textbf{$M$} & \\textbf{$B$}",
    "& \\textbf{P1} & \\textbf{P2} & \\textbf{P3} & \\textbf{P4} \\\\\n\\hline\n")

# Row order = exact order of Table 1 in the paper
rows <- list(
  c("1","Baseline NCI",          "$A$",                   "$\\mathbf{1}\\mathbf{1}^\\top - I$",              "\\checkmark","\\checkmark","\\checkmark","\\checkmark"),
  c("2","Density-adjusted",      "$A$",                   "$\\delta(A)(\\mathbf{1}\\mathbf{1}^\\top-I)$",    "\\checkmark","$\\circ$",   "\\checkmark","\\checkmark"),
  c("3","Null-model NCI",        "$A$",                   "$\\mathbb{E}[A\\mid\\mathcal{N}]$",               "\\checkmark","$\\circ$",   "$\\circ$",   "\\checkmark"),
  c("4","Degree-constrained NCI","$A$",                   "$\\max_{B\\in\\mathcal{G}(d)}B$",                 "\\checkmark","$\\circ$",   "$\\circ$",   "$\\circ$"),
  c("5","Weighted NCI",          "$W=[\\gamma_{ij}]$",    "$\\mathbf{1}\\mathbf{1}^\\top - I$",              "\\checkmark","\\checkmark","\\checkmark","$\\circ$"),
  c("6","Transformed-data NCI",  "$A^{(\\mathcal{T})}$",  "$\\mathbf{1}\\mathbf{1}^\\top - I$",              "\\checkmark","\\checkmark","\\checkmark","$\\circ$"),
  c("7","Multi-layer NCI",       "$A^{(\\alpha)}$",       "$\\mathbf{1}\\mathbf{1}^\\top - I$",              "\\checkmark","\\checkmark","\\checkmark","$\\circ$")
)
for (r in rows)
  cat(sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
              r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8]))
cat("\\hline\n")
cat("\\multicolumn{8}{l}{\\footnotesize \\ck = property holds; $\\circ$ = does not hold.}\n")
cat("\\end{tabular}\n\\end{table}\n")

cat("\n=== Completed ===\n")
