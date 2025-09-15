######### Visualization #######



##### revision gpt5
# =======================
# Pretty plots with boxplots + mean lines + adjusted brackets
# Requirements: resp_delta, post1 objects already created (from your preprocessing pipeline)
# Pattern order: Majority, Minority, Diffusion
# Brackets: only shown for significant pairs (adjusted p < .05)
# Default multiple-comparison: Bonferroni (set p_adj_method = "holm" to switch)
# =======================
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------------------------
# Global aesthetics and adjustable parameters
# -------------------------------------------------
pattern_levels <- c("Majority","Minority","Diffusion")
pattern_fill   <- c("Majority"="#e3f2fd", "Minority"="#ffe0b2", "Diffusion"="#e0f2f1")
pattern_color  <- c("Majority"="#1f77b4","Minority"="#ff7f0e","Diffusion"="#2ca02c")

# All geoms use the same dodge width
dodge_w <- 0.65

# Exact dodge offsets so brackets align with dodged boxes
get_dodge_offsets <- function(levels, dodge_width) {
  n <- length(levels)
  offs <- ((seq_len(n) - 1) - (n - 1) / 2) * (dodge_width / n)
  names(offs) <- levels
  offs
}

# -------------------------------------------------
# Spacing tweaks (top margin reduced; brackets slightly higher than before)
# -------------------------------------------------
behavior_top_expand_mult <- 0.06  # smaller top space for behavior plots

# Bracket positions (very slightly above the data)
behavior_bracket_top_pad <- 2.3   # previously around ~2.0 → slight raise
behavior_bracket_step_u  <- 2.0
self_top_pad             <- 0.35
self_bracket_top_pad     <- 0.30  # slight raise
self_bracket_step_u      <- 0.30
agent_top_pad            <- 0.45
agent_bracket_top_pad    <- 0.30  # slight raise
agent_bracket_step_u     <- 0.30

# Bracket rendering (manual; robust against missing groups)
# - tip height (downward tick) and star offset
behavior_tip_h     <- 0.8
behavior_star_off  <- 0.35
likert_tip_h       <- 0.08
likert_star_off    <- 0.06
bracket_size       <- 0.4
bracket_color      <- "black"
star_size          <- 6.5

# -------------------------------------------------
# Multiple-comparison correction
# -------------------------------------------------
p_adj_method <- "bonferroni"  # or "holm"

# -------------------------------------------------
# Helpers
# -------------------------------------------------
p_to_star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

fit_lmm_safe <- function(formula, data) {
  data_name <- deparse(substitute(data))
  m <- try(lmer(formula, data = data, REML = FALSE), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- lmer(
      update(formula, . ~ . - (1 + time_num | participant_id) + (1 | participant_id)),
      data = data, REML = FALSE
    )
  }
  m@call$data <- as.name(data_name)
  m
}

# dynamic y-limits ensuring brackets are visible but keeping top space tight
beh_y_limits <- function(data, dv, top_expand_mult = behavior_top_expand_mult, bracket_top = NA_real_) {
  if (grepl("_abs$", dv)) {
    data_max <- suppressWarnings(max(data[[dv]], na.rm = TRUE))
    if (!is.finite(data_max)) data_max <- 0
    upper <- max(100, data_max) * (1 + top_expand_mult)
    if (is.finite(bracket_top)) upper <- max(upper, bracket_top + 0.5)
    c(0, upper)
  } else {
    rng <- range(data[[dv]], na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)
    r <- diff(rng); if (r <= 0) r <- 1
    upper <- rng[2] + r * top_expand_mult
    if (is.finite(bracket_top)) upper <- max(upper, bracket_top + 0.02 * r)
    c(rng[1] - r * 0.02, upper)
  }
}

# Build cell-wise dodge offsets that reflect actually present groups (robust when some patterns are missing)
cell_offsets_behavior <- function(data) {
  data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(task_type, timeF) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      offs <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(offs), offset = as.numeric(offs))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

cell_offsets_likert <- function(data) {
  data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(task_type) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      offs <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(offs), offset = as.numeric(offs))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

# -------------------------------------------------
# A) Behavioral (Δ from T0) - Opinion / Confidence
# -------------------------------------------------
make_sig_df_behavior <- function(model, data, dv = c("opinion_delta","conf_delta","opinion_delta_abs","conf_delta_abs"),
                                 p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  # panel stats for y placement
  panel_stats <- data %>%
    group_by(task_type, timeF) %>%
    summarise(
      panel_max = max(.data[[dv]], na.rm = TRUE),
      panel_min = min(.data[[dv]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      panel_max = ifelse(is.finite(panel_max), panel_max, 0),
      panel_min = ifelse(is.finite(panel_min), panel_min, 0)
    )
  
  # cell-wise offsets (robust to missing groups)
  offs <- cell_offsets_behavior(data)
  
  # numeric x index for each timeF level
  x_levels <- levels(droplevels(factor(data$timeF)))
  prs %>%
    mutate(
      timeF = factor(timeF, levels = x_levels),
      x_idx = as.numeric(timeF)
    ) %>%
    left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","timeF","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","timeF","g2")) %>%
    left_join(panel_stats, by = c("task_type","timeF")) %>%
    mutate(
      # y placement (slightly above data)
      y = panel_max + behavior_bracket_top_pad + (dplyr::row_number() - 1L) * behavior_bracket_step_u,
      xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
      xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
      label = p_to_star(p.value),
      text_y = y + behavior_star_off,
      tip_y  = y - behavior_tip_h
    ) %>%
    select(task_type, timeF, xmin, xmax, y, text_y, tip_y, label)
}

plot_behavior_pretty <- function(data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
                                 model = NULL, p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  
  # Fit model if not provided (fixed formula)
  if (is.null(model)) {
    fml <- as.formula(paste0(
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"
    ))
    model <- fit_lmm_safe(fml, data)
  }
  
  sig_df <- make_sig_df_behavior(model, data, dv = dv, p_adjust = p_adjust)
  
  # y-limits (ensure stars/brackets visible but keep top margin tight)
  bracket_top_needed <- if (nrow(sig_df)) max(sig_df$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult, bracket_top_needed)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Opinion Delta| (Tk − T0)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Confidence Delta (Tk − T0)"
  )
  
  p <- ggplot(data, aes(x = timeF, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line", aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ task_type, nrow = 1) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.02, 0))) +
    labs(x = "Time (relative to T0)", y = ylab) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  
  if (nrow(sig_df)) {
    # Manual brackets + centered stars (robust and exact)
    p <- p +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_df,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# B) Self-reported (Compliance / Conversion)
# -------------------------------------------------
make_sig_df_self <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                             p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  panel_stats <- data %>%
    group_by(task_type) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  offs <- cell_offsets_likert(data)
  
  prs %>%
    mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
    left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","g2")) %>%
    left_join(panel_stats, by = "task_type") %>%
    group_by(task_type) %>%
    arrange(desc(abs(offset_g1 - offset_g2)), .by_group = TRUE) %>%
    mutate(
      y = panel_max + self_bracket_top_pad + (row_number() - 1) * self_bracket_step_u,
      xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
      xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
      label = p_to_star(p.value),
      text_y = y + likert_star_off,
      tip_y  = y - likert_tip_h
    ) %>%
    ungroup() %>%
    select(task_type, xmin, xmax, y, text_y, tip_y, label)
}

plot_self_pretty <- function(data, dv = c("compliance_mean","conversion_mean"),
                             model = NULL, p_adjust = p_adj_method,
                             y_top_pad = self_top_pad) {
  dv <- match.arg(dv)
  data <- data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels),
           task_type = factor(task_type, levels = c("Normative","Informative")))
  
  if (is.null(model)) {
    fml <- as.formula(paste0(dv, " ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id)"))
    model <- lmer(fml, data = data, REML = TRUE)
  }
  
  sig_df <- make_sig_df_self(model, data, dv = dv, p_adjust = p_adjust)
  
  # ensure brackets and stars visible
  y_upper_needed <- if (nrow(sig_df)) max(sig_df$text_y, na.rm = TRUE) + 0.08 else (7 + y_top_pad)
  
  p <- ggplot(data, aes(x = task_type, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed), expand = expansion(mult = c(0.01, 0))) +
    labs(
      x = "Task Type",
      y = ifelse(dv == "compliance_mean", "Perceived Compliance", "Perceived Conversion")
    ) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  
  if (nrow(sig_df)) {
    p <- p +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_df,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# C) Agent Perception (7 measures)
# -------------------------------------------------
agent_dims <- c("competence","predictability","integrity","understanding","utility","affect","trust")
agent_vars <- paste0("agent_", agent_dims)

make_agent_long <- function(post1) {
  post1 %>%
    pivot_longer(cols = any_of(agent_vars), names_to = "measure", values_to = "score") %>%
    mutate(
      measure   = sub("^agent_", "", measure),
      measure   = factor(measure, levels = agent_dims),
      pattern   = factor(pattern, levels = pattern_levels),
      task_type = factor(task_type, levels = c("Normative","Informative"))
    )
}

make_sig_df_agent <- function(long_dat, p_adjust = p_adj_method) {
  out <- list()
  for (meas in agent_dims) {
    d <- long_dat %>% filter(measure == meas)
    if (!nrow(d)) next
    
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id), data = d, REML = TRUE)
    emm <- emmeans(m, ~ pattern | task_type, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    
    prs <- prs %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    panel_stats <- d %>%
      group_by(task_type) %>%
      summarise(panel_max = max(score, na.rm = TRUE), .groups = "drop") %>%
      mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
    
    offs <- cell_offsets_likert(d)
    
    sig_df <- prs %>%
      mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
      left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","g1")) %>%
      left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","g2")) %>%
      left_join(panel_stats, by = "task_type") %>%
      group_by(task_type) %>%
      arrange(desc(abs(offset_g1 - offset_g2)), .by_group = TRUE) %>%
      mutate(
        y = panel_max + agent_bracket_top_pad + (row_number() - 1) * agent_bracket_step_u,
        xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
        xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
        label = p_to_star(p.value),
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h
      ) %>%
      ungroup() %>%
      transmute(measure = meas, task_type, xmin, xmax, y, text_y, tip_y, label)
    
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

plot_agent_pretty <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long(post1)
  sig_df <- make_sig_df_agent(long_dat, p_adjust = p_adjust)
  
  y_upper_needed <- if (nrow(sig_df)) max(sig_df$text_y, na.rm = TRUE) + 0.08 else (7 + y_top_pad)
  
  p <- ggplot(long_dat, aes(x = task_type, y = score, fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ measure, nrow = 2, ncol = 4) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed), expand = expansion(mult = c(0.01, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  
  if (nrow(sig_df)) {
    p <- p +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_df,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# Run plots (데이터 전처리 후 실행)
# -------------------------------------------------
# Behavioral: |Δ Opinion|
fig_opinion_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_opinion_abs)

# Behavioral: |Δ Confidence|
fig_conf_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta_abs",
  model = if (exists("m_conf_abs")) m_conf_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_conf_abs)

# Self-reported: Compliance
fig_comp <- plot_self_pretty(
  data = post1, dv = "compliance_mean",
  model = if (exists("m_comp")) m_comp else NULL,
  p_adjust = p_adj_method,
  y_top_pad = self_top_pad
)
print(fig_comp)

# Self-reported: Conversion
fig_conv <- plot_self_pretty(
  data = post1, dv = "conversion_mean",
  model = if (exists("m_conv")) m_conv else NULL,
  p_adjust = p_adj_method,
  y_top_pad = self_top_pad
)
print(fig_conv)

# Agent Perception: 2행 4열
fig_agent <- plot_agent_pretty(
  post1 = post1,
  p_adjust = p_adj_method,
  y_top_pad = agent_top_pad
)
print(fig_agent)




# ======== # 대박 개선안
############################
# =======================
# Pretty plots with boxplots + mean lines + adjusted brackets
# Requirements: resp_delta, post1 already exist
# Pattern order: Majority, Minority, Diffusion
# Brackets: Bonferroni-adjusted p < .05
# New:
#  - Fix for the t1 rename error in time-wise comparisons
#  - "All-in-one" behavior plot: in a single figure, show
#       (a) pattern differences within same time & task,
#       (b) task differences (N vs I) within same time & pattern,
#       (c) time differences within same task & pattern (baseline/adjacent/all)
#  - Slightly more top margin and brackets placed just a bit higher
# =======================
library(tidyverse)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------------------------
# Global aesthetics and adjustable parameters
# -------------------------------------------------
pattern_levels <- c("Majority","Minority","Diffusion")
pattern_fill   <- c("Majority"="#e3f2fd", "Minority"="#ffe0b2", "Diffusion"="#e0f2f1")
pattern_color  <- c("Majority"="#1f77b4","Minority"="#ff7f0e","Diffusion"="#2ca02c")

dodge_w <- 0.65

# Exact dodge offsets so brackets align with dodged boxes
get_dodge_offsets <- function(levels, dodge_width) {
  n <- length(levels)
  offs <- ((seq_len(n) - 1) - (n - 1) / 2) * (dodge_width / n)
  names(offs) <- levels
  offs
}

# -------------------------------------------------
# Spacing tweaks (slightly more top space)
# -------------------------------------------------
behavior_top_expand_mult <- 0.08  # a bit more than before

# Bracket positions (slightly above the data)
behavior_bracket_top_pad <- 2.4
behavior_bracket_step_u  <- 2.0
behavior_between_group_gap <- 0.6

self_top_pad             <- 0.45
self_bracket_top_pad     <- 0.32
self_bracket_step_u      <- 0.32
self_between_group_gap   <- 0.12

agent_top_pad            <- 0.55
agent_bracket_top_pad    <- 0.32
agent_bracket_step_u     <- 0.32
agent_between_group_gap  <- 0.12

# Bracket rendering (manual)
behavior_tip_h     <- 0.9
behavior_star_off  <- 0.45
likert_tip_h       <- 0.09
likert_star_off    <- 0.12
bracket_size       <- 0.4
bracket_color      <- "black"
star_size          <- 4.0

# -------------------------------------------------
# Multiple-comparison correction
# -------------------------------------------------
p_adj_method <- "bonferroni"  # or "holm"

# -------------------------------------------------
# Helpers
# -------------------------------------------------
p_to_star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

fit_lmm_safe <- function(formula, data) {
  data_name <- deparse(substitute(data))
  m <- try(lmer(formula, data = data, REML = FALSE), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- lmer(
      update(formula, . ~ . - (1 + time_num | participant_id) + (1 | participant_id)),
      data = data, REML = FALSE
    )
  }
  m@call$data <- as.name(data_name)
  m
}

# Dynamic y-limits keeping top space tight but ensuring brackets visible
beh_y_limits <- function(data, dv, top_expand_mult = behavior_top_expand_mult, bracket_top = NA_real_) {
  if (grepl("_abs$", dv)) {
    data_max <- suppressWarnings(max(data[[dv]], na.rm = TRUE))
    if (!is.finite(data_max)) data_max <- 0
    upper <- max(100, data_max) * (1 + top_expand_mult)
    if (is.finite(bracket_top)) upper <- max(upper, bracket_top + 0.8)
    c(0, upper)
  } else {
    rng <- range(data[[dv]], na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)
    r <- diff(rng); if (r <= 0) r <- 1
    upper <- rng[2] + r * top_expand_mult
    if (is.finite(bracket_top)) upper <- max(upper, bracket_top + 0.03 * r)
    c(rng[1] - r * 0.02, upper)
  }
}

# cell-wise offsets (robust to missing groups in each cell)
cell_offsets_behavior <- function(data) {
  data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(task_type, timeF) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      offs <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(offs), offset = as.numeric(offs))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

cell_offsets_likert <- function(data) {
  data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(task_type) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      offs <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(offs), offset = as.numeric(offs))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

# -------------------------------------------------
# A) Behavioral (Δ from T0): Facet-by-task version
#     - pattern-wise within each time
#     - time-wise within each pattern (baseline/adjacent/all)
# Fix: rename bug inside time-wise function
# -------------------------------------------------
make_sig_df_behavior_pattern <- function(model, data, dv, p_adjust = p_adj_method) {
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  panel_stats <- data %>%
    group_by(task_type, timeF) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 0))
  
  offs <- cell_offsets_behavior(data)
  x_levels <- levels(droplevels(factor(data$timeF)))
  
  sig_df <- prs %>%
    mutate(
      timeF = factor(timeF, levels = x_levels),
      x_idx = as.numeric(timeF)
    ) %>%
    left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","timeF","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","timeF","g2")) %>%
    left_join(panel_stats, by = c("task_type","timeF")) %>%
    group_by(task_type, timeF) %>%
    arrange(desc(abs((x_idx + offset_g2) - (x_idx + offset_g1)))) %>%
    mutate(
      tier = row_number() - 1L,
      y = panel_max + behavior_bracket_top_pad + tier * behavior_bracket_step_u,
      xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
      xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
      label = p_to_star(p.value),
      text_y = y + behavior_star_off,
      tip_y  = y - behavior_tip_h,
      type   = "pattern"
    ) %>%
    ungroup() %>%
    select(type, task_type, timeF, xmin, xmax, y, text_y, tip_y, label, tier) %>%
    tidyr::drop_na(xmin, xmax)
  
  sig_df
}

# FIXED: rename() usage when joining offsets for t1/t2
make_sig_df_behavior_time_raw <- function(model, data, dv, p_adjust = p_adj_method,
                                          time_contrast = c("baseline","adjacent","all")) {
  time_contrast <- match.arg(time_contrast)
  emm <- emmeans(model, ~ timeF | task_type * pattern, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("t1","t2"), sep = " - ")
  
  x_levels <- levels(emm@grid$timeF)
  if (is.null(x_levels)) x_levels <- levels(droplevels(factor(data$timeF)))
  
  if (time_contrast == "baseline") {
    base <- x_levels[1]
    prs <- prs %>% filter(p.value < 0.05, t1 == base | t2 == base)
  } else if (time_contrast == "adjacent") {
    prs <- prs %>%
      filter(p.value < 0.05) %>%
      mutate(i1 = match(t1, x_levels), i2 = match(t2, x_levels)) %>%
      filter(abs(i1 - i2) == 1)
  } else {
    prs <- prs %>% filter(p.value < 0.05)
  }
  if (!nrow(prs)) return(tibble())
  
  offs <- cell_offsets_behavior(data)
  
  idx_tbl <- tibble(timeF = factor(x_levels, levels = x_levels),
                    x_idx = seq_along(x_levels))
  
  raw <- prs %>%
    mutate(t1 = factor(t1, levels = x_levels),
           t2 = factor(t2, levels = x_levels)) %>%
    left_join(idx_tbl %>% rename(t1 = timeF, x_idx1 = x_idx), by = "t1") %>%
    left_join(idx_tbl %>% rename(t2 = timeF, x_idx2 = x_idx), by = "t2") %>%
    left_join(offs %>% rename(t1 = timeF, offset1 = offset), by = c("task_type","pattern","t1")) %>%
    left_join(offs %>% rename(t2 = timeF, offset2 = offset), by = c("task_type","pattern","t2")) %>%
    mutate(
      xmin = pmin(x_idx1 + offset1, x_idx2 + offset2),
      xmax = pmax(x_idx1 + offset1, x_idx2 + offset2),
      label = p_to_star(p.value)
    ) %>%
    select(task_type, pattern, t1, t2, xmin, xmax, label) %>%
    tidyr::drop_na(xmin, xmax)
  
  raw
}

plot_behavior_pretty <- function(
    data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
    model = NULL, p_adjust = p_adj_method,
    add_time_contrasts = TRUE,
    time_contrast = c("baseline","adjacent","all")
) {
  dv <- match.arg(dv)
  time_contrast <- match.arg(time_contrast)
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  
  if (is.null(model)) {
    fml <- as.formula(paste0(
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"
    ))
    model <- fit_lmm_safe(fml, data)
  }
  
  # 1) pattern-wise within time
  sig_pat <- make_sig_df_behavior_pattern(model, data, dv, p_adjust = p_adjust)
  
  # 2) time-wise within pattern (raw)
  sig_time_raw <- if (add_time_contrasts) {
    make_sig_df_behavior_time_raw(model, data, dv, p_adjust = p_adjust, time_contrast = time_contrast)
  } else {
    tibble()
  }
  
  facet_stats <- data %>%
    group_by(task_type) %>%
    summarise(facet_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(facet_max = ifelse(is.finite(facet_max), facet_max, 0))
  
  n_pat_stack_by_facet <- sig_pat %>%
    group_by(task_type, timeF) %>%
    summarise(n_pat = n(), .groups = "drop") %>%
    group_by(task_type) %>%
    summarise(n_pat_max = ifelse(n(), max(n_pat), 0), .groups = "drop")
  
  sig_time <- if (nrow(sig_time_raw)) {
    sig_time_raw %>%
      group_by(task_type) %>%
      arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
      mutate(tier = row_number() - 1L) %>%
      ungroup() %>%
      left_join(facet_stats, by = "task_type") %>%
      left_join(n_pat_stack_by_facet, by = "task_type") %>%
      mutate(
        n_pat_max = tidyr::replace_na(n_pat_max, 0),
        base_y = facet_max + behavior_bracket_top_pad + n_pat_max * behavior_bracket_step_u + behavior_between_group_gap,
        y = base_y + tier * behavior_bracket_step_u,
        text_y = y + behavior_star_off,
        tip_y  = y - behavior_tip_h,
        type   = "time"
      ) %>%
      select(type, task_type, timeF = NULL, xmin, xmax, y, text_y, tip_y, label, tier)
  } else {
    tibble()
  }
  
  sig_all <- bind_rows(sig_pat %>% select(type, task_type, timeF, xmin, xmax, y, text_y, tip_y, label),
                       sig_time)
  bracket_top_needed <- if (nrow(sig_all)) max(sig_all$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult, bracket_top_needed)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Δ Opinion| (0–100)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Δ Confidence (Tk − T0)"
  )
  
  p <- ggplot(data, aes(x = timeF, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line", aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ task_type, nrow = 1) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.02, 0))) +
    labs(x = "Time (relative to T0)", y = ylab) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  if (nrow(sig_all)) {
    p <- p +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_all,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# NEW: Behavior "all-in-one" plot (no facets)
#   x = interaction(timeF, task_type) → e.g., T1-N, T1-I, T2-N, T2-I, …
#   Shows:
#     - pattern-wise within same time & task
#     - task-wise (N vs I) within same time & pattern
#     - time-wise within same task & pattern (baseline/adjacent/all)
# -------------------------------------------------
make_x_tt <- function(data) {
  # order: T1-N, T1-I, T2-N, T2-I, ...
  time_lv <- levels(droplevels(factor(data$timeF)))
  tt_lv <- c("Normative","Informative")
  grid <- expand_grid(timeF = factor(time_lv, levels = time_lv),
                      task_type = factor(tt_lv, levels = tt_lv)) %>%
    mutate(x_tt = factor(paste(timeF, substr(task_type, 1, 1), sep = "-"),
                         levels = paste(rep(time_lv, each = 2), c("N","I"), sep = "-")))
  data %>%
    mutate(
      timeF = factor(timeF, levels = time_lv),
      task_type = factor(task_type, levels = tt_lv)
    ) %>%
    left_join(grid, by = c("timeF","task_type"))
}

offsets_by_xtt <- function(data_xtt) {
  data_xtt %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(x_tt) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      o <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(o), offset = as.numeric(o))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

# pattern-wise within (timeF, task_type)
make_sig_allinone_pattern <- function(model, data_xtt, dv, p_adjust = p_adj_method) {
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data_xtt)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>% tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>% filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  x_map <- data_xtt %>% distinct(timeF, task_type, x_tt) %>% mutate(x_idx = as.numeric(x_tt))
  offs <- offsets_by_xtt(data_xtt)
  x_max <- data_xtt %>% group_by(x_tt) %>% summarise(x_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop")
  
  prs %>%
    left_join(x_map, by = c("timeF","task_type")) %>%
    left_join(offs %>% rename(g1 = pattern, off1 = offset), by = c("x_tt","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, off2 = offset), by = c("x_tt","g2")) %>%
    left_join(x_max, by = "x_tt") %>%
    group_by(x_tt) %>%
    arrange(desc(abs((x_idx + off2) - (x_idx + off1)))) %>%
    mutate(
      tier = row_number() - 1L,
      xmin = pmin(x_idx + off1, x_idx + off2),
      xmax = pmax(x_idx + off1, x_idx + off2),
      base_y = x_max + behavior_bracket_top_pad,
      y = base_y + tier * behavior_bracket_step_u,
      text_y = y + behavior_star_off,
      tip_y = y - behavior_tip_h,
      label = p_to_star(p.value),
      type = "pattern"
    ) %>%
    ungroup() %>%
    select(type, x_tt, xmin, xmax, y, text_y, tip_y, label, tier) %>%
    tidyr::drop_na(xmin, xmax)
}

# task-wise (N vs I) within same (timeF, pattern)
make_sig_allinone_task <- function(model, data_xtt, dv, p_adjust = p_adj_method,
                                   top_from_pattern) {
  emm <- emmeans(model, ~ task_type | timeF * pattern, data = data_xtt)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>% tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>% filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  x_map <- data_xtt %>% distinct(timeF, task_type, x_tt) %>% mutate(x_idx = as.numeric(x_tt))
  offs <- offsets_by_xtt(data_xtt)
  
  # top of pattern stacks per x_tt (if none existed, use data max + pad)
  top_pat_by_x <- top_from_pattern %>%
    group_by(x_tt) %>%
    summarise(top_pat = max(y, na.rm = TRUE), .groups = "drop")
  
  x_max <- data_xtt %>% group_by(x_tt) %>% summarise(x_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop")
  base_by_x <- x_max %>%
    left_join(top_pat_by_x, by = "x_tt") %>%
    mutate(base_x = pmax(x_max + behavior_bracket_top_pad, coalesce(top_pat, -Inf)))
  
  prs %>%
    mutate(g1 = factor(g1, levels = c("Normative","Informative")),
           g2 = factor(g2, levels = c("Normative","Informative"))) %>%
    left_join(x_map %>% rename(x_tt1 = x_tt, x1 = x_idx), by = c("timeF","g1" = "task_type")) %>%
    left_join(x_map %>% rename(x_tt2 = x_tt, x2 = x_idx), by = c("timeF","g2" = "task_type")) %>%
    left_join(offs %>% rename(off1 = offset), by = c("x_tt1" = "x_tt","pattern")) %>%
    left_join(offs %>% rename(off2 = offset), by = c("x_tt2" = "x_tt","pattern")) %>%
    left_join(base_by_x %>% rename(base1 = base_x), by = c("x_tt1" = "x_tt")) %>%
    left_join(base_by_x %>% rename(base2 = base_x), by = c("x_tt2" = "x_tt")) %>%
    mutate(
      xmin = pmin(x1 + off1, x2 + off2),
      xmax = pmax(x1 + off1, x2 + off2),
      local_base = pmax(base1, base2) + behavior_between_group_gap
    ) %>%
    group_by(timeF) %>%
    arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
    mutate(
      tier = row_number() - 1L,
      y = local_base + tier * behavior_bracket_step_u,
      text_y = y + behavior_star_off,
      tip_y  = y - behavior_tip_h,
      label = p_to_star(p.value),
      type = "task"
    ) %>%
    ungroup() %>%
    select(type, x_tt1, x_tt2, xmin, xmax, y, text_y, tip_y, label, tier) %>%
    tidyr::drop_na(xmin, xmax)
}

# time-wise within same (task_type, pattern)
make_sig_allinone_time <- function(model, data_xtt, dv, p_adjust = p_adj_method,
                                   time_contrast = c("baseline","adjacent","all"),
                                   top_from_pattern, top_from_task) {
  time_contrast <- match.arg(time_contrast)
  emm <- emmeans(model, ~ timeF | task_type * pattern, data = data_xtt)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>% tidyr::separate(contrast, into = c("t1","t2"), sep = " - ")
  tlv <- levels(emm@grid$timeF)
  if (is.null(tlv)) tlv <- levels(droplevels(factor(data_xtt$timeF)))
  
  if (time_contrast == "baseline") {
    base <- tlv[1]
    prs <- prs %>% filter(p.value < 0.05, t1 == base | t2 == base)
  } else if (time_contrast == "adjacent") {
    prs <- prs %>%
      filter(p.value < 0.05) %>%
      mutate(i1 = match(t1, tlv), i2 = match(t2, tlv)) %>%
      filter(abs(i1 - i2) == 1)
  } else {
    prs <- prs %>% filter(p.value < 0.05)
  }
  if (!nrow(prs)) return(tibble())
  
  x_map <- data_xtt %>% distinct(timeF, task_type, x_tt) %>% mutate(x_idx = as.numeric(x_tt))
  offs <- offsets_by_xtt(data_xtt)
  
  # top bases
  top_pat_by_x <- top_from_pattern %>%
    group_by(x_tt) %>%
    summarise(top_pat = max(y, na.rm = TRUE), .groups = "drop")
  top_task_by_tp <- top_from_task %>%
    mutate(tp = paste(task_type, pattern)) %>%
    group_by(tp) %>%
    summarise(top_task = max(y, na.rm = TRUE), .groups = "drop")
  
  x_max <- data_xtt %>% group_by(x_tt) %>% summarise(x_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop")
  
  base_by_x <- x_max %>%
    left_join(top_pat_by_x, by = "x_tt") %>%
    mutate(base_x = pmax(x_max + behavior_bracket_top_pad, coalesce(top_pat, -Inf)))
  
  prs %>%
    mutate(
      t1 = factor(t1, levels = tlv),
      t2 = factor(t2, levels = tlv),
      tp = paste(task_type, pattern)
    ) %>%
    left_join(x_map %>% rename(x_tt1 = x_tt, x1 = x_idx), by = c("t1" = "timeF","task_type")) %>%
    left_join(x_map %>% rename(x_tt2 = x_tt, x2 = x_idx), by = c("t2" = "timeF","task_type")) %>%
    left_join(offs %>% rename(off1 = offset), by = c("x_tt1" = "x_tt","pattern")) %>%
    left_join(offs %>% rename(off2 = offset), by = c("x_tt2" = "x_tt","pattern")) %>%
    left_join(base_by_x %>% rename(base1 = base_x), by = c("x_tt1" = "x_tt")) %>%
    left_join(base_by_x %>% rename(base2 = base_x), by = c("x_tt2" = "x_tt")) %>%
    left_join(top_task_by_tp, by = "tp") %>%
    mutate(
      xmin = pmin(x1 + off1, x2 + off2),
      xmax = pmax(x1 + off1, x2 + off2),
      # ensure above both pattern and task brackets
      local_base = pmax(base1, base2, coalesce(top_task, -Inf)) + behavior_between_group_gap
    ) %>%
    group_by(task_type, pattern) %>%
    arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
    mutate(
      tier = row_number() - 1L,
      y = local_base + tier * behavior_bracket_step_u,
      text_y = y + behavior_star_off,
      tip_y  = y - behavior_tip_h,
      label  = p_to_star(p.value),
      type   = "time"
    ) %>%
    ungroup() %>%
    select(type, x_tt1, x_tt2, xmin, xmax, y, text_y, tip_y, label, tier) %>%
    tidyr::drop_na(xmin, xmax)
}

plot_behavior_allinone <- function(
    data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
    model = NULL, p_adjust = p_adj_method,
    time_contrast = c("baseline","adjacent","all")
) {
  dv <- match.arg(dv)
  time_contrast <- match.arg(time_contrast)
  
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  data_xtt <- make_x_tt(data)
  
  if (is.null(model)) {
    fml <- as.formula(paste0(
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"
    ))
    model <- fit_lmm_safe(fml, data_xtt)
  }
  
  # 1) pattern-wise (within x_tt)
  sig_pat <- make_sig_allinone_pattern(model, data_xtt, dv, p_adjust = p_adjust)
  
  # 2) task-wise (N vs I) within (time, pattern)
  sig_task <- make_sig_allinone_task(model, data_xtt, dv, p_adjust = p_adjust,
                                     top_from_pattern = sig_pat)
  
  # 3) time-wise (within task & pattern)
  sig_time <- make_sig_allinone_time(model, data_xtt, dv, p_adjust = p_adjust,
                                     time_contrast = time_contrast,
                                     top_from_pattern = sig_pat,
                                     top_from_task = sig_task)
  
  sig_all <- bind_rows(
    sig_pat %>% mutate(grouping = x_tt),
    sig_task %>% mutate(grouping = NA_character_),
    sig_time %>% mutate(grouping = NA_character_)
  )
  
  bracket_top_needed <- if (nrow(sig_all)) max(sig_all$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data_xtt, dv, behavior_top_expand_mult, bracket_top_needed)
  
  # nicer x labels: "T1-N", "T1-I", ...
  x_labels_tbl <- data_xtt %>%
    distinct(x_tt, timeF, task_type) %>%
    mutate(lbl = paste0(as.character(timeF), "-", substr(as.character(task_type), 1, 1)))
  x_lvs <- levels(data_xtt$x_tt)
  x_labs <- x_labels_tbl %>% arrange(match(x_tt, x_lvs)) %>% pull(lbl)
  names(x_labs) <- x_lvs
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Δ Opinion| (0–100)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Δ Confidence (Tk − T0)"
  )
  
  p <- ggplot(data_xtt, aes(x = x_tt, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = interaction(pattern, task_type)),
                 position = position_dodge(width = dodge_w), size = 0.9, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_x_discrete(labels = x_labs) +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.02, 0))) +
    labs(x = "Time-Task (e.g., T1-N, T1-I)", y = ylab) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  if (nrow(sig_all)) {
    p <- p +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_all,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# B) Self-reported (Compliance / Conversion)
#    - pattern-wise within task_type
#    - task-wise (N vs I) within pattern
# -------------------------------------------------
make_sig_df_self_pattern <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                                     p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  panel_stats <- data %>%
    group_by(task_type) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  offs <- cell_offsets_likert(data)
  
  prs %>%
    mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
    left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","g2")) %>%
    left_join(panel_stats, by = "task_type") %>%
    group_by(task_type) %>%
    arrange(desc(abs(offset_g1 - offset_g2)), .by_group = TRUE) %>%
    mutate(
      tier = row_number() - 1L,
      y = panel_max + self_bracket_top_pad + tier * self_bracket_step_u,
      xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
      xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
      label = p_to_star(p.value),
      text_y = y + likert_star_off,
      tip_y  = y - likert_tip_h,
      type   = "pattern"
    ) %>%
    ungroup() %>%
    select(type, task_type, xmin, xmax, y, text_y, tip_y, label, tier) %>%
    tidyr::drop_na(xmin, xmax)
}

make_sig_df_self_task_raw <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                                      p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ task_type | pattern, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"), levels = c("Normative","Informative")),
                      x_idx = c(1,2))
  offs <- cell_offsets_likert(data)
  
  prs %>%
    mutate(g1 = factor(g1, levels = c("Normative","Informative")),
           g2 = factor(g2, levels = c("Normative","Informative"))) %>%
    left_join(x_idx_tbl %>% rename(x1 = x_idx), by = c("g1" = "task_type")) %>%
    left_join(x_idx_tbl %>% rename(x2 = x_idx), by = c("g2" = "task_type")) %>%
    left_join(offs %>% rename(offset1 = offset), by = c("g1" = "task_type","pattern")) %>%
    left_join(offs %>% rename(offset2 = offset), by = c("g2" = "task_type","pattern")) %>%
    mutate(
      xmin = pmin(x1 + offset1, x2 + offset2),
      xmax = pmax(x1 + offset1, x2 + offset2),
      label = p_to_star(p.value)
    ) %>%
    select(pattern, xmin, xmax, label) %>%
    tidyr::drop_na(xmin, xmax)
}

plot_self_pretty <- function(data, dv = c("compliance_mean","conversion_mean"),
                             model = NULL, p_adjust = p_adj_method,
                             y_top_pad = self_top_pad) {
  dv <- match.arg(dv)
  data <- data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels),
           task_type = factor(task_type, levels = c("Normative","Informative")))
  
  if (is.null(model)) {
    fml <- as.formula(paste0(dv, " ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id)"))
    model <- lmer(fml, data = data, REML = TRUE)
  }
  
  sig_pat  <- make_sig_df_self_pattern(model, data, dv = dv, p_adjust = p_adjust)
  sig_task_raw <- make_sig_df_self_task_raw(model, data, dv = dv, p_adjust = p_adjust)
  
  panel_max_by_task <- data %>%
    group_by(task_type) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  pat_stack_by_task <- sig_pat %>%
    group_by(task_type) %>%
    summarise(n_pat = n(), .groups = "drop")
  
  base_task <- panel_max_by_task %>%
    left_join(pat_stack_by_task, by = "task_type") %>%
    mutate(n_pat = tidyr::replace_na(n_pat, 0),
           base_y = panel_max + self_bracket_top_pad + n_pat * self_bracket_step_u)
  
  max_base_y <- max(base_task$base_y, na.rm = TRUE)
  
  sig_task <- if (nrow(sig_task_raw)) {
    sig_task_raw %>%
      arrange(desc(abs(xmax - xmin))) %>%
      mutate(
        tier   = row_number() - 1L,
        y      = max_base_y + self_between_group_gap + tier * self_bracket_step_u,
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h,
        type   = "task"
      ) %>%
      select(type, xmin, xmax, y, text_y, tip_y, label, tier)
  } else tibble()
  
  sig_all <- bind_rows(sig_pat %>% select(type, xmin, xmax, y, text_y, tip_y, label, tier),
                       sig_task)
  y_upper_needed <- if (nrow(sig_all)) max(sig_all$text_y, na.rm = TRUE) + 0.12 else (7 + y_top_pad)
  
  p <- ggplot(data, aes(x = task_type, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed), expand = expansion(mult = c(0.01, 0))) +
    labs(
      x = "Task Type",
      y = ifelse(dv == "compliance_mean", "Perceived Compliance", "Perceived Conversion")
    ) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  if (nrow(sig_all)) {
    p <- p +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_all,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# C) Agent Perception (7 measures)
# -------------------------------------------------
agent_dims <- c("competence","predictability","integrity","understanding","utility","affect","trust")
agent_vars <- paste0("agent_", agent_dims)

make_agent_long <- function(post1) {
  post1 %>%
    pivot_longer(cols = any_of(agent_vars), names_to = "measure", values_to = "score") %>%
    mutate(
      measure   = sub("^agent_", "", measure),
      measure   = factor(measure, levels = agent_dims),
      pattern   = factor(pattern, levels = pattern_levels),
      task_type = factor(task_type, levels = c("Normative","Informative"))
    )
}

make_sig_df_agent_pattern <- function(long_dat, p_adjust = p_adj_method) {
  out <- list()
  for (meas in agent_dims) {
    d <- long_dat %>% filter(measure == meas)
    if (!nrow(d)) next
    
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id), data = d, REML = TRUE)
    emm <- emmeans(m, ~ pattern | task_type, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    
    prs <- prs %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    panel_stats <- d %>%
      group_by(task_type) %>%
      summarise(panel_max = max(score, na.rm = TRUE), .groups = "drop") %>%
      mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
    
    offs <- cell_offsets_likert(d)
    
    sig_df <- prs %>%
      mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
      left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","g1")) %>%
      left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","g2")) %>%
      left_join(panel_stats, by = "task_type") %>%
      group_by(task_type) %>%
      arrange(desc(abs(offset_g1 - offset_g2)), .by_group = TRUE) %>%
      mutate(
        tier = row_number() - 1L,
        y = panel_max + agent_bracket_top_pad + tier * agent_bracket_step_u,
        xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
        xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
        label = p_to_star(p.value),
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h,
        type   = "pattern",
        measure = meas
      ) %>%
      ungroup() %>%
      select(measure, type, task_type, xmin, xmax, y, text_y, tip_y, label, tier) %>%
      tidyr::drop_na(xmin, xmax)
    
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

make_sig_df_agent_task_raw <- function(long_dat, p_adjust = p_adj_method) {
  out <- list()
  for (meas in agent_dims) {
    d <- long_dat %>% filter(measure == meas)
    if (!nrow(d)) next
    
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id), data = d, REML = TRUE)
    emm <- emmeans(m, ~ task_type | pattern, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    
    prs <- prs %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"), levels = c("Normative","Informative")),
                        x_idx = c(1,2))
    offs <- cell_offsets_likert(d)
    
    sig_df <- prs %>%
      mutate(g1 = factor(g1, levels = c("Normative","Informative")),
             g2 = factor(g2, levels = c("Normative","Informative"))) %>%
      left_join(x_idx_tbl %>% rename(x1 = x_idx), by = c("g1" = "task_type")) %>%
      left_join(x_idx_tbl %>% rename(x2 = x_idx), by = c("g2" = "task_type")) %>%
      left_join(offs %>% rename(offset1 = offset), by = c("g1" = "task_type","pattern")) %>%
      left_join(offs %>% rename(offset2 = offset), by = c("g2" = "task_type","pattern")) %>%
      mutate(
        xmin = pmin(x1 + offset1, x2 + offset2),
        xmax = pmax(x1 + offset1, x2 + offset2),
        label = p_to_star(p.value),
        measure = meas
      ) %>%
      select(measure, pattern, xmin, xmax, label) %>%
      tidyr::drop_na(xmin, xmax)
    
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

plot_agent_pretty <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long(post1)
  
  sig_pat      <- make_sig_df_agent_pattern(long_dat, p_adjust = p_adjust)
  sig_task_raw <- make_sig_df_agent_task_raw(long_dat, p_adjust = p_adjust)
  
  panel_max_by_task_by_meas <- long_dat %>%
    group_by(measure, task_type) %>%
    summarise(panel_max = max(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  pat_stack_by_task_by_meas <- sig_pat %>%
    group_by(measure, task_type) %>%
    summarise(n_pat = n(), .groups = "drop")
  
  base_task <- panel_max_by_task_by_meas %>%
    left_join(pat_stack_by_task_by_meas, by = c("measure","task_type")) %>%
    mutate(n_pat = tidyr::replace_na(n_pat, 0),
           base_y = panel_max + agent_bracket_top_pad + n_pat * agent_bracket_step_u)
  
  base_task_max <- base_task %>%
    group_by(measure) %>%
    summarise(max_base_y = max(base_y, na.rm = TRUE), .groups = "drop")
  
  sig_task <- if (nrow(sig_task_raw)) {
    sig_task_raw %>%
      group_by(measure) %>%
      arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
      mutate(tier = row_number() - 1L) %>%
      ungroup() %>%
      left_join(base_task_max, by = "measure") %>%
      mutate(
        y      = max_base_y + agent_between_group_gap + tier * agent_bracket_step_u,
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h,
        type   = "task"
      ) %>%
      select(measure, type, xmin, xmax, y, text_y, tip_y, label, tier)
  } else tibble()
  
  sig_all <- bind_rows(sig_pat %>% select(measure, type, xmin, xmax, y, text_y, tip_y, label, tier),
                       sig_task)
  
  y_upper_needed <- if (nrow(sig_all)) max(sig_all$text_y, na.rm = TRUE) + 0.12 else (7 + y_top_pad)
  
  p <- ggplot(long_dat, aes(x = task_type, y = score, fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ measure, nrow = 2, ncol = 4) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed), expand = expansion(mult = c(0.01, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  if (nrow(sig_all)) {
    p <- p +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_segment(data = sig_all,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color) +
      geom_text(data = sig_all,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold")
  }
  
  p
}

# -------------------------------------------------
# Run plots (데이터 전처리 후 실행)
# -------------------------------------------------
# Agent Perception
fig_agent <- plot_agent_pretty(
  post1 = post1,
  p_adjust = p_adj_method,
  y_top_pad = agent_top_pad
)
print(fig_agent)

# Behavioral: Facet-by-task (pattern-wise + time-wise)
fig_opinion_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_opinion_abs)

fig_conf_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_conf_abs)

# Behavioral: All-in-one (pattern, task, time – in one panel)
fig_opinion_all <- plot_behavior_allinone(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method,
  time_contrast = "baseline"   # "adjacent" or "all" also available (more brackets)
)
print(fig_opinion_all)

fig_conf_all <- plot_behavior_allinone(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf else NULL,
  p_adjust = p_adj_method,
  time_contrast = "baseline"
)
print(fig_conf_all)



















#### 대박 개선안 이게 진짜 최종
# =======================
# Pretty plots with boxplots + mean lines + colored brackets by pattern
# What changed (simple tweaks):
#  - Bracket color: within-pattern comparisons (time-wise or task-wise) now use pattern colors
#      Majority = blue, Minority = orange, Diffusion = green
#    Between-pattern comparisons remain black (since they don’t belong to a single pattern)
#  - Stars a bit bigger
#  - Bracket vertical spacing slightly increased (less cramped)
# =======================

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# -------------------------------------------------
# Global aesthetics and adjustable parameters
# -------------------------------------------------
pattern_levels <- c("Majority","Minority","Diffusion")
pattern_fill   <- c("Majority"="#e3f2fd", "Minority"="#ffe0b2", "Diffusion"="#e0f2f1")
pattern_color  <- c("Majority"="#1f77b4","Minority"="#ff7f0e","Diffusion"="#2ca02c")

# All geoms use the same dodge width
dodge_w <- 0.65

# Exact dodge offsets so brackets align with dodged boxes (same as ggplot's position_dodge)
get_dodge_offsets <- function(levels, dodge_width) {
  n <- length(levels)
  offs <- ((seq_len(n) - 1) - (n - 1) / 2) * (dodge_width / n)
  names(offs) <- levels
  offs
}

# -------------------------------------------------
# Spacing tweaks: slightly more space between brackets, bigger stars
# -------------------------------------------------
behavior_top_expand_mult <- 0.08
behavior_bracket_top_pad <- 2.4
behavior_bracket_step_u  <- 2.30   # was 2.0 → slightly larger
behavior_between_group_gap <- 0.75 # was 0.6 → slightly larger

self_top_pad             <- 0.45
self_bracket_top_pad     <- 0.32
self_bracket_step_u      <- 0.36   # was 0.32 → slightly larger
self_between_group_gap   <- 0.14   # was 0.12 → slightly larger

agent_top_pad            <- 0.55
agent_bracket_top_pad    <- 0.32
agent_bracket_step_u     <- 0.36   # was 0.32 → slightly larger
agent_between_group_gap  <- 0.14   # was 0.12 → slightly larger

# Bracket rendering
behavior_tip_h     <- 0.9
behavior_star_off  <- 0.45
likert_tip_h       <- 0.09
likert_star_off    <- 0.13  # slightly larger top space for star text
bracket_size       <- 0.45  # a touch thicker helps perceived size with colors
bracket_color_default <- "black"
star_size          <- 6.0   # bigger stars (text size used below is star_size/3)

# -------------------------------------------------
# Multiple-comparison correction
# -------------------------------------------------
p_adj_method <- "holm"  # change to "bonferroni" or "fdr" if you prefer

# -------------------------------------------------
# Helpers
# -------------------------------------------------
p_to_star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

fit_lmm_safe <- function(formula, data) {
  m <- try(lmer(formula, data = data, REML = FALSE), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- lmer(
      update(formula, . ~ . - (1 + time_num | participant_id) + (1 | participant_id)),
      data = data, REML = FALSE
    )
  }
  m
}

# dynamic y-limits ensuring brackets are visible while keeping the top tight
beh_y_limits <- function(data, dv, top_expand_mult = behavior_top_expand_mult, bracket_top = NA_real_) {
  if (grepl("_abs$", dv)) {
    data_max <- suppressWarnings(max(data[[dv]], na.rm = TRUE))
    if (!is.finite(data_max)) data_max <- 0
    upper <- max(100, data_max) * (1 + top_expand_mult)
    if (is.finite(bracket_top)) upper <- max(upper, bracket_top + 0.8)
    c(0, upper)
  } else {
    rng <- range(data[[dv]], na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)
    r <- diff(rng); if (r <= 0) r <- 1
    upper <- rng[2] + r * top_expand_mult
    if (is.finite(bracket_top)) upper <- max(upper, bracket_top + 0.03 * r)
    c(rng[1] - r * 0.02, upper)
  }
}

# cell-wise offsets (robust to missing groups in each cell)
cell_offsets_behavior <- function(data) {
  data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(task_type, timeF) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      offs <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(offs), offset = as.numeric(offs))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

cell_offsets_likert <- function(data) {
  data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(task_type) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      offs <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(offs), offset = as.numeric(offs))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
}

# -------------------------------------------------
# A) Behavioral (Δ from T0) - Opinion / Confidence
#     - pattern-wise within each time  → brackets = black
#     - time-wise within each pattern  → brackets = pattern color
# -------------------------------------------------

# 1) Pattern-wise within time (black brackets)
make_sig_df_behavior_pattern <- function(model, data, dv, p_adjust = p_adj_method) {
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  panel_stats <- data %>%
    group_by(task_type, timeF) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 0))
  
  offs <- cell_offsets_behavior(data)
  x_levels <- levels(droplevels(factor(data$timeF)))
  
  sig_df <- prs %>%
    mutate(
      timeF = factor(timeF, levels = x_levels),
      x_idx = as.numeric(timeF)
    ) %>%
    left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","timeF","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","timeF","g2")) %>%
    left_join(panel_stats, by = c("task_type","timeF")) %>%
    group_by(task_type, timeF) %>%
    arrange(desc(abs((x_idx + offset_g2) - (x_idx + offset_g1)))) %>%
    mutate(
      tier = row_number() - 1L,
      y = panel_max + behavior_bracket_top_pad + tier * behavior_bracket_step_u,
      xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
      xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
      label = p_to_star(p.value),
      text_y = y + behavior_star_off,
      tip_y  = y - behavior_tip_h,
      type   = "pattern"  # between-pattern
    ) %>%
    ungroup() %>%
    select(type, task_type, timeF, xmin, xmax, y, text_y, tip_y, label, tier)
  
  sig_df
}

# 2) Time-wise within pattern (colored by pattern)
make_sig_df_behavior_time_raw <- function(model, data, dv, p_adjust = p_adj_method,
                                          time_contrast = c("baseline","adjacent","all")) {
  time_contrast <- match.arg(time_contrast)
  emm <- emmeans(model, ~ timeF | task_type * pattern, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- prs %>%
    tidyr::separate(contrast, into = c("t1","t2"), sep = " - ")
  
  x_levels <- levels(emm@grid$timeF)
  if (is.null(x_levels)) x_levels <- levels(droplevels(factor(data$timeF)))
  
  if (time_contrast == "baseline") {
    base <- x_levels[1]
    prs <- prs %>% filter(p.value < 0.05) %>% filter(t1 == base | t2 == base)
  } else if (time_contrast == "adjacent") {
    prs <- prs %>%
      filter(p.value < 0.05) %>%
      mutate(i1 = match(t1, x_levels), i2 = match(t2, x_levels)) %>%
      filter(abs(i1 - i2) == 1)
  } else {
    prs <- prs %>% filter(p.value < 0.05)
  }
  if (!nrow(prs)) return(tibble())
  
  offs <- cell_offsets_behavior(data)
  idx_tbl <- tibble(timeF = factor(x_levels, levels = x_levels),
                    x_idx = seq_along(x_levels))
  
  raw <- prs %>%
    mutate(t1 = factor(t1, levels = x_levels),
           t2 = factor(t2, levels = x_levels)) %>%
    left_join(idx_tbl %>% rename(t1 = timeF, x_idx1 = x_idx), by = "t1") %>%
    left_join(idx_tbl %>% rename(t2 = timeF, x_idx2 = x_idx), by = "t2") %>%
    left_join(offs %>% rename(offset1 = offset), by = c("task_type","pattern","t1" = "timeF")) %>%
    left_join(offs %>% rename(offset2 = offset), by = c("task_type","pattern","t2" = "timeF")) %>%
    mutate(
      xmin = pmin(x_idx1 + offset1, x_idx2 + offset2),
      xmax = pmax(x_idx1 + offset1, x_idx2 + offset2),
      label = p_to_star(p.value)
    ) %>%
    select(task_type, pattern, t1, t2, xmin, xmax, label)
  
  raw
}

# -------------------------------------------------
# Main Behavior plot (facet by task_type): pattern-within-time + time-within-pattern
#    - pattern-within-time (black)
#    - time-within-pattern  (colored by pattern)
# -------------------------------------------------
plot_behavior_pretty <- function(
    data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
    model = NULL, p_adjust = p_adj_method,
    add_time_contrasts = TRUE,
    time_contrast = c("baseline","adjacent","all")
) {
  dv <- match.arg(dv)
  time_contrast <- match.arg(time_contrast)
  
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  
  if (is.null(model)) {
    fml <- as.formula(paste0(
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"
    ))
    model <- fit_lmm_safe(fml, data)
  }
  
  # 1) pattern-wise within time (black)
  sig_pat <- make_sig_df_behavior_pattern(model, data, dv, p_adjust = p_adjust)
  
  # 2) time-wise within pattern (colored)
  sig_time_raw <- if (add_time_contrasts) {
    make_sig_df_behavior_time_raw(model, data, dv, p_adjust = p_adjust, time_contrast = time_contrast)
  } else {
    tibble()
  }
  
  facet_stats <- data %>%
    group_by(task_type) %>%
    summarise(facet_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(facet_max = ifelse(is.finite(facet_max), facet_max, 0))
  
  # Max stack height of pattern-wise brackets per facet
  n_pat_stack_by_facet <- sig_pat %>%
    group_by(task_type, timeF) %>%
    summarise(n_pat = n(), .groups = "drop") %>%
    group_by(task_type) %>%
    summarise(n_pat_max = ifelse(n(), max(n_pat), 0), .groups = "drop")
  
  # Finalize y for time-wise brackets (keep pattern to color)
  sig_time <- if (nrow(sig_time_raw)) {
    sig_time_raw %>%
      group_by(task_type) %>%
      arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
      mutate(tier = row_number() - 1L) %>%
      ungroup() %>%
      left_join(facet_stats, by = "task_type") %>%
      left_join(n_pat_stack_by_facet, by = "task_type") %>%
      mutate(
        n_pat_max = tidyr::replace_na(n_pat_max, 0),
        base_y = facet_max + behavior_bracket_top_pad + n_pat_max * behavior_bracket_step_u + behavior_between_group_gap,
        y = base_y + tier * behavior_bracket_step_u,
        text_y = y + behavior_star_off,
        tip_y  = y - behavior_tip_h,
        type   = "time"
      ) %>%
      select(type, task_type, pattern, xmin, xmax, y, text_y, tip_y, label, tier)
  } else {
    tibble()
  }
  
  # For y-limit calculation
  sig_all_for_ylim <- bind_rows(
    sig_pat %>% select(y, text_y),
    sig_time %>% select(y, text_y)
  )
  bracket_top_needed <- if (nrow(sig_all_for_ylim)) max(sig_all_for_ylim$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult, bracket_top_needed)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Opinion Delta| (Tk − T0)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Confidence Delta (Tk − T0)"
  )
  
  p <- ggplot(data, aes(x = timeF, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line", aes(group = pattern),
                 position = position_dodge(width = dodge_w), linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ task_type, nrow = 1) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.02, 0))) +
    labs(x = "Time (relative to T0)", y = ylab) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  # --- Brackets: draw in 4 small layers to allow per-pattern colors without changing global color scale ---
  # 1) Between-pattern brackets (black)
  if (nrow(sig_pat)) {
    p <- p +
      geom_segment(data = sig_pat,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_segment(data = sig_pat,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_segment(data = sig_pat,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_text(data = sig_pat,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = bracket_color_default)
  }
  # 2) Within-pattern (time-wise) brackets colored by pattern
  if (nrow(sig_time)) {
    for (pat in pattern_levels) {
      df <- sig_time %>% filter(pattern == pat)
      if (!nrow(df)) next
      col <- unname(pattern_color[[pat]])
      p <- p +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmax, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_text(data = df,
                  aes(x = (xmin + xmax)/2, y = text_y, label = label),
                  inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = col)
    }
  }
  
  p
}

# -------------------------------------------------
# Companion Behavior plot: between-task (N vs I) within each pattern/time
#    - within-pattern brackets colored by pattern
# -------------------------------------------------
make_sig_df_behavior_task <- function(model, data, dv, p_adjust = p_adj_method) {
  emm <- emmeans(model, ~ task_type | pattern * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"), levels = c("Normative","Informative")),
                      x_idx = c(1,2))
  
  offs <- data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    group_by(timeF, task_type) %>%
    summarise(present = list(intersect(pattern_levels, unique(as.character(pattern)))), .groups = "drop") %>%
    mutate(offset_tbl = purrr::map(present, ~{
      o <- get_dodge_offsets(.x, dodge_w)
      tibble(pattern = names(o), offset = as.numeric(o))
    })) %>%
    select(-present) %>%
    tidyr::unnest(offset_tbl)
  
  panel_stats <- data %>%
    group_by(timeF) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 0))
  
  sig_df <- prs %>%
    mutate(g1 = factor(g1, levels = c("Normative","Informative")),
           g2 = factor(g2, levels = c("Normative","Informative"))) %>%
    left_join(x_idx_tbl %>% rename(x1 = x_idx), by = c("g1" = "task_type")) %>%
    left_join(x_idx_tbl %>% rename(x2 = x_idx), by = c("g2" = "task_type")) %>%
    left_join(offs %>% rename(offset1 = offset), by = c("timeF","g1" = "task_type","pattern")) %>%
    left_join(offs %>% rename(offset2 = offset), by = c("timeF","g2" = "task_type","pattern")) %>%
    left_join(panel_stats, by = "timeF") %>%
    mutate(
      x_l = x1 + offset1,
      x_r = x2 + offset2,
      xmin = pmin(x_l, x_r),
      xmax = pmax(x_l, x_r),
      label = p_to_star(p.value)
    ) %>%
    group_by(timeF, pattern) %>%
    arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
    mutate(
      tier = row_number() - 1L,
      y = panel_max + behavior_bracket_top_pad + tier * behavior_bracket_step_u,
      text_y = y + behavior_star_off,
      tip_y  = y - behavior_tip_h
    ) %>%
    ungroup() %>%
    select(pattern, timeF, xmin, xmax, y, text_y, tip_y, label)
  
  sig_df
}

plot_behavior_taskdiff_pretty <- function(
    data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
    model = NULL, p_adjust = p_adj_method
) {
  dv <- match.arg(dv)
  
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  
  if (is.null(model)) {
    fml <- as.formula(paste0(
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"
    ))
    model <- fit_lmm_safe(fml, data)
  }
  
  sig_df <- make_sig_df_behavior_task(model, data, dv, p_adjust = p_adjust)
  bracket_top_needed <- if (nrow(sig_df)) max(sig_df$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult, bracket_top_needed)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Opinion Delta| (Tk − T0)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Confidence Delta (Tk − T0)"
  )
  
  p <- ggplot(data, aes(x = task_type, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ timeF, nrow = 1) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.02, 0))) +
    labs(x = "Task Type", y = ylab) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  # Colored brackets by pattern
  if (nrow(sig_df)) {
    for (pat in pattern_levels) {
      df <- sig_df %>% filter(pattern == pat)
      if (!nrow(df)) next
      col <- unname(pattern_color[[pat]])
      p <- p +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmax, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_text(data = df,
                  aes(x = (xmin + xmax)/2, y = text_y, label = label),
                  inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = col)
    }
  }
  
  p
}

# -------------------------------------------------
# B) Self-reported (Compliance / Conversion)
#    - pattern-wise within task_type (black)
#    - task-wise (N vs I) within pattern (colored by pattern)
# -------------------------------------------------
make_sig_df_self_pattern <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                                     p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  panel_stats <- data %>%
    group_by(task_type) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  offs <- cell_offsets_likert(data)
  
  prs %>%
    mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
    left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","g1")) %>%
    left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","g2")) %>%
    left_join(panel_stats, by = "task_type") %>%
    group_by(task_type) %>%
    arrange(desc(abs(offset_g1 - offset_g2)), .by_group = TRUE) %>%
    mutate(
      tier = row_number() - 1L,
      y = panel_max + self_bracket_top_pad + tier * self_bracket_step_u,
      xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
      xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
      label = p_to_star(p.value),
      text_y = y + likert_star_off,
      tip_y  = y - likert_tip_h,
      type   = "pattern"
    ) %>%
    ungroup() %>%
    select(type, task_type, xmin, xmax, y, text_y, tip_y, label, tier)
}

make_sig_df_self_task_raw <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                                      p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ task_type | pattern, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"), levels = c("Normative","Informative")),
                      x_idx = c(1,2))
  offs <- cell_offsets_likert(data)
  
  prs %>%
    mutate(g1 = factor(g1, levels = c("Normative","Informative")),
           g2 = factor(g2, levels = c("Normative","Informative"))) %>%
    left_join(x_idx_tbl %>% rename(x1 = x_idx), by = c("g1" = "task_type")) %>%
    left_join(x_idx_tbl %>% rename(x2 = x_idx), by = c("g2" = "task_type")) %>%
    left_join(offs %>% rename(offset1 = offset), by = c("g1" = "task_type","pattern")) %>%
    left_join(offs %>% rename(offset2 = offset), by = c("g2" = "task_type","pattern")) %>%
    mutate(
      xmin = pmin(x1 + offset1, x2 + offset2),
      xmax = pmax(x1 + offset1, x2 + offset2),
      label = p_to_star(p.value)
    ) %>%
    select(pattern, xmin, xmax, label)
}

plot_self_pretty <- function(data, dv = c("compliance_mean","conversion_mean"),
                             model = NULL, p_adjust = p_adj_method,
                             y_top_pad = self_top_pad) {
  dv <- match.arg(dv)
  
  data <- data %>%
    mutate(pattern = factor(pattern, levels = pattern_levels),
           task_type = factor(task_type, levels = c("Normative","Informative")))
  
  if (is.null(model)) {
    fml <- as.formula(paste0(dv, " ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id)"))
    model <- lmer(fml, data = data, REML = TRUE)
  }
  
  sig_pat  <- make_sig_df_self_pattern(model, data, dv = dv, p_adjust = p_adjust)
  sig_task_raw <- make_sig_df_self_task_raw(model, data, dv = dv, p_adjust = p_adjust)
  
  panel_max_by_task <- data %>%
    group_by(task_type) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  pat_stack_by_task <- sig_pat %>%
    group_by(task_type) %>%
    summarise(n_pat = n(), .groups = "drop")
  
  base_task <- panel_max_by_task %>%
    left_join(pat_stack_by_task, by = "task_type") %>%
    mutate(n_pat = tidyr::replace_na(n_pat, 0),
           base_y = panel_max + self_bracket_top_pad + n_pat * self_bracket_step_u)
  
  max_base_y <- max(base_task$base_y, na.rm = TRUE)
  
  sig_task <- if (nrow(sig_task_raw)) {
    sig_task_raw %>%
      arrange(desc(abs(xmax - xmin))) %>%
      mutate(
        tier   = row_number() - 1L,
        y      = max_base_y + self_between_group_gap + tier * self_bracket_step_u,
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h,
        type   = "task"
      ) %>%
      select(type, pattern, xmin, xmax, y, text_y, tip_y, label, tier)
  } else tibble()
  
  sig_all_for_ylim <- bind_rows(
    sig_pat %>% select(y, text_y),
    sig_task %>% select(y, text_y)
  )
  y_upper_needed <- if (nrow(sig_all_for_ylim)) max(sig_all_for_ylim$text_y, na.rm = TRUE) + 0.12 else (7 + y_top_pad)
  
  p <- ggplot(data, aes(x = task_type, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed), expand = expansion(mult = c(0.01, 0))) +
    labs(
      x = "Task Type",
      y = ifelse(dv == "compliance_mean", "Perceived Compliance", "Perceived Conversion")
    ) +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  # Brackets
  if (nrow(sig_pat)) {
    p <- p +
      geom_segment(data = sig_pat,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_segment(data = sig_pat,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_segment(data = sig_pat,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_text(data = sig_pat,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = bracket_color_default)
  }
  if (nrow(sig_task)) {
    for (pat in pattern_levels) {
      df <- sig_task %>% filter(pattern == pat)
      if (!nrow(df)) next
      col <- unname(pattern_color[[pat]])
      p <- p +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmax, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_text(data = df,
                  aes(x = (xmin + xmax)/2, y = text_y, label = label),
                  inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = col)
    }
  }
  
  p
}

# -------------------------------------------------
# C) Agent Perception (7 measures)
#    - pattern-wise within task_type (black)
#    - task-wise (N vs I) within pattern (colored by pattern)
# -------------------------------------------------
agent_dims <- c("competence","predictability","integrity","understanding","utility","affect","trust")
agent_vars <- paste0("agent_", agent_dims)

make_agent_long <- function(post1) {
  post1 %>%
    pivot_longer(cols = any_of(agent_vars), names_to = "measure", values_to = "score") %>%
    mutate(
      measure   = sub("^agent_", "", measure),
      measure   = factor(measure, levels = agent_dims),
      pattern   = factor(pattern, levels = pattern_levels),
      task_type = factor(task_type, levels = c("Normative","Informative"))
    )
}

make_sig_df_agent_pattern <- function(long_dat, p_adjust = p_adj_method) {
  out <- list()
  for (meas in agent_dims) {
    d <- long_dat %>% filter(measure == meas)
    if (!nrow(d)) next
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id), data = d, REML = TRUE)
    emm <- emmeans(m, ~ pattern | task_type, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    prs <- prs %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    panel_stats <- d %>%
      group_by(task_type) %>%
      summarise(panel_max = max(score, na.rm = TRUE), .groups = "drop") %>%
      mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
    
    offs <- cell_offsets_likert(d)
    
    sig_df <- prs %>%
      mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
      left_join(offs %>% rename(g1 = pattern, offset_g1 = offset), by = c("task_type","g1")) %>%
      left_join(offs %>% rename(g2 = pattern, offset_g2 = offset), by = c("task_type","g2")) %>%
      left_join(panel_stats, by = "task_type") %>%
      group_by(task_type) %>%
      arrange(desc(abs(offset_g1 - offset_g2)), .by_group = TRUE) %>%
      mutate(
        tier = row_number() - 1L,
        y = panel_max + agent_bracket_top_pad + tier * agent_bracket_step_u,
        xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
        xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
        label = p_to_star(p.value),
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h,
        type   = "pattern",
        measure = meas
      ) %>%
      ungroup() %>%
      select(measure, type, task_type, xmin, xmax, y, text_y, tip_y, label, tier)
    
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

make_sig_df_agent_task_raw <- function(long_dat, p_adjust = p_adj_method) {
  out <- list()
  for (meas in agent_dims) {
    d <- long_dat %>% filter(measure == meas)
    if (!nrow(d)) next
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id), data = d, REML = TRUE)
    emm <- emmeans(m, ~ task_type | pattern, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    prs <- prs %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"), levels = c("Normative","Informative")),
                        x_idx = c(1,2))
    offs <- cell_offsets_likert(d)
    
    sig_df <- prs %>%
      mutate(g1 = factor(g1, levels = c("Normative","Informative")),
             g2 = factor(g2, levels = c("Normative","Informative"))) %>%
      left_join(x_idx_tbl %>% rename(x1 = x_idx), by = c("g1" = "task_type")) %>%
      left_join(x_idx_tbl %>% rename(x2 = x_idx), by = c("g2" = "task_type")) %>%
      left_join(offs %>% rename(offset1 = offset), by = c("g1" = "task_type","pattern")) %>%
      left_join(offs %>% rename(offset2 = offset), by = c("g2" = "task_type","pattern")) %>%
      mutate(
        xmin = pmin(x1 + offset1, x2 + offset2),
        xmax = pmax(x1 + offset1, x2 + offset2),
        label = p_to_star(p.value),
        measure = meas
      ) %>%
      select(measure, pattern, xmin, xmax, label)
    
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

plot_agent_pretty <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long(post1)
  
  sig_pat      <- make_sig_df_agent_pattern(long_dat, p_adjust = p_adjust)
  sig_task_raw <- make_sig_df_agent_task_raw(long_dat, p_adjust = p_adjust)
  
  panel_max_by_task_by_meas <- long_dat %>%
    group_by(measure, task_type) %>%
    summarise(panel_max = max(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  pat_stack_by_task_by_meas <- sig_pat %>%
    group_by(measure, task_type) %>%
    summarise(n_pat = n(), .groups = "drop")
  
  base_task <- panel_max_by_task_by_meas %>%
    left_join(pat_stack_by_task_by_meas, by = c("measure","task_type")) %>%
    mutate(n_pat = tidyr::replace_na(n_pat, 0),
           base_y = panel_max + agent_bracket_top_pad + n_pat * agent_bracket_step_u)
  
  base_task_max <- base_task %>%
    group_by(measure) %>%
    summarise(max_base_y = max(base_y, na.rm = TRUE), .groups = "drop")
  
  sig_task <- if (nrow(sig_task_raw)) {
    sig_task_raw %>%
      group_by(measure) %>%
      arrange(desc(abs(xmax - xmin)), .by_group = TRUE) %>%
      mutate(tier = row_number() - 1L) %>%
      ungroup() %>%
      left_join(base_task_max, by = "measure") %>%
      mutate(
        y      = max_base_y + agent_between_group_gap + tier * agent_bracket_step_u,
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h,
        type   = "task"
      ) %>%
      select(measure, type, pattern, xmin, xmax, y, text_y, tip_y, label, tier)
  } else tibble()
  
  # --- 안전한 y-limit 계산: 비어있는 경우도 처리 ---
  tops <- c()
  if (!is.null(sig_pat)  && nrow(sig_pat))  tops <- c(tops,  max(sig_pat$text_y,  na.rm = TRUE))
  if (!is.null(sig_task) && nrow(sig_task)) tops <- c(tops, max(sig_task$text_y, na.rm = TRUE))
  y_upper_needed <- if (length(tops)) max(tops) + 0.12 else (7 + y_top_pad)
  
  p <- ggplot(long_dat, aes(x = task_type, y = score, fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ measure, nrow = 2, ncol = 4) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed), expand = expansion(mult = c(0.01, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  
  # Brackets: pattern-wise (black)
  if (!is.null(sig_pat) && nrow(sig_pat)) {
    p <- p +
      geom_segment(data = sig_pat,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_segment(data = sig_pat,
                   aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_segment(data = sig_pat,
                   aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                   inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color_default) +
      geom_text(data = sig_pat,
                aes(x = (xmin + xmax)/2, y = text_y, label = label),
                inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = bracket_color_default)
  }
  
  # Brackets: task-wise (colored by pattern)
  if (!is.null(sig_task) && nrow(sig_task)) {
    for (pat in pattern_levels) {
      df <- sig_task %>% filter(pattern == pat)
      if (!nrow(df)) next
      col <- unname(pattern_color[[pat]])
      p <- p +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmax, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmin, xend = xmin, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_segment(data = df,
                     aes(x = xmax, xend = xmax, y = y, yend = tip_y),
                     inherit.aes = FALSE, linewidth = bracket_size, color = col) +
        geom_text(data = df,
                  aes(x = (xmin + xmax)/2, y = text_y, label = label),
                  inherit.aes = FALSE, size = star_size/3, fontface = "bold", color = col)
    }
  }
  
  p
}
# -------------------------------------------------
# Run plots (데이터 전처리 후 실행)
# -------------------------------------------------
# Behavioral main: includes pattern-wise (black) and time-wise (colored) contrasts
fig_opinion_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_opinion_abs)

fig_opinion_signed <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta",
  model = if (exists("m_opinion")) m_opinion else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_opinion_signed)

fig_conf_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta_abs",
  model = if (exists("m_conf_abs")) m_conf_abs else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_conf_abs)

fig_conf_signed <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_conf_signed)

# Behavioral companion: between-task (N vs I) within each pattern and time (colored)
fig_opinion_abs_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_opinion_abs_taskdiff)

fig_conf_abs_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_conf_abs_taskdiff)

fig_conf_abs_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "conf_delta_abs",
  model = if (exists("m_conf_abs")) m_conf_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_conf_abs_taskdiff)

# Self-reported: includes pattern-wise (black) and task-wise (colored) brackets
fig_comp <- plot_self_pretty(
  data = post1, dv = "compliance_mean",
  model = if (exists("m_comp")) m_comp else NULL,
  p_adjust = p_adj_method,
  y_top_pad = self_top_pad
)
print(fig_comp)

fig_conv <- plot_self_pretty(
  data = post1, dv = "conversion_mean",
  model = if (exists("m_conv")) m_conv else NULL,
  p_adjust = p_adj_method,
  y_top_pad = self_top_pad
)
print(fig_conv)

# Agent Perception: both pattern-wise (black) and task-wise (colored) brackets
fig_agent <- plot_agent_pretty(
  post1 = post1,
  p_adjust = p_adj_method,
  y_top_pad = agent_top_pad
)
print(fig_agent)

