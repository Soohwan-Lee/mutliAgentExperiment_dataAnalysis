#### 최종 시각화 코드: T0 기준 sign-flip + T0~T4 raw / Δ / |Δ|에 대한 Bonferroni 꺾쇠 표시
#### (빈 비교 결과에서도 에러 없이 동작하도록 개선한 버전)

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
star_size          <- 6.0   # text size used below is star_size/3

# -------------------------------------------------
# Multiple-comparison correction
# -------------------------------------------------
p_adj_method <- "bonferroni"  # Bonferroni correction

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

# 분석 스크립트와 일관되게 time_num_c 랜덤슬롭 사용
fit_lmm_safe <- function(formula, data) {
  m <- try(lmer(formula, data = data, REML = FALSE), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- lmer(
      update(formula, . ~ . - (1 + time_num_c | participant_id) + (1 | participant_id)),
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
# A) Behavioral (raw / Δ / |Δ|) - Opinion / Confidence
#  - pattern-wise within each time  → black brackets
#  - time-wise within each pattern → pattern-colored brackets
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
# Main Behavior plot (facet by task_type): raw / Δ / |Δ|
#    - pattern-within-time (black)
#    - time-within-pattern  (colored by pattern)
# -------------------------------------------------
plot_behavior_pretty <- function(
    data,
    dv = c("opinion_delta_abs","conf_delta_abs",
           "opinion_delta","conf_delta",
           "opinion","confidence"),
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
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c | participant_id)"
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
  n_pat_stack_by_facet <- if (nrow(sig_pat)) {
    sig_pat %>%
      group_by(task_type, timeF) %>%
      summarise(n_pat = n(), .groups = "drop") %>%
      group_by(task_type) %>%
      summarise(n_pat_max = ifelse(n(), max(n_pat), 0), .groups = "drop")
  } else {
    # sig_pat이 완전히 비어 있는 경우: 각 task_type별로 0으로 설정
    facet_stats %>%
      transmute(task_type, n_pat_max = 0L)
  }
  
  # Finalize y for time-wise brackets
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
  
  # For y-limit calculation (safe even if no sigs)
  sig_all_for_ylim <- bind_rows(
    if (nrow(sig_pat))  sig_pat  %>% select(y, text_y)  else NULL,
    if (nrow(sig_time)) sig_time %>% select(y, text_y) else NULL
  )
  bracket_top_needed <- if (nrow(sig_all_for_ylim)) max(sig_all_for_ylim$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult, bracket_top_needed)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Δ Opinion| (Tk − T0)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    dv == "conf_delta"        ~ "Δ Confidence (Tk − T0)",
    dv == "opinion"           ~ "Opinion (sign-flipped, −100 to 100)",
    dv == "confidence"        ~ "Confidence (0–100)"
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
  
  # --- Brackets ---
  
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
  
  x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"),
                                         levels = c("Normative","Informative")),
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
    data,
    dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
    model = NULL, p_adjust = p_adj_method
) {
  dv <- match.arg(dv)
  
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  
  if (is.null(model)) {
    fml <- as.formula(paste0(
      dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c | participant_id)"
    ))
    model <- fit_lmm_safe(fml, data)
  }
  
  sig_df <- make_sig_df_behavior_task(model, data, dv, p_adjust = p_adjust)
  bracket_top_needed <- if (nrow(sig_df)) max(sig_df$text_y, na.rm = TRUE) else NA_real_
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult, bracket_top_needed)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Δ Opinion| (Tk − T0)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Δ Confidence (Tk − T0)"
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
  
  x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"),
                                         levels = c("Normative","Informative")),
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
  
  sig_pat      <- make_sig_df_self_pattern(model, data, dv = dv, p_adjust = p_adjust)
  sig_task_raw <- make_sig_df_self_task_raw(model, data, dv = dv, p_adjust = p_adjust)
  
  panel_max_by_task <- data %>%
    group_by(task_type) %>%
    summarise(panel_max = max(.data[[dv]], na.rm = TRUE), .groups = "drop") %>%
    mutate(panel_max = ifelse(is.finite(panel_max), panel_max, 1))
  
  pat_stack_by_task <- if (nrow(sig_pat)) {
    sig_pat %>%
      group_by(task_type) %>%
      summarise(n_pat = n(), .groups = "drop")
  } else {
    # 패턴 간 유의차가 전혀 없을 때: 각 task_type별 n_pat = 0
    panel_max_by_task %>%
      transmute(task_type, n_pat = 0L)
  }
  
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
  
  # 안전한 y-limit 계산
  tops <- c()
  if (nrow(sig_pat))  tops <- c(tops,  max(sig_pat$text_y,  na.rm = TRUE))
  if (nrow(sig_task)) tops <- c(tops, max(sig_task$text_y, na.rm = TRUE))
  y_upper_needed <- if (length(tops)) max(tops) + 0.12 else (7 + y_top_pad)
  
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
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed),
                       expand = expansion(mult = c(0.01, 0))) +
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
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
              data = d, REML = TRUE)
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
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
              data = d, REML = TRUE)
    emm <- emmeans(m, ~ task_type | pattern, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    prs <- prs %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    x_idx_tbl <- tibble(task_type = factor(c("Normative","Informative"),
                                           levels = c("Normative","Informative")),
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
  
  pat_stack_by_task_by_meas <- if (nrow(sig_pat)) {
    sig_pat %>%
      group_by(measure, task_type) %>%
      summarise(n_pat = n(), .groups = "drop")
  } else {
    # 어떤 measure에서도 패턴 간 유의차가 없을 때: 각 셀마다 n_pat = 0
    panel_max_by_task_by_meas %>%
      transmute(measure, task_type, n_pat = 0L)
  }
  
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
  if (nrow(sig_pat))  tops <- c(tops,  max(sig_pat$text_y,  na.rm = TRUE))
  if (nrow(sig_task)) tops <- c(tops, max(sig_task$text_y, na.rm = TRUE))
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
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper_needed),
                       expand = expansion(mult = c(0.01, 0))) +
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
  
  # Brackets: task-wise (colored by pattern)
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





##### 개인용 꺾은선 플롯:
plot_behavior_indiv <- function(
    data,
    dv = c("opinion", "confidence"),
    facet = c("task_only", "task_by_pattern")
) {
  dv    <- match.arg(dv)
  facet <- match.arg(facet)
  
  data <- data %>%
    mutate(
      pattern   = factor(pattern, levels = pattern_levels),
      task_type = factor(task_type, levels = c("Normative", "Informative"))
    ) %>%
    arrange(participant_id, task_type, time_num)
  
  # 공통 ggplot 객체
  p <- ggplot(
    data,
    aes(
      x     = time_num,
      y     = .data[[dv]],
      group = participant_id,      # 참가자별 꺾은선
      color = pattern              # 패턴 색으로 구분
    )
  ) +
    geom_line(alpha = 0.25, linewidth = 0.4) +
    geom_point(alpha = 0.35, size = 0.8) +
    scale_x_continuous(
      breaks = 0:4,
      labels = paste0("T", 0:4),
      minor_breaks = NULL
    ) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = "Pattern") +
    labs(
      x = "Time (relative to T0)",
      y = if (dv == "opinion") "Opinion (sign-flipped, −100 to 100)" else "Confidence (0–100)"
    ) +
    theme_gray(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # facet 설정
  if (facet == "task_only") {
    p <- p + facet_wrap(~ task_type, nrow = 1)
  } else { # task_by_pattern
    p <- p + facet_grid(pattern ~ task_type)
  }
  
  # (선 위에 집단 평균 꺾은선도 보고 싶으면 아래 주석 해제)
  # p <- p +
  #   stat_summary(
  #     fun = mean,
  #     geom = "line",
  #     aes(group = interaction(pattern, task_type)),
  #     linewidth = 1.0,
  #     color = "black"
  #   )
  
  p
}

add_indiv_lines_to_behavior <- function(
    p,
    data,
    dv = c("opinion", "confidence"),
    dodge_width = 0.3,   # boxplot에서 쓰는 dodge width와 비슷하게
    alpha_line  = 0.12,
    alpha_point = 0.20,
    sample_n    = NULL   # 너무 복잡하면 일부 참가자만 보기
) {
  dv <- match.arg(dv)
  
  # pattern 레벨 맞추기
  if (exists("pattern_levels")) {
    pat_lvls <- pattern_levels
  } else {
    pat_lvls <- levels(factor(data$pattern))
  }
  
  data <- data %>%
    dplyr::mutate(
      pattern   = factor(pattern, levels = pat_lvls),
      task_type = factor(task_type, levels = c("Normative", "Informative")),
      time_num  = as.numeric(timeF)
    )
  
  # 원하면 일부 참가자만 샘플링
  if (!is.null(sample_n)) {
    ids_all  <- unique(data$participant_id)
    sample_n <- min(sample_n, length(ids_all))
    set.seed(123)
    ids_sub  <- sample(ids_all, sample_n)
    data     <- dplyr::filter(data, participant_id %in% ids_sub)
  }
  
  # 패턴별 x offset 계산
  n_pat <- length(pat_lvls)
  # [-dodge_width/2, +dodge_width/2] 사이에 고르게 배치
  offset_vals <- seq(-dodge_width/2, dodge_width/2, length.out = n_pat)
  names(offset_vals) <- pat_lvls
  
  data <- data %>%
    dplyr::mutate(
      x_indiv = time_num + offset_vals[as.character(pattern)]
    )
  
  p +
    # 참가자 × pattern별 꺾은선
    geom_line(
      data = data,
      aes(
        x     = x_indiv,
        y     = .data[[dv]],
        group = interaction(participant_id, pattern),  # 참가자×조건별 라인
        color = pattern
      ),
      alpha       = alpha_line,
      linewidth   = 0.35,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    # 참가자 × pattern별 점
    geom_point(
      data = data,
      aes(
        x     = x_indiv,
        y     = .data[[dv]],
        group = interaction(participant_id, pattern),
        color = pattern
      ),
      alpha       = alpha_point,
      size        = 0.7,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
}



# -------------------------------------------------
# Run plots (데이터 전처리 후 실행)
# resp_behavior: T0~T4 raw (sign-flip), resp_delta: Δ / |Δ|
# -------------------------------------------------

## 1) Raw T0~T4 (Opinion / Confidence)
fig_opinion_raw <- plot_behavior_pretty(
  data = resp_behavior, dv = "opinion",
  model = if (exists("m_opinion_raw")) m_opinion_raw else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"  # T0 vs T1~T4
)
print(fig_opinion_raw)

fig_conf_raw <- plot_behavior_pretty(
  data = resp_behavior, dv = "confidence",
  model = if (exists("m_conf_raw")) m_conf_raw else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_conf_raw)


# 1-1) 참가자별 꺾은선 추가 -> RAW Data 시각화는 이거 사용
## Opinion: 박스플롯 + 참가자별 꺾은선
fig_opinion_raw_indiv <- plot_behavior_pretty(
  data = resp_behavior, dv = "opinion",
  model = if (exists("m_opinion_raw")) m_opinion_raw else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)

fig_opinion_raw_indiv <- add_indiv_lines_to_behavior(
  p           = fig_opinion_raw_indiv,
  data        = resp_behavior,
  dv          = "opinion",
  dodge_width = 0.43   # geom_boxplot()의 position_dodge(width=?)와 비슷하게
  # sample_n  = 50
)

fig_opinion_raw_indiv <- fig_opinion_raw_indiv +
  labs(
    x = "Time",      # 원하는 x축 라벨
    y = "Opinion"
  )

print(fig_opinion_raw_indiv)


## Confidence: 박스플롯 + 참가자별 꺾은선
fig_conf_raw_indiv <- plot_behavior_pretty(
  data = resp_behavior, dv = "confidence",
  model = if (exists("m_conf_raw")) m_conf_raw else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)

fig_conf_raw_indiv <- add_indiv_lines_to_behavior(
  p           = fig_conf_raw_indiv,
  data        = resp_behavior,
  dv          = "confidence",
  dodge_width = 0.43
)

fig_conf_raw_indiv <- fig_conf_raw_indiv +
  labs(
    x = "Time",      # 원하는 x축 라벨
    y = "Confidence"
  )

print(fig_conf_raw_indiv)






## 2) Δ / |Δ| (T1~T4, baseline = T0)
# Opinion Delta - Abs
fig_opinion_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"  # baseline = T1 (Δ는 이미 T0 차)
)
print(fig_opinion_abs)

# Opinion Delta - Signed
fig_opinion_signed <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta",
  model = if (exists("m_opinion_signed")) m_opinion_signed else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_opinion_signed)

# Confidence Delta - Abs
fig_conf_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta_abs",
  model = if (exists("m_conf_abs")) m_conf_abs else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_conf_abs)

# Confidence Delta - Signed
fig_conf_signed <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf_signed")) m_conf_signed else NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE, time_contrast = "baseline"
)
print(fig_conf_signed)


## 3) Behavioral companion: between-task (N vs I) within each pattern/time
fig_opinion_abs_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_opinion_abs_taskdiff)

fig_opinion_signed_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "opinion_delta",
  model = if (exists("m_opinion_signed")) m_opinion_signed else NULL,
  p_adjust = p_adj_method
)
print(fig_opinion_signed_taskdiff)

fig_conf_abs_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "conf_delta_abs",
  model = if (exists("m_conf_abs")) m_conf_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_conf_abs_taskdiff)

fig_conf_signed_taskdiff <- plot_behavior_taskdiff_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf_signed")) m_conf_signed else NULL,
  p_adjust = p_adj_method
)
print(fig_conf_signed_taskdiff)

## 4) Self-reported: includes pattern-wise (black) and task-wise (colored) brackets
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

## 5) Agent Perception: both pattern-wise (black) and task-wise (colored) brackets
fig_agent <- plot_agent_pretty(
  post1 = post1,
  p_adjust = p_adj_method,
  y_top_pad = agent_top_pad
)
print(fig_agent)


