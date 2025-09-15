# ============================================
# Agent Perception (split Minority/Diffusion by agent1 vs agent3)
# Requirements:
# - Libraries already loaded: tidyverse, lme4, lmerTest, emmeans, ggsignif
# - Helpers from your script exist: p_to_star, dodge_w, agent_top_pad, agent_bracket_step_u, p_adj_method
# - Data 'post1' includes agent_*, agent1_*, agent3_* columns
# ============================================

# Fallbacks (safe defaults) if not defined upstream
if (!exists("dodge_w")) dodge_w <- 0.65
if (!exists("agent_top_pad")) agent_top_pad <- 0.40
if (!exists("agent_bracket_step_u")) agent_bracket_step_u <- 0.18
if (!exists("p_adj_method")) p_adj_method <- "bonferroni"
if (!exists("p_to_star")) {
  p_to_star <- function(p) {
    dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  }
}

# 1) helper: compute dodge centers for an arbitrary number of groups
compute_dodge_offsets <- function(levels_vec, width = dodge_w) {
  k <- length(levels_vec)
  idx <- seq_len(k)
  offsets <- (idx - (k + 1) / 2) * (width / k)
  setNames(offsets, levels_vec)
}

# 2) long data builder that keeps Majority as single agent, splits Minority/Diffusion into A1/A3
make_agent_long_split <- function(post1) {
  agent_dims  <- c("competence","predictability","integrity","understanding","utility","affect","trust")
  agent_vars  <- paste0("agent_",  agent_dims)
  agent1_vars <- paste0("agent1_", agent_dims)
  agent3_vars <- paste0("agent3_", agent_dims)
  
  if (!any(names(post1) %in% agent1_vars) || !any(names(post1) %in% agent3_vars)) {
    warning("agent1_* 또는 agent3_* 컬럼이 없습니다. Minority/Diffusion에서 A1/A3 분리가 불가합니다.")
  }
  
  post1 %>%
    pivot_longer(
      cols = any_of(c(agent_vars, agent1_vars, agent3_vars)),
      names_to = c("source_tag","measure"),
      names_pattern = "^(agent(?:1|3)?)_(.+)$",
      values_to = "score",
      values_drop_na = TRUE
    ) %>%
    mutate(
      measure   = factor(measure, levels = agent_dims),
      pattern   = factor(pattern, levels = c("Majority","Minority","Diffusion")),
      task_type = factor(task_type, levels = c("Normative","Informative"))
    ) %>%
    # Majority -> agent, Minority/Diffusion -> agent1 or agent3
    filter(
      (pattern == "Majority"  & source_tag == "agent") |
        (pattern %in% c("Minority","Diffusion") & source_tag %in% c("agent1","agent3"))
    ) %>%
    mutate(
      target = dplyr::recode(source_tag, agent1 = "A1", agent3 = "A3", agent = "Agg"),
      pattern_group = dplyr::case_when(
        pattern == "Majority"                   ~ "Majority",
        pattern == "Minority"  & target == "A1" ~ "Minority A1",
        pattern == "Minority"  & target == "A3" ~ "Minority A3",
        pattern == "Diffusion" & target == "A1" ~ "Diffusion A1",
        pattern == "Diffusion" & target == "A3" ~ "Diffusion A3"
      ),
      pattern_group = factor(
        pattern_group,
        levels = c("Majority","Minority A1","Minority A3","Diffusion A1","Diffusion A3")
      )
    )
}

# 3) significance df: only A1 vs A3 within Minority and within Diffusion
make_sig_df_agent_split <- function(long_dat, p_adjust = p_adj_method, y_upper = 7 + agent_top_pad) {
  out <- list()
  # Dynamic offsets from the actual groups present
  grp_levels  <- levels(droplevels(long_dat$pattern_group))
  grp_offsets <- compute_dodge_offsets(grp_levels, width = dodge_w)
  
  agent_dims <- levels(long_dat$measure)
  for (meas in agent_dims) {
    d <- long_dat %>% dplyr::filter(measure == meas)
    if (!nrow(d)) next
    
    m <- lmer(score ~ pattern_group * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
              data = d, REML = TRUE)
    emm <- emmeans(m, ~ pattern_group | task_type)
    prs_all <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs_all)) next
    
    prs <- prs_all %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
      # Keep A1 vs A3 only within each pattern family
      dplyr::filter(
        (g1 %in% c("Minority A1","Minority A3") & g2 %in% c("Minority A1","Minority A3")) |
          (g1 %in% c("Diffusion A1","Diffusion A3") & g2 %in% c("Diffusion A1","Diffusion A3"))
      ) %>%
      dplyr::filter(p.value < 0.05) %>%
      dplyr::mutate(
        x_idx   = match(task_type, c("Normative","Informative")),
        x_left  = x_idx + grp_offsets[g1],
        x_right = x_idx + grp_offsets[g2],
        xmin    = pmin(x_left, x_right),
        xmax    = pmax(x_left, x_right)
      ) %>%
      dplyr::group_by(task_type) %>%
      dplyr::arrange(dplyr::desc(xmax - xmin), .by_group = TRUE) %>%
      dplyr::mutate(
        y     = y_upper - (dplyr::row_number() - 1) * agent_bracket_step_u,
        label = p_to_star(p.value)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(measure = meas, task_type, xmin, xmax, y, label)
    
    if (nrow(prs)) out[[meas]] <- prs
  }
  dplyr::bind_rows(out)
}

# 4) plotting function (5 groups = Majority, Minority A1/A3, Diffusion A1/A3)
plot_agent_pretty_split <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long_split(post1)
  
  # color/fill for 5 groups
  pattern5_levels <- c("Majority","Minority A1","Minority A3","Diffusion A1","Diffusion A3")
  pattern5_fill <- c(
    "Majority"     = "#e3f2fd",
    "Minority A1"  = "#ffe0b2",
    "Minority A3"  = "#ffcc80",
    "Diffusion A1" = "#e0f2f1",
    "Diffusion A3" = "#b2dfdb"
  )
  pattern5_color <- c(
    "Majority"     = "#1f77b4",
    "Minority A1"  = "#ff7f0e",
    "Minority A3"  = "#d95f0e",
    "Diffusion A1" = "#2ca02c",
    "Diffusion A3" = "#1b7f1b"
  )
  
  y_upper <- 7 + y_top_pad
  sig_df  <- make_sig_df_agent_split(long_dat, p_adjust = p_adjust, y_upper = y_upper)
  
  p <- ggplot(long_dat, aes(x = task_type, y = score,
                            fill = pattern_group, color = pattern_group)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern_group),
                 position = position_dodge(width = dodge_w), linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8) +
    facet_wrap(~ measure, nrow = 2, ncol = 4) +
    scale_fill_manual(values = pattern5_fill, breaks = pattern5_levels, name = NULL) +
    scale_color_manual(values = pattern5_color, breaks = pattern5_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(0.95, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    coord_cartesian(clip = "off") +
    theme_gray(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 3.6, color = "black", vjust = 0.5
    )
  }
  p
}

# 5) Run
fig_agent_split <- plot_agent_pretty_split(
  post1 = post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad
)
print(fig_agent_split)






### revision
# ============================================================
# Agent Perception (split Minority/Diffusion by agent1 vs agent3)
# - Fix emmeans pairs -> contrast(pairwise)
# - More top/bottom margin; robust no-clipping brackets
# - Show both:
#   (A) A1 vs A3 within each family (Minority, Diffusion) for each task_type
#   (B) Normative vs Informative within each pattern_group (incl. Majority)
# - Returns a ggplot with attributes "sig_df" and "sig_summary"
# ============================================================

# ---------- Fallback params if missing ----------
if (!exists("dodge_w")) dodge_w <- 0.65
if (!exists("p_adj_method")) p_adj_method <- "bonferroni"
if (!exists("p_to_star")) {
  p_to_star <- function(p) dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

# Visual params
likert_tip_h       <- 0.08   # small downward tick length
likert_star_off    <- 0.07   # stars slightly above bracket line
bracket_size       <- 0.40
bracket_color      <- "black"
star_size          <- 6.5

# Default spacing controls (tunable)
agent_bracket_step_u <- 0.28 # distance between stacked brackets (data units)
default_y_top_pad    <- 1.20 # base top margin above 7 (7 + this)
default_y_bottom_pad <- 0.16 # base bottom margin below 1 (1 - this)
default_bracket_gap  <- 0.35 # gap from panel max data to first bracket (data units)
default_extra_top    <- 0.40 # additional headroom above highest star (data units)

# 1) compute dodge centers for arbitrary number of groups
compute_dodge_offsets <- function(levels_vec, width = dodge_w) {
  k <- length(levels_vec)
  offs <- (seq_len(k) - (k + 1) / 2) * (width / k)
  stats::setNames(offs, levels_vec)
}

# 2) long data builder: Majority single agent, Minority/Diffusion split into A1/A3
make_agent_long_split <- function(post1) {
  agent_dims  <- c("competence","predictability","integrity","understanding","utility","affect","trust")
  agent_vars  <- paste0("agent_",  agent_dims)
  agent1_vars <- paste0("agent1_", agent_dims)
  agent3_vars <- paste0("agent3_", agent_dims)
  
  post1 %>%
    tidyr::pivot_longer(
      cols = tidyselect::any_of(c(agent_vars, agent1_vars, agent3_vars)),
      names_to = c("source_tag","measure"),
      names_pattern = "^(agent(?:1|3)?)_(.+)$",
      values_to = "score",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(
      measure   = factor(measure, levels = agent_dims),
      pattern   = factor(pattern, levels = c("Majority","Minority","Diffusion")),
      task_type = factor(task_type, levels = c("Normative","Informative"))
    ) %>%
    dplyr::filter(
      (pattern == "Majority"  & source_tag == "agent") |
        (pattern %in% c("Minority","Diffusion") & source_tag %in% c("agent1","agent3"))
    ) %>%
    dplyr::mutate(
      target = dplyr::recode(source_tag, agent1 = "A1", agent3 = "A3", agent = "Agg"),
      pattern_group = dplyr::case_when(
        pattern == "Majority"                   ~ "Majority",
        pattern == "Minority"  & target == "A1" ~ "Minority A1",
        pattern == "Minority"  & target == "A3" ~ "Minority A3",
        pattern == "Diffusion" & target == "A1" ~ "Diffusion A1",
        pattern == "Diffusion" & target == "A3" ~ "Diffusion A3"
      ),
      pattern_group = factor(
        pattern_group,
        levels = c("Majority","Minority A1","Minority A3","Diffusion A1","Diffusion A3")
      )
    )
}

# 3) cell-wise dodge offsets per task_type (robust when some groups are missing)
cell_offsets_by_task_split <- function(long_dat) {
  pg_levels <- levels(long_dat$pattern_group)
  long_dat %>%
    dplyr::mutate(pattern_group = factor(pattern_group, levels = pg_levels)) %>%
    dplyr::group_by(task_type) %>%
    dplyr::summarise(
      present = list(intersect(pg_levels, unique(as.character(pattern_group)))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(offset_tbl = purrr::map(present, ~{
      offs <- compute_dodge_offsets(.x, width = dodge_w)
      tibble::tibble(pattern_group = names(offs), offset = as.numeric(offs))
    })) %>%
    dplyr::select(-present) %>%
    tidyr::unnest(offset_tbl)
}

# 4) significance df: A1 vs A3 (within family) + Norm vs Inf (within pattern_group)
#    Bracket y-positions are anchored to data maxima (panel-wise), not to the top limit.
make_sig_df_agent_split3 <- function(long_dat, p_adjust = p_adj_method,
                                     bracket_gap = default_bracket_gap) {
  out <- list()
  agent_dims <- levels(long_dat$measure)
  offs_by_task <- cell_offsets_by_task_split(long_dat)
  
  # panel max by measure and task_type for y placement
  panel_max <- long_dat %>%
    dplyr::group_by(measure, task_type) %>%
    dplyr::summarise(ymax = max(score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(ymax = dplyr::if_else(is.finite(ymax), ymax, 1))
  
  for (meas in agent_dims) {
    d <- long_dat %>% dplyr::filter(measure == meas)
    if (!nrow(d)) next
    
    m <- lme4::lmer(score ~ pattern_group * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
                    data = d, REML = TRUE)
    
    # (A) A1 vs A3 within each family (Minority, Diffusion) and task_type
    emm_A <- emmeans::emmeans(m, ~ pattern_group | task_type)
    prs_A_all <- as.data.frame(emmeans::contrast(emm_A, method = "pairwise", adjust = p_adjust))
    sig_A <- tibble::tibble()
    if (nrow(prs_A_all)) {
      sig_A <- prs_A_all %>%
        tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
        dplyr::filter(
          (g1 %in% c("Minority A1","Minority A3") & g2 %in% c("Minority A1","Minority A3")) |
            (g1 %in% c("Diffusion A1","Diffusion A3") & g2 %in% c("Diffusion A1","Diffusion A3"))
        ) %>%
        dplyr::filter(p.value < 0.05) %>%
        dplyr::mutate(x_idx = match(task_type, c("Normative","Informative"))) %>%
        dplyr::left_join(
          offs_by_task %>% dplyr::rename(g1 = pattern_group, offset_g1 = offset),
          by = c("task_type","g1")
        ) %>%
        dplyr::left_join(
          offs_by_task %>% dplyr::rename(g2 = pattern_group, offset_g2 = offset),
          by = c("task_type","g2")
        ) %>%
        dplyr::left_join(panel_max %>% dplyr::filter(measure == meas),
                         by = "task_type") %>%
        dplyr::transmute(
          measure = meas, type = "A1_vs_A3",
          xmin = pmin(x_idx + offset_g1, x_idx + offset_g2),
          xmax = pmax(x_idx + offset_g1, x_idx + offset_g2),
          base_y = ymax + bracket_gap,
          label = p_to_star(p.value)
        )
    }
    
    # (B) Normative vs Informative within each pattern_group
    emm_B <- emmeans::emmeans(m, ~ task_type | pattern_group)
    prs_B_all <- as.data.frame(emmeans::contrast(emm_B, method = "pairwise", adjust = p_adjust))
    sig_B <- tibble::tibble()
    if (nrow(prs_B_all)) {
      offs_N <- offs_by_task %>% dplyr::filter(task_type == "Normative")   %>% dplyr::rename(offset_N = offset)
      offs_I <- offs_by_task %>% dplyr::filter(task_type == "Informative") %>% dplyr::rename(offset_I = offset)
      # panel max per pattern_group (max of both task panels)
      pg_max <- d %>%
        dplyr::group_by(pattern_group) %>%
        dplyr::summarise(ymax_pg = max(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(ymax_pg = dplyr::if_else(is.finite(ymax_pg), ymax_pg, 1))
      
      sig_B <- prs_B_all %>%
        dplyr::filter(p.value < 0.05) %>%
        dplyr::left_join(offs_N, by = "pattern_group") %>%
        dplyr::left_join(offs_I, by = "pattern_group") %>%
        dplyr::left_join(pg_max, by = "pattern_group") %>%
        dplyr::mutate(
          x_norm = 1 + offset_N,
          x_info = 2 + offset_I
        ) %>%
        dplyr::filter(is.finite(x_norm), is.finite(x_info)) %>%
        dplyr::transmute(
          measure = meas, type = "N_vs_I",
          xmin = pmin(x_norm, x_info),
          xmax = pmax(x_norm, x_info),
          base_y = ymax_pg + bracket_gap,
          label = p_to_star(p.value)
        )
    }
    
    sig_meas <- dplyr::bind_rows(sig_A, sig_B)
    if (!nrow(sig_meas)) next
    
    # Stack brackets within the measure: longer span higher
    sig_meas <- sig_meas %>%
      dplyr::mutate(width = xmax - xmin) %>%
      dplyr::arrange(dplyr::desc(width)) %>%
      dplyr::mutate(
        y      = base_y + (dplyr::row_number() - 1) * agent_bracket_step_u,
        text_y = y + likert_star_off,
        tip_y  = y - likert_tip_h
      ) %>%
      dplyr::select(measure, type, xmin, xmax, y, text_y, tip_y, label)
    
    out[[meas]] <- sig_meas
  }
  dplyr::bind_rows(out)
}

# 5) plotting function (5 groups: Majority, Minority A1/A3, Diffusion A1/A3)
plot_agent_pretty_split <- function(post1, p_adjust = p_adj_method,
                                    y_top_pad = default_y_top_pad,
                                    y_bottom_pad = default_y_bottom_pad,
                                    bracket_gap = default_bracket_gap,
                                    extra_top_space = default_extra_top,
                                    verbose = TRUE) {
  long_dat <- make_agent_long_split(post1)
  
  # colors for 5 groups
  pattern5_levels <- c("Majority","Minority A1","Minority A3","Diffusion A1","Diffusion A3")
  pattern5_fill <- c(
    "Majority"     = "#e3f2fd",
    "Minority A1"  = "#ffe0b2",
    "Minority A3"  = "#ffcc80",
    "Diffusion A1" = "#e0f2f1",
    "Diffusion A3" = "#b2dfdb"
  )
  pattern5_color <- c(
    "Majority"     = "#1f77b4",
    "Minority A1"  = "#ff7f0e",
    "Minority A3"  = "#d95f0e",
    "Diffusion A1" = "#2ca02c",
    "Diffusion A3" = "#1b7f1b"
  )
  
  # Significance (anchored to data maxima)
  sig_df  <- make_sig_df_agent_split3(
    long_dat, p_adjust = p_adjust, bracket_gap = bracket_gap
  )
  
  # Compute y-limits:
  # - base limits from Likert endpoints with pads
  # - ensure all stars fit: top = max(text_y) + extra_top_space
  y_upper_base   <- 7 + y_top_pad
  y_lower        <- 1 - y_bottom_pad
  y_upper_needed <- if (nrow(sig_df)) max(sig_df$text_y, na.rm = TRUE) + extra_top_space else y_upper_base
  y_upper_final  <- max(y_upper_base, y_upper_needed)
  
  if (verbose) {
    by_meas <- if (nrow(sig_df)) sig_df %>% dplyr::count(measure, name = "n_brackets") else tibble::tibble()
    msg <- paste0(
      "[Agent split] total brackets: ", if (nrow(sig_df)) nrow(sig_df) else 0,
      if (nrow(by_meas)) paste0("\n  by measure:\n", paste0("   - ", by_meas$measure, ": ", by_meas$n_brackets, collapse = "\n")) else ""
    )
    message(msg)
  }
  
  p <- ggplot2::ggplot(long_dat, ggplot2::aes(x = task_type, y = score,
                                              fill = pattern_group, color = pattern_group)) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = dodge_w), width = 0.6,
      outlier.shape = 1, outlier.size = 1.8, alpha = 0.9
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "line",
      ggplot2::aes(group = pattern_group),
      position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.9
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "point",
      position = ggplot2::position_dodge(width = dodge_w), size = 1.8
    ) +
    ggplot2::facet_wrap(~ measure, nrow = 2, ncol = 4) +
    ggplot2::scale_fill_manual(values = pattern5_fill, breaks = pattern5_levels, name = NULL) +
    ggplot2::scale_color_manual(values = pattern5_color, breaks = pattern5_levels, name = NULL) +
    ggplot2::scale_y_continuous(
      breaks = 1:7,
      limits = c(y_lower, y_upper_final),
      expand = ggplot2::expansion(mult = c(0.00, 0.00))
    ) +
    ggplot2::labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_gray(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      # give the figure itself some outer room (prevents device cropping)
      plot.margin = ggplot2::margin(t = 16, r = 10, b = 14, l = 10)
    )
  
  if (nrow(sig_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = sig_df,
        ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
        inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color
      ) +
      ggplot2::geom_segment(
        data = sig_df,
        ggplot2::aes(x = xmin, xend = xmin, y = y, yend = tip_y),
        inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color
      ) +
      ggplot2::geom_segment(
        data = sig_df,
        ggplot2::aes(x = xmax, xend = xmax, y = y, yend = tip_y),
        inherit.aes = FALSE, linewidth = bracket_size, color = bracket_color
      ) +
      ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(x = (xmin + xmax)/2, y = text_y, label = label),
        inherit.aes = FALSE, size = star_size / 3, fontface = "bold"
      )
  }
  
  # Attach sig info for inspection
  attr(p, "sig_df") <- sig_df
  if (nrow(sig_df)) {
    attr(p, "sig_summary") <- list(
      total = nrow(sig_df),
      by_measure = sig_df %>% dplyr::count(measure, name = "n_brackets"),
      by_type    = sig_df %>% dplyr::count(type,    name = "n_brackets")
    )
  } else {
    attr(p, "sig_summary") <- list(total = 0L)
  }
  p
}

# 6) Run with generous headroom (adjust these if still clipped)
fig_agent_split <- plot_agent_pretty_split(
  post1 = post1,
  p_adjust = p_adj_method,
  y_top_pad = 1.20,        # 상단 기본 여백 (7 -> 8.20)
  y_bottom_pad = 0.16,     # 하단 기본 여백 (1 -> 0.84)
  bracket_gap = 0.40,      # 데이터 위로 브라켓 시작 간격
  extra_top_space = 0.50,  # 가장 높은 별 위에 추가 여백 (잘림 방지에 가장 효과적)
  verbose = TRUE
)
print(fig_agent_split)

# 유의 브라켓 개수 확인
attr(fig_agent_split, "sig_summary")
