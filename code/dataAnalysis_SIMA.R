### anotehr real working code (+ mean sd 추가)
# =======================
# Multi-Agent Experiment: Analysis Script (fixed & robust)
# =======================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(effectsize)
})

# Type-III ANOVA 권장 대비(contr.sum)
options(contrasts = c("contr.sum", "contr.poly"))

# -------------------------------------------------------
# 0) 경로
# -------------------------------------------------------
base_dir <- "./Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data/"

# -------------------------------------------------------
# 1) 데이터 로드
# -------------------------------------------------------
bg   <- readr::read_csv(file.path(base_dir, "FINAL_background_surveys_rows_refined.csv"))
post <- readr::read_csv(file.path(base_dir, "FINAL_post_self_surveys_rows_refined.csv"))
resp <- readr::read_csv(file.path(base_dir, "FINAL_responses_rows_refined_experiment_only.csv"))

# -------------------------------------------------------
# 2) ID 타입 통일
# -------------------------------------------------------
bg$participant_id   <- as.character(bg$participant_id)
post$participant_id <- as.character(post$participant_id)
resp$participant_id <- as.character(resp$participant_id)

# -------------------------------------------------------
# 3) 공변량: 문항 평균 → z-표준화
# -------------------------------------------------------
get_mean_scale <- function(df, pattern_regex) {
  cols <- grep(pattern_regex, names(df), value = TRUE)
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  out <- rowMeans(dplyr::select(df, all_of(cols)), na.rm = TRUE)
  out[is.nan(out)] <- NA_real_  # 모든 항목 NA인 경우 NaN -> NA
  out
}

bg_cov <- bg %>%
  mutate(
    SII_raw   = get_mean_scale(., "^sii_\\d+$"),
    NFC_raw   = get_mean_scale(., "^nfc_\\d+$"),
    AIacc_raw = get_mean_scale(., "^ai_acceptance_\\d+$")
  ) %>%
  mutate(
    SII_z   = as.numeric(scale(SII_raw)),
    NFC_z   = as.numeric(scale(NFC_raw)),
    AIacc_z = as.numeric(scale(AIacc_raw))
  ) %>%
  select(participant_id, SII_z, NFC_z, AIacc_z)

# -------------------------------------------------------
# 4) 행동 데이터: 패턴/과업 지정 + 시간 정렬(t_ord=0..4)
# -------------------------------------------------------
resp1 <- resp %>%
  mutate(
    pattern = dplyr::recode(condition_type,
                            "majority" = "Majority",
                            "minority" = "Minority",
                            "minorityDiffusion" = "Diffusion"),
    task_type = dplyr::case_when(
      session_key %in% c("normative", "Normative")     ~ "Normative",
      session_key %in% c("informative", "Informative") ~ "Informative",
      TRUE ~ NA_character_
    ),
    created_at_parsed = suppressWarnings(lubridate::ymd_hms(created_at, quiet = TRUE)),
    rt_log = log1p(rt_ms)
  ) %>%
  filter(!is.na(pattern), !is.na(task_type)) %>%
  arrange(participant_id, task_type, created_at_parsed, response_index) %>%
  group_by(participant_id, task_type) %>%
  mutate(t_ord = row_number() - 1L) %>%
  ungroup() %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

# -------------------------------------------------------
# 5) 베이스라인(T0) 대비 델타(T1~T4) + 절대값 델타
# -------------------------------------------------------
resp_delta <- resp1 %>%
  group_by(participant_id, task_type) %>%
  mutate(
    opinion_T0 = opinion[t_ord == 0][1],
    conf_T0    = confidence[t_ord == 0][1],
    rt_T0      = rt_log[t_ord == 0][1],
    # signed delta
    opinion_delta = opinion - opinion_T0,
    conf_delta    = confidence - conf_T0,
    rt_log_delta  = rt_log - rt_T0
  ) %>%
  ungroup() %>%
  filter(t_ord > 0) %>%                                   # T1~T4만
  mutate(
    time_num = as.integer(t_ord),                         # 1..4
    timeF    = factor(time_num, levels = 1:4, labels = c("T1","T2","T3","T4")),
    # 절대값 델타
    opinion_delta_abs = abs(opinion_delta),
    conf_delta_abs    = abs(conf_delta)
  ) %>%
  left_join(bg_cov, by = "participant_id")

# -------------------------------------------------------
# 6) Sign flip 지표(참고; 출력은 하지 않음)
# -------------------------------------------------------
resp_flip <- resp1 %>%
  arrange(participant_id, task_type, t_ord) %>%
  group_by(participant_id, pattern, task_type) %>%
  summarise(
    T0 = opinion[which.min(t_ord)],  # t_ord==0
    flip_support_to_opp = as.integer(!is.na(T0) && T0 > 0 && any(opinion[t_ord > 0] < 0, na.rm = TRUE)),
    ever_changed_any    = as.integer(any(opinion[t_ord > 0] != T0, na.rm = TRUE)),
    flip_count_any = {
      s <- sign(opinion)
      s_prev <- dplyr::lag(s)
      sum(!is.na(s) & !is.na(s_prev) & s != s_prev)
    },
    flip_count_strict = {
      s <- sign(opinion)
      s_nz <- s[!is.na(s) & s != 0]
      if (length(s_nz) <= 1) 0L else sum(diff(s_nz) != 0)
    },
    .groups = "drop"
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

# -------------------------------------------------------
# 7) Self-report: compliance/ conversion 평균 + agent_* 채우기
# -------------------------------------------------------
post1 <- post %>%
  mutate(
    pattern = dplyr::recode(condition_type,
                            "majority" = "Majority",
                            "minority" = "Minority",
                            "minorityDiffusion" = "Diffusion"),
    task_type = dplyr::case_when(
      task_type %in% c("normative","Normative")     ~ "Normative",
      task_type %in% c("informative","Informative") ~ "Informative",
      TRUE ~ task_type
    )
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

comp_cols <- grep("^perceived_compliance_\\d+$", names(post1), value = TRUE)
conv_cols <- grep("^perceived_conversion_\\d+$", names(post1), value = TRUE)
if (length(comp_cols) > 0) {
  post1$compliance_mean <- rowMeans(post1[, comp_cols], na.rm = TRUE)
  post1$compliance_mean[is.nan(post1$compliance_mean)] <- NA_real_
}
if (length(conv_cols) > 0) {
  post1$conversion_mean <- rowMeans(post1[, conv_cols], na.rm = TRUE)
  post1$conversion_mean[is.nan(post1$conversion_mean)] <- NA_real_
}

agent_dims <- c("competence","predictability","integrity","understanding","utility","affect","trust")
for (d in agent_dims) {
  a_col  <- paste0("agent_", d)
  a1_col <- paste0("agent1_", d)
  a3_col <- paste0("agent3_", d)
  if (!(a_col %in% names(post1))) post1[[a_col]] <- NA_real_
  if (all(c(a1_col, a3_col) %in% names(post1))) {
    idx <- which(post1$pattern %in% c("Minority","Diffusion") & is.na(post1[[a_col]]))
    if (length(idx) > 0) {
      post1[[a_col]][idx] <- rowMeans(cbind(post1[[a1_col]][idx], post1[[a3_col]][idx]), na.rm = TRUE)
    }
  }
}

# -------------------------------------------------------
# 8) LMM/GLMM: 델타(T1~T4), self-report, agent 인식
# -------------------------------------------------------

# 공용 step-down LMM 적합기(수렴 안정화 + emmeans 호환 call 정리)
fit_lmm_step_general <- function(form_full, form_uncorr, form_int, data, reml = FALSE) {
  ctrl1 <- lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore")
  ctrl2 <- lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE, check.conv.singular = "ignore")
  
  m <- try(lmer(form_full, data = data, REML = reml, control = ctrl1), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- try(lmer(form_full, data = data, REML = reml, control = ctrl2), silent = TRUE)
  }
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- try(lmer(form_uncorr, data = data, REML = reml, control = ctrl1), silent = TRUE)
    if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
      m <- try(lmer(form_uncorr, data = data, REML = reml, control = ctrl2), silent = TRUE)
    }
  }
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- lmer(form_int, data = data, REML = reml, control = ctrl1)
  }
  
  # emmeans가 model@call$data를 통해 실제 data를 재참조할 수 있도록 정리
  data_name <- deparse(substitute(data))
  m@call$formula <- formula(m)
  m@call$data    <- as.name(data_name)
  if ("control" %in% names(m@call)) m@call$control <- NULL
  m
}

# time 중심화 + 그룹변수 factor 보장
resp_delta <- resp_delta %>%
  mutate(time_num_c = as.numeric(scale(time_num, scale = FALSE)))
resp_delta$participant_id <- factor(resp_delta$participant_id)
post1$participant_id      <- factor(post1$participant_id)

# ---------- 8-1) 행동 델타 LMM (opinion/confidence: 절대값, rt_log: signed delta) ----------
form_op_abs_full   <- opinion_delta_abs ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_op_abs_uncorr <- opinion_delta_abs ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_op_abs_int    <- opinion_delta_abs ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_abs_full   <- conf_delta ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_cf_abs_uncorr <- conf_delta ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_cf_abs_int    <- conf_delta ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 | participant_id)


# 적합
m_opinion_abs <- fit_lmm_step_general(form_op_abs_full, form_op_abs_uncorr, form_op_abs_int, data = resp_delta, reml = FALSE)
m_conf_abs    <- fit_lmm_step_general(form_cf_abs_full, form_cf_abs_uncorr, form_cf_abs_int, data = resp_delta, reml = FALSE)

cat("\n=== |ΔOpinion| (Tk−T0) ===\n");  print(summary(m_opinion_abs)); print(anova(m_opinion_abs, type=3))
cat("\n=== |ΔConfidence| (Tk−T0) ===\n");print(summary(m_conf_abs));    print(anova(m_conf_abs, type=3))

# 사후검정(Bonferroni) + 효과크기
emm_op_abs <- emmeans(m_opinion_abs, specs = "pattern", by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_op_abs, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_abs, sigma = sigma(m_opinion_abs), edf = df.residual(m_opinion_abs)))

emm_cf_abs <- emmeans(m_conf_abs, specs = "pattern", by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_cf_abs, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_abs, sigma = sigma(m_conf_abs), edf = df.residual(m_conf_abs)))

# ---------- 8-3) Self-report: compliance/ conversion LMM ----------
m_comp <- lmer(
  compliance_mean ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = post1, REML = TRUE
)
m_conv <- lmer(
  conversion_mean ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = post1, REML = TRUE
)

cat("\n=== Self-report: Compliance ===\n"); print(summary(m_comp)); print(anova(m_comp, type=3))
emm_comp <- emmeans(m_comp, specs = "pattern", by = "task_type", data = post1)
print(pairs(emm_comp, adjust="bonferroni"))
print(emmeans::eff_size(emm_comp, sigma = sigma(m_comp), edf = df.residual(m_comp)))

cat("\n=== Self-report: Conversion ===\n"); print(summary(m_conv)); print(anova(m_conv, type=3))
emm_conv <- emmeans(m_conv, specs = "pattern", by = "task_type", data = post1)
print(pairs(emm_conv, adjust="bonferroni"))
print(emmeans::eff_size(emm_conv, sigma = sigma(m_conv), edf = df.residual(m_conv)))





# ---------- 8-4) Agent 인식(각 척도별 LMM + Bonferroni 사후비교) ----------
fit_agent_agg <- function(var) {
  # 안전장치: 컬럼과 유효 데이터 확인
  if (!var %in% names(post1)) {
    return(list(var = var, error = sprintf("Column '%s' not found in post1", var)))
  }
  if (all(is.na(post1[[var]]))) {
    return(list(var = var, error = sprintf("Column '%s' is all NA; cannot fit model", var)))
  }
  
  frm <- as.formula(paste0(var, " ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id)"))
  res <- try({
    m <- lmer(frm, data = post1, REML = TRUE)
    emm <- emmeans(m, specs = "pattern", by = "task_type", data = post1)
    list(
      var = var,
      model = m,
      anova = anova(m, type = 3),
      pairs = pairs(emm, adjust = "bonferroni"),
      eff = emmeans::eff_size(emm, sigma = sigma(m), edf = df.residual(m))
    )
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    return(list(var = var, error = as.character(res)))
  } else {
    return(res)
  }
}

agent_models_agg <- lapply(paste0("agent_", agent_dims), fit_agent_agg)

# 결과 출력
for (res in agent_models_agg) {
  cat("\n=== Agent perception: ", res$var, " ===\n", sep = "")
  if (!is.null(res$error)) {
    cat("Error: ", res$error, "\n", sep = "")
    next
  }
  print(summary(res$model))
  print(res$anova)
  print(res$pairs)  # Bonferroni 보정 사후비교
  print(res$eff)
}

# -------------------------------------------------------
# 9) 기술통계 표: (Informative/Normative) x (Majority/Minority/Diffusion)
#    + 축별 "All" 포함 → 3 x 4 = 12개 셀의 mean, sd
#    measures: opinion_delta_abs, confidence_delta(conf_delta), 
#              compliance_mean, conversion_mean, agent_* (7개)
#    필요 객체: resp_delta, post1, agent_dims
# -------------------------------------------------------

# 안전장치: 필수 객체 확인
stopifnot(exists("resp_delta"), exists("post1"))

# agent_dims 미정의 시 기본 정의
if (!exists("agent_dims")) {
  agent_dims <- c("competence","predictability","integrity","understanding","utility","affect","trust")
}

# NaN -> NA 유틸
nan_to_na <- function(x) { x[is.nan(x)] <- NA_real_; x }

# 사용할 측정 변수 목록(존재하는 것만 사용)
meas_resp <- intersect(c("opinion_delta_abs", "conf_delta"), names(resp_delta))
meas_post <- intersect(c("compliance_mean", "conversion_mean", paste0("agent_", agent_dims)), names(post1))

if (length(meas_resp) == 0) stop("resp_delta에 요약할 측정변수가 없습니다.")
if (length(meas_post) == 0) stop("post1에 요약할 측정변수가 없습니다.")

# 9-1) resp_delta: 참가자×패턴×과업별(시간평균) 요약값
resp_delta_pp <- resp_delta %>%
  group_by(participant_id, pattern, task_type) %>%
  summarise(across(all_of(meas_resp), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 9-2) post1: 참가자×패턴×과업별 요약값(개인평균)
post1_pp <- post1 %>%
  group_by(participant_id, pattern, task_type) %>%
  summarise(across(all_of(meas_post), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 9-3) long 포맷으로 결합
resp_long <- resp_delta_pp %>%
  tidyr::pivot_longer(cols = all_of(meas_resp), names_to = "measure", values_to = "value") %>%
  mutate(measure = dplyr::recode(measure, conf_delta = "confidence_delta"))

post_long <- post1_pp %>%
  tidyr::pivot_longer(cols = all_of(meas_post), names_to = "measure", values_to = "value")

long_pp <- bind_rows(resp_long, post_long)

# 레벨 고정(결측 셀도 생성할 수 있게)
task_levels    <- c("Normative","Informative","All")
pattern_levels <- c("Majority","Minority","Diffusion","All")

# 9-4) 6개 기본 셀(Informative/Normative x Majority/Minority/Diffusion)
base6 <- long_pp %>%
  group_by(measure, task_type, pattern) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean = nan_to_na(mean), sd = nan_to_na(sd))

# 9-5) 축별 All + 전체 All-All
# "패턴 All"(task_type 별): 참가자×과업×측정항목으로 패턴 평균 → 집단 통계
all_pattern <- long_pp %>%
  group_by(participant_id, task_type, measure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(measure, task_type) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pattern = "All", mean = nan_to_na(mean), sd = nan_to_na(sd))

# "과업 All"(pattern 별): 참가자×패턴×측정항목으로 과업 평균 → 집단 통계
all_task <- long_pp %>%
  group_by(participant_id, pattern, measure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(measure, pattern) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(task_type = "All", mean = nan_to_na(mean), sd = nan_to_na(sd))

# "전체 All-All": 참가자×측정항목 전체 평균 → 집단 통계
all_all <- long_pp %>%
  group_by(participant_id, measure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(measure) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(task_type = "All", pattern = "All", mean = nan_to_na(mean), sd = nan_to_na(sd))

# 9-6) 12개 셀 결합 + 모든 셀 강제 생성(결측이면 NA)
cells12 <- bind_rows(
  base6,
  all_pattern %>% select(measure, task_type, pattern, n, mean, sd),
  all_task    %>% select(measure, task_type, pattern, n, mean, sd),
  all_all     %>% select(measure, task_type, pattern, n, mean, sd)
) %>%
  mutate(
    task_type = factor(task_type, levels = task_levels),
    pattern   = factor(pattern,   levels = pattern_levels)
  ) %>%
  tidyr::complete(
    measure,
    task_type = factor(task_levels, levels = task_levels),
    pattern   = factor(pattern_levels, levels = pattern_levels),
    fill = list(n = 0, mean = NA_real_, sd = NA_real_)
  ) %>%
  arrange(measure, task_type, pattern)

# 9-7) Wide 포맷 1: "Mean (SD)" 문자열로 12개 셀
table_wide_str12 <- cells12 %>%
  mutate(
    mean_round = round(mean, 2),
    sd_round   = round(sd, 2),
    mean_sd    = ifelse(is.na(mean), NA_character_, sprintf("%.2f (%.2f)", mean_round, sd_round)),
    cell       = paste(task_type, pattern, sep = " - ")
  ) %>%
  select(measure, cell, mean_sd) %>%
  tidyr::pivot_wider(names_from = cell, values_from = mean_sd) %>%
  arrange(measure)

cat("\n=== Descriptive: Mean (SD) for 12 cells (Task x Pattern incl. Alls) ===\n")
print(as.data.frame(table_wide_str12), row.names = FALSE)

# 9-8) Wide 포맷 2: 각 셀당 두 컬럼(mu, sigma) = 총 24 컬럼
table_wide_split12 <- cells12 %>%
  mutate(cell = paste(task_type, pattern, sep = " - ")) %>%
  transmute(measure, cell, mu = round(mean, 2), sigma = round(sd, 2)) %>%
  tidyr::pivot_longer(cols = c(mu, sigma), names_to = "stat", values_to = "val") %>%
  tidyr::pivot_wider(
    names_from = c(cell, stat),
    values_from = val
  ) %>%
  arrange(measure)

cat("\n=== Descriptive: Separate mu and sigma columns (12 cells x 2) ===\n")
print(as.data.frame(table_wide_split12), row.names = FALSE)

# 9-9) CSV 저장(옵션)
out_dir <- if (exists("base_dir")) file.path(base_dir, "analysis_outputs") else "analysis_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(cells12,            file.path(out_dir, "descriptives_cells12_tidy.csv"))
readr::write_csv(table_wide_str12,   file.path(out_dir, "descriptives_cells12_wide_mean_sd_str.csv"))
readr::write_csv(table_wide_split12, file.path(out_dir, "descriptives_cells12_wide_mu_sigma.csv"))

cat("\nDescriptive tables saved to:\n", out_dir, "\n")







##### this is real working code
# =======================
# Multi-Agent Experiment: Analysis Script (fixed & robust)
# =======================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(effectsize)
})

# Type-III ANOVA 권장 대비(contr.sum)
options(contrasts = c("contr.sum", "contr.poly"))

# -------------------------------------------------------
# 0) 경로
# -------------------------------------------------------
base_dir <- "./Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data/"

# -------------------------------------------------------
# 1) 데이터 로드
# -------------------------------------------------------
bg   <- readr::read_csv(file.path(base_dir, "FINAL_background_surveys_rows_refined.csv"))
post <- readr::read_csv(file.path(base_dir, "FINAL_post_self_surveys_rows_refined.csv"))
resp <- readr::read_csv(file.path(base_dir, "FINAL_responses_rows_refined_experiment_only.csv"))

# -------------------------------------------------------
# 2) ID 타입 통일
# -------------------------------------------------------
bg$participant_id   <- as.character(bg$participant_id)
post$participant_id <- as.character(post$participant_id)
resp$participant_id <- as.character(resp$participant_id)

# -------------------------------------------------------
# 3) 공변량: 문항 평균 → z-표준화
# -------------------------------------------------------
get_mean_scale <- function(df, pattern_regex) {
  cols <- grep(pattern_regex, names(df), value = TRUE)
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  out <- rowMeans(dplyr::select(df, all_of(cols)), na.rm = TRUE)
  out[is.nan(out)] <- NA_real_  # 모든 항목 NA인 경우 NaN -> NA
  out
}

bg_cov <- bg %>%
  mutate(
    SII_raw   = get_mean_scale(., "^sii_\\d+$"),
    NFC_raw   = get_mean_scale(., "^nfc_\\d+$"),
    AIacc_raw = get_mean_scale(., "^ai_acceptance_\\d+$")
  ) %>%
  mutate(
    SII_z   = as.numeric(scale(SII_raw)),
    NFC_z   = as.numeric(scale(NFC_raw)),
    AIacc_z = as.numeric(scale(AIacc_raw))
  ) %>%
  select(participant_id, SII_z, NFC_z, AIacc_z)

# -------------------------------------------------------
# 4) 행동 데이터: 패턴/과업 지정 + 시간 정렬(t_ord=0..4)
# -------------------------------------------------------
resp1 <- resp %>%
  mutate(
    pattern = dplyr::recode(condition_type,
                            "majority" = "Majority",
                            "minority" = "Minority",
                            "minorityDiffusion" = "Diffusion"),
    task_type = dplyr::case_when(
      session_key %in% c("normative", "Normative")     ~ "Normative",
      session_key %in% c("informative", "Informative") ~ "Informative",
      TRUE ~ NA_character_
    ),
    created_at_parsed = suppressWarnings(lubridate::ymd_hms(created_at, quiet = TRUE)),
    rt_log = log1p(rt_ms)
  ) %>%
  filter(!is.na(pattern), !is.na(task_type)) %>%
  arrange(participant_id, task_type, created_at_parsed, response_index) %>%
  group_by(participant_id, task_type) %>%
  mutate(t_ord = row_number() - 1L) %>%
  ungroup() %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

# -------------------------------------------------------
# 5) 베이스라인(T0) 대비 델타(T1~T4) + 절대값 델타
# -------------------------------------------------------
resp_delta <- resp1 %>%
  group_by(participant_id, task_type) %>%
  mutate(
    opinion_T0 = opinion[t_ord == 0][1],
    conf_T0    = confidence[t_ord == 0][1],
    rt_T0      = rt_log[t_ord == 0][1],
    # signed delta
    opinion_delta = opinion - opinion_T0,
    conf_delta    = confidence - conf_T0,
    rt_log_delta  = rt_log - rt_T0
  ) %>%
  ungroup() %>%
  filter(t_ord > 0) %>%                                   # T1~T4만
  mutate(
    time_num = as.integer(t_ord),                         # 1..4
    timeF    = factor(time_num, levels = 1:4, labels = c("T1","T2","T3","T4")),
    # 절대값 델타
    opinion_delta_abs = abs(opinion_delta),
    conf_delta_abs    = abs(conf_delta)
  ) %>%
  left_join(bg_cov, by = "participant_id")

# -------------------------------------------------------
# 6) Sign flip 지표(참고; 출력은 하지 않음)
# -------------------------------------------------------
resp_flip <- resp1 %>%
  arrange(participant_id, task_type, t_ord) %>%
  group_by(participant_id, pattern, task_type) %>%
  summarise(
    T0 = opinion[which.min(t_ord)],  # t_ord==0
    flip_support_to_opp = as.integer(!is.na(T0) && T0 > 0 && any(opinion[t_ord > 0] < 0, na.rm = TRUE)),
    ever_changed_any    = as.integer(any(opinion[t_ord > 0] != T0, na.rm = TRUE)),
    flip_count_any = {
      s <- sign(opinion)
      s_prev <- dplyr::lag(s)
      sum(!is.na(s) & !is.na(s_prev) & s != s_prev)
    },
    flip_count_strict = {
      s <- sign(opinion)
      s_nz <- s[!is.na(s) & s != 0]
      if (length(s_nz) <= 1) 0L else sum(diff(s_nz) != 0)
    },
    .groups = "drop"
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

# -------------------------------------------------------
# 7) Self-report: compliance/ conversion 평균 + agent_* 채우기
# -------------------------------------------------------
post1 <- post %>%
  mutate(
    pattern = dplyr::recode(condition_type,
                            "majority" = "Majority",
                            "minority" = "Minority",
                            "minorityDiffusion" = "Diffusion"),
    task_type = dplyr::case_when(
      task_type %in% c("normative","Normative")     ~ "Normative",
      task_type %in% c("informative","Informative") ~ "Informative",
      TRUE ~ task_type
    )
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

comp_cols <- grep("^perceived_compliance_\\d+$", names(post1), value = TRUE)
conv_cols <- grep("^perceived_conversion_\\d+$", names(post1), value = TRUE)
if (length(comp_cols) > 0) {
  post1$compliance_mean <- rowMeans(post1[, comp_cols], na.rm = TRUE)
  post1$compliance_mean[is.nan(post1$compliance_mean)] <- NA_real_
}
if (length(conv_cols) > 0) {
  post1$conversion_mean <- rowMeans(post1[, conv_cols], na.rm = TRUE)
  post1$conversion_mean[is.nan(post1$conversion_mean)] <- NA_real_
}

agent_dims <- c("competence","predictability","integrity","understanding","utility","affect","trust")
for (d in agent_dims) {
  a_col  <- paste0("agent_", d)
  a1_col <- paste0("agent1_", d)
  a3_col <- paste0("agent3_", d)
  if (!(a_col %in% names(post1))) post1[[a_col]] <- NA_real_
  if (all(c(a1_col, a3_col) %in% names(post1))) {
    idx <- which(post1$pattern %in% c("Minority","Diffusion") & is.na(post1[[a_col]]))
    if (length(idx) > 0) {
      post1[[a_col]][idx] <- rowMeans(cbind(post1[[a1_col]][idx], post1[[a3_col]][idx]), na.rm = TRUE)
    }
  }
}

# -------------------------------------------------------
# 8) LMM/GLMM: 델타(T1~T4), self-report, agent 인식
# -------------------------------------------------------

# 공용 step-down LMM 적합기(수렴 안정화 + emmeans 호환 call 정리)
fit_lmm_step_general <- function(form_full, form_uncorr, form_int, data, reml = FALSE) {
  ctrl1 <- lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore")
  ctrl2 <- lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE, check.conv.singular = "ignore")
  
  m <- try(lmer(form_full, data = data, REML = reml, control = ctrl1), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- try(lmer(form_full, data = data, REML = reml, control = ctrl2), silent = TRUE)
  }
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- try(lmer(form_uncorr, data = data, REML = reml, control = ctrl1), silent = TRUE)
    if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
      m <- try(lmer(form_uncorr, data = data, REML = reml, control = ctrl2), silent = TRUE)
    }
  }
  if (inherits(m, "try-error") || isSingular(m, tol = 1e-4)) {
    m <- lmer(form_int, data = data, REML = reml, control = ctrl1)
  }
  
  # emmeans가 model@call$data를 통해 실제 data를 재참조할 수 있도록 정리
  data_name <- deparse(substitute(data))
  m@call$formula <- formula(m)
  m@call$data    <- as.name(data_name)
  if ("control" %in% names(m@call)) m@call$control <- NULL
  m
}

# time 중심화 + 그룹변수 factor 보장
resp_delta <- resp_delta %>%
  mutate(time_num_c = as.numeric(scale(time_num, scale = FALSE)))
resp_delta$participant_id <- factor(resp_delta$participant_id)
post1$participant_id      <- factor(post1$participant_id)

# ---------- 8-1) 행동 델타 LMM (opinion/confidence: 절대값, rt_log: signed delta) ----------
# 교호작용은 pattern * task_type * timeF
form_op_abs_full   <- opinion_delta_abs ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_op_abs_uncorr <- opinion_delta_abs ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_op_abs_int    <- opinion_delta_abs ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_abs_full   <- conf_delta ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_cf_abs_uncorr <- conf_delta ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_cf_abs_int    <- conf_delta_ ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 | participant_id)


# 적합
m_opinion_abs <- fit_lmm_step_general(form_op_abs_full, form_op_abs_uncorr, form_op_abs_int, data = resp_delta, reml = FALSE)
m_conf_abs    <- fit_lmm_step_general(form_cf_abs_full, form_cf_abs_uncorr, form_cf_abs_int, data = resp_delta, reml = FALSE)

cat("\n=== |ΔOpinion| (Tk−T0) ===\n");  print(summary(m_opinion_abs)); print(anova(m_opinion_abs, type=3))
cat("\n=== |ΔConfidence| (Tk−T0) ===\n");print(summary(m_conf_abs));    print(anova(m_conf_abs, type=3))

# 사후검정(Bonferroni) + 효과크기
emm_op_abs <- emmeans(m_opinion_abs, specs = "pattern", by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_op_abs, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_abs, sigma = sigma(m_opinion_abs), edf = df.residual(m_opinion_abs)))

emm_cf_abs <- emmeans(m_conf_abs, specs = "pattern", by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_cf_abs, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_abs, sigma = sigma(m_conf_abs), edf = df.residual(m_conf_abs)))

# ---------- 8-3) Self-report: compliance/ conversion LMM ----------
m_comp <- lmer(
  compliance_mean ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = post1, REML = TRUE
)
m_conv <- lmer(
  conversion_mean ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = post1, REML = TRUE
)

cat("\n=== Self-report: Compliance ===\n"); print(summary(m_comp)); print(anova(m_comp, type=3))
emm_comp <- emmeans(m_comp, specs = "pattern", by = "task_type", data = post1)
print(pairs(emm_comp, adjust="bonferroni"))
print(emmeans::eff_size(emm_comp, sigma = sigma(m_comp), edf = df.residual(m_comp)))

cat("\n=== Self-report: Conversion ===\n"); print(summary(m_conv)); print(anova(m_conv, type=3))
emm_conv <- emmeans(m_conv, specs = "pattern", by = "task_type", data = post1)
print(pairs(emm_conv, adjust="bonferroni"))
print(emmeans::eff_size(emm_conv, sigma = sigma(m_conv), edf = df.residual(m_conv)))

# ---------- 8-4) Agent 인식(각 척도별 LMM + Bonferroni 사후비교) ----------
fit_agent_agg <- function(var) {
  # 안전장치: 컬럼과 유효 데이터 확인
  if (!var %in% names(post1)) {
    return(list(var = var, error = sprintf("Column '%s' not found in post1", var)))
  }
  if (all(is.na(post1[[var]]))) {
    return(list(var = var, error = sprintf("Column '%s' is all NA; cannot fit model", var)))
  }
  
  frm <- as.formula(paste0(var, " ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id)"))
  res <- try({
    m <- lmer(frm, data = post1, REML = TRUE)
    emm <- emmeans(m, specs = "pattern", by = "task_type", data = post1)
    list(
      var = var,
      model = m,
      anova = anova(m, type = 3),
      pairs = pairs(emm, adjust = "bonferroni"),
      eff = emmeans::eff_size(emm, sigma = sigma(m), edf = df.residual(m))
    )
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    return(list(var = var, error = as.character(res)))
  } else {
    return(res)
  }
}

agent_models_agg <- lapply(paste0("agent_", agent_dims), fit_agent_agg)

# 결과 출력
for (res in agent_models_agg) {
  cat("\n=== Agent perception: ", res$var, " ===\n", sep = "")
  if (!is.null(res$error)) {
    cat("Error: ", res$error, "\n", sep = "")
    next
  }
  print(summary(res$model))
  print(res$anova)
  print(res$pairs)  # Bonferroni 보정 사후비교
  print(res$eff)
}

cat("\nAll analyses completed.\n")

























######### Visualization ####### -> 이것도 옛날 코드

# =======================
# Pretty plots with boxplots + mean lines + Holm/Bonnferroni-adjusted brackets
# Requirements: resp_delta, post1 objects already created (from your preprocessing pipeline)
# Pattern order: Majority, Minority, Diffusion
# Brackets: only shown for significant pairs (adjusted p < .05)
# Default multiple-comparison: Holm (set p_adjust = "bonferroni" to switch)
# =======================
# =======================
# Pretty plots with boxplots + mean lines + Holm/Bonnferroni-adjusted brackets
# =======================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)

# -------------------------------------------------
# Global aesthetics and adjustable parameters
# -------------------------------------------------
pattern_levels <- c("Majority","Minority","Diffusion")
pattern_fill   <- c("Majority"="#e3f2fd", "Minority"="#ffe0b2", "Diffusion"="#e0f2f1")
pattern_color  <- c("Majority"="#1f77b4","Minority"="#ff7f0e","Diffusion"="#2ca02c")

dodge_w <- 0.65
pattern_offsets <- c(Majority = -0.22, Minority = 0.00, Diffusion = 0.22)

## ✨ 최종 파라미터: 별 잘림 및 겹침 문제 해결을 위해 여백과 간격을 크게 설정
# ------------------------------------------------------------------------------------
# 1. y축 상단 전체 여백을 더 많이 확보
behavior_top_expand_mult <- 0.25  # 행동 지표 플롯 높이 25% 증가
self_top_pad             <- 0.8   # Self-report 플롯 상단 여백
agent_top_pad            <- 1.3   # Agent 플롯 상단 여백 (겹침 방지를 위해 대폭 증가)

# 2. 브라켓 시작 위치를 더 내리고, 브라켓 간 간격도 넓게 설정
behavior_bracket_top_pad <- 12.0  # 첫 브라켓 시작 위치 (y축 상단으로부터의 거리)
behavior_bracket_step_u  <- 5.0   # 브라켓 간 수직 간격

self_bracket_top_pad     <- 0.4   
self_bracket_step_u      <- 0.4   

agent_bracket_top_pad    <- 0.4   # 첫 브라켓 시작 위치
agent_bracket_step_u     <- 0.35  # 브라켓 간 수직 간격 (겹침 방지를 위해 대폭 증가)
# ------------------------------------------------------------------------------------

# 다중비교 보정
p_adj_method <- "bonferroni"


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

# -------------------------------------------------
# A) Behavioral (Δ from T0) - Opinion / Confidence
# -------------------------------------------------
beh_y_limits <- function(data, dv, top_expand_mult = behavior_top_expand_mult) {
  if (grepl("_abs$", dv)) {
    ymax <- 100 + 100 * top_expand_mult
    return(c(0, ymax))
  } else {
    rng <- range(data[[dv]], na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)
    r <- diff(rng); if (r <= 0) r <- 1
    return(c(rng[1] - r*0.02, rng[2] + r*top_expand_mult))
  }
}

make_sig_df_behavior <- function(model, data, dv = c("opinion_delta","conf_delta","opinion_delta_abs"),
                                 p_adjust = p_adj_method, y_limits = NULL) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  
  if (!nrow(prs)) return(tibble())
  
  y_upper <- y_limits[2]
  
  sig_df <- prs %>%
    mutate(
      x_idx = as.numeric(factor(timeF, levels = levels(data$timeF))),
      xmin = x_idx + pattern_offsets[g1],
      xmax = x_idx + pattern_offsets[g2]
    ) %>%
    group_by(task_type, timeF) %>%
    arrange(desc(abs(pattern_offsets[g1] - pattern_offsets[g2]))) %>%
    mutate(
      y = y_upper - behavior_bracket_top_pad - (row_number() - 1) * behavior_bracket_step_u,
      label = p_to_star(p.value)
    ) %>%
    ungroup() %>%
    select(task_type, timeF, xmin, xmax, y, label)
  
  sig_df
}

plot_behavior_pretty <- function(data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
                                 model = NULL, p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  data <- data %>% mutate(
    pattern = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  if (is.null(model)) {
    fml <- as.formula(paste0(dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"))
    model <- try(lmer(fml, data = data, REML = FALSE), silent = TRUE)
    if (inherits(model, "try-error") || isSingular(model, tol=1e-4)) {
      model <- lmer(update(fml, . ~ . - (1 + time_num | participant_id) + (1 | participant_id)), data = data, REML = FALSE)
    }
  }
  
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult)
  sig_df <- make_sig_df_behavior(model, data, dv = dv, p_adjust = p_adjust, y_limits = ylims)
  
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
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.0, 0.0))) +
    labs(x = "Time (relative to T0)", y = ylab) +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 4,
      vjust = -0.15, color = "black"
    )
  }
  p
}

# -------------------------------------------------
# B) Self-reported (Compliance / Conversion)
# -------------------------------------------------
make_sig_df_self <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                             p_adjust = p_adj_method, y_upper = 7 + self_top_pad) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05)
  
  if (!nrow(prs)) return(tibble())
  
  sig_df <- prs %>%
    mutate(
      x_idx = match(task_type, c("Normative","Informative")),
      xmin = x_idx + pattern_offsets[g1],
      xmax = x_idx + pattern_offsets[g2]
    ) %>%
    group_by(task_type) %>%
    arrange(desc(abs(pattern_offsets[g1] - pattern_offsets[g2]))) %>%
    mutate(
      y = y_upper - self_bracket_top_pad - (row_number() - 1) * self_bracket_step_u,
      label = p_to_star(p.value)
    ) %>%
    ungroup() %>%
    select(task_type, xmin, xmax, y, label)
  
  sig_df
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
  
  y_upper <- 7 + y_top_pad
  sig_df <- make_sig_df_self(model, data, dv = dv, p_adjust = p_adjust, y_upper = y_upper)
  
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
    scale_y_continuous(breaks = 1:7, limits = c(1, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(
      x = "Task Type",
      y = ifelse(dv == "compliance_mean", "Perceived Compliance", "Perceived Conversion")
    ) +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 4,
      vjust = -0.15, color = "black"
    )
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

make_sig_df_agent <- function(long_dat, p_adjust = p_adj_method, y_upper = 7 + agent_top_pad) {
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
    
    sig_df <- prs %>%
      mutate(
        x_idx = match(task_type, c("Normative","Informative")),
        xmin = x_idx + pattern_offsets[g1],
        xmax = x_idx + pattern_offsets[g2]
      ) %>%
      group_by(task_type) %>%
      arrange(desc(abs(pattern_offsets[g1] - pattern_offsets[g2]))) %>%
      mutate(
        y = y_upper - agent_bracket_top_pad - (row_number() - 1) * agent_bracket_step_u,
        label = p_to_star(p.value)
      ) %>%
      ungroup() %>%
      transmute(measure = meas, task_type, xmin, xmax, y, label)
    
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

plot_agent_pretty <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long(post1)
  y_upper <- 7 + y_top_pad
  sig_df <- make_sig_df_agent(long_dat, p_adjust = p_adjust, y_upper = y_upper)
  
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
    scale_y_continuous(breaks = 1:7, limits = c(1, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 4,
      vjust = -0.15, color = "black"
    )
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
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf_abs else NULL,
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





#### GPT5 - this is real !!!!!!!!!!!!!!!!!!!!!
###############################
######### Visualization #######
# =======================
# Pretty plots with boxplots + mean lines + Holm/Bonferroni-adjusted brackets
# Requirements: resp_delta, post1 objects already created
# Pattern order: Majority, Minority, Diffusion
# Brackets: only shown for significant pairs (adjusted p < .05)
# Default multiple-comparison: Holm (set p_adjust = "bonferroni" to switch)
# =======================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)

# -------------------------------------------------
# Global aesthetics and adjustable parameters
# -------------------------------------------------
pattern_levels <- c("Majority","Minority","Diffusion")
pattern_fill   <- c("Majority"="#e3f2fd", "Minority"="#ffe0b2", "Diffusion"="#e0f2f1")
pattern_color  <- c("Majority"="#1f77b4","Minority"="#ff7f0e","Diffusion"="#2ca02c")

dodge_w <- 0.65

# ggplot2 position_dodge와 동일한 오프셋 계산식 (n개 그룹이면 step = width/n)
compute_dodge_offsets <- function(levels_chr, width) {
  n <- length(levels_chr)
  offs <- ( (seq_len(n) - (n + 1)/2) * (width / n) )   # n=3 -> (-1,0,1)*(width/3)
  setNames(offs, levels_chr)
}
pattern_offsets <- compute_dodge_offsets(pattern_levels, dodge_w)

# ===== 여백 파라미터 (상단 과도 여백 ↓, 하단 소량 여백 ↑) =====
behavior_top_expand_mult <- 0.12  # 행동 지표 상단 확장(데이터 범위 대비)
behavior_bottom_pad      <- 2.0   # |Δ| 지표의 하단 여백 (단위: y값)

self_top_pad   <- 0.6   # Self-report 상단 여백
agent_top_pad  <- 0.6   # Agent 상단 여백

# 브라켓(별) 수직 위치
behavior_bracket_top_pad <- 8.0
behavior_bracket_step_u  <- 4.0
self_bracket_top_pad     <- 0.35
self_bracket_step_u      <- 0.30
agent_bracket_top_pad    <- 0.35
agent_bracket_step_u     <- 0.30

# 다중비교 보정
p_adj_method <- "bonferroni"

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

# -------------------------------------------------
# A) Behavioral (Δ from T0) - Opinion / Confidence
# -------------------------------------------------
beh_y_limits <- function(data, dv, top_expand_mult = behavior_top_expand_mult) {
  if (grepl("_abs$", dv)) {
    ymax <- 100 + 100 * top_expand_mult
    ymin <- 0 - behavior_bottom_pad
    return(c(ymin, ymax))
  } else {
    rng <- range(data[[dv]], na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)
    r <- diff(rng); if (r <= 0) r <- 1
    return(c(rng[1] - r*0.04, rng[2] + r*top_expand_mult))
  }
}

make_sig_df_behavior <- function(model, data, dv = c("opinion_delta","conf_delta","opinion_delta_abs","conf_delta_abs"),
                                 p_adjust = p_adj_method, y_limits = NULL) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
    filter(p.value < 0.05) %>%
    mutate(
      g1 = as.character(factor(g1, levels = pattern_levels)),
      g2 = as.character(factor(g2, levels = pattern_levels))
    )
  if (!nrow(prs)) return(tibble())
  
  y_upper <- y_limits[2]
  time_levels <- levels(data$timeF)
  
  prs %>%
    mutate(
      x_idx = match(timeF, time_levels),
      xmin  = x_idx + as.numeric(pattern_offsets[g1]),
      xmax  = x_idx + as.numeric(pattern_offsets[g2])
    ) %>%
    group_by(task_type, timeF) %>%
    arrange(desc(abs(as.numeric(pattern_offsets[g1]) - as.numeric(pattern_offsets[g2])))) %>%
    mutate(
      y     = y_upper - behavior_bracket_top_pad - (row_number() - 1) * behavior_bracket_step_u,
      label = p_to_star(p.value)
    ) %>%
    ungroup() %>%
    select(task_type, timeF, xmin, xmax, y, label)
}

plot_behavior_pretty <- function(data, dv = c("opinion_delta_abs","conf_delta","opinion_delta","conf_delta_abs"),
                                 model = NULL, p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  data <- data %>% mutate(
    pattern   = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  
  if (is.null(model)) {
    fml <- as.formula(paste0(dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"))
    model <- try(lmer(fml, data = data, REML = FALSE), silent = TRUE)
    if (inherits(model, "try-error") || isSingular(model, tol=1e-4)) {
      model <- lmer(update(fml, . ~ . - (1 + time_num | participant_id) + (1 | participant_id)), data = data, REML = FALSE)
    }
  }
  
  ylims  <- beh_y_limits(data, dv, behavior_top_expand_mult)
  sig_df <- make_sig_df_behavior(model, data, dv = dv, p_adjust = p_adjust, y_limits = ylims)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Δ Opinion| (0–100)",
    dv == "conf_delta_abs"    ~ "|Δ Confidence| (0–100)",
    dv == "opinion_delta"     ~ "Δ Opinion (Tk − T0)",
    TRUE                      ~ "Δ Confidence (Tk − T0)"
  )
  
  p <- ggplot(data, aes(x = timeF, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "line", aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8, na.rm = TRUE) +
    facet_wrap(~ task_type, nrow = 1) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0))) +
    labs(x = "Time (relative to T0)", y = ylab) +
    theme_gray(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 11),
      axis.title.x     = element_text(margin = margin(t = 8)),
      axis.title.y     = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
    ) +
    coord_cartesian(clip = "off")
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 4,
      vjust = -0.15, color = "black"
    )
  }
  p
}

# -------------------------------------------------
# B) Self-reported (Compliance / Conversion)
# -------------------------------------------------
make_sig_df_self <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                             p_adjust = p_adj_method, y_upper = 7 + self_top_pad) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  
  prs <- prs %>%
    tidyr::separate(contrast, into = c("g1","g2"), sep = " - ") %>%
    filter(p.value < 0.05) %>%
    mutate(
      g1 = as.character(factor(g1, levels = pattern_levels)),
      g2 = as.character(factor(g2, levels = pattern_levels))
    )
  if (!nrow(prs)) return(tibble())
  
  prs %>%
    mutate(
      x_idx = match(task_type, c("Normative","Informative")),
      xmin  = x_idx + as.numeric(pattern_offsets[g1]),
      xmax  = x_idx + as.numeric(pattern_offsets[g2])
    ) %>%
    group_by(task_type) %>%
    arrange(desc(abs(as.numeric(pattern_offsets[g1]) - as.numeric(pattern_offsets[g2])))) %>%
    mutate(
      y     = y_upper - self_bracket_top_pad - (row_number() - 1) * self_bracket_step_u,
      label = p_to_star(p.value)
    ) %>%
    ungroup() %>%
    select(task_type, xmin, xmax, y, label)
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
  
  y_upper <- 7 + y_top_pad
  sig_df  <- make_sig_df_self(model, data, dv = dv, p_adjust = p_adj_method, y_upper = y_upper)
  
  p <- ggplot(data, aes(x = task_type, y = .data[[dv]], fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8, na.rm = TRUE) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(1, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(
      x = "Task Type",
      y = ifelse(dv == "compliance_mean", "Perceived Compliance", "Perceived Conversion")
    ) +
    theme_gray(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 11),
      axis.title.x     = element_text(margin = margin(t = 8)),
      axis.title.y     = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
    ) +
    coord_cartesian(clip = "off")
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 4,
      vjust = -0.15, color = "black"
    )
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

make_sig_df_agent <- function(long_dat, p_adjust = p_adj_method, y_upper = 7 + agent_top_pad) {
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
      filter(p.value < 0.05) %>%
      mutate(
        g1 = as.character(factor(g1, levels = pattern_levels)),
        g2 = as.character(factor(g2, levels = pattern_levels))
      )
    if (!nrow(prs)) next
    
    out[[meas]] <- prs %>%
      mutate(
        x_idx = match(task_type, c("Normative","Informative")),
        xmin  = x_idx + as.numeric(pattern_offsets[g1]),
        xmax  = x_idx + as.numeric(pattern_offsets[g2])
      ) %>%
      group_by(task_type) %>%
      arrange(desc(abs(as.numeric(pattern_offsets[g1]) - as.numeric(pattern_offsets[g2])))) %>%
      mutate(
        y     = y_upper - agent_bracket_top_pad - (row_number() - 1) * agent_bracket_step_u,
        label = p_to_star(p.value)
      ) %>%
      ungroup() %>%
      transmute(measure = meas, task_type, xmin, xmax, y, label)
  }
  dplyr::bind_rows(out)
}

plot_agent_pretty <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long(post1)
  y_upper  <- 7 + y_top_pad
  sig_df   <- make_sig_df_agent(long_dat, p_adjust = p_adjust, y_upper = y_upper)
  
  p <- ggplot(long_dat, aes(x = task_type, y = score, fill = pattern, color = pattern)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.6,
                 outlier.shape = 1, outlier.size = 1.8, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = pattern),
                 position = position_dodge(width = dodge_w), size = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = dodge_w), size = 1.8, na.rm = TRUE) +
    facet_wrap(~ measure, nrow = 2, ncol = 4) +
    scale_fill_manual(values = pattern_fill, breaks = pattern_levels, name = NULL) +
    scale_color_manual(values = pattern_color, breaks = pattern_levels, name = NULL) +
    scale_y_continuous(breaks = 1:7, limits = c(1, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    theme_gray(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
    ) +
    coord_cartesian(clip = "off")
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 4,
      vjust = -0.15, color = "black"
    )
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

# Behavioral: Δ Confidence (signed)
fig_conf <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf else NULL,
  p_adjust = p_adj_method
)
print(fig_conf)

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
























###############################
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)

# -------------------------------------------------
# Global aesthetics and adjustable parameters
# -------------------------------------------------
pattern_levels <- c("Majority","Minority","Diffusion")
pattern_fill   <- c("Majority"="#e3f2fd", "Minority"="#ffe0b2", "Diffusion"="#e0f2f1")
pattern_color  <- c("Majority"="#1f77b4","Minority"="#ff7f0e","Diffusion"="#2ca02c")

dodge_w <- 0.65
pattern_offsets <- c(Majority = -0.22, Minority = 0.00, Diffusion = 0.22) # bracket anchor at box centers

# 브라켓/상단 여백 컨트롤(Top-anchored brackets)
behavior_top_expand_mult <- 0.10  # |Δ|는 0~100이므로, 위로 10% (= +10) 여백
self_top_pad             <- 0.45  # 1~7 스케일에서 7 위로 +0.45
agent_top_pad            <- 0.65  # 1~7 스케일에서 7 위로 +0.65

# 브라켓을 "패널 상단"에 고정 배치하기 위한 패딩/간격(데이터 단위)
behavior_bracket_top_pad <- 1.2   # 패널 상단(y_upper)에서 아래로 내리는 시작 위치
behavior_bracket_step_u  <- 1.1   # 브라켓 간 수직 간격
self_bracket_top_pad     <- 0.08  # 7 위 패딩(리커트)
self_bracket_step_u      <- 0.08
agent_bracket_top_pad    <- 0.08
agent_bracket_step_u     <- 0.08

# 별(★)을 브라켓 바로 위에 두기 위한 오프셋(데이터 단위)
behavior_star_offset_u <- 0.25
self_star_offset_u     <- 0.04
agent_star_offset_u    <- 0.04

# 다중비교 보정
#p_adj_method <- "holm"   # "holm" 또는 "bonferroni"
p_adj_method <- "bonferroni"   # "holm" 또는 "bonferroni"


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

# -------------------------------------------------
# A) Behavioral (Δ from T0) - Opinion / Confidence
# -------------------------------------------------
# 절대값이면 0~100(+여백) 강제
beh_y_limits <- function(data, dv, top_expand_mult = behavior_top_expand_mult) {
  if (grepl("_abs$", dv)) {
    ymax <- 100 + 100 * top_expand_mult
    return(c(0, ymax))
  } else {
    rng <- range(data[[dv]], na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)
    r <- diff(rng); if (r <= 0) r <- 1
    return(c(rng[1] - r*0.02, rng[2] + r*top_expand_mult))
  }
}

# (1) 행동 델타 브라켓 좌표 생성 함수 교체
make_sig_df_behavior <- function(model, data, dv = c("opinion_delta","conf_delta","opinion_delta_abs"),
                                 p_adjust = p_adj_method, y_limits = NULL) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type * timeF, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- tidyr::separate(prs, contrast, into = c("g1","g2"), sep = " - ")
  prs <- prs %>% filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  # facet-wise stats
  sum_df <- data %>%
    group_by(pattern, task_type, timeF) %>%
    summarise(m = mean(.data[[dv]], na.rm = TRUE), .groups = "drop")
  
  # 글로벌 y-limits 확보(패널 밖으로 못 나가게만 clamp)
  if (is.null(y_limits)) {
    rng <- range(data[[dv]], na.rm = TRUE); if (!all(is.finite(rng))) rng <- c(-1,1)
    r <- diff(rng); if (r <= 0) r <- 1
    y_limits <- c(rng[1] - 0.02*r, rng[2] + 0.25*r)  # 여기의 0.25는 이전 설정과 일치시켜 주세요
  }
  y_upper <- y_limits[2]
  y_range <- diff(y_limits)
  
  ybase <- sum_df %>%
    group_by(task_type, timeF) %>%
    summarise(y_max = max(m, na.rm = TRUE), y_min = min(m, na.rm = TRUE), .groups = "drop") %>%
    mutate(margin = pmax((y_max - y_min)*0.18, 0.6))
  
  sig_df <- prs %>%
    mutate(
      x_idx = as.numeric(factor(timeF, levels = levels(data$timeF))),
      xmin = x_idx + pattern_offsets[g1],
      xmax = x_idx + pattern_offsets[g2]
    ) %>%
    left_join(ybase, by = c("task_type","timeF")) %>%
    group_by(task_type, timeF) %>%
    arrange(p.value) %>%
    mutate(
      y_desired = y_max + margin + (row_number()-1)*margin*behavior_bracket_step,
      # 오프셋만 추가(여백은 그대로)
      y_shifted = y_desired + behavior_bracket_bump,
      # 패널 밖으로 나가지 않게 살짝만 클램프
      y = pmin(y_shifted, y_upper - y_range*0.01),
      label = p_to_star(p.value)
    ) %>%
    ungroup() %>%
    select(task_type, timeF, xmin, xmax, y, label)
  sig_df
}

# emmeans 대비 브라켓 좌표 함수는 기존 것을 사용해도 OK
# 단, plot 내부의 수식에 * 사용
plot_behavior_pretty <- function(data, dv = c("opinion_delta_abs","conf_delta_abs","opinion_delta","conf_delta"),
                                 model = NULL, p_adjust = p_adj_method) {
  dv <- match.arg(dv)
  data <- data %>% mutate(
    pattern = factor(pattern, levels = pattern_levels),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
  if (is.null(model)) {
    fml <- as.formula(paste0(dv, " ~ pattern * task_type * timeF + SII_z + NFC_z + AIacc_z + (1 + time_num | participant_id)"))
    # 간단 버전 안전적합(원 코드의 fit_lmm_safe 써도 OK)
    model <- try(lmer(fml, data = data, REML = FALSE), silent = TRUE)
    if (inherits(model, "try-error") || isSingular(model, tol=1e-4)) {
      model <- lmer(update(fml, . ~ . - (1 + time_num | participant_id) + (1 | participant_id)), data = data, REML = FALSE)
    }
  }
  ylims <- beh_y_limits(data, dv, behavior_top_expand_mult)
  sig_df <- make_sig_df_behavior(model, data, dv = dv, p_adjust = p_adjust, y_limits = ylims)
  
  ylab <- dplyr::case_when(
    dv == "opinion_delta_abs" ~ "|Δ Opinion| (|Tk − T0|)",
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
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0.0, 0.0))) +
    labs(x = "Time (relative to T0)", y = ylab) +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 3.8, vjust = 0.1, color = "black"
    )
  }
  p
}

# -------------------------------------------------
# B) Self-reported (Compliance / Conversion)
# -------------------------------------------------
# (2) self-reported 브라켓 좌표 생성 함수 교체
make_sig_df_self <- function(model, data, dv = c("compliance_mean","conversion_mean"),
                             p_adjust = p_adj_method, y_upper = 7 + self_top_pad) {
  dv <- match.arg(dv)
  emm <- emmeans(model, ~ pattern | task_type, data = data)
  prs <- as.data.frame(pairs(emm, adjust = p_adjust))
  if (!nrow(prs)) return(tibble())
  prs <- tidyr::separate(prs, contrast, into = c("g1","g2"), sep = " - ")
  prs <- prs %>% filter(p.value < 0.05)
  if (!nrow(prs)) return(tibble())
  
  sum_df <- data %>%
    group_by(pattern, task_type) %>%
    summarise(m = mean(.data[[dv]], na.rm = TRUE), .groups = "drop")
  ybase <- sum_df %>%
    group_by(task_type) %>%
    summarise(y_max = max(m, na.rm = TRUE), y_min = min(m, na.rm = TRUE), .groups = "drop") %>%
    mutate(margin = pmax((y_max - y_min)*0.20, 0.25))
  
  # 1~7 스케일에서 패널 상단 살짝 아래까지 허용
  clamp_top <- y_upper - 0.05
  
  sig_df <- prs %>%
    mutate(
      x_idx = match(task_type, c("Normative","Informative")),
      xmin = x_idx + pattern_offsets[g1],
      xmax = x_idx + pattern_offsets[g2]
    ) %>%
    left_join(ybase, by = "task_type") %>%
    group_by(task_type) %>%
    arrange(p.value) %>%
    mutate(
      y_desired = y_max + margin + (row_number()-1)*margin*self_bracket_step,
      y_shifted = y_desired + self_bracket_bump,   # 오프셋 추가
      y = pmin(y_shifted, clamp_top),
      label = p_to_star(p.value)
    ) %>%
    ungroup() %>%
    select(task_type, xmin, xmax, y, label)
  sig_df
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
  
  y_upper <- 7 + y_top_pad
  sig_df <- make_sig_df_self(model, data, dv = dv, p_adjust = p_adjust, y_upper = y_upper)
  
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
    scale_y_continuous(breaks = 1:7, limits = c(1, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(
      x = "Task Type",
      y = ifelse(dv == "compliance_mean", "Perceived Compliance", "Perceived Conversion")
    ) +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 3.8, vjust = 0.1, color = "black"
    )
  }
  p
}

# -------------------------------------------------
# C) Agent Perception (7 measures) - 2행 4열 facet
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

# (3) agent 브라켓 좌표 생성 함수 교체
make_sig_df_agent <- function(long_dat, p_adjust = p_adj_method, y_upper = 7 + agent_top_pad) {
  out <- list()
  for (meas in agent_dims) {
    d <- long_dat %>% filter(measure == meas)
    if (!nrow(d)) next
    m <- lmer(score ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id), data = d, REML = TRUE)
    emm <- emmeans(m, ~ pattern | task_type, data = d)
    prs <- as.data.frame(pairs(emm, adjust = p_adjust))
    if (!nrow(prs)) next
    prs <- tidyr::separate(prs, contrast, into = c("g1","g2"), sep = " - ")
    prs <- prs %>% filter(p.value < 0.05)
    if (!nrow(prs)) next
    
    sum_df <- d %>%
      group_by(pattern, task_type) %>%
      summarise(m = mean(score, na.rm = TRUE), .groups = "drop")
    ybase <- sum_df %>%
      group_by(task_type) %>%
      summarise(y_max = max(m, na.rm = TRUE), y_min = min(m, na.rm = TRUE), .groups = "drop") %>%
      mutate(margin = pmax((y_max - y_min)*0.20, 0.25))
    
    clamp_top <- y_upper - 0.05
    
    sig_df <- prs %>%
      mutate(
        x_idx = match(task_type, c("Normative","Informative")),
        xmin = x_idx + pattern_offsets[g1],
        xmax = x_idx + pattern_offsets[g2]
      ) %>%
      left_join(ybase, by = "task_type") %>%
      group_by(task_type) %>%
      arrange(p.value) %>%
      mutate(
        y_desired = y_max + margin + (row_number()-1)*margin*agent_bracket_step,
        y_shifted = y_desired + agent_bracket_bump,  # 오프셋 추가
        y = pmin(y_shifted, clamp_top),
        label = p_to_star(p.value)
      ) %>%
      ungroup() %>%
      transmute(measure = meas, task_type, xmin, xmax, y, label)
    out[[meas]] <- sig_df
  }
  dplyr::bind_rows(out)
}

plot_agent_pretty <- function(post1, p_adjust = p_adj_method, y_top_pad = agent_top_pad) {
  long_dat <- make_agent_long(post1)
  y_upper <- 7 + y_top_pad
  sig_df <- make_sig_df_agent(long_dat, p_adjust = p_adjust, y_upper = y_upper)
  
  # 박스플롯 + 평균 라인/점
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
    scale_y_continuous(breaks = 1:7, limits = c(1, y_upper), expand = expansion(mult = c(0, 0))) +
    labs(x = "Task Type", y = "Agent Perception (Likert 1–7)") +
    theme_gray(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(sig_df)) {
    p <- p + ggsignif::geom_signif(
      data = sig_df,
      aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y),
      manual = TRUE, inherit.aes = FALSE,
      tip_length = 0.01, textsize = 3.8, vjust = 0.1, color = "black"
    )
  }
  p
}

# -------------------------------------------------
# Run plots (after your preprocessing created resp_delta, post1)
# -------------------------------------------------

# 브라켓/상단 여백 컨트롤
behavior_top_expand_mult <- 0.25   # 행동 플롯: y 상단 여백(데이터 범위의 25%)
behavior_bracket_step    <- 0.90   # 행동 플롯: 브라켓 간 수직 간격(여백의 60%)
self_top_pad             <- 0.45   # self-reported: y 상단 여백(7 위로 0.45)
self_bracket_step        <- 0.70   # self-reported: 브라켓 간 수직 간격(여백의 70%)
agent_top_pad            <- 0.65   # agent: y 상단 여백(7 위로 0.45)
agent_bracket_step       <- 0.70   # agent: 브라켓 간 수직 간격(여백의 70%)

# 상단 여백(절대값은 0~100이라 10%면 110까지)
behavior_top_expand_mult <- 0.10

# Behavioral: |Δ Opinion|
fig_opinion_abs <- plot_behavior_pretty(
  data = resp_delta, dv = "opinion_delta_abs",
  model = if (exists("m_opinion_abs")) m_opinion_abs else NULL,
  p_adjust = p_adj_method
)
print(fig_opinion_abs)

# Behavioral: Δ Confidence
fig_conf <- plot_behavior_pretty(
  data = resp_delta, dv = "conf_delta",
  model = if (exists("m_conf")) m_conf else NULL,
  p_adjust = p_adj_method
)
print(fig_conf)

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

