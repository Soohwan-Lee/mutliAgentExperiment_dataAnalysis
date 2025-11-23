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
# resp 원본 보존하고 싶으면 (선택 사항)
resp <- resp %>%
  mutate(opinion_raw = opinion)   # 원본 opinion 따로 저장(선택)

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
  mutate(
    t_ord = row_number() - 1L,
    # ---- 여기부터 새로 추가: T0가 음수면 0~4 모두 부호 뒤집기 ----
    op_sign_T0 = sign(opinion[t_ord == 0][1]),
    op_sign_T0 = ifelse(is.na(op_sign_T0), 0, op_sign_T0),
    opinion    = dplyr::if_else(op_sign_T0 < 0, -opinion, opinion)
  ) %>%
  ungroup() %>%
  select(-op_sign_T0) %>%
  mutate(
    pattern   = factor(pattern,   levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )
# -------------------------------------------------------
# 4-1) 행동 데이터: 원값(T0~T4) + 공변량 조인
# -------------------------------------------------------
resp_behavior <- resp1 %>%
  mutate(
    time_num = as.integer(t_ord),                               # 0..4
    timeF    = factor(time_num, levels = 0:4,
                      labels = c("T0","T1","T2","T3","T4"))
  ) %>%
  left_join(bg_cov, by = "participant_id")


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


# time 중심화 + 그룹변수 factor 보장
resp_behavior <- resp_behavior %>%
  mutate(time_num_c = as.numeric(scale(time_num, scale = FALSE)))
resp_behavior$participant_id <- factor(resp_behavior$participant_id)

resp_delta <- resp_delta %>%
  mutate(time_num_c = as.numeric(scale(time_num, scale = FALSE)))
resp_delta$participant_id <- factor(resp_delta$participant_id)

post1$participant_id <- factor(post1$participant_id)


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

# ---------- 8-1) 행동 원값 LMM (opinion / confidence) ----------
form_op_raw_full   <- opinion ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_op_raw_uncorr <- opinion ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_op_raw_int    <- opinion ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_raw_full   <- confidence ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_cf_raw_uncorr <- confidence ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_cf_raw_int    <- confidence ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

m_opinion_raw <- fit_lmm_step_general(
  form_op_raw_full, form_op_raw_uncorr, form_op_raw_int,
  data = resp_behavior, reml = FALSE
)
m_conf_raw <- fit_lmm_step_general(
  form_cf_raw_full, form_cf_raw_uncorr, form_cf_raw_int,
  data = resp_behavior, reml = FALSE
)


### Opinion Raw - Signed
cat("\n=== Opinion (raw T0~T4) ===\n")
print(summary(m_opinion_raw));print(anova(m_opinion_raw, type = 3))

# Post-hoc - Pattern
emm_op_raw_pattern <- emmeans(m_opinion_raw,
                      specs = "pattern",
                      by   = c("task_type", "timeF"))
pairs(emm_op_raw_pattern, adjust = "bonferroni")
emmeans::eff_size(emm_op_raw_pattern, sigma = sigma(m_opinion_raw),
                  edf = df.residual(m_opinion_raw))
# Post-hoc - Task_type
emm_op_raw_task_type <- emmeans(m_opinion_raw,
                      specs = "task_type",
                      by   = c("pattern", "timeF"))
pairs(emm_op_raw_task_type, adjust = "bonferroni")
emmeans::eff_size(emm_op_raw_task_type, sigma = sigma(m_opinion_raw),
                  edf = df.residual(m_opinion_raw))
# Post-hoc - timeF
emm_op_raw_timeF <- emmeans(m_opinion_raw,
                      specs = "timeF",
                      by   = c("pattern", "task_type"))
pairs(emm_op_raw_timeF, adjust = "bonferroni")
emmeans::eff_size(emm_op_raw_timeF, sigma = sigma(m_opinion_raw),
                  edf = df.residual(m_opinion_raw))


### Confidence Raw - Signed
cat("\n=== Confidence (raw T0~T4) ===\n")
print(summary(m_conf_raw));    print(anova(m_conf_raw, type = 3))

# Post-hoc - Pattern
emm_cf_raw_pattern <- emmeans(m_conf_raw,
                      specs = "pattern",
                      by   = c("task_type", "timeF"))
pairs(emm_cf_raw_pattern, adjust = "bonferroni")
emmeans::eff_size(emm_cf_raw_pattern, sigma = sigma(m_conf_raw),
                  edf = df.residual(m_conf_raw))
# Post-hoc - task_type
emm_cf_raw_task_type <- emmeans(m_conf_raw,
                      specs = "task_type",
                      by   = c("pattern", "timeF"))
pairs(emm_cf_raw_task_type, adjust = "bonferroni")
emmeans::eff_size(emm_cf_raw_task_type, sigma = sigma(m_conf_raw),
                  edf = df.residual(m_conf_raw))
# Post-hoc - timeF
emm_cf_raw_timeF <- emmeans(m_conf_raw,
                      specs = "timeF",
                      by   = c("pattern", "task_type"))
pairs(emm_cf_raw_timeF, adjust = "bonferroni")
emmeans::eff_size(emm_cf_raw_timeF, sigma = sigma(m_conf_raw),
                  edf = df.residual(m_conf_raw))



# ---------- 8-2) 행동 델타 LMM (signed delta: ΔOpinion, ΔConfidence) ----------
form_op_signed_full   <- opinion_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_op_signed_uncorr <- opinion_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_op_signed_int    <- opinion_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_signed_full   <- conf_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_cf_signed_uncorr <- conf_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_cf_signed_int    <- conf_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

m_opinion_signed <- fit_lmm_step_general(
  form_op_signed_full, form_op_signed_uncorr, form_op_signed_int,
  data = resp_delta, reml = FALSE
)
m_conf_signed <- fit_lmm_step_general(
  form_cf_signed_full, form_cf_signed_uncorr, form_cf_signed_int,
  data = resp_delta, reml = FALSE
)


### Opinion Delta - Signed
cat("\n=== ΔOpinion (signed, Tk−T0) ===\n")
print(summary(m_opinion_signed)); print(anova(m_opinion_signed, type = 3))

# Post-hoc - Pattern
emm_op_signed_pattern <- emmeans(m_opinion_signed, specs = "pattern",
                         by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_op_signed_pattern, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_signed_pattern, sigma = sigma(m_opinion_signed),
                        edf = df.residual(m_opinion_signed)))
# Post-hoc - Task_type
emm_op_signed_task_type <- emmeans(m_opinion_signed, specs = "task_type",
                         by = c("pattern","timeF"), data = resp_delta)
print(pairs(emm_op_signed_task_type, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_signed_task_type, sigma = sigma(m_opinion_signed),
                        edf = df.residual(m_opinion_signed)))
# Post-hoc - timeF
emm_op_signed_timeF <- emmeans(m_opinion_signed, specs = "timeF",
                         by = c("pattern","task_type"), data = resp_delta)
print(pairs(emm_op_signed_timeF, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_signed_timeF, sigma = sigma(m_opinion_signed),
                        edf = df.residual(m_opinion_signed)))


### Confidence Delta- Signed
cat("\n=== ΔConfidence (signed, Tk−T0) ===\n")
print(summary(m_conf_signed));    print(anova(m_conf_signed, type = 3))

# Post-hoc - Pattern
emm_cf_signed_pattern <- emmeans(m_conf_signed, specs = "pattern",
                      by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_cf_signed_pattern, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_signed_pattern, sigma = sigma(m_conf_signed),
                        edf = df.residual(m_conf_signed)))
# Post-hoc - task_type
emm_cf_signed_task_type <- emmeans(m_conf_signed, specs = "task_type",
                         by = c("pattern","timeF"), data = resp_delta)
print(pairs(emm_cf_signed_task_type, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_signed_task_type, sigma = sigma(m_conf_signed),
                        edf = df.residual(m_conf_signed)))
# Post-hoc - timeF
emm_cf_signed_timeF <- emmeans(m_conf_signed, specs = "timeF",
                         by = c("task_type","pattern"), data = resp_delta)
print(pairs(emm_cf_signed_timeF, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_signed_timeF, sigma = sigma(m_conf_signed),
                        edf = df.residual(m_conf_signed)))



# ---------- 8-3) 행동 델타 LMM (absolute delta: |ΔOpinion|, |ΔConfidence|) ----------
form_op_abs_full   <- opinion_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_op_abs_uncorr <- opinion_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_op_abs_int    <- opinion_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_abs_full   <- conf_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c |  participant_id)
form_cf_abs_uncorr <- conf_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 + time_num_c || participant_id)
form_cf_abs_int    <- conf_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

m_opinion_abs <- fit_lmm_step_general(
  form_op_abs_full, form_op_abs_uncorr, form_op_abs_int,
  data = resp_delta, reml = FALSE
)
m_conf_abs <- fit_lmm_step_general(
  form_cf_abs_full, form_cf_abs_uncorr, form_cf_abs_int,
  data = resp_delta, reml = FALSE
)

### Opinion Delta - Abs
cat("\n=== |ΔOpinion| (Tk−T0) ===\n")
print(summary(m_opinion_abs)); print(anova(m_opinion_abs, type = 3))

# Post-hoc - Pattern
emm_op_abs_pattern <- emmeans(m_opinion_abs, specs = "pattern",
                      by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_op_abs_pattern, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_abs_pattern, sigma = sigma(m_opinion_abs),
                        edf = df.residual(m_opinion_abs)))
# Post-hoc - task_type
emm_op_abs_task_type <- emmeans(m_opinion_abs, specs = "task_type",
                      by = c("pattern","timeF"), data = resp_delta)
print(pairs(emm_op_abs_task_type, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_abs_task_type, sigma = sigma(m_opinion_abs),
                        edf = df.residual(m_opinion_abs)))
# Post-hoc - timeF
emm_op_abs_timeF <- emmeans(m_opinion_abs, specs = "timeF",
                      by = c("task_type","pattern"), data = resp_delta)
print(pairs(emm_op_abs_timeF, adjust = "bonferroni"))
print(emmeans::eff_size(emm_op_abs_timeF, sigma = sigma(m_opinion_abs),
                        edf = df.residual(m_opinion_abs)))




### Confidence Delta - Abs
cat("\n=== |ΔConfidence| (Tk−T0) ===\n")
print(summary(m_conf_abs));    print(anova(m_conf_abs, type = 3))

# Post-hoc - Pattern
emm_cf_abs_pattern <- emmeans(m_conf_abs, specs = "pattern",
                      by = c("task_type","timeF"), data = resp_delta)
print(pairs(emm_cf_abs_pattern, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_abs_pattern, sigma = sigma(m_conf_abs),
                        edf = df.residual(m_conf_abs)))
# Post-hoc - task_type
emm_cf_abs_task_type <- emmeans(m_conf_abs, specs = "task_type",
                      by = c("pattern","timeF"), data = resp_delta)
print(pairs(emm_cf_abs_task_type, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_abs_task_type, sigma = sigma(m_conf_abs),
                        edf = df.residual(m_conf_abs)))
# Post-hoc - timeF
emm_cf_abs_timeF <- emmeans(m_conf_abs, specs = "timeF",
                      by = c("task_type","pattern"), data = resp_delta)
print(pairs(emm_cf_abs_timeF, adjust = "bonferroni"))
print(emmeans::eff_size(emm_cf_abs_timeF, sigma = sigma(m_conf_abs),
                        edf = df.residual(m_conf_abs)))




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

