###############################################
# Multi-Agent Experiment
# Trajectory-based Opinion / Confidence Analysis
#  - 데이터 로드부터 Trajectory LMM 분석까지 하나의 스크립트
#  - contrasts: contr.sum / contr.poly
#  - 랜덤효과: (1 | participant_id) 로 통일
#  - EMMs + Bonferroni post-hoc + 효과크기 출력
###############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(effectsize)
})

## -------------------------------------------------------
## 0) 옵션 및 경로 설정
## -------------------------------------------------------

## Type-III ANOVA + sum / polynomial contrasts
options(contrasts = c("contr.sum", "contr.poly"))

# base_dir <- "./Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data/"   # MAC
base_dir <- "C:\\Users\\USER\\Desktop\\github\\mutliAgentExperiment_dataAnalysis\\refined_data\\"   # Windows

out_dir <- base_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


## -------------------------------------------------------
## 1) 데이터 로드
##    - 배경설문(bg), 행동(resp)만 사용 (post는 여기선 필요X)
## -------------------------------------------------------

bg   <- readr::read_csv(file.path(base_dir, "FINAL_background_surveys_rows_refined.csv"))
resp <- readr::read_csv(file.path(base_dir, "FINAL_responses_rows_refined_experiment_only.csv"))

## ID를 문자형으로 통일
bg$participant_id   <- as.character(bg$participant_id)
resp$participant_id <- as.character(resp$participant_id)


## -------------------------------------------------------
## 2) 공변량 생성: 문항 평균 후 z-표준화 (SII, NFC, AI 수용성)
## -------------------------------------------------------

get_mean_scale <- function(df, pattern_regex) {
  cols <- grep(pattern_regex, names(df), value = TRUE)
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  out <- rowMeans(dplyr::select(df, all_of(cols)), na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
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


## -------------------------------------------------------
## 3) 행동 데이터 전처리
##    - 패턴/과업 지정, 시간 순서(t_ord=0..4), opinion sign-flip
## -------------------------------------------------------

resp <- resp %>%
  mutate(opinion_raw = opinion)   # 원본 opinion 저장(선택)

resp1 <- resp %>%
  mutate(
    pattern = dplyr::recode(
      condition_type,
      "majority"          = "Majority",
      "minority"          = "Minority",
      "minorityDiffusion" = "Diffusion"
    ),
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
    t_ord = row_number() - 1L,  # 0..4
    # T0가 음수면 0~4 opinion 부호 뒤집기
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


## -------------------------------------------------------
## 4) 행동 데이터: 원값(T0~T4) + 공변량 조인 → resp_behavior
## -------------------------------------------------------

resp_behavior <- resp1 %>%
  mutate(
    time_num = as.integer(t_ord),                               # 0..4
    timeF    = factor(time_num, levels = 0:4,
                      labels = c("T0","T1","T2","T3","T4"))
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    participant_id = factor(participant_id)
  )


## -------------------------------------------------------
## 5) 베이스라인(T0) 대비 델타(T1~T4) + 절대값 델타 → resp_delta
## -------------------------------------------------------

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
    timeF    = factor(time_num, levels = 1:4,
                      labels = c("T1","T2","T3","T4")),
    # 절대값 델타
    opinion_delta_abs = abs(opinion_delta),
    conf_delta_abs    = abs(conf_delta)
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    participant_id = factor(participant_id)
  )


## -------------------------------------------------------
## 6) Trajectory 분류 (opinion / confidence)
##    - 단위: participant_id × task_type
##    - 그룹: always_pos / always_neg / mixed
## -------------------------------------------------------

cat("\n\n==============================\n")
cat("Trajectory-based grouping (opinion / confidence)\n")
cat("==============================\n")

# 6-1) Trajectory 분류 함수
#  - x: T0~T4 값 벡터
#  - diff(x)를 계산한 후:
#     * 0은 무시
#     * non-zero diff가 모두 >0  → always_pos
#     * non-zero diff가 모두 <0  → always_neg
#     * 양/음 섞이거나, 모든 diff가 0 → mixed
classify_traj <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_character_)
  d <- diff(x)
  d_nz <- d[d != 0]
  if (length(d_nz) == 0L) {
    return("mixed")
  }
  if (all(d_nz > 0)) return("always_pos")
  if (all(d_nz < 0)) return("always_neg")
  return("mixed")
}

time_levels <- paste0("T", 0:4)


## 6-2) Opinion trajectory 분류 (resp_behavior$opinion 사용)
##      - opinion은 이미 sign-flip된 값이라고 가정

traj_opinion <- resp_behavior %>%
  filter(timeF %in% time_levels) %>%
  mutate(timeF = factor(timeF, levels = time_levels)) %>%
  arrange(participant_id, task_type, timeF) %>%
  group_by(participant_id, task_type) %>%
  summarise(
    pattern     = dplyr::first(pattern),
    n_t_unique  = n_distinct(timeF),
    opinion_seq = list(opinion),
    .groups     = "drop"
  )

# 시점이 5개가 아닌 경우 경고 출력
traj_opinion_problem <- traj_opinion %>% filter(n_t_unique != 5L)
if (nrow(traj_opinion_problem) > 0) {
  cat("\n[경고] Opinion: T0~T4 5시점이 아닌 (participant_id, task_type) 조합이 있습니다:\n")
  print(traj_opinion_problem %>% select(participant_id, task_type, pattern, n_t_unique))
}

traj_opinion <- traj_opinion %>%
  filter(n_t_unique == 5L) %>%
  mutate(
    opinion_traj_group = purrr::map_chr(opinion_seq, classify_traj),
    opinion_traj_group = factor(
      opinion_traj_group,
      levels = c("always_pos", "always_neg", "mixed")
    )
  )


## 6-3) Confidence trajectory 분류 (resp_behavior$confidence 사용)

traj_conf <- resp_behavior %>%
  filter(timeF %in% time_levels) %>%
  mutate(timeF = factor(timeF, levels = time_levels)) %>%
  arrange(participant_id, task_type, timeF) %>%
  group_by(participant_id, task_type) %>%
  summarise(
    pattern    = dplyr::first(pattern),
    n_t_unique = n_distinct(timeF),
    conf_seq   = list(confidence),
    .groups    = "drop"
  )

traj_conf_problem <- traj_conf %>% filter(n_t_unique != 5L)
if (nrow(traj_conf_problem) > 0) {
  cat("\n[경고] Confidence: T0~T4 5시점이 아닌 (participant_id, task_type) 조합이 있습니다:\n")
  print(traj_conf_problem %>% select(participant_id, task_type, pattern, n_t_unique))
}

traj_conf <- traj_conf %>%
  filter(n_t_unique == 5L) %>%
  mutate(
    conf_traj_group = purrr::map_chr(conf_seq, classify_traj),
    conf_traj_group = factor(
      conf_traj_group,
      levels = c("always_pos", "always_neg", "mixed")
    )
  )


## 6-4) Trajectory 그룹 요약 (count + percentage 출력)

cat("\n=== Opinion trajectory groups (단위: participant × task_type, 즉 trajectory 수) ===\n")

opinion_group_overall <- traj_opinion %>%
  count(opinion_traj_group, name = "n_traj") %>%
  mutate(
    prop = n_traj / sum(n_traj),
    pct  = round(100 * prop, 1)
  )
print(opinion_group_overall)

cat("\n--- Opinion trajectory groups by task_type (각 task 안에서 참가자 비율) ---\n")
opinion_group_by_task <- traj_opinion %>%
  count(task_type, opinion_traj_group, name = "n_participants") %>%
  group_by(task_type) %>%
  mutate(
    prop = n_participants / sum(n_participants),
    pct  = round(100 * prop, 1)
  ) %>%
  ungroup()
print(opinion_group_by_task)

cat("\n--- Opinion trajectory groups: pattern × task_type (trajectory 기준) ---\n")
opinion_group_by_pattern_task <- traj_opinion %>%
  count(pattern, task_type, opinion_traj_group, name = "n_traj") %>%
  group_by(pattern, task_type) %>%
  mutate(
    prop = n_traj / sum(n_traj),
    pct  = round(100 * prop, 1)
  ) %>%
  ungroup()
print(opinion_group_by_pattern_task)


cat("\n=== Confidence trajectory groups (단위: participant × task_type, 즉 trajectory 수) ===\n")

conf_group_overall <- traj_conf %>%
  count(conf_traj_group, name = "n_traj") %>%
  mutate(
    prop = n_traj / sum(n_traj),
    pct  = round(100 * prop, 1)
  )
print(conf_group_overall)

cat("\n--- Confidence trajectory groups by task_type (각 task 안에서 참가자 비율) ---\n")
conf_group_by_task <- traj_conf %>%
  count(task_type, conf_traj_group, name = "n_participants") %>%
  group_by(task_type) %>%
  mutate(
    prop = n_participants / sum(n_participants),
    pct  = round(100 * prop, 1)
  ) %>%
  ungroup()
print(conf_group_by_task)

cat("\n--- Confidence trajectory groups: pattern × task_type (trajectory 기준) ---\n")
conf_group_by_pattern_task <- traj_conf %>%
  count(pattern, task_type, conf_traj_group, name = "n_traj") %>%
  group_by(pattern, task_type) %>%
  mutate(
    prop = n_traj / sum(n_traj),
    pct  = round(100 * prop, 1)
  ) %>%
  ungroup()
print(conf_group_by_pattern_task)


## 6-5) participant_no 목록 CSV (resp_behavior에 존재할 때만)

if ("participant_no" %in% names(resp_behavior)) {
  opinion_ids_long <- traj_opinion %>%
    left_join(
      resp_behavior %>% distinct(participant_id, participant_no),
      by = "participant_id"
    ) %>%
    distinct(task_type, opinion_traj_group, participant_no) %>%
    arrange(task_type, opinion_traj_group, participant_no)
  
  opinion_csv_path <- file.path(out_dir, "FINAL_3Trajectory_Opinion-ParticipantNo.csv")
  write.csv(opinion_ids_long, opinion_csv_path,
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n[파일 저장 완료] Opinion ID 목록 저장 위치:\n", opinion_csv_path, "\n")
  
  conf_ids_long <- traj_conf %>%
    left_join(
      resp_behavior %>% distinct(participant_id, participant_no),
      by = "participant_id"
    ) %>%
    distinct(task_type, conf_traj_group, participant_no) %>%
    arrange(task_type, conf_traj_group, participant_no)
  
  conf_csv_path <- file.path(out_dir, "FINAL_3Trajectory_Confidence-ParticipantNo.csv")
  write.csv(conf_ids_long, conf_csv_path,
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n[파일 저장 완료] Confidence ID 목록 저장 위치:\n", conf_csv_path, "\n")
} else {
  cat("\n[알림] resp_behavior에 participant_no 변수가 없어 ID CSV는 저장하지 않습니다.\n")
}


## -------------------------------------------------------
## 7) LMM 포뮬러 정의
##    - 랜덤효과는 (1 | participant_id) 로 통일
##    - 예시 스크립트 구조와 동일
## -------------------------------------------------------

## Opinion (raw / signed delta / abs delta)
form_op_raw <- opinion ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_op_signed <- opinion_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_op_abs <- opinion_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

## Confidence (raw / signed delta / abs delta)
form_cf_raw <- confidence ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_signed <- conf_delta ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)

form_cf_abs <- conf_delta_abs ~ pattern * task_type * timeF +
  SII_z + NFC_z + AIacc_z + (1 | participant_id)









###########################################################
# 8) OPINION trajectory 그룹별 분석
#    - 각 그룹("always_pos","always_neg","mixed")를
#      차례로 아래 블록들을 실행하면 됨
###########################################################


###########################################################
# 8-1) OPINION: 그룹 = "always_pos"
###########################################################

## [STEP 0] 분석할 Opinion trajectory 그룹 선택
op_group <- "always_pos"

## [STEP 1] 이 그룹에 속하는 (participant_id, task_type) 추출
ids_op <- traj_opinion %>%
  dplyr::filter(opinion_traj_group == op_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_op$participant_id)), "\n")

## [STEP 2] 이 그룹에 해당하는 데이터 만들기
dat_op_raw <- resp_behavior %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_op_delta <- resp_delta %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")

## [STEP 3] 모델 1: Opinion raw (T0~T4)
cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")

m_op_raw <- lmer(
  form_op_raw,
  data = dat_op_raw,
  REML = FALSE
)

cat("\n[요약]\n")
print(summary(m_op_raw))

cat("\n[Type-III ANOVA]\n")
print(anova(m_op_raw, type = 3))

cat("\n[EMMs] Opinion raw: pattern × task_type × timeF\n")
emm_op_raw <- emmeans(m_op_raw, ~ pattern * task_type * timeF)
print(emm_op_raw)

## Post-hoc 1: pattern
cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
emm_op_raw_pattern <- emmeans(
  m_op_raw,
  specs = "pattern",
  by    = c("task_type", "timeF"),
  data  = dat_op_raw
)
print(emm_op_raw_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_raw_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(
  emm_op_raw_pattern,
  sigma = sigma(m_op_raw),
  edf   = df.residual(m_op_raw)
))

## Post-hoc 2: task_type
cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
emm_op_raw_task <- emmeans(
  m_op_raw,
  specs = "task_type",
  by    = c("pattern", "timeF"),
  data  = dat_op_raw
)
print(emm_op_raw_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_raw_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(
  emm_op_raw_task,
  sigma = sigma(m_op_raw),
  edf   = df.residual(m_op_raw)
))

## Post-hoc 3: timeF
cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
emm_op_raw_time <- emmeans(
  m_op_raw,
  specs = "timeF",
  by    = c("pattern", "task_type"),
  data  = dat_op_raw
)
print(emm_op_raw_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_raw_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(
  emm_op_raw_time,
  sigma = sigma(m_op_raw),
  edf   = df.residual(m_op_raw)
))

## [STEP 4] 모델 2: ΔOpinion (signed, Tk−T0)
cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")

m_op_signed <- lmer(
  form_op_signed,
  data = dat_op_delta,
  REML = FALSE
)

cat("\n[요약]\n")
print(summary(m_op_signed))

cat("\n[Type-III ANOVA]\n")
print(anova(m_op_signed, type = 3))

cat("\n[EMMs] ΔOpinion signed: pattern × task_type × timeF\n")
emm_op_signed <- emmeans(m_op_signed, ~ pattern * task_type * timeF)
print(emm_op_signed)

## Post-hoc 1: pattern
cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
emm_op_signed_pattern <- emmeans(
  m_op_signed,
  specs = "pattern",
  by    = c("task_type", "timeF"),
  data  = dat_op_delta
)
print(emm_op_signed_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_signed_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(
  emm_op_signed_pattern,
  sigma = sigma(m_op_signed),
  edf   = df.residual(m_op_signed)
))

## Post-hoc 2: task_type
cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
emm_op_signed_task <- emmeans(
  m_op_signed,
  specs = "task_type",
  by    = c("pattern", "timeF"),
  data  = dat_op_delta
)
print(emm_op_signed_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_signed_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(
  emm_op_signed_task,
  sigma = sigma(m_op_signed),
  edf   = df.residual(m_op_signed)
))

## Post-hoc 3: timeF
cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
emm_op_signed_time <- emmeans(
  m_op_signed,
  specs = "timeF",
  by    = c("pattern", "task_type"),
  data  = dat_op_delta
)
print(emm_op_signed_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_signed_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(
  emm_op_signed_time,
  sigma = sigma(m_op_signed),
  edf   = df.residual(m_op_signed)
))

## [STEP 5] 모델 3: |ΔOpinion| (absolute delta)
cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")

m_op_abs <- lmer(
  form_op_abs,
  data = dat_op_delta,
  REML = FALSE
)

cat("\n[요약]\n")
print(summary(m_op_abs))

cat("\n[Type-III ANOVA]\n")
print(anova(m_op_abs, type = 3))

cat("\n[EMMs] |ΔOpinion|: pattern × task_type × timeF\n")
emm_op_abs <- emmeans(m_op_abs, ~ pattern * task_type * timeF)
print(emm_op_abs)

## Post-hoc 1: pattern
cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
emm_op_abs_pattern <- emmeans(
  m_op_abs,
  specs = "pattern",
  by    = c("task_type", "timeF"),
  data  = dat_op_delta
)
print(emm_op_abs_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_abs_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(
  emm_op_abs_pattern,
  sigma = sigma(m_op_abs),
  edf   = df.residual(m_op_abs)
))

## Post-hoc 2: task_type
cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
emm_op_abs_task <- emmeans(
  m_op_abs,
  specs = "task_type",
  by    = c("pattern", "timeF"),
  data  = dat_op_delta
)
print(emm_op_abs_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_abs_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(
  emm_op_abs_task,
  sigma = sigma(m_op_abs),
  edf   = df.residual(m_op_abs)
))

## Post-hoc 3: timeF
cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
emm_op_abs_time <- emmeans(
  m_op_abs,
  specs = "timeF",
  by    = c("pattern", "task_type"),
  data  = dat_op_delta
)
print(emm_op_abs_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_abs_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(
  emm_op_abs_time,
  sigma = sigma(m_op_abs),
  edf   = df.residual(m_op_abs)
))










###########################################################
# 8-2) OPINION: 그룹 = "always_neg"
#      (위와 동일한 코드 흐름, op_group만 변경)
###########################################################

op_group <- "always_neg"

ids_op <- traj_opinion %>%
  dplyr::filter(opinion_traj_group == op_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_op$participant_id)), "\n")

dat_op_raw <- resp_behavior %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_op_delta <- resp_delta %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")

cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")
m_op_raw <- lmer(
  form_op_raw,
  data = dat_op_raw,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_op_raw))
cat("\n[Type-III ANOVA]\n"); print(anova(m_op_raw, type = 3))
cat("\n[EMMs] Opinion raw: pattern × task_type × timeF\n")
emm_op_raw <- emmeans(m_op_raw, ~ pattern * task_type * timeF); print(emm_op_raw)

cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
emm_op_raw_pattern <- emmeans(m_op_raw, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_op_raw)
print(emm_op_raw_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_raw_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_op_raw_pattern, sigma = sigma(m_op_raw),
                        edf = df.residual(m_op_raw)))

cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
emm_op_raw_task <- emmeans(m_op_raw, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_op_raw)
print(emm_op_raw_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_raw_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_op_raw_task, sigma = sigma(m_op_raw),
                        edf = df.residual(m_op_raw)))

cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
emm_op_raw_time <- emmeans(m_op_raw, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_op_raw)
print(emm_op_raw_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_raw_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_op_raw_time, sigma = sigma(m_op_raw),
                        edf = df.residual(m_op_raw)))

cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")
m_op_signed <- lmer(
  form_op_signed,
  data = dat_op_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_op_signed))
cat("\n[Type-III ANOVA]\n"); print(anova(m_op_signed, type = 3))
cat("\n[EMMs] ΔOpinion signed: pattern × task_type × timeF\n")
emm_op_signed <- emmeans(m_op_signed, ~ pattern * task_type * timeF); print(emm_op_signed)

cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
emm_op_signed_pattern <- emmeans(m_op_signed, specs = "pattern",
                                 by = c("task_type", "timeF"), data = dat_op_delta)
print(emm_op_signed_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_signed_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_op_signed_pattern, sigma = sigma(m_op_signed),
                        edf = df.residual(m_op_signed)))

cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
emm_op_signed_task <- emmeans(m_op_signed, specs = "task_type",
                              by = c("pattern", "timeF"), data = dat_op_delta)
print(emm_op_signed_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_signed_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_op_signed_task, sigma = sigma(m_op_signed),
                        edf = df.residual(m_op_signed)))

cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
emm_op_signed_time <- emmeans(m_op_signed, specs = "timeF",
                              by = c("pattern", "task_type"), data = dat_op_delta)
print(emm_op_signed_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_signed_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_op_signed_time, sigma = sigma(m_op_signed),
                        edf = df.residual(m_op_signed)))

cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")
m_op_abs <- lmer(
  form_op_abs,
  data = dat_op_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_op_abs))
cat("\n[Type-III ANOVA]\n"); print(anova(m_op_abs, type = 3))
cat("\n[EMMs] |ΔOpinion|: pattern × task_type × timeF\n")
emm_op_abs <- emmeans(m_op_abs, ~ pattern * task_type * timeF); print(emm_op_abs)

cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
emm_op_abs_pattern <- emmeans(m_op_abs, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_op_delta)
print(emm_op_abs_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_abs_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_op_abs_pattern, sigma = sigma(m_op_abs),
                        edf = df.residual(m_op_abs)))

cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
emm_op_abs_task <- emmeans(m_op_abs, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_op_delta)
print(emm_op_abs_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_abs_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_op_abs_task, sigma = sigma(m_op_abs),
                        edf = df.residual(m_op_abs)))

cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
emm_op_abs_time <- emmeans(m_op_abs, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_op_delta)
print(emm_op_abs_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_abs_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_op_abs_time, sigma = sigma(m_op_abs),
                        edf = df.residual(m_op_abs)))










###########################################################
# 8-3) OPINION: 그룹 = "mixed"
###########################################################

op_group <- "mixed"

ids_op <- traj_opinion %>%
  dplyr::filter(opinion_traj_group == op_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_op$participant_id)), "\n")

dat_op_raw <- resp_behavior %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_op_delta <- resp_delta %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")

cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")
m_op_raw <- lmer(
  form_op_raw,
  data = dat_op_raw,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_op_raw))
cat("\n[Type-III ANOVA]\n"); print(anova(m_op_raw, type = 3))
cat("\n[EMMs] Opinion raw: pattern × task_type × timeF\n")
emm_op_raw <- emmeans(m_op_raw, ~ pattern * task_type * timeF); print(emm_op_raw)

cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
emm_op_raw_pattern <- emmeans(m_op_raw, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_op_raw)
print(emm_op_raw_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_raw_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_op_raw_pattern, sigma = sigma(m_op_raw),
                        edf = df.residual(m_op_raw)))

cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
emm_op_raw_task <- emmeans(m_op_raw, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_op_raw)
print(emm_op_raw_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_raw_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_op_raw_task, sigma = sigma(m_op_raw),
                        edf = df.residual(m_op_raw)))

cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
emm_op_raw_time <- emmeans(m_op_raw, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_op_raw)
print(emm_op_raw_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_raw_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_op_raw_time, sigma = sigma(m_op_raw),
                        edf = df.residual(m_op_raw)))

cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")
m_op_signed <- lmer(
  form_op_signed,
  data = dat_op_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_op_signed))
cat("\n[Type-III ANOVA]\n"); print(anova(m_op_signed, type = 3))
cat("\n[EMMs] ΔOpinion signed: pattern × task_type × timeF\n")
emm_op_signed <- emmeans(m_op_signed, ~ pattern * task_type * timeF); print(emm_op_signed)

cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
emm_op_signed_pattern <- emmeans(m_op_signed, specs = "pattern",
                                 by = c("task_type", "timeF"), data = dat_op_delta)
print(emm_op_signed_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_signed_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_op_signed_pattern, sigma = sigma(m_op_signed),
                        edf = df.residual(m_op_signed)))

cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
emm_op_signed_task <- emmeans(m_op_signed, specs = "task_type",
                              by = c("pattern", "timeF"), data = dat_op_delta)
print(emm_op_signed_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_signed_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_op_signed_task, sigma = sigma(m_op_signed),
                        edf = df.residual(m_op_signed)))

cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
emm_op_signed_time <- emmeans(m_op_signed, specs = "timeF",
                              by = c("pattern", "task_type"), data = dat_op_delta)
print(emm_op_signed_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_signed_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_op_signed_time, sigma = sigma(m_op_signed),
                        edf = df.residual(m_op_signed)))

cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")
m_op_abs <- lmer(
  form_op_abs,
  data = dat_op_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_op_abs))
cat("\n[Type-III ANOVA]\n"); print(anova(m_op_abs, type = 3))
cat("\n[EMMs] |ΔOpinion|: pattern × task_type × timeF\n")
emm_op_abs <- emmeans(m_op_abs, ~ pattern * task_type * timeF); print(emm_op_abs)

cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
emm_op_abs_pattern <- emmeans(m_op_abs, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_op_delta)
print(emm_op_abs_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_op_abs_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_op_abs_pattern, sigma = sigma(m_op_abs),
                        edf = df.residual(m_op_abs)))

cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
emm_op_abs_task <- emmeans(m_op_abs, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_op_delta)
print(emm_op_abs_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_op_abs_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_op_abs_task, sigma = sigma(m_op_abs),
                        edf = df.residual(m_op_abs)))

cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
emm_op_abs_time <- emmeans(m_op_abs, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_op_delta)
print(emm_op_abs_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_op_abs_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_op_abs_time, sigma = sigma(m_op_abs),
                        edf = df.residual(m_op_abs)))









###########################################################
# 9) CONFIDENCE trajectory 그룹별 분석
#    - 각 그룹("always_pos","always_neg","mixed")
###########################################################


###########################################################
# 9-1) CONFIDENCE: 그룹 = "always_pos"
###########################################################

cf_group <- "always_pos"

ids_cf <- traj_conf %>%
  dplyr::filter(conf_traj_group == cf_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_cf$participant_id)), "\n")

dat_cf_raw <- resp_behavior %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_cf_delta <- resp_delta %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[CONFIDENCE] dat_cf_raw 행수 :", nrow(dat_cf_raw), "\n")
cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")

cat("\n[CONFIDENCE] --- 모델 1: Confidence raw (T0~T4) ---\n")
m_cf_raw <- lmer(
  form_cf_raw,
  data = dat_cf_raw,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_raw))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_raw, type = 3))
cat("\n[EMMs] Confidence raw: pattern × task_type × timeF\n")
emm_cf_raw <- emmeans(m_cf_raw, ~ pattern * task_type * timeF); print(emm_cf_raw)

cat("\n[Post-hoc] Confidence raw － pattern (by task_type × timeF)\n")
emm_cf_raw_pattern <- emmeans(m_cf_raw, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_cf_raw)
print(emm_cf_raw_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_raw_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_raw_pattern, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[Post-hoc] Confidence raw － task_type (by pattern × timeF)\n")
emm_cf_raw_task <- emmeans(m_cf_raw, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_cf_raw)
print(emm_cf_raw_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_raw_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_raw_task, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[Post-hoc] Confidence raw － timeF (by pattern × task_type)\n")
emm_cf_raw_time <- emmeans(m_cf_raw, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_cf_raw)
print(emm_cf_raw_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_raw_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_raw_time, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[CONFIDENCE] --- 모델 2: ΔConfidence (signed, Tk−T0) ---\n")
m_cf_signed <- lmer(
  form_cf_signed,
  data = dat_cf_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_signed))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_signed, type = 3))
cat("\n[EMMs] ΔConfidence signed: pattern × task_type × timeF\n")
emm_cf_signed <- emmeans(m_cf_signed, ~ pattern * task_type * timeF); print(emm_cf_signed)

cat("\n[Post-hoc] ΔConfidence signed － pattern (by task_type × timeF)\n")
emm_cf_signed_pattern <- emmeans(m_cf_signed, specs = "pattern",
                                 by = c("task_type", "timeF"), data = dat_cf_delta)
print(emm_cf_signed_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_signed_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_signed_pattern, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[Post-hoc] ΔConfidence signed － task_type (by pattern × timeF)\n")
emm_cf_signed_task <- emmeans(m_cf_signed, specs = "task_type",
                              by = c("pattern", "timeF"), data = dat_cf_delta)
print(emm_cf_signed_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_signed_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_signed_task, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[Post-hoc] ΔConfidence signed － timeF (by pattern × task_type)\n")
emm_cf_signed_time <- emmeans(m_cf_signed, specs = "timeF",
                              by = c("pattern", "task_type"), data = dat_cf_delta)
print(emm_cf_signed_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_signed_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_signed_time, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[CONFIDENCE] --- 모델 3: |ΔConfidence| (absolute delta) ---\n")
m_cf_abs <- lmer(
  form_cf_abs,
  data = dat_cf_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_abs))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_abs, type = 3))
cat("\n[EMMs] |ΔConfidence|: pattern × task_type × timeF\n")
emm_cf_abs <- emmeans(m_cf_abs, ~ pattern * task_type * timeF); print(emm_cf_abs)

cat("\n[Post-hoc] |ΔConfidence| － pattern (by task_type × timeF)\n")
emm_cf_abs_pattern <- emmeans(m_cf_abs, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_cf_delta)
print(emm_cf_abs_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_abs_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_abs_pattern, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))

cat("\n[Post-hoc] |ΔConfidence| － task_type (by pattern × timeF)\n")
emm_cf_abs_task <- emmeans(m_cf_abs, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_cf_delta)
print(emm_cf_abs_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_abs_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_abs_task, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))

cat("\n[Post-hoc] |ΔConfidence| － timeF (by pattern × task_type)\n")
emm_cf_abs_time <- emmeans(m_cf_abs, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_cf_delta)
print(emm_cf_abs_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_abs_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_abs_time, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))









###########################################################
# 9-2) CONFIDENCE: 그룹 = "always_neg"
###########################################################

cf_group <- "always_neg"

ids_cf <- traj_conf %>%
  dplyr::filter(conf_traj_group == cf_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_cf$participant_id)), "\n")

dat_cf_raw <- resp_behavior %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_cf_delta <- resp_delta %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[CONFIDENCE] dat_cf_raw 행수 :", nrow(dat_cf_raw), "\n")
cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")

cat("\n[CONFIDENCE] --- 모델 1: Confidence raw (T0~T4) ---\n")
m_cf_raw <- lmer(
  form_cf_raw,
  data = dat_cf_raw,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_raw))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_raw, type = 3))
cat("\n[EMMs] Confidence raw: pattern × task_type × timeF\n")
emm_cf_raw <- emmeans(m_cf_raw, ~ pattern * task_type * timeF); print(emm_cf_raw)

cat("\n[Post-hoc] Confidence raw － pattern (by task_type × timeF)\n")
emm_cf_raw_pattern <- emmeans(m_cf_raw, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_cf_raw)
print(emm_cf_raw_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_raw_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_raw_pattern, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[Post-hoc] Confidence raw － task_type (by pattern × timeF)\n")
emm_cf_raw_task <- emmeans(m_cf_raw, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_cf_raw)
print(emm_cf_raw_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_raw_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_raw_task, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[Post-hoc] Confidence raw － timeF (by pattern × task_type)\n")
emm_cf_raw_time <- emmeans(m_cf_raw, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_cf_raw)
print(emm_cf_raw_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_raw_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_raw_time, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[CONFIDENCE] --- 모델 2: ΔConfidence (signed, Tk−T0) ---\n")
m_cf_signed <- lmer(
  form_cf_signed,
  data = dat_cf_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_signed))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_signed, type = 3))
cat("\n[EMMs] ΔConfidence signed: pattern × task_type × timeF\n")
emm_cf_signed <- emmeans(m_cf_signed, ~ pattern * task_type * timeF); print(emm_cf_signed)

cat("\n[Post-hoc] ΔConfidence signed － pattern (by task_type × timeF)\n")
emm_cf_signed_pattern <- emmeans(m_cf_signed, specs = "pattern",
                                 by = c("task_type", "timeF"), data = dat_cf_delta)
print(emm_cf_signed_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_signed_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_signed_pattern, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[Post-hoc] ΔConfidence signed － task_type (by pattern × timeF)\n")
emm_cf_signed_task <- emmeans(m_cf_signed, specs = "task_type",
                              by = c("pattern", "timeF"), data = dat_cf_delta)
print(emm_cf_signed_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_signed_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_signed_task, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[Post-hoc] ΔConfidence signed － timeF (by pattern × task_type)\n")
emm_cf_signed_time <- emmeans(m_cf_signed, specs = "timeF",
                              by = c("pattern", "task_type"), data = dat_cf_delta)
print(emm_cf_signed_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_signed_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_signed_time, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[CONFIDENCE] --- 모델 3: |ΔConfidence| (absolute delta) ---\n")
m_cf_abs <- lmer(
  form_cf_abs,
  data = dat_cf_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_abs))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_abs, type = 3))
cat("\n[EMMs] |ΔConfidence|: pattern × task_type × timeF\n")
emm_cf_abs <- emmeans(m_cf_abs, ~ pattern * task_type * timeF); print(emm_cf_abs)

cat("\n[Post-hoc] |ΔConfidence| － pattern (by task_type × timeF)\n")
emm_cf_abs_pattern <- emmeans(m_cf_abs, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_cf_delta)
print(emm_cf_abs_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_abs_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_abs_pattern, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))

cat("\n[Post-hoc] |ΔConfidence| － task_type (by pattern × timeF)\n")
emm_cf_abs_task <- emmeans(m_cf_abs, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_cf_delta)
print(emm_cf_abs_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_abs_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_abs_task, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))

cat("\n[Post-hoc] |ΔConfidence| － timeF (by pattern × task_type)\n")
emm_cf_abs_time <- emmeans(m_cf_abs, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_cf_delta)
print(emm_cf_abs_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_abs_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_abs_time, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))









###########################################################
# 9-3) CONFIDENCE: 그룹 = "mixed"
###########################################################

cf_group <- "mixed"

ids_cf <- traj_conf %>%
  dplyr::filter(conf_traj_group == cf_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_cf$participant_id)), "\n")

dat_cf_raw <- resp_behavior %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_cf_delta <- resp_delta %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[CONFIDENCE] dat_cf_raw 행수 :", nrow(dat_cf_raw), "\n")
cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")

cat("\n[CONFIDENCE] --- 모델 1: Confidence raw (T0~T4) ---\n")
m_cf_raw <- lmer(
  form_cf_raw,
  data = dat_cf_raw,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_raw))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_raw, type = 3))
cat("\n[EMMs] Confidence raw: pattern × task_type × timeF\n")
emm_cf_raw <- emmeans(m_cf_raw, ~ pattern * task_type * timeF); print(emm_cf_raw)

cat("\n[Post-hoc] Confidence raw － pattern (by task_type × timeF)\n")
emm_cf_raw_pattern <- emmeans(m_cf_raw, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_cf_raw)
print(emm_cf_raw_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_raw_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_raw_pattern, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[Post-hoc] Confidence raw － task_type (by pattern × timeF)\n")
emm_cf_raw_task <- emmeans(m_cf_raw, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_cf_raw)
print(emm_cf_raw_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_raw_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_raw_task, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[Post-hoc] Confidence raw － timeF (by pattern × task_type)\n")
emm_cf_raw_time <- emmeans(m_cf_raw, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_cf_raw)
print(emm_cf_raw_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_raw_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_raw_time, sigma = sigma(m_cf_raw),
                        edf = df.residual(m_cf_raw)))

cat("\n[CONFIDENCE] --- 모델 2: ΔConfidence (signed, Tk−T0) ---\n")
m_cf_signed <- lmer(
  form_cf_signed,
  data = dat_cf_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_signed))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_signed, type = 3))
cat("\n[EMMs] ΔConfidence signed: pattern × task_type × timeF\n")
emm_cf_signed <- emmeans(m_cf_signed, ~ pattern * task_type * timeF); print(emm_cf_signed)

cat("\n[Post-hoc] ΔConfidence signed － pattern (by task_type × timeF)\n")
emm_cf_signed_pattern <- emmeans(m_cf_signed, specs = "pattern",
                                 by = c("task_type", "timeF"), data = dat_cf_delta)
print(emm_cf_signed_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_signed_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_signed_pattern, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[Post-hoc] ΔConfidence signed － task_type (by pattern × timeF)\n")
emm_cf_signed_task <- emmeans(m_cf_signed, specs = "task_type",
                              by = c("pattern", "timeF"), data = dat_cf_delta)
print(emm_cf_signed_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_signed_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_signed_task, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[Post-hoc] ΔConfidence signed － timeF (by pattern × task_type)\n")
emm_cf_signed_time <- emmeans(m_cf_signed, specs = "timeF",
                              by = c("pattern", "task_type"), data = dat_cf_delta)
print(emm_cf_signed_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_signed_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_signed_time, sigma = sigma(m_cf_signed),
                        edf = df.residual(m_cf_signed)))

cat("\n[CONFIDENCE] --- 모델 3: |ΔConfidence| (absolute delta) ---\n")
m_cf_abs <- lmer(
  form_cf_abs,
  data = dat_cf_delta,
  REML = FALSE
)
cat("\n[요약]\n"); print(summary(m_cf_abs))
cat("\n[Type-III ANOVA]\n"); print(anova(m_cf_abs, type = 3))
cat("\n[EMMs] |ΔConfidence|: pattern × task_type × timeF\n")
emm_cf_abs <- emmeans(m_cf_abs, ~ pattern * task_type * timeF); print(emm_cf_abs)

cat("\n[Post-hoc] |ΔConfidence| － pattern (by task_type × timeF)\n")
emm_cf_abs_pattern <- emmeans(m_cf_abs, specs = "pattern",
                              by = c("task_type", "timeF"), data = dat_cf_delta)
print(emm_cf_abs_pattern)
cat("\nBonferroni (pattern within task_type × timeF):\n")
print(pairs(emm_cf_abs_pattern, adjust = "bonferroni"))
cat("\nEffect sizes (pattern within task_type × timeF):\n")
print(emmeans::eff_size(emm_cf_abs_pattern, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))

cat("\n[Post-hoc] |ΔConfidence| － task_type (by pattern × timeF)\n")
emm_cf_abs_task <- emmeans(m_cf_abs, specs = "task_type",
                           by = c("pattern", "timeF"), data = dat_cf_delta)
print(emm_cf_abs_task)
cat("\nBonferroni (task_type within pattern × timeF):\n")
print(pairs(emm_cf_abs_task, adjust = "bonferroni"))
cat("\nEffect sizes (task_type within pattern × timeF):\n")
print(emmeans::eff_size(emm_cf_abs_task, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))

cat("\n[Post-hoc] |ΔConfidence| － timeF (by pattern × task_type)\n")
emm_cf_abs_time <- emmeans(m_cf_abs, specs = "timeF",
                           by = c("pattern", "task_type"), data = dat_cf_delta)
print(emm_cf_abs_time)
cat("\nBonferroni (timeF within pattern × task_type):\n")
print(pairs(emm_cf_abs_time, adjust = "bonferroni"))
cat("\nEffect sizes (timeF within pattern × task_type):\n")
print(emmeans::eff_size(emm_cf_abs_time, sigma = sigma(m_cf_abs),
                        edf = df.residual(m_cf_abs)))





cat("\n\n모든 trajectory 기반 분석이 완료되었습니다.\n")



# ## 제대로 실행이 안될 시 dataAnalysis_SIMA_FINAL.R 먼저 실행.
# 
# ## ---------------------------------------------------------
# ## 9) Trajectory 분류 (opinion / confidence)
# ##    - 단위: participant_id × task_type
# ##    - 그룹: always_pos / always_neg / mixed
# ## ---------------------------------------------------------
# 
# library(dplyr)
# library(purrr)
# library(emmeans)
# 
# cat("\n\n==============================\n")
# cat("9) Trajectory-based grouping (opinion / confidence)\n")
# cat("==============================\n")
# 
# # 9-1) 분류 함수
# #  - x: T0~T4 값 벡터
# #  - diff(x)를 계산한 후:
# #     * 0은 무시
# #     * non-zero diff가 모두 >0  → always_pos
# #     * non-zero diff가 모두 <0  → always_neg
# #     * 양/음 섞이거나, 모든 diff가 0 → mixed
# classify_traj <- function(x) {
#   x <- x[!is.na(x)]
#   if (length(x) < 2L) return(NA_character_)
#   d <- diff(x)
#   d_nz <- d[d != 0]
#   if (length(d_nz) == 0L) {
#     # 모든 diff가 0 인 경우 → mixed로 분류
#     return("mixed")
#   }
#   if (all(d_nz > 0)) return("always_pos")
#   if (all(d_nz < 0)) return("always_neg")
#   return("mixed")
# }
# 
# time_levels <- paste0("T", 0:4)
# 
# ## 9-2) Opinion trajectory 분류 (resp_behavior$opinion 사용)
# ##  - opinion은 이미 sign-flip된 값이라고 가정
# traj_opinion <- resp_behavior %>%
#   filter(timeF %in% time_levels) %>%
#   mutate(timeF = factor(timeF, levels = time_levels)) %>%
#   arrange(participant_id, task_type, timeF) %>%
#   group_by(participant_id, task_type) %>%
#   summarise(
#     pattern     = dplyr::first(pattern),
#     n_t_unique  = n_distinct(timeF),
#     opinion_seq = list(opinion),
#     .groups     = "drop"
#   )
# 
# # 시점이 5개가 아닌 경우 경고 출력
# traj_opinion_problem <- traj_opinion %>% filter(n_t_unique != 5L)
# if (nrow(traj_opinion_problem) > 0) {
#   cat("\n[경고] Opinion: T0~T4 5시점이 아닌 (participant_id, task_type) 조합이 있습니다:\n")
#   print(traj_opinion_problem %>% select(participant_id, task_type, pattern, n_t_unique))
# }
# 
# traj_opinion <- traj_opinion %>%
#   filter(n_t_unique == 5L) %>%
#   mutate(
#     opinion_traj_group = map_chr(opinion_seq, classify_traj),
#     opinion_traj_group = factor(opinion_traj_group,
#                                 levels = c("always_pos", "always_neg", "mixed"))
#   )
# 
# ## 9-3) Confidence trajectory 분류 (resp_behavior$confidence 사용)
# traj_conf <- resp_behavior %>%
#   filter(timeF %in% time_levels) %>%
#   mutate(timeF = factor(timeF, levels = time_levels)) %>%
#   arrange(participant_id, task_type, timeF) %>%
#   group_by(participant_id, task_type) %>%
#   summarise(
#     pattern    = dplyr::first(pattern),
#     n_t_unique = n_distinct(timeF),
#     conf_seq   = list(confidence),
#     .groups    = "drop"
#   )
# 
# traj_conf_problem <- traj_conf %>% filter(n_t_unique != 5L)
# if (nrow(traj_conf_problem) > 0) {
#   cat("\n[경고] Confidence: T0~T4 5시점이 아닌 (participant_id, task_type) 조합이 있습니다:\n")
#   print(traj_conf_problem %>% select(participant_id, task_type, pattern, n_t_unique))
# }
# 
# traj_conf <- traj_conf %>%
#   filter(n_t_unique == 5L) %>%
#   mutate(
#     conf_traj_group = map_chr(conf_seq, classify_traj),
#     conf_traj_group = factor(conf_traj_group,
#                              levels = c("always_pos", "always_neg", "mixed"))
#   )
# 
# 
# ## 9-4) Trajectory 그룹 요약 (count + percentage 출력)
# ##  - 주의: traj_*의 한 행 = (participant_id, task_type) 1개의 trajectory
# ##          → 전체 합계는 "trajectory 개수" (참가자 수 × task 수)
# ##  - task_type별 요약을 보면, 각 task 안에서의 참가자 수로 해석 가능
# 
# cat("\n=== Opinion trajectory groups (단위: participant × task_type, 즉 trajectory 수) ===\n")
# 
# opinion_group_overall <- traj_opinion %>%
#   count(opinion_traj_group, name = "n_traj") %>%      # 그룹별 trajectory 개수
#   mutate(
#     prop = n_traj / sum(n_traj),                      # 비율 (0~1)
#     pct  = round(100 * prop, 1)                       # 퍼센트
#   )
# print(opinion_group_overall)
# 
# cat("\n--- Opinion trajectory groups by task_type (각 task 안에서 참가자 비율) ---\n")
# opinion_group_by_task <- traj_opinion %>%
#   count(task_type, opinion_traj_group, name = "n_participants") %>%  # task 내 참가자 수
#   group_by(task_type) %>%
#   mutate(
#     prop = n_participants / sum(n_participants),
#     pct  = round(100 * prop, 1)
#   ) %>%
#   ungroup()
# print(opinion_group_by_task)
# 
# ## === Opinion: task_type × opinion_traj_group × participant_no 목록 만들기 ===
# 
# opinion_ids_long <- traj_opinion %>%
#   # participant_id -> participant_no 매핑 붙이기
#   dplyr::left_join(
#     resp_behavior %>%
#       dplyr::distinct(participant_id, participant_no),
#     by = "participant_id"
#   ) %>%
#   # 각 참가자를 한 줄씩 (task_type × opinion_traj_group × participant_no)
#   dplyr::distinct(task_type, opinion_traj_group, participant_no) %>%
#   dplyr::arrange(task_type, opinion_traj_group, participant_no)
# 
# ## 경로 지정 (원하는 대로 바꾸세요)
# # out_dir <- "./Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data/"  # MAC
# our_dir <- "C:\\Users\\USER\\Desktop\\github\\mutliAgentExperiment_dataAnalysis\\refined_data\\"   #  Window
# 
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 
# opinion_csv_path <- file.path(out_dir, "FINAL_3Trajectory_Opinion-ParticipantNo.csv")
# 
# write.csv(
#   opinion_ids_long,
#   opinion_csv_path,
#   row.names   = FALSE,
#   fileEncoding = "UTF-8"
# )
# 
# cat("\n[파일 저장 완료] Opinion ID 목록 저장 위치:\n", opinion_csv_path, "\n")
# 
# ## --- 추가 코드 끝 ---
# 
# cat("\n--- Opinion trajectory groups: pattern × task_type (trajectory 기준) ---\n")
# opinion_group_by_pattern_task <- traj_opinion %>%
#   count(pattern, task_type, opinion_traj_group, name = "n_traj") %>%
#   group_by(pattern, task_type) %>%
#   mutate(
#     prop = n_traj / sum(n_traj),
#     pct  = round(100 * prop, 1)
#   ) %>%
#   ungroup()
# print(opinion_group_by_pattern_task)
# 
# 
# cat("\n=== Confidence trajectory groups (단위: participant × task_type, 즉 trajectory 수) ===\n")
# 
# conf_group_overall <- traj_conf %>%
#   count(conf_traj_group, name = "n_traj") %>%
#   mutate(
#     prop = n_traj / sum(n_traj),
#     pct  = round(100 * prop, 1)
#   )
# print(conf_group_overall)
# 
# cat("\n--- Confidence trajectory groups by task_type (각 task 안에서 참가자 비율) ---\n")
# conf_group_by_task <- traj_conf %>%
#   count(task_type, conf_traj_group, name = "n_participants") %>%
#   group_by(task_type) %>%
#   mutate(
#     prop = n_participants / sum(n_participants),
#     pct  = round(100 * prop, 1)
#   ) %>%
#   ungroup()
# print(conf_group_by_task)
# 
# ## === Confidence: task_type × conf_traj_group × participant_no 목록 만들기 ===
# conf_ids_long <- traj_conf %>%
#   dplyr::left_join(
#     resp_behavior %>%
#       dplyr::distinct(participant_id, participant_no),
#     by = "participant_id"
#   ) %>%
#   dplyr::distinct(task_type, conf_traj_group, participant_no) %>%
#   dplyr::arrange(task_type, conf_traj_group, participant_no)
# 
# ## 위에서 쓴 out_dir 그대로 사용하거나, 새로 지정해도 됩니다.
# conf_csv_path <- file.path(out_dir, "FINAL_3Trajectory_Confidence_ParticipantNo.csv")
# 
# write.csv(
#   conf_ids_long,
#   conf_csv_path,
#   row.names   = FALSE,
#   fileEncoding = "UTF-8"
# )
# 
# cat("\n[파일 저장 완료] Confidence ID 목록 저장 위치:\n", conf_csv_path, "\n")
# ## --- 추가 코드 끝 ---
# 
# cat("\n--- Confidence trajectory groups: pattern × task_type (trajectory 기준) ---\n")
# conf_group_by_pattern_task <- traj_conf %>%
#   count(pattern, task_type, conf_traj_group, name = "n_traj") %>%
#   group_by(pattern, task_type) %>%
#   mutate(
#     prop = n_traj / sum(n_traj),
#     pct  = round(100 * prop, 1)
#   ) %>%
#   ungroup()
# print(conf_group_by_pattern_task)
# 
# 
# 
# 
# 
# 
# 
# ## =========================================================
# ## OPINION trajectory 그룹별 분석 (함수/반복문 없이)
# ## ---------------------------------------------------------
# ## 전제:
# ##  - traj_opinion : 9번에서 만든 trajectory 요약 (participant_id, task_type, opinion_traj_group 등)
# ##  - resp_behavior: 원래 T0~T4 반복측정 데이터 (opinion, confidence 포함)
# ##  - resp_delta   : Tk−T0 delta 데이터
# ##  - fit_lmm_step_general(), form_op_* 들이 이미 정의되어 있다고 가정
# ## =========================================================
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Opinion trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ##   - 아래 줄에서 원하는 그룹 이름으로 바꾼 뒤,
# ##     STEP 1 ~ STEP 5 를 순서대로 실행
# ## ---------------------------------------------------------
# op_group <- "always_pos"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 다시 실행 가능
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_op <- traj_opinion %>%
#   dplyr::filter(opinion_traj_group == op_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_op$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_op_raw   : T0~T4 원자료 (opinion값)
# ##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta 등)
# ## ---------------------------------------------------------
# 
# dat_op_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_op_delta <- resp_delta %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
# cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
# cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Opinion raw (T0~T4) 모델 적합 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_raw  : LMM 모델
# ##       * emm_op_raw_pattern, emm_op_raw_task, emm_op_raw_time : emmeans 결과
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")
# 
# m_op_raw <- fit_lmm_step_general(
#   form_op_raw_full,    # 풀 모형
#   form_op_raw_uncorr,  # 상호작용 축소 모형 (있다면)
#   form_op_raw_int,     # 최소 모형 (있다면)
#   data = dat_op_raw,
#   reml = FALSE
# )
# 
# summary(m_op_raw)
# anova(m_op_raw, type = 3)
# 
# ## Post-hoc 1: 패턴(pattern) 비교 (task_type × timeF 별)
# cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
# 
# emm_op_raw_pattern <- emmeans::emmeans(
#   m_op_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_raw     # 안전하게 data 명시
# )
# pairs(emm_op_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_pattern,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 2: task_type 비교 (pattern × timeF 별)
# cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
# 
# emm_op_raw_task <- emmeans::emmeans(
#   m_op_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_task,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 3: timeF 비교 (pattern × task_type 별)
# cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
# 
# emm_op_raw_time <- emmeans::emmeans(
#   m_op_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_time,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔOpinion (signed, Tk−T0) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_signed
# ##       * emm_op_signed_pattern, emm_op_signed_task, emm_op_signed_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")
# 
# m_op_signed <- fit_lmm_step_general(
#   form_op_signed_full,
#   form_op_signed_uncorr,
#   form_op_signed_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_signed)
# anova(m_op_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
# 
# emm_op_signed_pattern <- emmeans::emmeans(
#   m_op_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_pattern,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
# 
# emm_op_signed_task <- emmeans::emmeans(
#   m_op_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_task,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
# 
# emm_op_signed_time <- emmeans::emmeans(
#   m_op_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_time,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 5] |ΔOpinion| (absolute delta) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_abs
# ##       * emm_op_abs_pattern, emm_op_abs_task, emm_op_abs_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")
# 
# m_op_abs <- fit_lmm_step_general(
#   form_op_abs_full,
#   form_op_abs_uncorr,
#   form_op_abs_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_abs)
# anova(m_op_abs, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
# 
# emm_op_abs_pattern <- emmeans::emmeans(
#   m_op_abs,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_pattern,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
# 
# emm_op_abs_task <- emmeans::emmeans(
#   m_op_abs,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_task,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
# 
# emm_op_abs_time <- emmeans::emmeans(
#   m_op_abs,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_time,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# 
# 
# 
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Opinion trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ##   - 아래 줄에서 원하는 그룹 이름으로 바꾼 뒤,
# ##     STEP 1 ~ STEP 5 를 순서대로 실행
# ## ---------------------------------------------------------
# op_group <- "always_neg"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 다시 실행 가능
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_op <- traj_opinion %>%
#   dplyr::filter(opinion_traj_group == op_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_op$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_op_raw   : T0~T4 원자료 (opinion값)
# ##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta 등)
# ## ---------------------------------------------------------
# 
# dat_op_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_op_delta <- resp_delta %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
# cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
# cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Opinion raw (T0~T4) 모델 적합 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_raw  : LMM 모델
# ##       * emm_op_raw_pattern, emm_op_raw_task, emm_op_raw_time : emmeans 결과
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")
# 
# m_op_raw <- fit_lmm_step_general(
#   form_op_raw_full,    # 풀 모형
#   form_op_raw_uncorr,  # 상호작용 축소 모형 (있다면)
#   form_op_raw_int,     # 최소 모형 (있다면)
#   data = dat_op_raw,
#   reml = FALSE
# )
# 
# summary(m_op_raw)
# anova(m_op_raw, type = 3)
# 
# ## Post-hoc 1: 패턴(pattern) 비교 (task_type × timeF 별)
# cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
# 
# emm_op_raw_pattern <- emmeans::emmeans(
#   m_op_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_raw     # 안전하게 data 명시
# )
# pairs(emm_op_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_pattern,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 2: task_type 비교 (pattern × timeF 별)
# cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
# 
# emm_op_raw_task <- emmeans::emmeans(
#   m_op_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_task,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 3: timeF 비교 (pattern × task_type 별)
# cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
# 
# emm_op_raw_time <- emmeans::emmeans(
#   m_op_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_time,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔOpinion (signed, Tk−T0) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_signed
# ##       * emm_op_signed_pattern, emm_op_signed_task, emm_op_signed_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")
# 
# m_op_signed <- fit_lmm_step_general(
#   form_op_signed_full,
#   form_op_signed_uncorr,
#   form_op_signed_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_signed)
# anova(m_op_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
# 
# emm_op_signed_pattern <- emmeans::emmeans(
#   m_op_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_pattern,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
# 
# emm_op_signed_task <- emmeans::emmeans(
#   m_op_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_task,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
# 
# emm_op_signed_time <- emmeans::emmeans(
#   m_op_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_time,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 5] |ΔOpinion| (absolute delta) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_abs
# ##       * emm_op_abs_pattern, emm_op_abs_task, emm_op_abs_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")
# 
# m_op_abs <- fit_lmm_step_general(
#   form_op_abs_full,
#   form_op_abs_uncorr,
#   form_op_abs_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_abs)
# anova(m_op_abs, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
# 
# emm_op_abs_pattern <- emmeans::emmeans(
#   m_op_abs,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_pattern,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
# 
# emm_op_abs_task <- emmeans::emmeans(
#   m_op_abs,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_task,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
# 
# emm_op_abs_time <- emmeans::emmeans(
#   m_op_abs,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_time,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# 
# 
# 
# 
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Opinion trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ##   - 아래 줄에서 원하는 그룹 이름으로 바꾼 뒤,
# ##     STEP 1 ~ STEP 5 를 순서대로 실행
# ## ---------------------------------------------------------
# op_group <- "always_neg"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 다시 실행 가능
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_op <- traj_opinion %>%
#   dplyr::filter(opinion_traj_group == op_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_op$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_op_raw   : T0~T4 원자료 (opinion값)
# ##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta 등)
# ## ---------------------------------------------------------
# 
# dat_op_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_op_delta <- resp_delta %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
# cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
# cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Opinion raw (T0~T4) 모델 적합 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_raw  : LMM 모델
# ##       * emm_op_raw_pattern, emm_op_raw_task, emm_op_raw_time : emmeans 결과
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")
# 
# m_op_raw <- fit_lmm_step_general(
#   form_op_raw_full,    # 풀 모형
#   form_op_raw_uncorr,  # 상호작용 축소 모형 (있다면)
#   form_op_raw_int,     # 최소 모형 (있다면)
#   data = dat_op_raw,
#   reml = FALSE
# )
# 
# summary(m_op_raw)
# anova(m_op_raw, type = 3)
# 
# ## Post-hoc 1: 패턴(pattern) 비교 (task_type × timeF 별)
# cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
# 
# emm_op_raw_pattern <- emmeans::emmeans(
#   m_op_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_raw     # 안전하게 data 명시
# )
# pairs(emm_op_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_pattern,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 2: task_type 비교 (pattern × timeF 별)
# cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
# 
# emm_op_raw_task <- emmeans::emmeans(
#   m_op_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_task,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 3: timeF 비교 (pattern × task_type 별)
# cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
# 
# emm_op_raw_time <- emmeans::emmeans(
#   m_op_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_time,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔOpinion (signed, Tk−T0) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_signed
# ##       * emm_op_signed_pattern, emm_op_signed_task, emm_op_signed_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")
# 
# m_op_signed <- fit_lmm_step_general(
#   form_op_signed_full,
#   form_op_signed_uncorr,
#   form_op_signed_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_signed)
# anova(m_op_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
# 
# emm_op_signed_pattern <- emmeans::emmeans(
#   m_op_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_pattern,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
# 
# emm_op_signed_task <- emmeans::emmeans(
#   m_op_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_task,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
# 
# emm_op_signed_time <- emmeans::emmeans(
#   m_op_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_time,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 5] |ΔOpinion| (absolute delta) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_abs
# ##       * emm_op_abs_pattern, emm_op_abs_task, emm_op_abs_time
# ## ---------------------------------------------------------
# # 
# # cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")
# # 
# # m_op_abs <- fit_lmm_step_general(
# #   form_op_abs_full,
# #   form_op_abs_uncorr,
# #   form_op_abs_int,
# #   data = dat_op_delta,
# #   reml = FALSE
# # )
# # 
# # summary(m_op_abs)
# # anova(m_op_abs, type = 3)
# # 
# # ## Post-hoc 1: 패턴
# # cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
# # 
# # emm_op_abs_pattern <- emmeans::emmeans(
# #   m_op_abs,
# #   specs = "pattern",
# #   by    = c("task_type", "timeF"),
# #   data  = dat_op_delta
# # )
# # pairs(emm_op_abs_pattern, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_op_abs_pattern,
# #   sigma = sigma(m_op_abs),
# #   edf   = df.residual(m_op_abs)
# # )
# # 
# # ## Post-hoc 2: task_type
# # cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
# # 
# # emm_op_abs_task <- emmeans::emmeans(
# #   m_op_abs,
# #   specs = "task_type",
# #   by    = c("pattern", "timeF"),
# #   data  = dat_op_delta
# # )
# # pairs(emm_op_abs_task, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_op_abs_task,
# #   sigma = sigma(m_op_abs),
# #   edf   = df.residual(m_op_abs)
# # )
# # 
# # ## Post-hoc 3: timeF
# # cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
# # 
# # emm_op_abs_time <- emmeans::emmeans(
# #   m_op_abs,
# #   specs = "timeF",
# #   by    = c("pattern", "task_type"),
# #   data  = dat_op_delta
# # )
# # pairs(emm_op_abs_time, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_op_abs_time,
# #   sigma = sigma(m_op_abs),
# #   edf   = df.residual(m_op_abs)
# # )
# 
# 
# 
# 
# 
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Opinion trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ##   - 아래 줄에서 원하는 그룹 이름으로 바꾼 뒤,
# ##     STEP 1 ~ STEP 5 를 순서대로 실행
# ## ---------------------------------------------------------
# op_group <- "mixed"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 다시 실행 가능
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_op <- traj_opinion %>%
#   dplyr::filter(opinion_traj_group == op_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_op$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_op_raw   : T0~T4 원자료 (opinion값)
# ##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta 등)
# ## ---------------------------------------------------------
# 
# dat_op_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_op_delta <- resp_delta %>%
#   dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[OPINION] dat_op_raw 행수 :", nrow(dat_op_raw), "\n")
# cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
# cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Opinion raw (T0~T4) 모델 적합 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_raw  : LMM 모델
# ##       * emm_op_raw_pattern, emm_op_raw_task, emm_op_raw_time : emmeans 결과
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 1: Opinion raw (T0~T4) ---\n")
# 
# m_op_raw <- fit_lmm_step_general(
#   form_op_raw_full,    # 풀 모형
#   form_op_raw_uncorr,  # 상호작용 축소 모형 (있다면)
#   form_op_raw_int,     # 최소 모형 (있다면)
#   data = dat_op_raw,
#   reml = FALSE
# )
# 
# summary(m_op_raw)
# anova(m_op_raw, type = 3)
# 
# ## Post-hoc 1: 패턴(pattern) 비교 (task_type × timeF 별)
# cat("\n[Post-hoc] Opinion raw － pattern (by task_type × timeF)\n")
# 
# emm_op_raw_pattern <- emmeans::emmeans(
#   m_op_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_raw     # 안전하게 data 명시
# )
# pairs(emm_op_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_pattern,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 2: task_type 비교 (pattern × timeF 별)
# cat("\n[Post-hoc] Opinion raw － task_type (by pattern × timeF)\n")
# 
# emm_op_raw_task <- emmeans::emmeans(
#   m_op_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_task,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# ## Post-hoc 3: timeF 비교 (pattern × task_type 별)
# cat("\n[Post-hoc] Opinion raw － timeF (by pattern × task_type)\n")
# 
# emm_op_raw_time <- emmeans::emmeans(
#   m_op_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_raw
# )
# pairs(emm_op_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_raw_time,
#   sigma = sigma(m_op_raw),
#   edf   = df.residual(m_op_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔOpinion (signed, Tk−T0) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_signed
# ##       * emm_op_signed_pattern, emm_op_signed_task, emm_op_signed_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 2: ΔOpinion (signed, Tk−T0) ---\n")
# 
# m_op_signed <- fit_lmm_step_general(
#   form_op_signed_full,
#   form_op_signed_uncorr,
#   form_op_signed_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_signed)
# anova(m_op_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔOpinion signed － pattern (by task_type × timeF)\n")
# 
# emm_op_signed_pattern <- emmeans::emmeans(
#   m_op_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_pattern,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔOpinion signed － task_type (by pattern × timeF)\n")
# 
# emm_op_signed_task <- emmeans::emmeans(
#   m_op_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_task,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔOpinion signed － timeF (by pattern × task_type)\n")
# 
# emm_op_signed_time <- emmeans::emmeans(
#   m_op_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_signed_time,
#   sigma = sigma(m_op_signed),
#   edf   = df.residual(m_op_signed)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 5] |ΔOpinion| (absolute delta) 모델 & Post-hoc
# ##   - 결과 객체:
# ##       * m_op_abs
# ##       * emm_op_abs_pattern, emm_op_abs_task, emm_op_abs_time
# ## ---------------------------------------------------------
# 
# cat("\n[OPINION] --- 모델 3: |ΔOpinion| (absolute delta) ---\n")
# 
# m_op_abs <- fit_lmm_step_general(
#   form_op_abs_full,
#   form_op_abs_uncorr,
#   form_op_abs_int,
#   data = dat_op_delta,
#   reml = FALSE
# )
# 
# summary(m_op_abs)
# anova(m_op_abs, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] |ΔOpinion| － pattern (by task_type × timeF)\n")
# 
# emm_op_abs_pattern <- emmeans::emmeans(
#   m_op_abs,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_pattern,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] |ΔOpinion| － task_type (by pattern × timeF)\n")
# 
# emm_op_abs_task <- emmeans::emmeans(
#   m_op_abs,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_task,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] |ΔOpinion| － timeF (by pattern × task_type)\n")
# 
# emm_op_abs_time <- emmeans::emmeans(
#   m_op_abs,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_op_delta
# )
# pairs(emm_op_abs_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_op_abs_time,
#   sigma = sigma(m_op_abs),
#   edf   = df.residual(m_op_abs)
# )
# 
# 
# 
# 
# 
# 
# 
# 
# ## =========================================================
# ## CONFIDENCE trajectory 그룹별 분석 (함수/반복문 없이)
# ## ---------------------------------------------------------
# ## 전제:
# ##  - traj_conf   : 9번에서 만든 trajectory 요약 (conf_traj_group 포함)
# ##  - resp_behavior, resp_delta 는 동일
# ##  - fit_lmm_step_general(), form_cf_* 들이 이미 정의되어 있다고 가정
# ## =========================================================
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Confidence trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ## ---------------------------------------------------------
# cf_group <- "always_pos"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 사용
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_cf <- traj_conf %>%
#   dplyr::filter(conf_traj_group == cf_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_cf$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_cf_raw   : T0~T4 원자료 (confidence값)
# ##   - dat_cf_delta : Tk−T0 delta 자료 (confidence_delta 등)
# ## ---------------------------------------------------------
# 
# dat_cf_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_cf_delta <- resp_delta %>%
#   dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[CONFIDENCE] dat_cf_raw 행수 :", nrow(dat_cf_raw), "\n")
# cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
# cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Confidence raw (T0~T4) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 1: Confidence raw (T0~T4) ---\n")
# 
# m_cf_raw <- fit_lmm_step_general(
#   form_cf_raw_full,
#   form_cf_raw_uncorr,
#   form_cf_raw_int,
#   data = dat_cf_raw,
#   reml = FALSE
# )
# 
# summary(m_cf_raw)
# anova(m_cf_raw, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] Confidence raw － pattern (by task_type × timeF)\n")
# 
# emm_cf_raw_pattern <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_pattern,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] Confidence raw － task_type (by pattern × timeF)\n")
# 
# emm_cf_raw_task <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_task,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] Confidence raw － timeF (by pattern × task_type)\n")
# 
# emm_cf_raw_time <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_time,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔConfidence (signed, Tk−T0) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 2: ΔConfidence (signed, Tk−T0) ---\n")
# 
# m_cf_signed <- fit_lmm_step_general(
#   form_cf_signed_full,
#   form_cf_signed_uncorr,
#   form_cf_signed_int,
#   data = dat_cf_delta,
#   reml = FALSE
# )
# 
# summary(m_cf_signed)
# anova(m_cf_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔConfidence signed － pattern (by task_type × timeF)\n")
# 
# emm_cf_signed_pattern <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_pattern,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔConfidence signed － task_type (by pattern × timeF)\n")
# 
# emm_cf_signed_task <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_task,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔConfidence signed － timeF (by pattern × task_type)\n")
# 
# emm_cf_signed_time <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_time,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# 
# # ## ---------------------------------------------------------
# # ## [STEP 5] |ΔConfidence| (absolute delta) 모델 & Post-hoc
# # ## ---------------------------------------------------------
# # 
# # cat("\n[CONFIDENCE] --- 모델 3: |ΔConfidence| (absolute delta) ---\n")
# # 
# # m_cf_abs <- fit_lmm_step_general(
# #   form_cf_abs_full,
# #   form_cf_abs_uncorr,
# #   form_cf_abs_int,
# #   data = dat_cf_delta,
# #   reml = FALSE
# # )
# # 
# # summary(m_cf_abs)
# # anova(m_cf_abs, type = 3)
# # 
# # ## Post-hoc 1: 패턴
# # cat("\n[Post-hoc] |ΔConfidence| － pattern (by task_type × timeF)\n")
# # 
# # emm_cf_abs_pattern <- emmeans::emmeans(
# #   m_cf_abs,
# #   specs = "pattern",
# #   by    = c("task_type", "timeF"),
# #   data  = dat_cf_delta
# # )
# # pairs(emm_cf_abs_pattern, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_cf_abs_pattern,
# #   sigma = sigma(m_cf_abs),
# #   edf   = df.residual(m_cf_abs)
# # )
# # 
# # ## Post-hoc 2: task_type
# # cat("\n[Post-hoc] |ΔConfidence| － task_type (by pattern × timeF)\n")
# # 
# # emm_cf_abs_task <- emmeans::emmeans(
# #   m_cf_abs,
# #   specs = "task_type",
# #   by    = c("pattern", "timeF"),
# #   data  = dat_cf_delta
# # )
# # pairs(emm_cf_abs_task, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_cf_abs_task,
# #   sigma = sigma(m_cf_abs),
# #   edf   = df.residual(m_cf_abs)
# # )
# # 
# # ## Post-hoc 3: timeF
# # cat("\n[Post-hoc] |ΔConfidence| － timeF (by pattern × task_type)\n")
# # 
# # emm_cf_abs_time <- emmeans::emmeans(
# #   m_cf_abs,
# #   specs = "timeF",
# #   by    = c("pattern", "task_type"),
# #   data  = dat_cf_delta
# # )
# # pairs(emm_cf_abs_time, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_cf_abs_time,
# #   sigma = sigma(m_cf_abs),
# #   edf   = df.residual(m_cf_abs)
# # )
# 
# 
# 
# 
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Confidence trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ## ---------------------------------------------------------
# cf_group <- "always_neg"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 사용
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_cf <- traj_conf %>%
#   dplyr::filter(conf_traj_group == cf_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_cf$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_cf_raw   : T0~T4 원자료 (confidence값)
# ##   - dat_cf_delta : Tk−T0 delta 자료 (confidence_delta 등)
# ## ---------------------------------------------------------
# 
# dat_cf_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_cf_delta <- resp_delta %>%
#   dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[CONFIDENCE] dat_cf_raw 행수 :", nrow(dat_cf_raw), "\n")
# cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
# cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Confidence raw (T0~T4) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 1: Confidence raw (T0~T4) ---\n")
# 
# m_cf_raw <- fit_lmm_step_general(
#   form_cf_raw_full,
#   form_cf_raw_uncorr,
#   form_cf_raw_int,
#   data = dat_cf_raw,
#   reml = FALSE
# )
# 
# summary(m_cf_raw)
# anova(m_cf_raw, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] Confidence raw － pattern (by task_type × timeF)\n")
# 
# emm_cf_raw_pattern <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_pattern,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] Confidence raw － task_type (by pattern × timeF)\n")
# 
# emm_cf_raw_task <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_task,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] Confidence raw － timeF (by pattern × task_type)\n")
# 
# emm_cf_raw_time <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_time,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔConfidence (signed, Tk−T0) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 2: ΔConfidence (signed, Tk−T0) ---\n")
# 
# m_cf_signed <- fit_lmm_step_general(
#   form_cf_signed_full,
#   form_cf_signed_uncorr,
#   form_cf_signed_int,
#   data = dat_cf_delta,
#   reml = FALSE
# )
# 
# summary(m_cf_signed)
# anova(m_cf_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔConfidence signed － pattern (by task_type × timeF)\n")
# 
# emm_cf_signed_pattern <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_pattern,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔConfidence signed － task_type (by pattern × timeF)\n")
# 
# emm_cf_signed_task <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_task,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔConfidence signed － timeF (by pattern × task_type)\n")
# 
# emm_cf_signed_time <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_time,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# # 
# # ## ---------------------------------------------------------
# # ## [STEP 5] |ΔConfidence| (absolute delta) 모델 & Post-hoc
# # ## ---------------------------------------------------------
# # 
# # cat("\n[CONFIDENCE] --- 모델 3: |ΔConfidence| (absolute delta) ---\n")
# # 
# # m_cf_abs <- fit_lmm_step_general(
# #   form_cf_abs_full,
# #   form_cf_abs_uncorr,
# #   form_cf_abs_int,
# #   data = dat_cf_delta,
# #   reml = FALSE
# # )
# # 
# # summary(m_cf_abs)
# # anova(m_cf_abs, type = 3)
# # 
# # ## Post-hoc 1: 패턴
# # cat("\n[Post-hoc] |ΔConfidence| － pattern (by task_type × timeF)\n")
# # 
# # emm_cf_abs_pattern <- emmeans::emmeans(
# #   m_cf_abs,
# #   specs = "pattern",
# #   by    = c("task_type", "timeF"),
# #   data  = dat_cf_delta
# # )
# # pairs(emm_cf_abs_pattern, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_cf_abs_pattern,
# #   sigma = sigma(m_cf_abs),
# #   edf   = df.residual(m_cf_abs)
# # )
# # 
# # ## Post-hoc 2: task_type
# # cat("\n[Post-hoc] |ΔConfidence| － task_type (by pattern × timeF)\n")
# # 
# # emm_cf_abs_task <- emmeans::emmeans(
# #   m_cf_abs,
# #   specs = "task_type",
# #   by    = c("pattern", "timeF"),
# #   data  = dat_cf_delta
# # )
# # pairs(emm_cf_abs_task, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_cf_abs_task,
# #   sigma = sigma(m_cf_abs),
# #   edf   = df.residual(m_cf_abs)
# # )
# # 
# # ## Post-hoc 3: timeF
# # cat("\n[Post-hoc] |ΔConfidence| － timeF (by pattern × task_type)\n")
# # 
# # emm_cf_abs_time <- emmeans::emmeans(
# #   m_cf_abs,
# #   specs = "timeF",
# #   by    = c("pattern", "task_type"),
# #   data  = dat_cf_delta
# # )
# # pairs(emm_cf_abs_time, adjust = "bonferroni")
# # emmeans::eff_size(
# #   emm_cf_abs_time,
# #   sigma = sigma(m_cf_abs),
# #   edf   = df.residual(m_cf_abs)
# # )
# # 
# # 
# 
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 0] 분석할 Confidence trajectory 그룹 선택
# ##   - 가능한 값: "always_pos", "always_neg", "mixed"
# ## ---------------------------------------------------------
# cf_group <- "mixed"   # <- 여기만 "always_pos", "always_neg" 또는 "mixed" 로 바꿔서 사용
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 1] 선택한 그룹에 속하는 (participant_id, task_type) 추출
# ## ---------------------------------------------------------
# 
# ids_cf <- traj_conf %>%
#   dplyr::filter(conf_traj_group == cf_group) %>%
#   dplyr::select(participant_id, task_type)
# 
# cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
# cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
# cat(" - 포함된 고유 참가자 수:",
#     length(unique(ids_cf$participant_id)), "\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 2] 이 그룹에 해당하는 데이터 만들기
# ##   - dat_cf_raw   : T0~T4 원자료 (confidence값)
# ##   - dat_cf_delta : Tk−T0 delta 자료 (confidence_delta 등)
# ## ---------------------------------------------------------
# 
# dat_cf_raw <- resp_behavior %>%
#   dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# dat_cf_delta <- resp_delta %>%
#   dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
#   droplevels()
# 
# cat("\n[CONFIDENCE] dat_cf_raw 행수 :", nrow(dat_cf_raw), "\n")
# cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
# cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 3] Confidence raw (T0~T4) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 1: Confidence raw (T0~T4) ---\n")
# 
# m_cf_raw <- fit_lmm_step_general(
#   form_cf_raw_full,
#   form_cf_raw_uncorr,
#   form_cf_raw_int,
#   data = dat_cf_raw,
#   reml = FALSE
# )
# 
# summary(m_cf_raw)
# anova(m_cf_raw, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] Confidence raw － pattern (by task_type × timeF)\n")
# 
# emm_cf_raw_pattern <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_pattern,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] Confidence raw － task_type (by pattern × timeF)\n")
# 
# emm_cf_raw_task <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_task,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] Confidence raw － timeF (by pattern × task_type)\n")
# 
# emm_cf_raw_time <- emmeans::emmeans(
#   m_cf_raw,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_raw
# )
# pairs(emm_cf_raw_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_raw_time,
#   sigma = sigma(m_cf_raw),
#   edf   = df.residual(m_cf_raw)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 4] ΔConfidence (signed, Tk−T0) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 2: ΔConfidence (signed, Tk−T0) ---\n")
# 
# m_cf_signed <- fit_lmm_step_general(
#   form_cf_signed_full,
#   form_cf_signed_uncorr,
#   form_cf_signed_int,
#   data = dat_cf_delta,
#   reml = FALSE
# )
# 
# summary(m_cf_signed)
# anova(m_cf_signed, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] ΔConfidence signed － pattern (by task_type × timeF)\n")
# 
# emm_cf_signed_pattern <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_pattern,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] ΔConfidence signed － task_type (by pattern × timeF)\n")
# 
# emm_cf_signed_task <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_task,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] ΔConfidence signed － timeF (by pattern × task_type)\n")
# 
# emm_cf_signed_time <- emmeans::emmeans(
#   m_cf_signed,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_signed_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_signed_time,
#   sigma = sigma(m_cf_signed),
#   edf   = df.residual(m_cf_signed)
# )
# 
# 
# ## ---------------------------------------------------------
# ## [STEP 5] |ΔConfidence| (absolute delta) 모델 & Post-hoc
# ## ---------------------------------------------------------
# 
# cat("\n[CONFIDENCE] --- 모델 3: |ΔConfidence| (absolute delta) ---\n")
# 
# m_cf_abs <- fit_lmm_step_general(
#   form_cf_abs_full,
#   form_cf_abs_uncorr,
#   form_cf_abs_int,
#   data = dat_cf_delta,
#   reml = FALSE
# )
# 
# summary(m_cf_abs)
# anova(m_cf_abs, type = 3)
# 
# ## Post-hoc 1: 패턴
# cat("\n[Post-hoc] |ΔConfidence| － pattern (by task_type × timeF)\n")
# 
# emm_cf_abs_pattern <- emmeans::emmeans(
#   m_cf_abs,
#   specs = "pattern",
#   by    = c("task_type", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_abs_pattern, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_abs_pattern,
#   sigma = sigma(m_cf_abs),
#   edf   = df.residual(m_cf_abs)
# )
# 
# ## Post-hoc 2: task_type
# cat("\n[Post-hoc] |ΔConfidence| － task_type (by pattern × timeF)\n")
# 
# emm_cf_abs_task <- emmeans::emmeans(
#   m_cf_abs,
#   specs = "task_type",
#   by    = c("pattern", "timeF"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_abs_task, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_abs_task,
#   sigma = sigma(m_cf_abs),
#   edf   = df.residual(m_cf_abs)
# )
# 
# ## Post-hoc 3: timeF
# cat("\n[Post-hoc] |ΔConfidence| － timeF (by pattern × task_type)\n")
# 
# emm_cf_abs_time <- emmeans::emmeans(
#   m_cf_abs,
#   specs = "timeF",
#   by    = c("pattern", "task_type"),
#   data  = dat_cf_delta
# )
# pairs(emm_cf_abs_time, adjust = "bonferroni")
# emmeans::eff_size(
#   emm_cf_abs_time,
#   sigma = sigma(m_cf_abs),
#   edf   = df.residual(m_cf_abs)
# )
# 
