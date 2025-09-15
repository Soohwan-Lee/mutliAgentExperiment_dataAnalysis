# 패키지
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)

# 경로
base_dir <- "./Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data/"

# 데이터 로드
bg <- readr::read_csv(file.path(base_dir, "FINAL_background_surveys_rows_refined.csv"))
resp <- readr::read_csv(file.path(base_dir, "FINAL_responses_rows_refined_experiment_only.csv"))

# 타입 통일
bg$participant_id <- as.character(bg$participant_id)
resp$participant_id <- as.character(resp$participant_id)

# 공변량: 문항 평균 → z
bg_cov <- bg %>%
  mutate(
    SII_raw = rowMeans(dplyr::select(., dplyr::matches("^sii_\d+")), na.rm = TRUE),
    NFC_raw = rowMeans(dplyr::select(., dplyr::matches("^nfc_\\d+")), na.rm = TRUE),
    AIacc_raw = rowMeans(dplyr::select(., dplyr::matches("^ai_acceptance_\d+$")), na.rm = TRUE),
    SII_z = as.numeric(scale(SII_raw)),
    NFC_z = as.numeric(scale(NFC_raw)),
    AIacc_z = as.numeric(scale(AIacc_raw))
  ) %>%
  select(participant_id, SII_z, NFC_z, AIacc_z)

# 패턴·과업 지정 + 원시 시간 인덱스 선택
resp1 <- resp %>%
  mutate(
    pattern = dplyr::recode(condition_type,
                            "majority" = "Majority",
                            "minority" = "Minority",
                            "minorityDiffusion" = "Diffusion"),
    task_type = dplyr::case_when(
      session_key %in% c("normative","Normative") ~ "Normative",
      session_key %in% c("informative","Informative") ~ "Informative",
      TRUE ~ NA_character_
    ),
    time_raw = dplyr::case_when(
      task_type == "Normative" ~ as.numeric(normative_task_index),
      task_type == "Informative" ~ as.numeric(informative_task_index),
      TRUE ~ NA_real_
    ),
    rt_log = log1p(rt_ms)
  ) %>%
  filter(!is.na(task_type), !is.na(pattern), !is.na(time_raw))

# 여기서 핵심 수정: 그룹 내 최소값을 0으로 평행이동하여 t_ord(0..4) 생성
resp1 <- resp1 %>%
  group_by(participant_id, task_type) %>%
  mutate(t_ord = as.integer(time_raw - min(time_raw, na.rm = TRUE))) %>%
  ungroup()

# T0/T1만 뽑아 wide로 피벗하여 delta 계산
resp_d1 <- resp1 %>%
  filter(t_ord %in% c(0, 1)) %>%
  select(participant_id, pattern, task_type, t_ord, opinion, confidence, rt_log) %>%
  tidyr::pivot_wider(
    names_from = t_ord,
    values_from = c(opinion, confidence, rt_log),
    names_prefix = "t"
  ) %>%
  
  
#컬럼 이름 정리
rename(
  opinion_T0 = opinion_t0, opinion_T1 = opinion_t1,
  conf_T0 = confidence_t0, conf_T1 = confidence_t1,
  rt_T0 = rt_log_t0, rt_T1 = rt_log_t1
) %>%
  

# T0/T1 둘 다 있는 케이스만
filter(!is.na(opinion_T0), !is.na(opinion_T1)) %>%
  mutate(
    delta_opinion_T1 = opinion_T1 - opinion_T0,
    delta_conf_T1 = conf_T1 - conf_T0,
    delta_rt_log_T1 = rt_T1 - rt_T0
  ) %>%
  left_join(bg_cov, by = "participant_id") %>%
  mutate(
    pattern = factor(pattern, levels = c("Majority","Minority","Diffusion")),
    task_type = factor(task_type, levels = c("Normative","Informative"))
  )

# sanity check
cat("resp_d1 rows:", nrow(resp_d1), "\n") # 기대: 참가자수*2(과업)
print(table(resp_d1$task_type, resp_d1$pattern, useNA="ifany"))
print(colSums(is.na(resp_d1[, c("delta_opinion_T1","delta_conf_T1","delta_rt_log_T1",
                                "SII_z","NFC_z","AIacc_z")])))

# 당신이 이미 갖고 있는 LMM 코드 그대로 실행
m_delta_op <- lmer(
  delta_opinion_T1 ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = resp_d1, REML = TRUE
)
m_delta_conf <- lmer(
  delta_conf_T1 ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = resp_d1, REML = TRUE
)
m_delta_rt <- lmer(
  delta_rt_log_T1 ~ pattern * task_type + SII_z + NFC_z + AIacc_z + (1 | participant_id),
  data = resp_d1, REML = TRUE
)

summary(m_delta_op); anova(m_delta_op, type=3)
summary(m_delta_conf); anova(m_delta_conf, type=3)
summary(m_delta_rt); anova(m_delta_rt, type=3)

# 사후검정 + 효과크기
emm_op <- emmeans(m_delta_op, ~ pattern | task_type)
pairs(emm_op, adjust="bonferroni")
eff_size(emm_op, sigma = sigma(m_delta_op), edf = df.residual(m_delta_op))

emm_cf <- emmeans(m_delta_conf, ~ pattern | task_type)
pairs(emm_cf, adjust="bonferroni")
eff_size(emm_cf, sigma = sigma(m_delta_conf), edf = df.residual(m_delta_conf))

emm_rt <- emmeans(m_delta_rt, ~ pattern | task_type)
pairs(emm_rt, adjust="bonferroni")
eff_size(emm_rt, sigma = sigma(m_delta_rt), edf = df.residual(m_delta_rt))