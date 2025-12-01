### dataVisualizaiton_SIMA_FINAL.R을 먼저 실행하고 모든 함수가 정의되어 있다고 가정.

## ---------------------------------------------------------
## [OPINION] 시각화하고 싶은 trajectory 그룹 선택
##   - 가능한 값: "always_pos", "always_neg", "mixed"
##   - 아래 op_group 값만 바꿔서,
##     밑에 코드 블럭들을 위에서부터 순서대로 실행하면 됨.
## ---------------------------------------------------------

#### Always Postiive - opinion

op_group <- "always_pos"   # 여기만 "always_neg" 또는 "mixed" 로 바꾸면 됨


## ---------------------------------------------------------
## [OPINION] 선택한 trajectory 그룹에 속하는 (participant_id, task_type) 추출
## ---------------------------------------------------------

ids_op <- traj_opinion %>%
  dplyr::filter(opinion_traj_group == op_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_op$participant_id)), "\n")


## ---------------------------------------------------------
## [OPINION] 이 그룹에 해당하는 raw / delta 데이터 만들기
##   - dat_op_raw   : T0~T4 원자료 (opinion)
##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta, opinion_delta_abs)
## ---------------------------------------------------------

dat_op_raw <- resp_behavior %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_op_delta <- resp_delta %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[OPINION] dat_op_raw 행수   :", nrow(dat_op_raw), "\n")
cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")


## ---------------------------------------------------------
## [OPINION] RAW (T0~T4) 플롯: pattern × time × task_type
##   - dv = "opinion"
##   - pattern-within-time (검은 꺾쇠)
##   - time-within-pattern (색깔 꺾쇠)
##   - 아래 plot_behavior_pretty() 는
##     이미 정의해 둔 최종 시각화 함수를 재사용하는 것.
## ---------------------------------------------------------

fig_opinion_raw_optraj <- plot_behavior_pretty(
  data = dat_op_raw,
  dv   = "opinion",
  model = NULL,                     # 여기서 새로 LMM 적합해서 꺾쇠 계산
  p_adjust = p_adj_method,          # 위에서 설정한 Bonferroni
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # T0 vs T1~T4
) +
#  ggtitle(paste0("Opinion (Raw Value) – Trajectory group: ", op_group)) +
  labs(
    y = "Opinion (Raw Value, -50 to 50)"   # <- 여기에 원하는 라벨 이름 넣기
  )  +
  theme(
    plot.title = element_text(
      hjust = 0.5,       # 가운데 정렬
    )
  )

print(fig_opinion_raw_optraj)


## ---------------------------------------------------------
## [선택] RAW (T0~T4) + 참가자별 꺾은선 (조금 복잡하지만 직관적인 그림)
## ---------------------------------------------------------

fig_opinion_raw_optraj_indiv <- plot_behavior_pretty(
  data = dat_op_raw,
  dv   = "opinion",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)

fig_opinion_raw_optraj_indiv <- add_indiv_lines_to_behavior(
  p           = fig_opinion_raw_optraj_indiv,
  data        = dat_op_raw,
  dv          = "opinion",
  dodge_width = 0.43   # 박스플롯 dodge와 비슷하게
  # sample_n  = 50     # 원하면 일부 참가자만 샘플링
)

fig_opinion_raw_optraj_indiv <- fig_opinion_raw_optraj_indiv +
  labs(
    x = "Time",
    y = "Opinion"
  ) +
  labs(
    y = "Opinion (Raw Value, -50 to 50)"   # <- 여기에 원하는 라벨 이름 넣기
  )  
  #ggtitle(paste0("Opinion (raw + individual lines) – Trajectory group: ", op_group))

print(fig_opinion_raw_optraj_indiv)




## ---------------------------------------------------------
## [OPINION] ΔOpinion (signed, Tk−T0) 플롯
##   - dv = "opinion_delta"
##   - resp_delta에서 T1~T4 데이터만 들어있기 때문에
##     time_contrast = "baseline" 은 T1 vs T2~T4 라는 의미가 됨
## ---------------------------------------------------------

fig_opinion_signed_optraj <- plot_behavior_pretty(
  data = dat_op_delta,
  dv   = "opinion_delta",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # baseline = T1
)  +
  labs(
    y = "Opinion Delta (Tk-T0)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
  #ggtitle(paste0("ΔOpinion (Tk − T0) – Trajectory group: ", op_group))

print(fig_opinion_signed_optraj)


# ## ---------------------------------------------------------
# ## [OPINION] |ΔOpinion| (absolute delta) 플롯
# ##   - dv = "opinion_delta_abs"
# ## ---------------------------------------------------------
# 
# fig_opinion_abs_optraj <- plot_behavior_pretty(
#   data = dat_op_delta,
#   dv   = "opinion_delta_abs",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) 
#   #ggtitle(paste0("|ΔOpinion| – Trajectory group: ", op_group))
# 
# print(fig_opinion_abs_optraj)







#### Always Negative - opinion

op_group <- "always_neg"   # 여기만 "always_neg" 또는 "mixed" 로 바꾸면 됨


## ---------------------------------------------------------
## [OPINION] 선택한 trajectory 그룹에 속하는 (participant_id, task_type) 추출
## ---------------------------------------------------------

ids_op <- traj_opinion %>%
  dplyr::filter(opinion_traj_group == op_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_op$participant_id)), "\n")


## ---------------------------------------------------------
## [OPINION] 이 그룹에 해당하는 raw / delta 데이터 만들기
##   - dat_op_raw   : T0~T4 원자료 (opinion)
##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta, opinion_delta_abs)
## ---------------------------------------------------------

dat_op_raw <- resp_behavior %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_op_delta <- resp_delta %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[OPINION] dat_op_raw 행수   :", nrow(dat_op_raw), "\n")
cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")


## ---------------------------------------------------------
## [OPINION] RAW (T0~T4) 플롯: pattern × time × task_type
##   - dv = "opinion"
##   - pattern-within-time (검은 꺾쇠)
##   - time-within-pattern (색깔 꺾쇠)
##   - 아래 plot_behavior_pretty() 는
##     이미 정의해 둔 최종 시각화 함수를 재사용하는 것.
## ---------------------------------------------------------

fig_opinion_raw_optraj <- plot_behavior_pretty(
  data = dat_op_raw,
  dv   = "opinion",
  model = NULL,                     # 여기서 새로 LMM 적합해서 꺾쇠 계산
  p_adjust = p_adj_method,          # 위에서 설정한 Bonferroni
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # T0 vs T1~T4
) +
  #  ggtitle(paste0("Opinion (Raw Value) – Trajectory group: ", op_group)) +
  labs(
    y = "Opinion (Raw Value, -50 to 50)"   # <- 여기에 원하는 라벨 이름 넣기
  )  +
  theme(
    plot.title = element_text(
      hjust = 0.5,       # 가운데 정렬
    )
  )

print(fig_opinion_raw_optraj)


## ---------------------------------------------------------
## [선택] RAW (T0~T4) + 참가자별 꺾은선 (조금 복잡하지만 직관적인 그림)
## ---------------------------------------------------------

fig_opinion_raw_optraj_indiv <- plot_behavior_pretty(
  data = dat_op_raw,
  dv   = "opinion",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)

fig_opinion_raw_optraj_indiv <- add_indiv_lines_to_behavior(
  p           = fig_opinion_raw_optraj_indiv,
  data        = dat_op_raw,
  dv          = "opinion",
  dodge_width = 0.43   # 박스플롯 dodge와 비슷하게
  # sample_n  = 50     # 원하면 일부 참가자만 샘플링
)

fig_opinion_raw_optraj_indiv <- fig_opinion_raw_optraj_indiv +
  labs(
    x = "Time",
    y = "Opinion"
  ) +
  labs(
    y = "Opinion (Raw Value, -50 to 50)"   # <- 여기에 원하는 라벨 이름 넣기
  )  
#ggtitle(paste0("Opinion (raw + individual lines) – Trajectory group: ", op_group))

print(fig_opinion_raw_optraj_indiv)




## ---------------------------------------------------------
## [OPINION] ΔOpinion (signed, Tk−T0) 플롯
##   - dv = "opinion_delta"
##   - resp_delta에서 T1~T4 데이터만 들어있기 때문에
##     time_contrast = "baseline" 은 T1 vs T2~T4 라는 의미가 됨
## ---------------------------------------------------------

fig_opinion_signed_optraj <- plot_behavior_pretty(
  data = dat_op_delta,
  dv   = "opinion_delta",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # baseline = T1
) +
  labs(
    y = "Opinion Delta (Tk-T0)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
#ggtitle(paste0("ΔOpinion (Tk − T0) – Trajectory group: ", op_group))

print(fig_opinion_signed_optraj)


# ## ---------------------------------------------------------
# ## [OPINION] |ΔOpinion| (absolute delta) 플롯
# ##   - dv = "opinion_delta_abs"
# ## ---------------------------------------------------------
# 
# fig_opinion_abs_optraj <- plot_behavior_pretty(
#   data = dat_op_delta,
#   dv   = "opinion_delta_abs",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) 
# #ggtitle(paste0("|ΔOpinion| – Trajectory group: ", op_group))
# 
# print(fig_opinion_abs_optraj)





#### Mixed - opinion

op_group <- "mixed"   # 여기만 "always_neg" 또는 "mixed" 로 바꾸면 됨


## ---------------------------------------------------------
## [OPINION] 선택한 trajectory 그룹에 속하는 (participant_id, task_type) 추출
## ---------------------------------------------------------

ids_op <- traj_opinion %>%
  dplyr::filter(opinion_traj_group == op_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[OPINION] 선택한 그룹:", op_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_op), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_op$participant_id)), "\n")


## ---------------------------------------------------------
## [OPINION] 이 그룹에 해당하는 raw / delta 데이터 만들기
##   - dat_op_raw   : T0~T4 원자료 (opinion)
##   - dat_op_delta : Tk−T0 delta 자료 (opinion_delta, opinion_delta_abs)
## ---------------------------------------------------------

dat_op_raw <- resp_behavior %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_op_delta <- resp_delta %>%
  dplyr::semi_join(ids_op, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[OPINION] dat_op_raw 행수   :", nrow(dat_op_raw), "\n")
cat("          dat_op_delta 행수 :", nrow(dat_op_delta), "\n")
cat("          (참가자 수:", length(unique(dat_op_raw$participant_id)), ")\n")


## ---------------------------------------------------------
## [OPINION] RAW (T0~T4) 플롯: pattern × time × task_type
##   - dv = "opinion"
##   - pattern-within-time (검은 꺾쇠)
##   - time-within-pattern (색깔 꺾쇠)
##   - 아래 plot_behavior_pretty() 는
##     이미 정의해 둔 최종 시각화 함수를 재사용하는 것.
## ---------------------------------------------------------

fig_opinion_raw_optraj <- plot_behavior_pretty(
  data = dat_op_raw,
  dv   = "opinion",
  model = NULL,                     # 여기서 새로 LMM 적합해서 꺾쇠 계산
  p_adjust = p_adj_method,          # 위에서 설정한 Bonferroni
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # T0 vs T1~T4
) +
  #  ggtitle(paste0("Opinion (Raw Value) – Trajectory group: ", op_group)) +
  labs(
    y = "Opinion (Raw Value, -50 to 50)"   # <- 여기에 원하는 라벨 이름 넣기
  )  +
  theme(
    plot.title = element_text(
      hjust = 0.5,       # 가운데 정렬
    )
  )

print(fig_opinion_raw_optraj)


## ---------------------------------------------------------
## [선택] RAW (T0~T4) + 참가자별 꺾은선 (조금 복잡하지만 직관적인 그림)
## ---------------------------------------------------------

fig_opinion_raw_optraj_indiv <- plot_behavior_pretty(
  data = dat_op_raw,
  dv   = "opinion",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)

fig_opinion_raw_optraj_indiv <- add_indiv_lines_to_behavior(
  p           = fig_opinion_raw_optraj_indiv,
  data        = dat_op_raw,
  dv          = "opinion",
  dodge_width = 0.43   # 박스플롯 dodge와 비슷하게
  # sample_n  = 50     # 원하면 일부 참가자만 샘플링
)

fig_opinion_raw_optraj_indiv <- fig_opinion_raw_optraj_indiv +
  labs(
    x = "Time",
    y = "Opinion"
  ) +
  labs(
    y = "Opinion (Raw Value, -50 to 50)"   # <- 여기에 원하는 라벨 이름 넣기
  )  
#ggtitle(paste0("Opinion (raw + individual lines) – Trajectory group: ", op_group))

print(fig_opinion_raw_optraj_indiv)




## ---------------------------------------------------------
## [OPINION] ΔOpinion (signed, Tk−T0) 플롯
##   - dv = "opinion_delta"
##   - resp_delta에서 T1~T4 데이터만 들어있기 때문에
##     time_contrast = "baseline" 은 T1 vs T2~T4 라는 의미가 됨
## ---------------------------------------------------------

fig_opinion_signed_optraj <- plot_behavior_pretty(
  data = dat_op_delta,
  dv   = "opinion_delta",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # baseline = T1
) +
  labs(
    y = "Opinion Delta (Tk-T0)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
#ggtitle(paste0("ΔOpinion (Tk − T0) – Trajectory group: ", op_group))

print(fig_opinion_signed_optraj)


## ---------------------------------------------------------
## [OPINION] |ΔOpinion| (absolute delta) 플롯
##   - dv = "opinion_delta_abs"
## ---------------------------------------------------------

fig_opinion_abs_optraj <- plot_behavior_pretty(
  data = dat_op_delta,
  dv   = "opinion_delta_abs",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
) +
  labs(
    y = "|Opinion Delta| (Tk-T0)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
#ggtitle(paste0("|ΔOpinion| – Trajectory group: ", op_group))

print(fig_opinion_abs_optraj)









###### 여기서부터는 CONFIDENCE  ######

## ---------------------------------------------------------
## [CONFIDENCE] 시각화하고 싶은 trajectory 그룹 선택
##   - 가능한 값: "always_pos", "always_neg", "mixed"
## ---------------------------------------------------------



# Always positive - confidence

cf_group <- "always_pos"   # "always_neg", "mixed" 로 바꿔 가면서 사용

## ---------------------------------------------------------
## [CONFIDENCE] 선택한 trajectory 그룹에 속하는 (participant_id, task_type) 추출
## ---------------------------------------------------------

ids_cf <- traj_conf %>%
  dplyr::filter(conf_traj_group == cf_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_cf$participant_id)), "\n")


## ---------------------------------------------------------
## [CONFIDENCE] 이 그룹에 해당하는 raw / delta 데이터 만들기
##   - dat_cf_raw   : T0~T4 원자료 (confidence)
##   - dat_cf_delta : Tk−T0 delta 자료 (conf_delta, conf_delta_abs)
## ---------------------------------------------------------

dat_cf_raw <- resp_behavior %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_cf_delta <- resp_delta %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[CONFIDENCE] dat_cf_raw 행수   :", nrow(dat_cf_raw), "\n")
cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")

# ## ---------------------------------------------------------
# ## [CONFIDENCE] RAW (T0~T4) 플롯
# ##   - dv = "confidence"
# ## ---------------------------------------------------------
# 
# fig_conf_raw_cftraj <- plot_behavior_pretty(
#   data = dat_cf_raw,
#   dv   = "confidence",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) +
#   labs(
#     y = "Confidence (Raw Value, 0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
#   ) 
# 
# #  ggtitle(paste0("Confidence (raw) – Trajectory group: ", cf_group))
# 
# print(fig_conf_raw_cftraj)



## ---------------------------------------------------------
## [선택] RAW + 참가자별 꺾은선
## ---------------------------------------------------------

fig_conf_raw_cftraj_indiv <- plot_behavior_pretty(
  data = dat_cf_raw,
  dv   = "confidence",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)

fig_conf_raw_cftraj_indiv <- add_indiv_lines_to_behavior(
  p           = fig_conf_raw_cftraj_indiv,
  data        = dat_cf_raw,
  dv          = "confidence",
  dodge_width = 0.43
)

fig_conf_raw_cftraj_indiv <- fig_conf_raw_cftraj_indiv +
  labs(
    x = "Time",
    y = "Confidence"
  ) +
  labs(
    y = "Confidence (Raw Value, 0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
  # ggtitle(paste0("Confidence (raw + individual lines) – Trajectory group: ", cf_group))

print(fig_conf_raw_cftraj_indiv)



## ---------------------------------------------------------
## [CONFIDENCE] ΔConfidence (signed, Tk−T0) 플롯
##   - dv = "conf_delta"
## ---------------------------------------------------------

fig_conf_signed_cftraj <- plot_behavior_pretty(
  data = dat_cf_delta,
  dv   = "conf_delta",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # baseline = T1
) +
  labs(
    y = "Confidence Delta (Tk-T0)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
 # ggtitle(paste0("ΔConfidence (Tk − T0) – Trajectory group: ", cf_group))

print(fig_conf_signed_cftraj)



# ## ---------------------------------------------------------
# ## [CONFIDENCE] |ΔConfidence| (absolute delta) 플롯
# ##   - dv = "conf_delta_abs"
# ## ---------------------------------------------------------
# 
# fig_conf_abs_cftraj <- plot_behavior_pretty(
#   data = dat_cf_delta,
#   dv   = "conf_delta_abs",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) 
#   #ggtitle(paste0("|ΔConfidence| – Trajectory group: ", cf_group))
# 
# print(fig_conf_abs_cftraj)






### Always Negative - confidence

cf_group <- "always_neg"   # "always_neg", "mixed" 로 바꿔 가면서 사용

## ---------------------------------------------------------
## [CONFIDENCE] 선택한 trajectory 그룹에 속하는 (participant_id, task_type) 추출
## ---------------------------------------------------------

ids_cf <- traj_conf %>%
  dplyr::filter(conf_traj_group == cf_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_cf$participant_id)), "\n")


## ---------------------------------------------------------
## [CONFIDENCE] 이 그룹에 해당하는 raw / delta 데이터 만들기
##   - dat_cf_raw   : T0~T4 원자료 (confidence)
##   - dat_cf_delta : Tk−T0 delta 자료 (conf_delta, conf_delta_abs)
## ---------------------------------------------------------

dat_cf_raw <- resp_behavior %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_cf_delta <- resp_delta %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[CONFIDENCE] dat_cf_raw 행수   :", nrow(dat_cf_raw), "\n")
cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")

# ## ---------------------------------------------------------
# ## [CONFIDENCE] RAW (T0~T4) 플롯
# ##   - dv = "confidence"
# ## ---------------------------------------------------------
# 
# fig_conf_raw_cftraj <- plot_behavior_pretty(
#   data = dat_cf_raw,
#   dv   = "confidence",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) +
#   labs(
#     y = "Confidence (Raw Value, 0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
#   ) 
# 
# #  ggtitle(paste0("Confidence (raw) – Trajectory group: ", cf_group))
# 
# print(fig_conf_raw_cftraj)



## ---------------------------------------------------------
## [선택] RAW + 참가자별 꺾은선
## ---------------------------------------------------------

fig_conf_raw_cftraj_indiv <- plot_behavior_pretty(
  data = dat_cf_raw,
  dv   = "confidence",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)

fig_conf_raw_cftraj_indiv <- add_indiv_lines_to_behavior(
  p           = fig_conf_raw_cftraj_indiv,
  data        = dat_cf_raw,
  dv          = "confidence",
  dodge_width = 0.43
)

fig_conf_raw_cftraj_indiv <- fig_conf_raw_cftraj_indiv +
  labs(
    x = "Time",
    y = "Confidence"
  ) +
  labs(
    y = "Confidence (Raw Value, 0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
# ggtitle(paste0("Confidence (raw + individual lines) – Trajectory group: ", cf_group))

print(fig_conf_raw_cftraj_indiv)



## ---------------------------------------------------------
## [CONFIDENCE] ΔConfidence (signed, Tk−T0) 플롯
##   - dv = "conf_delta"
## ---------------------------------------------------------

fig_conf_signed_cftraj <- plot_behavior_pretty(
  data = dat_cf_delta,
  dv   = "conf_delta",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # baseline = T1
) +
  labs(
    y = "Confidence Delta (Tk-T0)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
# ggtitle(paste0("ΔConfidence (Tk − T0) – Trajectory group: ", cf_group))

print(fig_conf_signed_cftraj)



# ## ---------------------------------------------------------
# ## [CONFIDENCE] |ΔConfidence| (absolute delta) 플롯
# ##   - dv = "conf_delta_abs"
# ## ---------------------------------------------------------
# 
# fig_conf_abs_cftraj <- plot_behavior_pretty(
#   data = dat_cf_delta,
#   dv   = "conf_delta_abs",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) 
#   #ggtitle(paste0("|ΔConfidence| – Trajectory group: ", cf_group))
# 
# print(fig_conf_abs_cftraj)





# Mixed - confidence

cf_group <- "mixed"   # "always_neg", "mixed" 로 바꿔 가면서 사용

## ---------------------------------------------------------
## [CONFIDENCE] 선택한 trajectory 그룹에 속하는 (participant_id, task_type) 추출
## ---------------------------------------------------------

ids_cf <- traj_conf %>%
  dplyr::filter(conf_traj_group == cf_group) %>%
  dplyr::select(participant_id, task_type)

cat("\n[CONFIDENCE] 선택한 그룹:", cf_group, "\n")
cat(" - trajectory 개수 (participant × task_type):", nrow(ids_cf), "\n")
cat(" - 포함된 고유 참가자 수:",
    length(unique(ids_cf$participant_id)), "\n")


## ---------------------------------------------------------
## [CONFIDENCE] 이 그룹에 해당하는 raw / delta 데이터 만들기
##   - dat_cf_raw   : T0~T4 원자료 (confidence)
##   - dat_cf_delta : Tk−T0 delta 자료 (conf_delta, conf_delta_abs)
## ---------------------------------------------------------

dat_cf_raw <- resp_behavior %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

dat_cf_delta <- resp_delta %>%
  dplyr::semi_join(ids_cf, by = c("participant_id", "task_type")) %>%
  droplevels()

cat("\n[CONFIDENCE] dat_cf_raw 행수   :", nrow(dat_cf_raw), "\n")
cat("             dat_cf_delta 행수 :", nrow(dat_cf_delta), "\n")
cat("             (참가자 수:", length(unique(dat_cf_raw$participant_id)), ")\n")

# ## ---------------------------------------------------------
# ## [CONFIDENCE] RAW (T0~T4) 플롯
# ##   - dv = "confidence"
# ## ---------------------------------------------------------
# 
# fig_conf_raw_cftraj <- plot_behavior_pretty(
#   data = dat_cf_raw,
#   dv   = "confidence",
#   model = NULL,
#   p_adjust = p_adj_method,
#   add_time_contrasts = TRUE,
#   time_contrast      = "baseline"
# ) +
#   labs(
#     y = "Confidence (Raw Value, 0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
#   ) 
# 
# #  ggtitle(paste0("Confidence (raw) – Trajectory group: ", cf_group))
# 
# print(fig_conf_raw_cftraj)



## ---------------------------------------------------------
## [선택] RAW + 참가자별 꺾은선
## ---------------------------------------------------------

fig_conf_raw_cftraj_indiv <- plot_behavior_pretty(
  data = dat_cf_raw,
  dv   = "confidence",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)

fig_conf_raw_cftraj_indiv <- add_indiv_lines_to_behavior(
  p           = fig_conf_raw_cftraj_indiv,
  data        = dat_cf_raw,
  dv          = "confidence",
  dodge_width = 0.43
)

fig_conf_raw_cftraj_indiv <- fig_conf_raw_cftraj_indiv +
  labs(
    x = "Time",
    y = "Confidence"
  ) +
  labs(
    y = "Confidence (Raw Value, 0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
# ggtitle(paste0("Confidence (raw + individual lines) – Trajectory group: ", cf_group))

print(fig_conf_raw_cftraj_indiv)



## ---------------------------------------------------------
## [CONFIDENCE] ΔConfidence (signed, Tk−T0) 플롯
##   - dv = "conf_delta"
## ---------------------------------------------------------

fig_conf_signed_cftraj <- plot_behavior_pretty(
  data = dat_cf_delta,
  dv   = "conf_delta",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"   # baseline = T1
)+
  labs(
    y = "Confidence Delta (-100 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
# ggtitle(paste0("ΔConfidence (Tk − T0) – Trajectory group: ", cf_group))

print(fig_conf_signed_cftraj)



## ---------------------------------------------------------
## [CONFIDENCE] |ΔConfidence| (absolute delta) 플롯
##   - dv = "conf_delta_abs"
## ---------------------------------------------------------

fig_conf_abs_cftraj <- plot_behavior_pretty(
  data = dat_cf_delta,
  dv   = "conf_delta_abs",
  model = NULL,
  p_adjust = p_adj_method,
  add_time_contrasts = TRUE,
  time_contrast      = "baseline"
)+
  labs(
    y = "|Confidence Delta| (0 to 100)"   # <- 여기에 원하는 라벨 이름 넣기
  ) 
  #ggtitle(paste0("|ΔConfidence| – Trajectory group: ", cf_group))

print(fig_conf_abs_cftraj)

