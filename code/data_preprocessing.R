# ============================================
# Refine participants_rows.csv by PROLIFIC APPROVED status
# Author: Soohwan's data wrangler script (R)
# ============================================

# Packages -------------------z--------------------------------------------
# (필요시 설치: install.packages(c("tidyverse","janitor")))
library(tidyverse)
library(janitor)
library(lubridate)




### Participant List Refining ###

# Paths (Windows에서도 슬래시는 / 로 써도 됩니다) -------------------------
raw_dir      <- "C:/Users/USER/Desktop/github/mutliAgentExperiment_dataAnalysis/raw_data"
refined_dir  <- "C:/Users/USER/Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data"

path_participants <- file.path(raw_dir, "participants_rows.csv")
# 파일명 오타 주의: 요청하신 경로 그대로 사용합니다 (demograhpic_1)
path_prolific_1   <- file.path(raw_dir, "prolific_export_demograhpic_1.csv")
path_prolific_2   <- file.path(raw_dir, "prolific_export_demographic_2.csv")

# Helper: 가장 먼저 매칭되는 컬럼명을 찾아주는 함수 -----------------------
pick_col <- function(df, candidates) {
  nm <- names(df)
  # 완전 일치 우선
  hit <- candidates[candidates %in% nm]
  if (length(hit) > 0) return(hit[1])
  # clean_names 기준(소문자+언더스코어)으로도 재시도
  df2 <- janitor::clean_names(df)
  nm2 <- names(df2)
  cand2 <- janitor::make_clean_names(candidates)
  hit2 <- cand2[cand2 %in% nm2]
  if (length(hit2) > 0) return(names(df)[match(hit2[1], nm2)])
  stop(sprintf(
    "Could not find any of columns: %s in [%s]",
    paste(candidates, collapse=", "),
    paste(names(df), collapse=", ")
  ))
}

# Read Prolific exports --------------------------------------------------
read_prolific <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- suppressMessages(readr::read_csv(path, guess_max = 1e5))
  df
}

p1 <- read_prolific(path_prolific_1) %>% clean_names()
p2 <- read_prolific(path_prolific_2) %>% clean_names()

# Prolific 파일 내 status / prolific_id 계열 컬럼 찾기 -------------------
# (Prolific export는 보통 'status'와 'participant_id' 또는 'prolific_id'가 있습니다)
status_candidates <- c("status", "submission_status", "participant_status")
id_candidates     <- c("prolific_id", "participant_id", "prolific_pid", "prolific id", "participant id", "pid")

get_approved_ids <- function(df) {
  status_col <- pick_col(df, status_candidates)
  id_col     <- pick_col(df, id_candidates)
  df %>%
    rename(.status = all_of(status_col),
           .pid    = all_of(id_col)) %>%
    mutate(.status = toupper(as.character(.status)),
           .pid    = as.character(.pid)) %>%
    filter(.status == "APPROVED", !is.na(.pid), .pid != "") %>%
    transmute(prolific_pid = .pid) %>%
    distinct()
}

approved_ids <- bind_rows(get_approved_ids(p1), get_approved_ids(p2)) %>% distinct()

# Read participants_rows -------------------------------------------------
participants <- suppressMessages(readr::read_csv(path_participants, guess_max = 1e5)) %>% clean_names()

# participants_rows에서 prolific id 컬럼 찾기 ----------------------------
part_id_col <- pick_col(participants, c("prolific_pid", "prolific_id", "participant_id", "pid"))
participants_norm <- participants %>%
  rename(prolific_pid = all_of(part_id_col)) %>%
  mutate(prolific_pid = as.character(prolific_pid))

# Filter: APPROVED ID만 남기기 ------------------------------------------
kept   <- participants_norm %>% semi_join(approved_ids, by = "prolific_pid")
dropped <- participants_norm %>% anti_join(approved_ids, by = "prolific_pid")

# 출력 디렉토리 생성 및 저장 -------------------------------------------
if (!dir.exists(refined_dir)) dir.create(refined_dir, recursive = TRUE)

# 엑셀 호환이 필요하면 write_excel_csv 사용(UTF-8 BOM 추가)
readr::write_csv(kept,    file.path(refined_dir, "participants_rows_refined.csv"))
readr::write_csv(dropped, file.path(refined_dir, "participants_rows_dropped_not_approved.csv"))
readr::write_csv(approved_ids, file.path(refined_dir, "approved_prolific_ids_from_exports.csv"))

# 로그 요약 --------------------------------------------------------------
cat("=== Refinement summary ===\n")
cat("Participants rows (original): ", nrow(participants_norm), "\n")
cat("Approved IDs (unique):       ", nrow(approved_ids), "\n")
cat("Kept rows:                    ", nrow(kept), "\n")
cat("Dropped rows (not approved):  ", nrow(dropped), "\n")
cat("Saved refined csv to:         ", file.path(refined_dir, "participants_rows_refined.csv"), "\n")


###############

# ----- Paths -----
refined_dir <- "C:/Users/USER/Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data"
infile      <- file.path(refined_dir, "participants_rows_refined.csv")
outfile     <- file.path(refined_dir, "participants_rows_refined_reindexed.csv")
mapfile     <- file.path(refined_dir, "participant_reindex_map.csv")
sumfile     <- file.path(refined_dir, "participants_summary_counts.csv")

# ----- Load -----
dat <- readr::read_csv(infile, guess_max = 1e5) %>% clean_names()

# 필수 컬럼 확인
stopifnot(all(c("participant_no","created_at") %in% names(dat)))

# created_at 파싱(문자열일 수 있음)
dat <- dat %>% mutate(
  created_at = suppressWarnings(ymd_hms(created_at, quiet = TRUE))
)

# ----- 중복 진단 -----
dup_counts <- dat %>%
  count(participant_no, name = "n") %>%
  filter(n > 1) %>%
  arrange(desc(n), participant_no)

n_dup_ids <- nrow(dup_counts)
dup_ids   <- dup_counts$participant_no

cat("중복 participant_no 개수:", n_dup_ids, "\n")
print(dup_counts)

# "later" 행 식별(같은 participant_no 내 created_at 순으로 2번째 이후)
dat_ranked <- dat %>%
  arrange(participant_no, created_at) %>%
  group_by(participant_no) %>%
  mutate(row_in_p = row_number()) %>%
  ungroup() %>%
  mutate(is_later_dup = participant_no %in% dup_ids & row_in_p > 1)

# later 행만 추출, 시간 순으로 정렬(재배정 순서)
later_rows <- dat_ranked %>% filter(is_later_dup) %>% arrange(created_at)

n_later <- nrow(later_rows)  # later 중복 행 수
cat("재배정할 later 행 수:", n_later, "\n")

# ----- 새 participant_no 할당 -----
# 옵션 1) 121..(120 + n_later) 로 확장(예: n_later=8이면 121..128)
# new_ids <- seq(121, 120 + n_later)

# 옵션 2) 55가 비어 있다면 55 먼저 쓰고, 나머지를 121..127
# (원한다면 아래 두 줄의 주석을 해제하고, 위 new_ids 라인을 주석 처리하세요)
 if (n_later == 8) {
   new_ids <- c(55, 121:127)
 }

# 충돌 방지 체크(이미 존재하는 id와 겹치지 않는지)
if (any(new_ids %in% dat$participant_no)) {
  stop("새로 배정할 participant_no 중 기존 값과 충돌이 있습니다: ",
       paste(intersect(new_ids, dat$participant_no), collapse=", "))
}

# 재배정(순서를 보장하기 위해 later_rows의 순서대로 맵핑)
reindex_map <- later_rows %>%
  mutate(new_participant_no = new_ids) %>%
  select(id, participant_no_old = participant_no, new_participant_no, created_at)

# dat_ranked에 적용
dat_reindexed <- dat_ranked %>%
  left_join(reindex_map %>% select(id, new_participant_no), by = "id") %>%
  mutate(participant_no = if_else(!is.na(new_participant_no), new_participant_no, participant_no)) %>%
  select(-new_participant_no)

# ----- 저장 -----
readr::write_csv(dat_reindexed %>% arrange(participant_no, created_at), outfile)
readr::write_csv(reindex_map, mapfile)

cat("재배정 완료. 저장 위치:\n")
cat(" - 데이터 :", outfile, "\n")
cat(" - 맵(로그):", mapfile, "\n")

# ----- 요약(분포 집계) -----
sum_condition <- dat_reindexed %>%
  count(condition_type, name = "n") %>% arrange(condition_type)
sum_order <- dat_reindexed %>%
  count(task_order, name = "n") %>% arrange(task_order)

# index 분포(0~5 고정)
count_0_5 <- function(x) {
  tibble(index = 0:5) %>%
    left_join(tibble(index = as.integer(x)) %>% count(index, name = "n"), by = "index") %>%
    mutate(n = replace_na(n, 0))
}
sum_info_idx <- count_0_5(dat_reindexed$informative_task_index) %>% mutate(which="informative_task_index")
sum_norm_idx <- count_0_5(dat_reindexed$normative_task_index)   %>% mutate(which="normative_task_index")

# 하나의 csv로 내보내기(섹션 구분을 위해 which/label 컬럼 포함)
summary_out <- bind_rows(
  sum_condition %>% mutate(which = "condition_type") %>% rename(label = condition_type),
  sum_order     %>% mutate(which = "task_order")      %>% rename(label = task_order),
  sum_info_idx  %>% rename(label = index),
  sum_norm_idx  %>% rename(label = index)
)

readr::write_csv(summary_out, sumfile)
cat("요약 통계 저장:", sumfile, "\n")



##################


# 필요한 패키지
# install.packages("tidyverse")
library(dplyr)
library(readr)
library(stringr)

# 경로 설정
base_raw <- "C:/Users/USER/Desktop/github/mutliAgentExperiment_dataAnalysis/raw_data"
base_out <- "C:/Users/USER/Desktop/github/mutliAgentExperiment_dataAnalysis/refined_data"
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

# 마스터 파일 탐색: raw_data 또는 refined_data 중 존재하는 곳 사용
master_candidates <- file.path(c(base_raw, base_out), "FINAL_participants_rows_refined_reindexed.csv")
master_candidates_exists <- file.exists(master_candidates)
if (!any(master_candidates_exists)) {
  stop("FINAL_participants_rows_refined_reindexed.csv 파일을 raw_data나 refined_data에서 찾을 수 없습니다. master_path를 직접 지정하세요.")
}
master_path <- master_candidates[which(master_candidates_exists)][1]

# 마스터 로드: id, participant_no 기준 테이블
master <- read_csv(master_path, show_col_types = FALSE) %>%
  select(id, participant_no) %>%
  distinct()

keep_ids <- master$id

# 입력 파일
in_files <- c(
  background = file.path(base_raw, "background_surveys_rows.csv"),
  post_self  = file.path(base_raw, "post_self_surveys_rows.csv"),
  post_open  = file.path(base_raw, "post_open_surveys_rows.csv")
)

# 출력 파일 (파일명은 동일, 위치만 refined_data로)
out_files <- file.path(base_out, basename(in_files))

# 정제 함수
refine_one <- function(in_path, out_path, master) {
  d <- read_csv(in_path, show_col_types = FALSE)
  n0 <- nrow(d)
  
  # 1) 마스터에 없는 participant_id 행 제거
  d_filtered <- d %>% semi_join(master, by = c("participant_id" = "id"))
  removed <- n0 - nrow(d_filtered)
  
  # 2) participant_no 불일치 개수(수정 전) 집계
  m_chk <- d_filtered %>%
    left_join(master %>% rename(right_no = participant_no), by = c("participant_id" = "id")) %>%
    mutate(mismatch = !is.na(participant_no) & !is.na(right_no) & participant_no != right_no)
  mismatches <- sum(m_chk$mismatch, na.rm = TRUE)
  
  # 3) participant_no를 마스터 기준으로 교정
  d_aligned <- d_filtered %>%
    left_join(master %>% rename(right_no = participant_no), by = c("participant_id" = "id")) %>%
    mutate(participant_no = coalesce(right_no, participant_no)) %>%
    select(-right_no) %>%
    arrange(participant_no)
  
  # 저장
  write_csv(d_aligned, out_path)
  
  list(
    file = basename(in_path),
    input_rows = n0,
    kept_rows = nrow(d_aligned),
    removed_not_in_master = removed,
    mismatches_fixed = mismatches,
    output_path = out_path
  )
}

# 세 파일 처리
summary_list <- Map(function(in_path, out_path) refine_one(in_path, out_path, master),
                    in_files, out_files)

# 요약 출력
summary <- bind_rows(lapply(summary_list, as.data.frame))
print(summary)

# 추가 점검: 저장된 파일의 participant_id가 모두 마스터에 존재하는지 확인
for (f in out_files) {
  d <- read_csv(f, show_col_types = FALSE)
  bad <- setdiff(d$participant_id, keep_ids)
  if (length(bad) > 0) {
    warning(sprintf("마스터에 없는 participant_id %d개 발견: %s", length(bad), basename(f)))
  }
}

