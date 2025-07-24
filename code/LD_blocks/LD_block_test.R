library(bigsnpr)
library(dplyr)
library(purrr)
library(writexl)

THR_R2 <- 0.01

# start_time <- Sys.time()
# bedfile <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr1_filtered_200k.bed"
# backingfile <- tempfile()  # Temporary file for bigSNP data
# obj.bigsnp <- snp_readBed2(
#   bedfile,
#   backingfile = backingfile
# ) %>% 
#   snp_attach()
# 
# obj.G<- obj.bigsnp$genotypes
# map <- obj.bigsnp$map
# POS <- obj.bigsnp$map$physical.pos
# 
# # genetic_map <- bigreadr::fread2(
# #    "/data2/morgante_lab/ukbiobank_projects/prs_effects/other/1000-genomes-genetic-maps/interpolated_OMNI/chr22.OMNI.interpolated_genetic_map.gz",
# #    col.names = c("id", "pos", "pos_cM"))
# # genetic_map$chr <- 22
# # 
# # POS2 = snp_asGeneticPos(map$chromosome, POS, 
# #                         dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/other/1000-genomes-genetic-maps/interpolated_OMNI/")
# 
# corr_matrix <- snp_cor(obj.G, infos.pos = POS, thr_r2 = THR_R2)
# 
# cost_time = (Sys.time() - start_time)
# 
# cat("Generating correlation matrix for Chromosome 1 needs", cost_time, "\n")
# 
# saveRDS(corr_matrix, file = paste0("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect/chr1_cormatrix_geneticmap.rds"))
# corr_matrix <- readRDS("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect/chr1_cormatrix_geneticmap.rds")
# Loop through R2 thresholds
# 
# for (r2_value in seq(0.3, 0.8, by = 0.1)) {
#   cat("\nRunning snp_ldsplit with max_r2 =", r2_value, "\n")
#   iter_start <- Sys.time()
#   
#   ld_blocks <- snp_ldsplit(
#     corr = corr_matrix,
#     thr_r2 = THR_R2,
#     min_size = 100,
#     max_size = 3000,
#     max_K = 300,
#     max_r2 = r2_value
#   )
#   
#   ld_blocks_fixed <- ld_blocks %>%
#     mutate(across(where(is.list), 
#                   ~ map_chr(., ~ paste(.x, collapse = ","))))
#   file_prefix <- sprintf("ld_blocks_chr1_r2_%.1f", r2_value)
#   
#   # Save as Excel file
#   
#   excel_file <- file.path("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect", 
#                           paste0(file_prefix, ".xlsx"))
#   write_xlsx(ld_blocks_fixed, path = excel_file)
#   
#   rds_file <- file.path("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect", 
#                         paste0(file_prefix, ".rds"))
#   saveRDS(ld_blocks, file = rds_file)
#   
#   iter_time <- Sys.time() - iter_start
#   cat("Completed in:", format(iter_time), "\n")
#   cat("Saved to:", excel_file, "\n")
#   
# }
# 
# 
# gc()
# rm(list = ls())
# 
# THR_R2 <- 0.01
# 
# start_time <- Sys.time()
# 
# bedfile <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr6_filtered_200k.bed"
# backingfile <- tempfile()  # Temporary file for bigSNP data
# obj.bigsnp <- snp_readBed2(
#   bedfile,
#   backingfile = backingfile
# ) %>% 
#   snp_attach()
# 
# obj.G<- obj.bigsnp$genotypes
# map <- obj.bigsnp$map
# POS <- obj.bigsnp$map$physical.pos

# genetic_map <- bigreadr::fread2(
#    "/data2/morgante_lab/ukbiobank_projects/prs_effects/other/1000-genomes-genetic-maps/interpolated_OMNI/chr22.OMNI.interpolated_genetic_map.gz",
#    col.names = c("id", "pos", "pos_cM"))
# genetic_map$chr <- 22
# 
# POS2 = snp_asGeneticPos(map$chromosome, POS, 
#                         dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/other/1000-genomes-genetic-maps/interpolated_OMNI/")

# corr_matrix <- snp_cor(obj.G, infos.pos = POS, thr_r2 = THR_R2)
# 
# cost_time = (Sys.time() - start_time)
# 
# cat("Generating correlation matrix for Chromosome 6 needs", cost_time, "\n")
# 
# saveRDS(corr_matrix, file = paste0("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect/chr6_cormatrix_geneticmap.rds"))

corr_matrix <- readRDS("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect/chr6_cormatrix_geneticmap.rds")
# Loop through R2 thresholds
for (r2_value in seq(0.1, 0.8, by = 0.1)) {
  cat("\nRunning snp_ldsplit with max_r2 =", r2_value, "\n")
  iter_start <- Sys.time()
  
  ld_blocks <- snp_ldsplit(
    corr_matrix,
    thr_r2 = THR_R2,
    min_size = 200,
    max_size = 5000,
    max_K = 300,
    max_r2 = r2_value
  )
  
  ld_blocks_fixed <- ld_blocks %>%
    mutate(across(where(is.list), 
                  ~ map_chr(., ~ paste(.x, collapse = ","))))
  file_prefix <- sprintf("ld_blocks_chr6_r2_%.1f", r2_value)
  
  # Save as Excel file
  
  excel_file <- file.path("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect", 
                          paste0(file_prefix, ".xlsx"))
  write_xlsx(ld_blocks_fixed, path = excel_file)
  
  rds_file <- file.path("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect", 
                        paste0(file_prefix, ".rds"))
  saveRDS(ld_blocks, file = rds_file)
  
  iter_time <- Sys.time() - iter_start
  cat("Completed in:", format(iter_time), "\n")
  cat("Saved to:", excel_file, "\n")
  
}

# 
# for (r2_value in seq(0.5, 0.8, by = 0.1)) {
#   cat("\nRunning snp_ldsplit with max_r2 =", r2_value, "\n")
#   iter_start <- Sys.time()
#   
#   ld_blocks <- snp_ldsplit(
#     corr_matrix,
#     thr_r2 = THR_R2,
#     min_size = 100,
#     max_size = 3000,
#     max_K = 300,
#     max_r2 = r2_value
#   )
#   
#   ld_blocks <- ld_blocks %>%
#     mutate(
#       block_starts_phys = map(all_last, ~ {
#         starts <- c(1, .x[-length(.x)] + 1)
#         map$physical.pos[starts]
#       }),
#       
#       block_ends_phys = map(all_last, ~ map$physical.pos[.x]),
#       
#       block_spans_phys = map2(block_ends_phys, block_starts_phys, ~ .x - .y),
#       
#       block_starts_str = map_chr(block_starts_phys, ~ paste(.x, collapse = ",")),
#       block_ends_str = map_chr(block_ends_phys, ~ paste(.x, collapse = ",")),
#       block_spans_str = map_chr(block_spans_phys, ~ paste(.x, collapse = ","))
#     )
#   
#   excel_file <- sprintf("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect/ld_blocks_withposition_r2_%.1f.xlsx", r2_value)
#   write_xlsx(ld_blocks, path = excel_file)
#   
#   rds_file <- file.path("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/com_ldetect", 
#                         paste0(file_prefix, ".rds"))
#   saveRDS(ld_blocks, file = rds_file)
#   
#   iter_time <- Sys.time() - iter_start
#   cat("Completed in:", format(iter_time), "\n")
#   cat("Saved to:", excel_file, "\n")
# 
# }


