# Interaction-values-positive-control-ratios-first-assay-boxplot
Interaction values/positive control ratios first assay boxplot

setwd("C:/Users/Gebruiker/Documents/PPH SLA")

split_luc_data <- read.table(
  file="split_luc_LBD16.txt",
  header = FALSE,
  sep = "\t"
)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

colnames(split_luc_data)[1] <- "Row"


colnames(split_luc_data)[2:13] <- as.character(1:12)

split_luc_data <- split_luc_data[, colSums(!is.na(split_luc_data)) > 0]

split_luc_data <- subset(split_luc_data, select = -V14)


data_long <- split_luc_data %>%
  pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%
  mutate(Well = paste0(Row, Column))

data_long <- data_long %>%
  filter(!Row %in% c("G", "H")) %>%
  mutate(Value = suppressWarnings(as.numeric(Value))) %>%  
  filter(!is.na(Value))  

data_long_log <- data_long %>%
  mutate(Value_log = log10(Value))

data_long_log <- data_long_log %>%
  mutate(Column = as.character(Column))


pos_ctrl_cols <- c("1", "5", "9")

pos_median <- data_long_log %>%
  filter(Column %in% pos_ctrl_cols) %>%
  summarise(median_pos = median(Value_log, na.rm = TRUE)) %>%
  pull(median_pos)



data_norm <- data_long_log %>% mutate(Norm_value = Value_log - pos_median)

summary_stats <- data_norm %>%
  group_by(Row) %>%
  summarise(
    mean_norm = mean(Norm_value, na.rm = TRUE),
    sd_norm   = sd(Norm_value, na.rm = TRUE),
    n         = n()
  )

print(summary_stats)

group_map <- data.frame(
  Column = as.character(1:12),
  Group = c(
    "nLUC_LBD16-cLUC_LBD16",   
    "nLUC_LBD16-cLUC_ZAT6",   
    "nLUC_LBD16-cLUC_LBD14",     
    "nLUC_LBD14-cLUC_LBD16",     
    "nLUC_LBD16-cLUC_LBD16",   
    "nLUC_LBD16-cLUC_ZAT6",   
    "nLUC_LBD16-cLUC_ARF7",    
    "nLUC_ARF7-cLUC_LBD16",    
    "nLUC_LBD16-cLUC_LBD16",   
    "nLUC_LBD16-cLUC_ZAT6",   
    "nLUC_LBD16-cLUC_ARF8",     
    "nLUC_SWI3B-cLUC_LBD16"      
  )
)

data_grouped <- data_norm %>%
  mutate(Column = as.character(Column)) %>%
  left_join(group_map, by = "Column")


data_grouped <- data_grouped %>%  
  mutate(
    Column = as.integer(Column),
    row_index   = match(Row, LETTERS) - 1,           
    leaf_in_row = ceiling(Column / 4),              
    leaf_number = row_index * 3 + leaf_in_row,       
    leaf_id     = paste0("Leaf ", leaf_number)
  ) %>%
  select(-row_index, -leaf_in_row, -leaf_number)


lighter_wells <- c("A1", "A2", "A3", "A4", "A9", "A10", "A11", "A12", "B9", "B10", "B11", "B12", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "D9", "D10", "D11", "D12", "E6", "E9", "E10", "E11", "E12", "F5", "F9", "F10", "F11", "F12") 

data_grouped <- data_grouped %>%
  mutate(LeafColor = ifelse(Well %in% lighter_wells, "Light", "Dark"))

library(dplyr)

control_label <- "nLUC_LBD16-cLUC_LBD16"

control_stats <- data_grouped %>%
  filter(Group == control_label) %>%
  group_by(LeafColor) %>%
  summarise(
    control_mean_log = mean(Value_log, na.rm = TRUE),
    control_sd_log   = sd(Value_log,   na.rm = TRUE),
    n_control        = n(),
    .groups = "drop"
  )

df_ratio <- data_grouped %>%
  left_join(control_stats, by = "LeafColor") %>%
  mutate(
    ratio_log10 = Value_log - control_mean_log,           
    ratio_fold  = 10^(ratio_log10)                        
  )

summary_tbl <- df_ratio %>%
  group_by(LeafColor, Group) %>%
  summarise(
    n           = n(),
    mean_log10  = mean(ratio_log10, na.rm = TRUE),
    median_log10= median(ratio_log10, na.rm = TRUE),
    sd_log10    = sd(ratio_log10, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)


df_ratio <- df_ratio %>%
  group_by(LeafColor, Group) %>%
  mutate(median_grp = median(ratio_log10, na.rm = TRUE)) %>%
  ungroup()


library(tidyverse)
library(ggbeeswarm) 
library(dplyr)

df_ratio <- df_ratio %>%
  mutate(Group_order = fct_reorder(Group, median_grp, .desc = TRUE))


control_stats <- data_grouped %>%
  dplyr::filter(Group == control_label) %>%
  dplyr::group_by(LeafColor) %>%
  dplyr::summarise(
    control_mean_log = mean(Value_log, na.rm = TRUE),
    .groups = "drop"
  )


df_ratio <- data_grouped %>%
  dplyr::left_join(control_stats, by = "LeafColor") %>%
  dplyr::mutate(
    log_ratio = Value_log - control_mean_log
  )


levels_wanted <- c("Light", "Dark")    
df_ld <- df_ratio %>%
  dplyr::filter(LeafColor %in% levels_wanted)

pos_dodge <- position_dodge(width = 0.6)


order_vec <- c("nLUC_LBD16-cLUC_LBD16", "nLUC_LBD16-cLUC_LBD14", "nLUC_LBD14-cLUC_LBD16", "nLUC_LBD16-cLUC_ARF7", "nLUC_ARF7-cLUC_LBD16", "nLUC_LBD16-cLUC_ARF8", "nLUC_SWI3B-cLUC_LBD16", "nLUC_LBD16-cLUC_ZAT6")  # <-- PAS AAN
df_ratio <- df_ratio %>%
  dplyr::mutate(Group = forcats::fct_relevel(Group, order_vec))


ggplot(df_ratio, aes(x = Group, y = log_ratio, color = Group)) +
  geom_boxplot(fill = NA, width = 0.6, outlier.shape = NA, size = 0.7) +
  ggbeeswarm::geom_quasirandom(width = 0.25, alpha = 0.9, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "tomato") +
  labs(
    title = "Log10-ratio Interaction / Positive control",
    x = "Interaction Pair",
    y = "log10-ratio (Interaction / Pos ctrl)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 55, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

ggplot(df_ratio, aes(x = Group, y = log_ratio)) +
   geom_boxplot(aes(color = Group), fill = NA, width = 0.6, outlier.shape = NA, size = 0.7, show.legend = FALSE) + 
   ggbeeswarm::geom_quasirandom( aes(color = LeafColor), width = 0.25, alpha = 0.9, size = 2 ) + 
   scale_color_manual( values = c("Dark" = "blue", "Light" = "red"), name = "Sample colour" ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "tomato") + 
  labs( title = "Log10-ratio Interaction / Positive control", 
        x = "Interaction Pair", 
        y = "log10-ratio (Interaction / Pos ctrl)" ) + 
  theme_bw(base_size = 12) + 
  theme( axis.text.x = element_text(angle = 55, hjust = 1), 
         plot.title = element_text(face = "bold") )
