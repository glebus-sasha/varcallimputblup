library(tidyverse)
library(vcfR)
library(plotly)
library(ggplot2)
library(DT)
library(htmlwidgets)
library(htmltools)
library(optparse)

# Устанавливаем параметры командной строки
option_list <- list(
  make_option(c("--bed_file"), type = "character", help = "Path to BED file"),
  make_option(c("--interpretation_file"), type = "character", help = "Path to interpretation file")
)

# Парсим параметры командной строки
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Загружаем файлы по переданным путям
bed_file <- opt$bed_file
interpretation_file <- opt$interpretation_file
vcf_folder <- getwd()

# Получаем список всех VCF-файлов в текущей папке
vcf_files <- list.files(path = vcf_folder, pattern = "\\.vcf$", full.names = TRUE)

# Функция для чтения VCF и преобразования в tidy-формат
read_vcf_as_df <- function(file) {
  vcf <- tryCatch(
    read.vcfR(file, verbose = FALSE), 
    error = function(e) return(NULL) # Пропускаем файл, если ошибка
  )

  vcf_df <- vcfR2tidy(vcf, single_frame = TRUE)$dat %>%
    mutate(sid = tools::file_path_sans_ext(basename(file))) %>%   # Добавляем имя файла
    select(sid, everything())
  
  return(vcf_df)
}

# Читаем все файлы и объединяем их в один data frame
vcf_data <- map(vcf_files, read_vcf_as_df) %>%
  compact() %>%  # Убираем NULL-значения (пустые файлы)
  bind_rows() %>% 
  mutate(gt_GT_alleles = str_replace_all(gt_GT_alleles, "\\|", "/"))

SNV_53_interpretation <- read_csv(interpretation_file)
bed <- read_table(bed_file, 
                  col_names = c("Chromosome", "start", "end", "name"))

SNV_53_interpretation_expanded <- SNV_53_interpretation %>% 
  separate_rows(ID, sep = "\\+") %>% 
  left_join(bed, by = c("ID" = "name"))

sid_data <- vcf_data %>%
  select(sid)

# Создаем таблицу с каждым сочетанием позиции и образца (sid)
cross_table <- SNV_53_interpretation_expanded %>%
  crossing(sid_data)

# Объединяем с vcf_data по соответствующим столбцам
status <- cross_table %>%
  group_by(sid, Number) %>%  # Группируем по sid и Number
  mutate(Number_unique = paste(Number, row_number(), sep = ".")) %>%  # Добавляем суффиксы внутри каждой группы sid
  ungroup() %>%  # Разгруппируем после добавления суффиксов
  left_join(vcf_data, by = c("Chromosome" = "CHROM", "end" = "POS", "sid")) %>% 
  mutate(
    Positive = str_split(Positive, "\\s+"),  
    Neutral = str_split(Neutral, "\\s+"),
    Negative = str_split(Negative, "\\s+")
  ) %>%
  group_by(sid, Number) %>%  # Группируем по сэмплу и номеру мутации
  mutate(status = case_when(
    is.na(gt_GT_alleles) ~ NA,
    gt_GT_alleles %in% unlist(Positive) ~ "Positive",
    gt_GT_alleles %in% unlist(Neutral) ~ "Neutral",
    gt_GT_alleles %in% unlist(Negative) ~ "Negative",
    TRUE ~ "undefined"
  )) %>% 
  ungroup() %>% 
  select(sid, Number_unique, ID = ID.x, Chromosome, POS = end, REF, ALT, gt_GT_alleles, Positive, Neutral, Negative, status, Interpretation)

# Строим тепловую карту с ggplot
p <- ggplot(status, aes(x = Number_unique, y = sid, fill = status, text = Interpretation)) +
  geom_tile(color = "black", size = 0.5) +  # Добавляем черные оконтовки с толщиной 0.5
  theme_minimal() +
  labs(title = "Heatmap of Status by Sample", x = "Number", y = "Sample", fill = "Status") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# Преобразуем ggplot в интерактивный график с ggplotly
interactive_plot <- ggplotly(p, tooltip = "text")

# Сохраняем график как HTML-файл
htmlwidgets::saveWidget(interactive_plot, "carpet.html")

# Сохраняем таблицу как интерактивный HTML-файл
datatable(status) %>%
  saveWidget(file = "table.html")

# Преобразуем все столбцы, которые содержат списки, в строки
status_cleaned <- status %>%
  mutate(across(where(is.list), ~sapply(., paste, collapse = ", ")))

# Сохраняем таблицу как CSV-файл
write.csv(status_cleaned, "table.csv", row.names = FALSE)