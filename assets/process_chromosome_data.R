

# Определение функции
process_chromosome_data <- function(sid, filename, reference_length_file, coverage_width_file, bcfstatsFile) {
  # Чтение таблицы из файла
  data <- read_table(filename)

  # Чтение значений reference_length и coverage_width из файлов
  reference_length <- as.numeric(read_lines(reference_length_file))
  coverage_width <- as.numeric(read_lines(coverage_width_file))

  # Фильтрация и суммирование для выбранных хромосом
  selected_chromosomes <- data %>%
      filter(chrom %in% c(as.character(1:29), 'X', 'Y', 'MT')) %>%
      summarise(
      length = sum(length),
      bases = sum(bases),
      mean = mean(mean),
      min = min(min),
      max = max(max)
      )

  # Фильтрация строки total из оригинальной таблицы
  original_total <- data %>%
      filter(chrom == 'total') %>%
      select(-chrom)

  # Объединение результатов в одну строку
  result <- bind_cols(
      original_total %>% rename_with(~ paste0("mosdepth:", ., "_all_chr")),
      selected_chromosomes %>% rename_with(~ paste0("mosdepth:", ., "_selected_chr"))
  )

  # Добавление колонки с процентом breadth
  result <- result %>%
      mutate(`breadth, %` = (coverage_width / reference_length) * 100) %>%
      mutate(sid = sid) %>%
      select(sid, everything())
  # Чтение данных из файла bcfstatsFile
  file_lines <- readLines(bcfstatsFile)

  # Найдите начало нужного отрывка с использованием регулярного выражения
  start_line <- grep("# SN\\t\\[2\\]id\\t\\[3\\]key\\t\\[4\\]value", file_lines)

  # Проверьте, найдена ли строка
  if (length(start_line) == 0) {
      stop("Не удалось найти начало нужного отрывка в файле.")
  }

  # Извлеките строки, начиная с найденной линии
  relevant_lines <- file_lines[(start_line + 1):(start_line + 10)]

  # Преобразуйте строки в датафрейм
  bcfstats_data <- read.table(text = relevant_lines, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  # Присвоим имена колонкам
  colnames(bcfstats_data) <- c("SN", "id", "key", "value")

  # Преобразуем данные в широкий формат
  bcfstats_wide <- bcfstats_data %>%
      select(key, value) %>%
      spread(key, value)

  # Добавление префикса к именам колонок
  colnames(bcfstats_wide) <- paste0("bcftools:", colnames(bcfstats_wide))

  combined_df <- bind_cols(result, bcfstats_wide)

  return(combined_df)
}