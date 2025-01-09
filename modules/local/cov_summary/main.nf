// Define the `COV_SUMMARY` process that generates a summary table for all samples
process COV_SUMMARY {
    container ''
    conda "${moduleDir}/environment.yml"
    tag 'summary'
    cpus 10
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(breadthFile), path(depthStatsFile), path(bcfstatsFile)

    output:
    path("combined_summary.txt")

    script:
    """
    #!/usr/bin/env Rscript

    # Определите функцию create_summary_table
    create_summary_table <- function(sid, depthStatsFile, bcfstatsFile, breadthFile) {
      # Чтение данных из файла depthStatsFile
      depth <- read_table(depthStatsFile)

      # Чтение данных из файла breadthFile
      breads <- read_csv(breadthFile, col_names = FALSE) %>% pull(1) %>% as.numeric()

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
      bcfstats_data <- read.table(text = relevant_lines, header = FALSE, sep = "\\t", stringsAsFactors = FALSE)

      # Присвоим имена колонкам
      colnames(bcfstats_data) <- c("SN", "id", "key", "value")

      # Преобразуем данные в широкий формат
      bcfstats_wide <- bcfstats_data %>%
        select(key, value) %>%
        spread(key, value)

      # Создайте датафрейм с sid
      sid_df <- data.frame(sid = sid)

      # Объедините все данные в один датафрейм
      combined_df <- bind_cols(sid_df, depth, breads = breads, bcfstats_wide)

      return(combined_df)
    }

    # Создайте сводную таблицу
    summary_table <- create_summary_table($sid, $depthStatsFile, $bcfstatsFile, $breadthFile)

    # Сохраните результат в файл
    write.csv(summary_table, file = outputFile, row.names = FALSE)
    """
}