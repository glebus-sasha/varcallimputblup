// Define the `COV_STATS` process that calculates the depth of coverage and generates a plot and summary statistics
process COV_STATS {
    container 'glebusasha/r_env_image:latest'
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_STATS"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(mosdepth_summary), path(coverage_width)
    path reference_length

    output:
    tuple val(sid), path("${sid}_stats.csv")        , emit: cov_stats
    tuple val(sid), path("${sid}_coverage_plot.png"), emit: coverage_plot

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(readr)
    library(ggplot2)
    library(gridExtra)

    # Определение функции
    process_chromosome_data <- function(sid, filename, reference_length_file, coverage_width_file) {
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

    # Запись результата в CSV файл
    write_csv(result, "${sid}_stats.csv")

    # Построение violin plot только для данных покрытия (столбец mean)
    data_long <- data %>%
        filter(chrom != 'total') %>%
        mutate(group = ifelse(chrom %in% c(as.character(1:29), 'X', 'Y', 'MT'), "selected", "all")) %>%
        select(chrom, mean, group) %>%
        pivot_longer(cols = mean, names_to = "metric", values_to = "value")

    # Добавление столбца для типа шкалы
    data_long <- data_long %>%
        mutate(scale_type = "Linear") %>%
        bind_rows(data_long %>% mutate(scale_type = "Log"))

    # Построение графиков
    plot <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
        geom_violin() +
        facet_wrap(~ scale_type, scales = "free_y") +
        scale_y_continuous(trans = "identity") +
        scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(title = "Violin Plot of Chromosome Coverage", x = "Group", y = "Coverage Mean") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))

    # Сохранение картинки
    ggsave(filename = paste0(sid, "_coverage_plot.png"), plot = plot, width = 10, height = 8)

    }
    process_chromosome_data('${sid}', '${mosdepth_summary}', '${reference_length}', '${coverage_width}')
    """
}
