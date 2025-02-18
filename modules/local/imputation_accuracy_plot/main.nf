// Define the `IMPUTATION_ACCURACY_PLOT` process that plots concordance accuracy
process IMPUTATION_ACCURACY_PLOT {
    conda "${moduleDir}/environment.yml"
    tag 'all_samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/IMPUTATION_ACCURACY_PLOT"
//	debug true
    errorStrategy 'ignore'
	
    input:
    tuple val(sid), path(errors_grp)

    output:
    tuple val(sid), path("${sid}_${chr}_accplot.png"), emit: accplot

    script:
    """
    #!/usr/bin/env Rscript    
    library(ggplot2)
    library(readr)
    library(dplyr)
    library(purrr)

    # Чтение данных и обработка
    df <- read_delim("${errors_grp}", delim = " ", col_names = FALSE, show_col_types = FALSE) %>%
        # Переименовываем колонки
        set_names(c("col1", "col2", "maf", "col4", "aggregate_r2")) %>%
        # Преобразуем столбцы в числовые значения
        mutate(
            maf = as.numeric(maf),
            aggregate_r2 = as.numeric(aggregate_r2)
        )

    # Определяем бины и подписи для оси X
    mybins <- c(0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.3, 0.4, 0.50, 0.60)
    labels <- c('0.02', '0.05', '0.1', '0.2', '0.5', '1', '2', '5', '10', '15', '20', '30', '40', '50', '60')

    # Строим график с использованием ggplot2
    p <- df %>%
        ggplot(aes(x = maf * 100, y = aggregate_r2)) +  # MAF в процентах
        geom_line(linetype = "dashed", linewidth = 1) +  # Изменили 'size' на 'linewidth'
        geom_point(size = 3, alpha = 0.8) +
        scale_x_log10(labels = labels) +  # Логарифмическая шкала с нужными подписями
        labs(
            x = "Minor allele frequency (%)",
            y = expression(r^2~"imputed vs true genotypes"),
            title = "Imputation accuracy"
        ) +
        theme_minimal(base_size = 20) +
        theme(
            panel.grid.major = element_line(linetype = "dashed"),
            legend.position = "bottom"
        ) +
        ylim(0.1, 1.05) +
        xlim(0.05, 60)  # X-ось тоже в процентах

    # Сохраняем график
    ggsave("${sid}_${chr}_accplot.png", p, width = 16, height = 12, dpi = 300)
    """
}