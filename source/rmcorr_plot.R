#' @param m1 Measure 1
#' @param m2 Measure 2
#' @param lab1 Label for measure 1
#' @param lab2 Label for measure 2
#' @return ggplot with repeated measures correlation

rmcorr_plot <- function(m1, m2){
  
  tmp <- dat |> 
    dplyr::select(subjectId, visitId, any_of(m1), any_of(m2)) |> 
    dplyr::filter(!is.na(.data[[m1]]) & !is.na(.data[[m2]]))
  
  # run rmcorr
  mod <- rmcorr::rmcorr(participant = subjectId,
                        measure1 = m1, measure2 = m2,
                        dataset = dat)
  
  lab1 <- res_lab |> dplyr::filter(measure1 == m1 & measure2 == m2) |> 
      dplyr::select(lab_m1, assaylab_m1)
  
  lab2 <- res_lab |> dplyr::filter(measure1 == m1 & measure2 == m2) |> 
      dplyr::select(lab_m2, assaylab_m2)
  
  # Rename labels
  labels <- data.frame(id = c(1, 2),
                       measures = c(m1, m2),
                       name = c(lab1$lab_m1, lab2$lab_m2),
                       assay = c(lab1$assaylab_m1, lab2$assaylab_m2))
  
  labels <- labels |> 
    dplyr::mutate(assay = paste0("<span style='color:",
                                 cols_for_assays[match(assay, names(cols_for_assays))],
                                 "'>", assay, "</span>"))
  
  r <- paste0("r<sub><i>rm</i></sub>", " = ",
              format(mod$r, digits = 2),
              " (p=", format(mod$p, digits = 2), ")")
  
  labels$label <- paste0("<b>", labels$assay, "<br>", labels$name, "</b><br>(scaled)")
  
  # plot
  plot <- tmp |> 
    dplyr::filter(!is.na(.data[[m1]] & !is.na(.data[[m2]]))) |> 
    dplyr::select(subjectId, visitId, .data[[m1]], .data[[m2]]) |> 
    ggplot(aes(x = .data[[m1]],
               y = .data[[m2]])) +
    geom_point(aes(colour = subjectId),
               size = 1.5,
               alpha = 0.9) +
    geom_line(aes(y = mod$model$fitted.values,
                  colour = subjectId),
              linetype = 1,
              alpha = 0.5) +
    geom_smooth(method = 'lm', se = F,
                linetype = 'dotted',
                colour = 'black') +
    scale_colour_manual(values = grDevices::colorRampPalette(cols[c(8, 1,2,3,5,4,6,7)])(156),
                        aesthetics = c('colour', 'fill')) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title = element_markdown()) +
    labs(x = labels$label[1],
         y = labels$label[2]) +
    annotate('richtext', label = r,
             Inf, -Inf,
             fill = NA, label.color = NA,
             hjust = 1,
             vjust = -5)
  
  return(plot)
}
