#' @name methylResponseResid
#' @description
#' A function to identify associated CpGs and calculate residuals
#' @param sampletype 
#' @param out output folder

methylResponseResid <- function(sampletype = 'blood',
                                out = here('1-analyses/response/1-out/'), 
                                return = F){
   
  cat('Loading data ...')
  load(here('data/data_raw.Rdata')) # Load colData
  data <- as.data.frame(colData(data))
  data <- data[data$interventionId!='S',]
  
  assign('beta', load(here(paste0('1-analyses/variance-analysis/1-out/reliable_', sampletype, '.Rdata'))))  # Load reliable matrix
  beta_matrix <- get(beta[[1]])
  pheno_df <- get(beta[[2]])
  
  load('~/Dropbox/data/tirolgesund/pheno/methylation_pheno_params.Rdata') # Load params
  params <- pheno
  pheno <- pheno_df[pheno_df$interventionId!='S',]
  beta <- beta_matrix[,pheno$basename]
  # identical(pheno$basename, colnames(beta))
  pheno <- pheno |> dplyr::left_join(data, by = c('interventionId', 'subjectId', 'visitId')) # Append information 
  ind <- match(pheno$basename, params$basename)
  params <- params[ind,]
  # identical(pheno$basename, params$basename)
  if(sampletype == 'blood'){
    pheno <- cbind(pheno, params[,grepl('hepidish_', colnames(params))])
  } else {
    pheno$ic <- params$ic
  }
  cat('done.\nComputing correlation...')
  
  # Basic correlation
  cor <- apply(beta, 1, function(x){
    cor(x, pheno$time, method = 'spearman')
  })
  saveRDS(cor, file.path(out, paste0('cg_cor_', sampletype, '.Rds')))

  cg_assoc <- cor[abs(cor)>0.05]
  cat('done:', length(cg_assoc), 'associated CpGs.\nComputing residuals...')
  
  beta_assoc <- beta[names(cg_assoc),] # Subset associated CpGs
  
  # Adjust for neutrophils
  beta_assoc_resid <- apply(as.matrix(beta_assoc), 1, function(x){
    if(sampletype == 'blood'){
      fit <- lm(as.numeric(x) ~ pheno$hepidish_Neutro)
    } else {
      fit <- lm(as.numeric(x) ~ pheno$ic)
    }
    residuals(fit)
  })
  
  beta_assoc_resid <- t(beta_assoc_resid)
  colnames(beta_assoc_resid) <- colnames(beta_assoc)
  saveRDS(beta_assoc_resid, file.path(out, paste0('assoc_resid_', sampletype, '.Rds')))
  saveRDS(beta_assoc, file.path(out, paste0('assoc_', sampletype, '.Rds')))
  saveRDS(pheno, file.path(out, paste0('pheno_', sampletype, '.Rds')))
  
  cat('done.\nPivoting and computing delta...')
  
  # Compute delta-beta via pivoting
  pheno <- pheno |> dplyr::arrange(subjectId, visitId) # arrange by subjectId and time
  beta_assoc <- beta_assoc[,pheno$basename]
  tmp <- cbind(pheno, t(beta_assoc)) |> 
    tidyr::pivot_longer(cols = rownames(beta_assoc),
                        names_to = 'cg',
                        values_to = 'beta')
  tmp_delta <- tmp |> 
    tidyr::pivot_wider(id_cols = c(subjectId, cg),
                       names_from = visitId,
                       values_from = beta) |> 
    dplyr::mutate(across(c(M2, M4, M6), ~ . - M0))
  
  # Pivot back to a matrix
  beta_delta <- tmp_delta |>
    tidyr::pivot_longer(M2:M6,
                        names_to = 'visitId',
                        values_to = 'beta') |>
    dplyr::filter(!is.na(beta)) |> 
    tidyr::pivot_wider(id_cols = c(cg),
                       names_from = c('subjectId', 'visitId'),
                       values_from = 'beta') |> 
    tibble::column_to_rownames('cg')
  
  # Subset delta pheno
  pheno_delta <- pheno |>
    dplyr::mutate(id = paste0(subjectId, "_", visitId)) |> 
    dplyr::slice(match(colnames(beta_delta), id))
  # identical(colnames(beta_delta), pheno_delta$id)
  cat('done.\nComputing residuals...')
  
  # Adjust for neutrophils
  beta_delta_resid <- apply(as.matrix(beta_delta), 1, function(x){
    if(sampletype == 'blood'){
      fit <- lm(as.numeric(x) ~ pheno_delta$hepidish_Neutro)
    } else {
      fit <- lm(as.numeric(x) ~ pheno_delta$ic)
    }
    residuals(fit)
  })
  
  beta_delta_resid <- t(beta_delta_resid)
  colnames(beta_delta_resid) <- colnames(beta_delta)
  cat('done.\nSaving outputs...')
  
  saveRDS(beta_delta_resid, file.path(out, paste0('delta_resid_', sampletype, '.Rds')))
  saveRDS(beta_delta, file.path(out, paste0('delta_', sampletype, '.Rds')))
  saveRDS(pheno_delta, file.path(out, paste0('pheno_delta_', sampletype, '.Rds')))
  
  if(return == T){
  cat('done.\nReturning files.')
  return(list(beta_delta_resid = beta_delta_resid,
              beta_delta = beta_delta,
              pheno_delta = pheno_delta))
  } else {
    cat('done.')
  }
  
  }