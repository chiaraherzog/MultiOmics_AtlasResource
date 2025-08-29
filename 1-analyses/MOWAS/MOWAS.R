#' @name MOWAS
#' @description MOWAS wrapper function
#' @param phenotype phenotype to be queried, must be one of age, bmi, smoking_current, smoking_ever, vo2max, menopause, activity, alcohol
#' @param out.folder output folder
#' @return Returns MOWAS results

MOWAS <- function(phenotype,
                   out.folder){
  
  cat("MOWAS for:", phenotype, '\n')

  # Load non-methylation data -----------------
  if(phenotype == 'vo2max'){
    load(here('data/data_raw.Rdata'))
    vo2 <- suppressMessages(suppressWarnings(as.data.frame(wideFormat(data['vo2max',data$visitId=='M0',], colDataCols = 'subjectId')) |> dplyr::filter(!is.na(subjectId))))
  }
  
  load(here('data/data_normalized.Rdata'))
  
  # Filter individuals
  data <- data[,data$visitId =='M0']  # Keep only M0 data
  if(phenotype == 'smoking_ever'){
    data <- data[,data$smoking_ever!='unknown',]
  } else if(phenotype == 'vo2max'){
    data <- data[,match(vo2$subjectId, data$subjectId)]
  } else if (phenotype == 'menopause'){
    data <- data[,!grepl("[*]", data$mpstatrs)]
  } else if (phenotype == 'activity'){
    data <- data[,!is.na(data$intactcurr)]
  } else if (phenotype == 'alcohol'){
    data <- data[,data$etoh_curr != 'unknown']
  } else if (phenotype == 'alcohol_units'){
    data <- data[,data$etohu_curr != 'unknown']
  }
  
  pheno <- as.data.frame(colData(data))
  
  if(phenotype == 'vo2max'){
    if(!identical(pheno$subjectId, vo2$subjectId)){stop("VO2max: names not identical, check")}
    pheno$vo2max <- vo2$Functional.sports.exam_vo2max
  } else if(phenotype == 'alcohol_units'){
    pheno <- pheno |> dplyr::mutate(etohu_curr_pseudonum = case_when(
      
                                                            etohu_curr == 'no alcohol' ~ 0,
                                                            etohu_curr == '1-3 units' ~ 2,
                                                            etohu_curr == '11-15 units' ~ 13,
                                                            etohu_curr == '16-20 units' ~ 18,
                                                            etohu_curr == '4-6 units' ~ 5,
                                                            etohu_curr == '7-10 units' ~ 8,
                                                            etohu_curr == 'less than 1 unit' ~ 0.5,
                                                            etohu_curr == 'no alcohol' ~ 0)
                                    )
    
  }
  
  # Get variables
  load(here("data/vars.Rdata"))
  
  # Prepare non-methylation data  -----------------
  suppressMessages(
    suppressWarnings({
      
      v <- unique(vars[vars$assay=='Blood haemogram',]$x)
      blood_haemogram <- t(as.matrix(data[v,,1][[1]]))
    
      v <- unique(vars[vars$assay=='Body composition',]$x)
      body <- t(as.matrix(data[v,,'Body composition'][[1]]))
    
      v <- unique(vars[vars$assay=='Skin histology and transepidermal water loss assay',]$x)
      skin <- t(as.matrix(data[v,,'Skin histology and transepidermal water loss assay'][[1]]))
    
      v <- unique(vars[vars$assay=='Functional sports exam',]$x)
      funct <- as.data.frame(data[v,,'Functional sports exam'][[1]]) |> janitor::remove_empty(which = 'rows') |> t()
    
      v <- unique(vars[vars$assay=='Vascular and body sonography',]$x)
      vasc <- t(as.matrix(data[v,,'Vascular and body sonography'][[1]]))
    
      v <- unique(vars[vars$assay=='Flow cytometry: T cell staining',]$x)
      flow_t <- t(as.matrix(data[v,,'Flow cytometry: T cell staining'][[1]]))
    
      v <- unique(vars[vars$assay=='Flow cytometry: white blood cell staining',]$x)
      flow_wb <- t(as.matrix(data[v,,'Flow cytometry: white blood cell staining'][[1]]))
    
      v <- unique(vars[vars$assay=='Flow cytometry: unstimulated T cells',]$x)
      flow_unstim <- t(as.matrix(data[v,,'Flow cytometry: unstimulated T cells'][[1]]))
      
      v <- unique(vars[vars$assay=='Flow cytometry: stimulated T cells',]$x)
      flow_stim <- t(as.matrix(data[v,,'Flow cytometry: stimulated T cells'][[1]]))
      
      v <- unique(vars[vars$assay=='Urine nuclear magnetic resonance: normalized',]$x)
      urine_metab <- t(as.matrix(data[v,,'Urine nuclear magnetic resonance: normalized'][[1]]))
    
      v <- unique(vars[vars$assay=='Saliva nuclear magnetic resonance: normalized',]$x)
      saliva_metab <- t(as.matrix(data[v,,'Saliva nuclear magnetic resonance: normalized'][[1]]))
    
      v <- unique(vars[vars$assay=='Stool microbiome: families_clr',]$x)
      faecal_fam <- t(as.matrix(data[v,,'Stool microbiome: families_clr'][[1]]))
    
      v <- unique(vars[vars$assay=='Stool microbiome: ASVs_clr',]$x)
      faecal <- t(as.matrix(data[v,,'Stool microbiome: ASVs_clr'][[1]]))
    
      v <- unique(vars[vars$assay=='Saliva microbiome: families_clr',]$x)
      saliva_microb_fam <- t(as.matrix(data[v,,'Saliva microbiome: families_clr'][[1]]))
      
      v <- unique(vars[vars$assay=='Saliva microbiome: ASVs_clr',]$x)
      saliva_microb <- t(as.matrix(data[v,,'Saliva microbiome: ASVs_clr'][[1]]))
      
      v <- unique(vars[vars$assay=='Composite methylation scores: blood',]$x)
      blood_comp <- t(as.matrix(data[v,,'Composite methylation scores: blood'][[1]]))
      
      v <- unique(vars[vars$assay=='Composite methylation scores: buccal',]$x)
      buccal_comp <- t(as.matrix(data[v,,'Composite methylation scores: buccal'][[1]]))
      
      v <- unique(vars[vars$assay=='Composite methylation scores: cervical',]$x)
      cervical_comp <- t(as.matrix(data[v,,'Composite methylation scores: cervical'][[1]]))
      
      })
  )
  
  for(x in c('blood_haemogram', 'body', 'skin', 'funct', 'vasc',
             'flow_t', 'flow_wb', 'urine_metab', 'saliva_metab',
             'faecal_fam', 'faecal', 'saliva_microb', 'saliva_microb_fam',
             'blood_comp', 'buccal_comp', 'cervical_comp',
             'flow_stim', 'flow_unstim')){
    mat <- get(x) # Retrieve the matrix object by its name
    rownames(mat) <- substr(rownames(mat), 1, 4) # Modify the rownames
    assign(x, mat) # Assign the modified matrix back to the variable
  }
  
  items <- list(blood = blood_haemogram,
                body = body, skin = skin, funct = funct, vasc = vasc,
                flow_t = flow_t,
                flow_wb = flow_wb,
                flow_stim = flow_stim,
                flow_unstim = flow_unstim,
                urine_metab = urine_metab,
                saliva_metab = saliva_metab,
                faecal = faecal, faecal_fam = faecal_fam,
                saliva_microb = saliva_microb,
                saliva_microb_fam = saliva_microb_fam,
                blood_comp = blood_comp,
                buccal_comp = buccal_comp,
                cervical_comp = cervical_comp)
  
  # Compute non-methylation associations  -----------------
  
  out <-  if(phenotype == 'age'){
    # AGE
    lapply(items, function(x){
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      
      apply(x, 2, function(var){
        fit <- lm(var ~ age) |> 
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |> 
        bind_rows(.id = 'var')
    })
    
  }  else if (phenotype == 'smoking_curr'){
    
    # Smoking current
    lapply(items, function(x){
    smk_curr <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$smoking_curr)
    age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
    
    apply(x, 2, function(var){
      
      tryCatch({
      fit <- lm(var ~ smk_curr + age) |> 
        broom::tidy() |> 
        dplyr::filter(term != '(Intercept)')
      }, error = function(e){
      tibble::tibble(term = NA_character_, estimate = NA_real_, std.error = NA_real_, 
                     statistic = NA_real_, p.value = NA_real_)
      }) 
      }) |> 
      bind_rows(.id = 'var')
    
    })
  
    } else if (phenotype == 'smoking_ever'){
      
    # Smoking ever
    lapply(items, function(x){
    smk_ever <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$smoking_ever)
    age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
    
    apply(x, 2, function(var){
      fit <- lm(var ~ smk_ever + age) |> 
        broom::tidy() |> 
        dplyr::filter(term != '(Intercept)')
      }) |> 
      bind_rows(.id = 'var')
    })
      
  } else if (phenotype == 'bmi'){
    # BMI
    lapply(items, function(x){
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      bmi <- pheno[match(rownames(x), pheno$subjectId),]$bmi_at_consent
      
      apply(x, 2, function(var){
        fit <- lm(var ~ bmi + age) |> 
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |> 
        bind_rows(.id = 'var')
    })
    
  }  else if (phenotype == 'vo2max'){
    
    # VO2Max 
    lapply(items, function(x){
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      vo2max <- pheno[match(rownames(x), pheno$subjectId),]$vo2max
      
      apply(x, 2, function(var){
        fit <- lm(var ~ vo2max + age) |> 
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |> 
        bind_rows(.id = 'var')
    })
    
  } else if (phenotype == 'menopause'){
    # menopause ever
    lapply(items, function(x){
      mp <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$mpstatrs)
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      
      apply(x, 2, function(var){
        fit <- glm(var ~ mp + age) |> 
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |> 
        bind_rows(.id = 'var')
    })
    
  } else if (phenotype == 'activity'){
    # Activity
    lapply(items, function(x){
      
      activity <- pheno[match(rownames(x), pheno$subjectId),]$intactcurr
      activity <- ifelse(activity == '0', 'no', 'yes')
      activity = factor(activity, levels = c('no', 'yes')) # code acitivity
      # smk <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$smoking_curr)
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      
      apply(x, 2, function(var) {
          fit <- lm(var ~ activity + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)') 
      }) |> 
        bind_rows(.id = 'var')
    })
    
  } else if (phenotype == 'alcohol'){
    
    # Alcohol
    lapply(items, function(x){
      
      alcohol <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$etoh_curr)
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      smoking <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$smoking_curr)
      
      apply(x, 2, function(var){
        
        tryCatch({
        fit <- glm(var ~ alcohol + age + smoking) |> 
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
        }, error = function(e){
          fit <- glm(var ~ alcohol + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        })
      }) |> 
        bind_rows(.id = 'var')
    })
  } else if(phenotype == 'alcohol_units'){
    
    
    # Alcohol
    lapply(items, function(x){
      
      alcohol <- pheno[match(rownames(x), pheno$subjectId),]$etohu_curr_pseudonum
      age <- pheno[match(rownames(x), pheno$subjectId),]$age_at_consent
      smoking <- as.factor(pheno[match(rownames(x), pheno$subjectId),]$smoking_curr)
      
      apply(x, 2, function(var){
        
        tryCatch({
          fit <- glm(var ~ alcohol + age + smoking) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }, error = function(e){
          fit <- glm(var ~ alcohol + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        })
      }) |> 
        bind_rows(.id = 'var')
    })
  }
  
  total <- out |> bind_rows(.id = 'type')
  save(total, file = file.path(out.folder, paste0("mowas_nonmethylation_", phenotype, ".Rdata")))
  
  # Prepare methylation data  -----------------
  cat("Non-methylation data complete; preparing methylation data\n")
  
  # Loading reliable subsets:
  load(here('1-analyses/variance-analysis/1-out/reliable_buccal.Rdata'))
  load(here('1-analyses/variance-analysis/1-out/reliable_blood.Rdata'))
  load(here('1-analyses/variance-analysis/1-out/reliable_cervical.Rdata'))
  
  # Getting relevant samples/basenames
  grabSamples <- function(beta, basename){
    methyl <- t(beta[,na.omit(pheno[[basename]])])
    rownames(methyl) <- pheno[match(rownames(methyl), pheno[[basename]]),]$subjectId
    return(methyl)
  }
  
  methyl_blood <- grabSamples(beta_reliable_blood, 'basename_blood')
  methyl_buccal <- grabSamples(beta_reliable_buccal, 'basename_buccal')
  methyl_cervical <- grabSamples(beta_reliable_cervical, 'basename_cervical')
  
  # Filter CpGs based on basic correlation with phenotype
  filterCorr <- function(beta, phenotype){
    pheno_tmp <- pheno[match(rownames(beta), pheno$subjectId),]
    
    if(phenotype == 'age'){
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$age_at_consent)})
    } else if(phenotype == 'bmi'){
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$bmi_at_consent)})
    } else if (phenotype == 'vo2max'){
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$vo2max)})
    } else if (phenotype == 'smoking_curr'){
      pheno_tmp$smoking_curr <- ifelse(pheno_tmp$smoking_curr=='no', 0, 1)
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$smoking_curr)})
    } else if (phenotype == 'smoking_ever'){
      pheno_tmp$smoking_ever <- ifelse(pheno_tmp$smoking_ever=='no', 0, 1)
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$smoking_ever)})
    } else if (phenotype == 'activity') {
      pheno_tmp$act <- ifelse(pheno_tmp$intactcurr=='0', 0, 1)
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$act)})
    } else if (phenotype == 'alcohol'){
      pheno_tmp$etoh <- ifelse(pheno_tmp$etoh_curr=='no', 0, 1)
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$etoh)})
    } else if (phenotype == 'menopause'){
      pheno_tmp$mp <- ifelse(pheno_tmp$mpstatrs=='no', 0, 1)
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$mp)})
    } else if (phenotype == 'alcohol_units'){
      corr <- apply(beta, 2, function(x){cor(x, pheno_tmp$etohu_curr_pseudonum)})
    } 
    
    beta <- beta[,abs(corr)>0.1]
    return(beta)
  }
  
  methyl_blood <- filterCorr(methyl_blood, phenotype)
  methyl_buccal <- filterCorr(methyl_buccal, phenotype)
  methyl_cervical <- filterCorr(methyl_cervical, phenotype)
  
  # Save file (interim)
  save(methyl_blood, methyl_buccal, methyl_cervical, file = file.path(out.folder, paste0("methyl_filtered_", phenotype, ".Rdata")))
  
  # Main pheno
  pheno_main <- pheno
  
  # Load in params
  load("~/Dropbox/data/tirolgesund/pheno/methylation_pheno_params.Rdata")
  params_m0 <- pheno[pheno$visitId=='M0',]
  blood_params <- pheno |> dplyr::filter(sampletype=='blood') |> 
    dplyr::slice(match(rownames(methyl_blood), subjectId)) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(gran = hepidish_Mono+hepidish_Neutro+hepidish_Eosino) |> 
    dplyr::ungroup()
  if(!identical(blood_params$subjectId, rownames(methyl_blood))){stop("Check methylation blood rownames and subjectId")}
  buccal_params <- pheno |> dplyr::filter(sampletype=='buccal') |> 
    dplyr::slice(match(rownames(methyl_buccal), subjectId))
  if(!identical(buccal_params$subjectId, rownames(methyl_buccal))){{stop("Check methylation buccal rownames and subjectId")}}
  cervical_params <- pheno |> dplyr::filter(sampletype=='cervical') |> 
    dplyr::slice(match(rownames(methyl_cervical), subjectId))
  if(!identical(cervical_params$subjectId, rownames(methyl_cervical))){{stop("Check methylation cervical rownames and subjectId")}}
  
  # Methylation data models -----------------
  
  lm_blood <- 
    
    if(phenotype == 'age'){
        age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
        apply(methyl_blood, 2, function(var){
          fit <- lm(var ~ age + blood_params$gran) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
          }) |>
          bind_rows(.id = 'var')

    } else if (phenotype == 'bmi'){
      
        age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
        bmi <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$bmi
        apply(methyl_blood, 2, function(var){
          fit <- lm(var ~ bmi + blood_params$gran + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')

    } else if (phenotype == 'vo2max'){
      
        age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
        vo2max <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$vo2max
        apply(methyl_blood, 2, function(var){
          fit <- lm(var ~ vo2max + blood_params$gran + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'smoking_curr') {
      
        age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
        smoking <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$smoking_curr)
        
        apply(methyl_blood, 2, function(var){
          fit <- lm(var ~ smoking + blood_params$gran + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'smoking_ever'){
      
        age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
        smoking <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$smoking_ever)
        
        apply(methyl_blood, 2, function(var){
          fit <- lm(var ~ smoking + blood_params$gran + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')

    } else if (phenotype == 'activity'){
      
      age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
      # smoking <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$smoking_ever)
      activity <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$intactcurr
      activity <- ifelse(activity == '0', 'no', 'yes')
      activity <- as.factor(activity)
      
      apply(methyl_blood, 2, function(var){
        fit <- lm(var ~ activity + blood_params$gran + age) |>
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'alcohol'){
      
      age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
      smoking <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$smoking_curr)
      alcohol <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$etoh_curr)
      
      apply(methyl_blood, 2, function(var){
        tryCatch({
        fit <- lm(var ~ alcohol + blood_params$gran + age + smoking) |>
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
        }, error = function(e){
          fit <- lm(var ~ alcohol + blood_params$gran + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        })
      }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'menopause'){
      
      age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
      mp <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$mpstatrs)
      
      apply(methyl_blood, 2, function(var){
        fit <- lm(var ~ blood_params$gran + mp + age) |>
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'alcohol_units'){
      
      age <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$age_at_consent
      smoking <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$smoking_curr)
      alcohol <- pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$etohu_curr_pseudonum
      
      apply(methyl_blood, 2, function(var){
        tryCatch({
          fit <- lm(var ~ alcohol + blood_params$gran + age + smoking) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }, error = function(e){
          fit <- lm(var ~ alcohol + blood_params$gran + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        })
      }) |>
        bind_rows(.id = 'var')
      
    } 
    
  save(lm_blood, file = file.path(out.folder, paste0("mowas_methylation_blood_", phenotype, ".Rdata")))

  lm_buccal <- 
    
    if(phenotype == 'age'){
        age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
        apply(methyl_buccal, 2, function(var){
          fit <- lm(var ~ age + buccal_params$ic) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')

    } else if (phenotype == 'bmi'){
      
        age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
        bmi <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$bmi
        apply(methyl_buccal, 2, function(var){
          fit <- lm(var ~ bmi + buccal_params$ic + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')

    } else if (phenotype == 'vo2max'){
      
      
        age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
        vo2max <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$vo2max
        apply(methyl_buccal, 2, function(var){
          fit <- lm(var ~ vo2max + buccal_params$ic + age) |> 
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')
        

    } else if (phenotype == 'smoking_curr') {
      
      
        age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
        smoking <- as.factor(pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$smoking_curr)
        
        apply(methyl_buccal, 2, function(var){
          fit <- lm(var ~ smoking + buccal_params$ic + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')

    } else if (phenotype == 'smoking_ever'){
      
        age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
        smoking <- as.factor(pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$smoking_ever)
        
        apply(methyl_buccal, 2, function(var){
          fit <- lm(var ~ smoking + buccal_params$ic + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
        bind_rows(.id = 'var')

    } else if (phenotype == 'activity'){
      
      age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
      # smoking <- as.factor(pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$smoking_ever)
      activity <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$intactcurr
      activity <- ifelse(activity == '0', 'no', 'yes')
      activity <- as.factor(activity)
      
      apply(methyl_buccal, 2, function(var){
        fit <- lm(var ~ activity + buccal_params$ic + age) |>
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'alcohol'){
      
      age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
      smoking <- as.factor(pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$smoking_curr)
      alcohol <- as.factor(pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$etoh_curr)
      
      apply(methyl_buccal, 2, function(var){
        
        tryCatch({
        fit <- lm(var ~ alcohol + buccal_params$ic + age + smoking) |>
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
        
        }, error = function(e){
          fit <- lm(var ~ alcohol + buccal_params$ic + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        })
      }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'menopause'){
      
      age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
      mp <- as.factor(pheno_main[match(rownames(methyl_blood), pheno_main$subjectId),]$mpstatrs)
      
      apply(methyl_buccal, 2, function(var){
        fit <- lm(var ~ buccal_params$ic + mp + age) |>
          broom::tidy() |> 
          dplyr::filter(term != '(Intercept)')
      }) |>
        bind_rows(.id = 'var')
      
    } else if (phenotype == 'alcohol_units'){
      
      age <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$age_at_consent
      smoking <- as.factor(pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$smoking_curr)
      alcohol <- pheno_main[match(rownames(methyl_buccal), pheno_main$subjectId),]$etohu_curr_pseudonum
      
      apply(methyl_buccal, 2, function(var){
        
        tryCatch({
          fit <- lm(var ~ alcohol + buccal_params$ic + age + smoking) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
          
        }, error = function(e){
          fit <- lm(var ~ alcohol + buccal_params$ic + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        })
      }) |>
        bind_rows(.id = 'var')
      
    } 
    
    save(lm_buccal, file = file.path(out.folder, paste0("mowas_methylation_buccal_", phenotype, ".Rdata")))
  
    
    lm_cervical <- 
      
      if(phenotype == 'age'){
          age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
          apply(methyl_cervical, 2, function(var){
            fit <- lm(var ~ age + cervical_params$ic) |> 
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          }) |>
          bind_rows(.id = 'var')

      } else if (phenotype == 'bmi'){
        
        
          age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
          bmi <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$bmi
          apply(methyl_cervical, 2, function(var){
            fit <- lm(var ~ bmi + cervical_params$ic + age) |> 
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          }) |>
          bind_rows(.id = 'var')

      } else if (phenotype == 'vo2max'){
        
        
          age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
          vo2max <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$vo2max
          apply(methyl_cervical, 2, function(var){
            fit <- lm(var ~ vo2max + cervical_params$ic + age) |> 
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          }) |>
          bind_rows(.id = 'var')

      } else if (phenotype == 'smoking_curr'){
        
          age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
          smoking <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$smoking_curr)
          
          apply(methyl_cervical, 2, function(var){
            fit <- lm(var ~ smoking + cervical_params$ic + age) |>
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          }) |>
          bind_rows(.id = 'var')

      } else if (phenotype == 'smoking_ever'){
        
          age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
          smoking <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$smoking_ever)
          
          apply(methyl_cervical, 2, function(var){
            fit <- lm(var ~ smoking + cervical_params$ic + age) |>
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          }) |>
          bind_rows(.id = 'var')

      } else if (phenotype == 'activity'){
        
        age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
        # smoking <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$smoking_ever)
        activity <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$intactcurr
        activity <- ifelse(activity == '0', 'no', 'yes')
        activity <- as.factor(activity)
        
        apply(methyl_cervical, 2, function(var){
          fit <- lm(var ~ activity + cervical_params$ic + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
          bind_rows(.id = 'var')
        
      } else if (phenotype == 'alcohol'){
        
        age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
        smoking <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$smoking_curr)
        alcohol <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$etoh_curr)
        
        apply(methyl_cervical, 2, function(var){
          
          tryCatch({
          fit <- lm(var ~ alcohol + cervical_params$ic + age + smoking) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
          }, error = function(e){
            fit <- lm(var ~ alcohol + cervical_params$ic + age) |>
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          })
        }) |>
          bind_rows(.id = 'var')
        
      } else if (phenotype == 'menopause'){
        
        age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
        mp <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$mpstatrs)
        
        apply(methyl_cervical, 2, function(var){
          fit <- lm(var ~ cervical_params$ic + mp + age) |>
            broom::tidy() |> 
            dplyr::filter(term != '(Intercept)')
        }) |>
          bind_rows(.id = 'var')
        
      } else if (phenotype == 'alcohol_units'){
        
        age <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$age_at_consent
        smoking <- as.factor(pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$smoking_curr)
        alcohol <- pheno_main[match(rownames(methyl_cervical), pheno_main$subjectId),]$etohu_curr_pseudonum
        
        apply(methyl_cervical, 2, function(var){
          
          tryCatch({
            fit <- lm(var ~ alcohol + cervical_params$ic + age + smoking) |>
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          }, error = function(e){
            fit <- lm(var ~ alcohol + cervical_params$ic + age) |>
              broom::tidy() |> 
              dplyr::filter(term != '(Intercept)')
          })
        }) |>
          bind_rows(.id = 'var')
        
      } 
      
      save(lm_cervical, file = file.path(out.folder, paste0("mowas_methylation_cervical_", phenotype, ".Rdata")))
     
      return(list(lm_cervical,
                  lm_buccal,
                  lm_blood,
                  total))
      
}
