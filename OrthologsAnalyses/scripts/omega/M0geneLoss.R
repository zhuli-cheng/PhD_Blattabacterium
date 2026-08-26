source("omega/OmegaMain.R")

#remove the gene (amiD) with omega > 1
M0.omega.COG <- M0.omega.COG %>%
  filter(omega < 1)

############# gene loss ~ statistics #############
#define a function
gene.loss <- function(statistic) {
  # Build formula and fit model
  formula_str <- paste0("cbind(number_family_with_loss, 7 - number_family_with_loss) ~ ", 
                        statistic, " + (1 | curated_COG_category)")
  model_formula <- as.formula(formula_str)  
  model.statistic <- glmer(model_formula, data = M0.omega.COG, family = binomial)
  
  newdata <- M0.omega.COG %>%
    group_by(curated_COG_category) %>%
    summarise(x = seq(min(.data[[statistic]], na.rm = TRUE), max(.data[[statistic]], na.rm = TRUE), length.out = 100), .groups = "drop") %>%
    tidyr::unnest(cols = c(x))
  
  colnames(newdata)[colnames(newdata) == "x"] <- statistic
  
  newdata$pred_prob <- predict(model.statistic, newdata = newdata, type = "response", re.form = NULL)
  
  newdata$pred_count <- newdata$pred_prob * 7
  
  pred <- predict(model.statistic, newdata = newdata, type = "link", se.fit = TRUE, re.form = NULL)
  
  newdata$fit_link <- pred$fit
  newdata$se_link  <- pred$se.fit
  
  newdata$upr <- plogis(newdata$fit_link + 1.96 * newdata$se_link)
  newdata$lwr <- plogis(newdata$fit_link - 1.96 * newdata$se_link)
  
  newdata$upr_count <- newdata$upr * 7
  newdata$lwr_count <- newdata$lwr * 7
  
  newdata$xvar <- newdata[[statistic]]
  
  p <- ggplot(M0.omega.COG, aes(x = .data[[statistic]], y = number_family_with_loss)) +
    geom_point() +
    geom_ribbon(data = newdata, aes(x = xvar, ymin = lwr_count, ymax = upr_count), alpha = 0.2, inherit.aes = FALSE) +
    geom_line(data = newdata, aes(x = xvar, y = pred_count), linewidth = 1.2, inherit.aes = FALSE) +
    facet_wrap(~curated_COG_category, labeller = label_wrap_gen(width = 30)) +
    theme_bw(base_size = 18)
  
  ggsave(paste0("../244/output/geneLoss/", statistic, ".pdf"), plot = p, width = 18, height = 12)
  
  #open the PDF
  system(paste0("open ../244/output/geneLoss/", statistic, ".pdf"))
  
  # Return model results and plot object
  return(model.statistic)
}


#GC
model.GC <- gene.loss("mean_GC_in_orthogroup")
summary(model.GC)
r.squaredGLMM(model.GC)

#dN
model.dN <- gene.loss("dN_standardized")
summary(model.dN)
r.squaredGLMM(model.dN)

#dS
model.dS <- gene.loss("dS_standardized")
summary(model.dS)
r.squaredGLMM(model.dS)

#omega
model.omega <- gene.loss("omega")
summary(model.omega)
r.squaredGLMM(model.omega)

#position
model.position <- gene.loss("mean_position_in_orthogroup")
summary(model.position)
r.squaredGLMM(model.position)

#only COG categories
model.COG <- glmer(cbind(number_family_with_loss, 7 - number_family_with_loss) ~ 
                      (1 | curated_COG_category), data = M0.omega.COG, family = binomial)
summary(model.COG)
r.squaredGLMM(model.COG)

#save model outputs
tab_model(model.GC, model.dN, model.dS, model.omega, model.position, model.COG, 
          transform = NULL, 
          show.re.var = TRUE,
          linebreak = TRUE,
          dv.labels = c("Model (GC and COG)", "Model (Standardized <i>d</i><sub>N</sub> and COG)", "Model (Standardized <i>d</i><sub>S</sub> and COG)", "Model (<i>d</i><sub>N</sub>/<i>d</i><sub>S</sub> and COG)", "Model (chromosome position and COG)", "Model (only COG)"),
          pred.labels = c("(Intercept)", "GC content", "Standardized <i>d</i><sub>N</sub>", "Standardized <i>d</i><sub>S</sub>", "<i>d</i><sub>N</sub>/<i>d</i><sub>S</sub>", "Chromosomal position"),
          file = "../244/output/geneLoss/gene_loss_models_summary.html")
