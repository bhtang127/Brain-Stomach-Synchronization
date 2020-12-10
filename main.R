source('utility.R')
library(grid)
library(gridExtra)
library(tidyverse)


#############################################################################
## Global parameter to control whether use amplitude weighted PLV
USE_AW_PLV = FALSE


#############################################################################
## Calculate mismatched PLV
plv_mri_egg = data.frame(session_i = NA, session_j = NA, 
                         run_i = NA, run_j = NA, 
                         component=NA, PLV = NA)
for(comp in 1:18){
  for(i in 1:19) { for(j in 1:19) {
    for(r1 in 1:2) { for(r2 in 1:2) {
      sig1 = readMRIEGG(i, r1) # EGG signal
      sig2 = readMRIEGG(j, r2, component=comp) # fMRI signal
      
      if(USE_AW_PLV){
        m1 = Mod(sig1); m2 = Mod(sig2)
        weight = m1 * m2 / sum(m1 *  m2)
      }
      else 
        weight = 1 / length(sig1)
      
      plv = Mod( sum( weight * exp(1i * (Arg(sig1) - Arg(sig2))) ) )
      plv_mri_egg = rbind(plv_mri_egg, c(i,j,r1,r2,comp,plv))
    }}
  }}
}
plv_mri_egg = na.omit(plv_mri_egg)
## Save PLV data
write.csv(plv_mri_egg, file = 'result/plv_mriegg_ds.csv', row.names = FALSE)


#############################################################################
## Get Summary table
mean_plv = c(); std_plv = c()
mean_mis_plv = c(); std_mis_plv = c()
raw_p = c(); fdr_p = c(); bf_p = c()
for(comp in 1:18){
  match = (
    plv_mri_egg %>% 
      dplyr::filter(component==comp & 
                    session_i==session_j & 
                    run_i==run_j)
    )$PLV
  mismatch = (
    plv_mri_egg %>% 
      dplyr::filter(component==comp & 
                      session_i!=session_j)
  )$PLV
  mean_plv[comp] = mean(match)
  std_plv[comp] = sqrt(var(match))
  mean_mis_plv[comp] = mean(mismatch)
  std_mis_plv[comp] = sqrt(var(mismatch))
  raw_p[comp] = wilcox.test(match, mismatch, alternative="greater")$p.value
}
fdr_p = p.adjust(raw_p, method="fdr")
bf_p = p.adjust(raw_p, method="bonferroni")
#### Save summary table
result_table = data.frame(
    id=1:18, region_name=region_name, 
    PLV_mean_empirical = mean_plv, 
    PLV_std_empirical = std_plv,
    PLV_mean_surrogate = mean_mis_plv,
    PLV_std_surrogate = std_mis_plv,
    raw_p_value = raw_p, 
    fdr_p_value = fdr_p,
    bonferroni_p_value = bf_p
)
write.csv(result_table, file = 'result/summary.csv', row.names = FALSE)


#############################################################################
## Function to plot PLV densities for certain region
## will show FDR corrected p-value
## Input: num: will plot summaries for region with name
##        region_name[region_rank[num]]
plot.mismatch = function(num){
  comp = region_rank[num]
  rname = region_name[comp]
  data = plv_mri_egg %>% 
    dplyr::filter(!(session_i == session_j & run_i != run_j) &
                    component == comp) %>%
    mutate(matched = session_i == session_j) 

  group = data$matched
  group[group == TRUE] = "Matched"
  group[group == FALSE] = "Surrogate"
  data$group = factor(group, levels = c("Surrogate","Matched"))
  p1 = data %>%
    ggplot(aes(x=PLV)) +
    geom_density(alpha = 0.15, aes(color=group, fill=group)) +
    xlim(c(0, max(plv_mri_egg$PLV))) + ylim(c(0, 6.3)) +
    ggtitle(paste("RSN:", rname, "\nadjusted p-value:", 
                  round(fdr_p[comp], 5))) +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black", size = 0.3), 
          legend.title=element_blank(),
          legend.position = c(0.8, 0.85))
  ggsave(paste0('image/RSN_',rname,'.png'), plot=p1, width = 6, height = 4)
  p1
}

## Combine all plots in a grid
gs = list()
for(comp in 1:18){
    out = plot.mismatch(comp)
    gs[[comp]] = out +
                 theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       legend.position = "none")
    if(comp == 18){
      gs[[comp]] = out +
                   theme(axis.title.x = element_blank(),
                         axis.title.y = element_blank())
    }
}

## Save images
png('image/combined_density_fdr.png', width = 600, height = 900)
pg = grid.arrange(grobs=gs, ncol=3, bottom="PLV", left = "density")
dev.off()
