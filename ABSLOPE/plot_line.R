library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())
library(gbm)


plot_line_signallevel = function(Mresults, pm, cr, plot_addr){
  
  ifelse(!dir.exists(file.path(plot_addr, paste0('ABSLOPE_pmiss',pm*100,'_corr',cr*10))),
         dir.create(file.path(plot_addr, paste0('ABSLOPE_pmiss',pm*100,'_corr',cr*10))),
         FALSE)
  
  Sresults2 <- Mresults %>% filter(method == 'ABSLOPE') %>% 
    filter(p.miss == pm & corr == cr) %>% dplyr::select(-p.miss, -corr) %>% 
    group_by(method, signallevel, nspr, crit) %>%
    summarise(mean_value = mean(value)) 
  #head(Sresults2)
  
  # Power
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_corr',cr*10,'/a_Power.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'pr') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(signallevel), group=factor(signallevel))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Power") +
    scale_colour_discrete(name  ="Signal strength") 
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  
  
  # FDR
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_corr',cr*10,'/b_FDR.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'fdr') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(signallevel), group=factor(signallevel))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("FDR") +
    scale_colour_discrete(name  ="Signal strength") +
    geom_hline(yintercept=0.1, linetype="dotted")
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  # bias_beta
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_corr',cr*10,'/c_bias_beta.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'bias_beta') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(signallevel), group=factor(signallevel))) +
    geom_line() + 
    geom_point() +
    #ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab(bquote("Relative MSE of" ~ beta )) +
    scale_colour_discrete(name  ="Signal strength")
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  # bias_sigma
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_corr',cr*10,'/d_bias_sigma.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'bias_sigma') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(signallevel), group=factor(signallevel))) +
    geom_line() + 
    geom_point() +
    #ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Bias of sigma") +
    scale_colour_discrete(name  ="Signal strength") 
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  # Prediction error (MSE)
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_corr',cr*10,'/e_pred_error.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'MSE') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(signallevel), group=factor(signallevel))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Prediction error") +
    scale_colour_discrete(name  ="Signal strength")
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
}




plot_line_percentna = function(Mresults, sl, cr, plot_addr){
  
  ifelse(!dir.exists(file.path(plot_addr, paste0('ABSLOPE_signallevel',sl,'_corr',cr*10))),
         dir.create(file.path(plot_addr, paste0('ABSLOPE_signallevel',sl,'_corr',cr*10))),
         FALSE)
  
  Sresults3 <- Mresults %>% filter(method == 'ABSLOPE') %>% 
    filter(signallevel == sl & corr == cr) %>% dplyr::select(-signallevel, -corr) %>% 
    group_by(method, p.miss, nspr, crit) %>%
    summarise(mean_value = mean(value)) 
  #head(Sresults3)
  
  # Power
  pdf(paste0(plot_addr,'ABSLOPE_signallevel',sl,'_corr',cr*10,'/a_Power.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults3 %>% filter(crit == 'pr') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(p.miss), group=factor(p.miss))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Power") +
    scale_colour_discrete(name  ="Percentage NA") 
    #ggtitle(paste0('ABSLOPE: signal strength=',sl,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  
  
  # FDR
  pdf(paste0(plot_addr,'ABSLOPE_signallevel',sl,'_corr',cr*10,'/b_FDR.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults3 %>% filter(crit == 'fdr') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(p.miss), group=factor(p.miss))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("FDR") +
    scale_colour_discrete(name  ="Percentage NA") +
    geom_hline(yintercept=0.1, linetype="dotted")
    #ggtitle(paste0('ABSLOPE: signal strength=',sl,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  # bias_beta
  pdf(paste0(plot_addr,'ABSLOPE_signallevel',sl,'_corr',cr*10,'/c_bias_beta.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults3 %>% filter(crit == 'bias_beta') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(p.miss), group=factor(p.miss))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab(bquote("Relative MSE of" ~ beta )) +
    scale_colour_discrete(name  ="Percentage NA")
    #ggtitle(paste0('ABSLOPE: signal strength=',sl,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  # bias_sigma
  pdf(paste0(plot_addr,'ABSLOPE_signallevel',sl,'_corr',cr*10,'/d_bias_sigma.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults3 %>% filter(crit == 'bias_sigma') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(p.miss), group=factor(p.miss))) +
    geom_line() + 
    geom_point() +
    #ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Bias of sigma") +
    scale_colour_discrete(name  ="Percentage NA") 
    #ggtitle(paste0('ABSLOPE: signal strength=',sl,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
  
  # Prediction error (MSE)
  pdf(paste0(plot_addr,'ABSLOPE_signallevel',sl,'_corr',cr*10,'/e_pred_error.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults3 %>% filter(crit == 'MSE') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(p.miss), group=factor(p.miss))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Prediction error") +
    scale_colour_discrete(name  ="Percentage NA") 
    #ggtitle(paste0('ABSLOPE: signal strength=',sl,', corr=',cr*10,', nb.simu=40'))
  print(p)
  dev.off()
}

plot_line_corr = function(Mresults, pm, sl, plot_addr){
  
  ifelse(!dir.exists(file.path(plot_addr, paste0('ABSLOPE_pmiss',pm*100,'_signallevel',sl))),
         dir.create(file.path(plot_addr, paste0('ABSLOPE_pmiss',pm*100,'_signallevel',sl))),
         FALSE)
  
  Sresults2 <- Mresults %>% filter(method == 'ABSLOPE') %>% 
    filter(p.miss == pm & signallevel == sl) %>% dplyr::select(-p.miss, -signallevel) %>% 
    group_by(method, corr, nspr, crit) %>%
    summarise(mean_value = mean(value)) 
  #head(Sresults2)
  
  # Power
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_signallevel',sl,'/a_Power.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'pr') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(corr), group=factor(corr))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Power") +
    scale_colour_discrete(name  ="Correlation") 
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', signallevel=',sl,', nb.simu=40'))
  print(p)
  dev.off()
  
  
  
  # FDR
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_signallevel',sl,'/b_FDR.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'fdr') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(corr), group=factor(corr))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("FDR") +
    scale_colour_discrete(name  ="Correlation") +
    geom_hline(yintercept=0.1, linetype="dotted")
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', signallevel=',sl,', nb.simu=40'))
  print(p)
  dev.off()
  
  # bias_beta
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_signallevel',sl,'/c_bias_beta.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'bias_beta') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(corr), group=factor(corr))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab(bquote("Relative MSE of" ~ beta ))
    scale_colour_discrete(name  ="Correlation")
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', signallevel=',sl,', nb.simu=40'))
  print(p)
  dev.off()
  
  # bias_sigma
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_signallevel',sl,'/d_bias_sigma.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'bias_sigma') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(corr), group=factor(corr))) +
    geom_line() + 
    geom_point() +
    #ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Bias of sigma") +
    scale_colour_discrete(name  ="Correlation") 
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', signallevel=',sl,', nb.simu=40'))
  print(p)
  dev.off()
  
  # Prediction error (MSE)
  pdf(paste0(plot_addr,'ABSLOPE_pmiss',pm*100,'_signallevel',sl,'/e_pred_error.pdf'), width = 6, height = 4,  paper='special')
  p <- Sresults2 %>% filter(crit == 'MSE') %>% rename(obj_value = mean_value) %>% 
    ggplot(aes(x=nspr, y=obj_value, colour=factor(corr), group=factor(corr))) +
    geom_line() + 
    geom_point() +
    ylim(0, 1) +
    xlab("Number of relevant features") +
    ylab("Prediction error") +
    scale_colour_discrete(name  ="Correlation") 
    #ggtitle(paste0('ABSLOPE: p.miss=',pm*100,', signallevel=',sl,', nb.simu=40'))
  print(p)
  dev.off()
}


