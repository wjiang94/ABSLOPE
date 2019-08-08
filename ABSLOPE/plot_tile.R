library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())
library(gbm)


plot_tile = function(Mresults, pm, cr, plot_addr){

  ifelse(!dir.exists(file.path(plot_addr, paste0('pmiss',pm*100,'_corr',cr*10))),
         dir.create(file.path(plot_addr, paste0('pmiss',pm*100,'_corr',cr*10))),
         FALSE)
  
  Sresults1 <- Mresults %>% filter(p.miss == pm & corr == cr)  %>% select(-p.miss, -corr) %>% 
    group_by(method, signallevel, nspr, crit) %>%
    summarise(mean_value = mean(value)) 
  # head(Sresults1)
  #create a new variable from incidence
  Sresults1$mean_factor <- cut(Sresults1$mean_value,
                            breaks = c(-1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,max(Sresults1$mean_value,na.rm=T)),
                            labels=c("<.1",".1 - .2",".2 - .3",".3 - .4",".4 - .5",".5 - .6",".6 - .7",".7 - .8","> .8"))
  Sresults1$mean_factor <- factor(as.character(Sresults1$mean_factor ), levels=rev(levels(Sresults1$mean_factor )))
  # Power
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/a_Power.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'pr') %>% rename(Power = mean_factor) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=Power)) +
    geom_tile() + 
    #scale_fill_gradient(low = "white", high = "#132B43")+
    #change the scale_fill_manual from previous code to below
    scale_fill_manual(values=rev(brewer.pal(9,"YlGnBu")),na.value="grey90")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength") 
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # FDR
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/b_FDR.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'fdr') %>% rename(FDR = mean_factor) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=FDR)) +
    geom_tile() + 
    #scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    facet_grid(method ~ .) +
    scale_fill_manual(values=rev(brewer.pal(9,"YlGnBu")),na.value="grey90")+
    xlab("number of relevant feautures") +
    ylab("signal strength") 
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # bias_beta
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/c_bias_beta.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'bias_beta') %>% rename(bias_beta = mean_factor) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=bias_beta)) +
    geom_tile() + 
    #scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    scale_fill_manual(values=rev(brewer.pal(9,"YlGnBu")),na.value="grey90")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength")
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # bias_sigma
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/d_bias_sigma.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'bias_sigma' & ! (method %in% c('MeanImp+LASSO', 'ncLASSO', 'MeanImp+adaLASSO','MeanImp+SLOPE'))) %>% rename(bias_sigma = mean_factor) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=bias_sigma)) +
    geom_tile() + 
    # scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    scale_fill_manual(values=rev(brewer.pal(9,"YlGnBu")),na.value="grey90")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength")
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # Prediction error (MSE)
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/e_pred_error.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'MSE') %>% rename(pred_error = mean_factor) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=pred_error)) +
    geom_tile() + 
    #scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    scale_fill_manual(values=rev(brewer.pal(9,"YlGnBu")),na.value="grey90")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength")
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
}

plot_tile2 = function(Mresults, pm, cr, plot_addr){
  
  ifelse(!dir.exists(file.path(plot_addr, paste0('pmiss',pm*100,'_corr',cr*10))),
         dir.create(file.path(plot_addr, paste0('pmiss',pm*100,'_corr',cr*10))),
         FALSE)
  
  Sresults1 <- Mresults %>% filter(p.miss == pm & corr == cr)  %>% select(-p.miss, -corr) %>% 
    group_by(method, signallevel, nspr, crit) %>%
    summarise(mean_value = mean(value)) 
  # head(Sresults1)

  # Power
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/a_Power.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'pr') %>% rename(Power = mean_value) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=Power)) +
    geom_tile() + 
    #scale_fill_gradient(low = "white", high = "#132B43")+
    #change the scale_fill_manual from previous code to below
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength")
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # FDR
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/b_FDR.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'fdr') %>% rename(FDR = mean_value) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=FDR)) +
    geom_tile() + 
    #scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    facet_grid(method ~ .) +
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")+
    xlab("number of relevant feautures") +
    ylab("signal strength") 
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # bias_beta
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/c_bias_beta.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'bias_beta') %>% rename(bias_beta = mean_value) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=bias_beta)) +
    geom_tile() + 
    #scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength") 
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # bias_sigma
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/d_bias_sigma.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'bias_sigma' & ! (method %in% c('MeanImp+LASSO', 'ncLASSO', 'MeanImp+adaLASSO','MeanImp+SLOPE'))) %>% rename(bias_sigma = mean_value) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=bias_sigma)) +
    geom_tile() + 
    # scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength") 
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
  
  # Prediction error (MSE)
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'/e_pred_error.pdf'), width = 7, height = 10, paper='special')
  p <- Sresults1 %>% filter(crit == 'MSE') %>% rename(pred_error = mean_value) %>% 
    ggplot(aes(x=nspr, y=signallevel, fill=pred_error)) +
    geom_tile() + 
    #scale_fill_gradient(low = "#56B1F7", high = "#132B43")+
    scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")+
    facet_grid(method ~ .) +
    xlab("number of relevant feautures") +
    ylab("signal strength") 
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', nb.simu = 200'))
  print(p)
  dev.off()
}



plot_line_sd = function(Mresults, pm, cr, sl, plot_addr){
  
  ifelse(!dir.exists(file.path(plot_addr, paste0('pmiss',pm*100,'_corr',cr*10,'_sl',sl))),
         dir.create(file.path(plot_addr, paste0('pmiss',pm*100,'_corr',cr*10,'_sl',sl))),
         FALSE)
  
  Sresults1 <- Mresults %>% filter(p.miss == pm & corr == cr & signallevel ==sl)  %>% select(-p.miss, -corr, -signallevel) %>% 
    group_by(method, nspr, crit) %>%
    summarise(mean_value = mean(value), sd_value = sd(value)) 
  
  # Power
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'_sl',sl,'/a_Power.pdf'), width = 7, height = 7, paper='special')
  p <- Sresults1 %>% filter(crit == 'pr') %>% 
    ggplot(aes(x=nspr, y=mean_value, group=method, color=method, linetype=method)) + #shape=method,
    #geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=0.2, 
    #                  position=position_dodge(0.2)) +
    geom_line() + 
    geom_point()+
    scale_color_brewer(palette="Paired")+
    xlab("number of relevant feautures") +
    ylab("Power") +
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', sl = ',sl,', nb.simu = 200'))+
    theme_minimal()
  print(p)
  dev.off()
  
  # FDR
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'_sl',sl,'/b_FDR.pdf'), width = 7, height = 7, paper='special')
  p <- Sresults1 %>% filter(crit == 'fdr') %>% 
    ggplot(aes(x=nspr, y=mean_value, group=method, color=method, linetype=method)) + #shape=method,
    #geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=0.2, 
    #              position=position_dodge(0.2)) +
    geom_line() + 
    geom_point()+
    scale_color_brewer(palette="Paired")+
    xlab("number of relevant feautures") +
    ylab("FDR") +
    geom_hline(yintercept=0.1, linetype="dotted") +
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', sl = ',sl,', nb.simu = 200'))+
    theme_minimal()
  print(p)
  dev.off()
  
  # bias_beta
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'_sl',sl,'/c_bias_beta.pdf'), width = 7, height = 7, paper='special')
  p <- Sresults1 %>% filter(crit == 'bias_beta') %>% 
    ggplot(aes(x=nspr, y=mean_value, group=method, color=method, linetype=method)) + #shape=method,
    #geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=0.2, 
    #              position=position_dodge(0.2)) +
    geom_line() + 
    geom_point()+
    scale_color_brewer(palette="Paired")+
    xlab("number of relevant feautures") +
    ylab("Bias of beta") +
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', sl = ',sl,', nb.simu = 200'))+
    theme_minimal()
  print(p)
  dev.off()
  
  # bias_sigma
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'_sl',sl,'/d_bias_sigma.pdf'), width = 7, height = 7, paper='special')
  p <- Sresults1 %>% filter(crit == 'bias_sigma' & ! (method %in% c('MeanImp+LASSO', 'ncLASSO', 'MeanImp+adaLASSO','MeanImp+SLOPE'))) %>% 
    ggplot(aes(x=nspr, y=mean_value, group=method, color=method, linetype=method)) + #shape=method,
    #geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=0.2, 
    #              position=position_dodge(0.2)) +
    geom_line() + 
    geom_point()+
    scale_color_brewer(palette="Paired")+
    xlab("number of relevant feautures") +
    ylab("Bias of sigma") +
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', sl = ',sl,', nb.simu = 200'))+
    theme_minimal()
  print(p)
  dev.off()
  
  # Prediction error (MSE)
  pdf(paste0(plot_addr,'pmiss',pm*100,'_corr',cr*10,'_sl',sl,'/e_pred_error.pdf'), width = 7, height = 7, paper='special')
  p <- Sresults1 %>% filter(crit == 'MSE') %>% 
    ggplot(aes(x=nspr, y=mean_value, group=method, color=method, linetype=method)) + #shape=method,
    #geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=0.2, 
    #              position=position_dodge(0.2)) +
    geom_line() + 
    geom_point()+
    scale_color_brewer(palette="Paired")+
    xlab("number of relevant feautures") +
    ylab("Prediction error") +
    #ggtitle(paste0('p.miss = ',pm*100,', corr = ',cr*10,', sl = ',sl,', nb.simu = 200'))+
    theme_minimal()
  print(p)
  dev.off()
}

