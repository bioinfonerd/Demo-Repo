#Purpose: To determine whether Binding Density is correlated to anything related to the imaging data

#Outline:
#Output: Scatterplot with R2 correlation plotted per pooled vs not group
#Filter any samples that were 'test' samples


# Libraries ---------------------------------------------------------------
Packages <- c('tidyverse','ggplot2','openxlsx','RColorBrewer',
              'rstatix','ggpubr','ggiraph','ggiraphExtra','plyr','jtools',
              'ggstance','broom.mixed','ggpmisc')
lapply(Packages, library, character.only = TRUE)

# Function ----------------------------------------------------------------
read_all_sheets = function(xlsxFile, ...) {
  sheet_names = openxlsx::getSheetNames(xlsxFile)
  sheet_list = as.list(rep(NA, length(sheet_names)))
  names(sheet_list) = sheet_names
  for (sn in sheet_names) {
    sheet_list[[sn]] = openxlsx::read.xlsx(xlsxFile, sheet=sn, ...)
  }
  return(sheet_list)
}

# Data --------------------------------------------------------------------

#metadata
metadata <- read.csv2('./data/2021.03.22 Nanostring_Metadata.csv',
                      sep=",",fileEncoding="latin1")

#image & binding information
image_bd <- read_all_sheets('./data/binding-density-vs-spheroids.xlsx')


# Images ------------------------------------------------------------------


# Trying Different Linear Models for Binding Density
fit  <- lm(BindingDensity~Total.Area,data=image_bd$`binding-density-vs-spheroids`)
fit2 <- lm(BindingDensity~Total.Area+Spheroid.Count+Live.Area+percent.Dead,
           data = image_bd$`binding-density-vs-spheroids`)
fit3 <- lm(BindingDensity~Total.Area+Spheroid.Count+percent.Dead,
           data = image_bd$`binding-density-vs-spheroids`)
fit4  <- lm(BindingDensity~Spheroid.Count,data=image_bd$`binding-density-vs-spheroids`)
fit5  <- lm(BindingDensity~Live.Area,data=image_bd$`binding-density-vs-spheroids`)
fit10  <- lm(BindingDensity~percent.Dead,data=image_bd$`binding-density-vs-spheroids`)


summary(fit)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)

pdf('./output/LinearModels_.pdf')
plot_summs(fit, fit2,fit3,fit4,fit5,fit10, scale = TRUE, plot.distributions = TRUE,
           model.names=c("Total.Area","All",'-Live.Area','Spheriod','Live.Area','percent.Dead'))
dev.off()

#
ggPredict(fit,se=TRUE,interactive=TRUE)
ggPredict(fit4,se=TRUE,interactive=TRUE)
ggPredict(fit5,se=TRUE,interactive=TRUE)


#try log scale (not useful)
image_bd$`binding-density-vs-spheroids`$log.Spheroid.Count=log(1+image_bd$`binding-density-vs-spheroids`$Spheroid.Count)
image_bd$`binding-density-vs-spheroids`$log.Total.Area=log(1+image_bd$`binding-density-vs-spheroids`$Total.Area)
image_bd$`binding-density-vs-spheroids`$log.Live.Area=log(1+image_bd$`binding-density-vs-spheroids`$Live.Area)
image_bd$`binding-density-vs-spheroids`$log.percent.Dead=log(1+image_bd$`binding-density-vs-spheroids`$percent.Dead)


fit6  <- lm(BindingDensity~log.Spheroid.Count,data=image_bd$`binding-density-vs-spheroids`)
summary(fit6)

fit7  <- lm(BindingDensity~log.Total.Area,data=image_bd$`binding-density-vs-spheroids`)
summary(fit7)

fit8  <- lm(BindingDensity~log.Live.Area,data=image_bd$`binding-density-vs-spheroids`)
summary(fit8)

fit9  <- lm(BindingDensity~log.percent.Dead,data=image_bd$`binding-density-vs-spheroids`)
summary(fit9)

pdf('./output/LinearModels_log_.pdf')
plot_summs(fit, fit2,fit5,fit6,fit7,fit8,fit9, scale = TRUE, plot.distributions = TRUE,
           model.names=c("Total.Area","All",'Live.Area','log.Spheroid.Count',
                         'log.Total.Area','log.Live.Area','log.percent.Dead'))


dev.off()


# Ggplot Image ------------------------------------------------------------

#including Kevin's way to calculate pooled BindingDensity~Total.Area
fit  <- lm(BindingDensity~Total.Area,data=image_bd$`binding-density-vs-spheroids`)
pdf('./output/Total_AreaxBindingDensity_pooled.pdf')
ggplot(image_bd$`binding-density-vs-spheroids`,
       aes(x = Total.Area, y = BindingDensity)) + 
       geom_point(aes(color = Pooled)) +
       theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
       geom_smooth(method = "lm") +
       labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                           "Intercept =",signif(fit$coef[[1]],5 ),
                           " Slope =",signif(fit$coef[[2]], 5),
                           " P =",signif(summary(fit)$coef[2,4], 5)))
dev.off()

#including Kevin's way to calculate pooled BindingDensity~Live.Area
fit  <- lm(BindingDensity~Live.Area,data=image_bd$`binding-density-vs-spheroids`)
pdf('./output/Live.AreaxBindingDensity_pooled.pdf')
ggplot(image_bd$`binding-density-vs-spheroids`,
       aes(x = Live.Area, y = BindingDensity)) + 
  geom_point(aes(color = Pooled)) +
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  geom_smooth(method = "lm") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
dev.off()


#remove pooled
data <- image_bd$`binding-density-vs-spheroids` %>% filter(is.na(Pooled))
fit  <- lm(BindingDensity~Total.Area,data=data)
pdf('./output/Total.AreaxBindingDensity_NOpooled.pdf')
ggplot(data, aes(x = Live.Area, y = BindingDensity)) + 
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
dev.off()


#remove pooled, percent.Dead
data <- image_bd$`binding-density-vs-spheroids` %>% filter(is.na(Pooled))
fit  <- lm(BindingDensity~percent.Dead,data=data)
pdf('./output/percent.DeadxBindingDensity_NOpooled.pdf')
ggplot(data, aes(x = percent.Dead, y = BindingDensity)) + 
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
dev.off()


#remove pooled, Spheroid.Count
data <- image_bd$`binding-density-vs-spheroids` %>% filter(is.na(Pooled))
fit  <- lm(BindingDensity~Spheroid.Count,data=data)
pdf('./output/Spheroid.CountxBindingDensity_NOpooled.pdf')
ggplot(data, aes(x = Spheroid.Count, y = BindingDensity)) + 
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
dev.off()


#remove pooled
data <- image_bd$`binding-density-vs-spheroids` %>% filter(is.na(Pooled))
fit  <- lm(BindingDensity~Total.Area,data=data)
pdf('./output/Total.AreaxBindingDensity_NOpooled.pdf')
ggplot(data, aes(x = Live.Area, y = BindingDensity)) + 
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
dev.off()

# Plot binding density of individual samples vs pooled
# Does 

pdf(file="cell_type.pdf")
ggplot(cell_type, aes(x = Cell_Type, y = Value, fill = pooled)) +    
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_flip()+
  ggtitle("Cell Type Over All Diseases") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()