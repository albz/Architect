#
#--- libraries ---#
library(ggplot2)
library(latex2exp)

#--- path ---#
path <- getwd()
#--- inputs ---#
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("\nAt least one argument must be supplied (input file).\n", call.=FALSE)
}
bunch_number <- as.integer(args[1])

#--- import data ---#
integrated_data <- as.data.frame( as.matrix( read.table(file.path(path,'out','integrated_diagnostics',sprintf('bunch_integrated_quantity_%d_dcut.dat',bunch_number))) ) )

columns_name <- c('z',"<X>","<Y>","<Z>","<Px>","<Py>","<Pz>","rmsX","rmsY","<rmsZ>","<rmsPx>","<rmsPy>","<rmsPz>",'EmittX','EmittY',"Gamma",'DGammasuGamma',"cov<xPx>","cov<yPy>","cov<zPz>")#,"n_over_ne")
colnames(integrated_data) <- columns_name

#--- physical conversions ---#
integrated_data$z             <- integrated_data$z             / 1e4 #now in [cm]
integrated_data$DGammasuGamma <- integrated_data$DGammasuGamma * 1e2 #now in [per cent]

#--- plot ---#
plotted <- ggplot(integrated_data, aes(x=z)) +
  geom_line(aes(y=rmsX,          colour = 's_x'),          size=.5) +
  geom_line(aes(y=rmsY,          colour = 's_y'),          size=.5) +
  geom_line(aes(y=EmittX,        colour = 'emit_x'),       size=.5) +
  geom_line(aes(y=EmittY,        colour = 'emit_y'),       size=.5) +
  geom_line(aes(y=DGammasuGamma, colour = 'sigma gamma/gamma'), size=.5) +
  geom_line(aes(y=Gamma*0.5/1e3, colour = 'gamma'),        size=.5) +
  theme(legend.position='bottom',legend.title=element_blank(), axis.text=element_text(size=5), axis.title = element_text(size = 6), legend.text=element_text(size=6)) +
  ylab(TeX('$\\mu m$')) +
  xlab('z [cm]') +
  #ylim(0, 2) +
  scale_color_discrete( labels=lapply(c('$\\epsilon_x$','$\\epsilon_y$','MeV/1000','$\\sigma_x','$\\sigma_y','$\\sigma_{\\gamma}/\\gamma$'), TeX) )

ggsave("integrated.pdf", width = 3.25, height = 3, units = 'in', useDingbats=FALSE)
