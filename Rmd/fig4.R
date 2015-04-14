library("ggplot2")
library("dplyr")
library("tidyr")
library("gsl")
library("gridExtra")

load("./data/ensemble.Rdata")

age.antitrans = function(x) 21 * lambert_W0(exp(21*x + 1)^(1/21)) - 1

y3 = data3$age.normed

cols = c(ppr="Regression",
         svmLinear="Kernel",
         pls="Regression",
         rvmLinear="Kernel",
         gbm="Boosting",
         rvmRadial="Kernel",
         pcr="Regression",
         svmRadialCost="Kernel",
         gaussprRadial="Kernel",
         rf="Boosting",
         ensemble="Ensemble",
         cubist="Boosting")

perf = cor(cbind(z3, y3.pred), y3) %>%
  as.data.frame %>%
  mutate(method=row.names(.)) %>%
  mutate(cor=V1) %>%
  mutate(family=cols[method]) %>%
  arrange(family, cor) %>%
  mutate(method=factor(method, levels=method))


performance = ggplot(perf, aes(method, V1, color=family)) + 
  geom_point(size=10) +
  theme_bw() +
  ggtitle("Performance on Testing Set")+
  coord_flip()+
  ylab("Correlation between observed and predicted age") +
  xlab("Prediction Method")

error = cbind(z3, y3.pred) %>%
  apply(2, function(x) abs(x-y3)) %>%
  as.data.frame %>%
  mutate(sample=1:nrow(.)) %>%
  gather(method, error, -sample) %>%
  mutate(method=factor(method, levels(perf$method))) %>%
  mutate(family=cols[as.character(method)]) 

performance2 = ggplot(error, aes(method, error, fill=family)) +
  geom_boxplot(outlier.colour=NA) +
  theme_bw() +
  coord_flip() +
  ylim(0, 2.5)+
  ggtitle("Absolute Error") +
  ylab("Error")

scatter.data = cbind(y3.pred, y3) %>%
  as.data.frame %>%
  mutate(age.pred=age.antitrans(ensemble)) %>%
  mutate(age=age.antitrans(y3))

scatter = ggplot(scatter.data, aes(age, age.pred)) +
  geom_point() +
  theme_bw() +
  xlab("Age") +
  ylab("Predicted Age") +
  ggtitle("Observed vs Predicted Age: Testing Set")

tiff("./figures/fig4.tiff", width=33, height=11, units="in", res=150)
grid.arrange(performance, performance2, scatter, ncol=3)
dev.off()