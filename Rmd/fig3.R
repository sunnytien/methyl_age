library("ggplot2")
library("ggthemes")

load("~/Projects/methyl_age/data/afr.Rdata")

afr %<>% 
  mutate(logP=-log10(pnorm(-abs(`t value`)))) %>%
  mutate(label=ifelse(logP>60 | Estimate > 0.5, gene, ""))

volcano = ggplot(afr, aes(Estimate, logP)) + 
  geom_point() + 
  geom_text(aes(label=label), vjust=2, size=3) +
  ylab("Log10 P-value") +
  xlab("Methylation Difference Between Africans and Europeans") + 
  xlim(-max(abs(afr$Estimate)), max(abs(afr$Estimate))) + 
  theme_bw() +
  ggtitle("Association between Ancestry and Methylation for All Promoters")

fmod.data = load("~/Projects/methyl_age/data/model_data/fmod.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

fmod = ggplot(fmod.data, aes(predicted.ancestry, M, color=predicted.ancestry)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-7, 5) +
  xlab("") +
  ylab("M value") +
  ggtitle("FMOD promoter")

or2.data = load("~/Projects/methyl_age/data/model_data/OR2L13.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

or2 = ggplot(or2.data, aes(predicted.ancestry, M, color=predicted.ancestry)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-7, 5) +
  xlab("") +
  ylab("M value") +
  ggtitle("OR2L13 promoter")

tiff("./figures/fig3.tiff", width=33, height=11, units="in", res=150)
grid.arrange(volcano, fmod, or2, ncol=2)
dev.off()