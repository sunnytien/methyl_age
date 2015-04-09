library("ggplot2")
library("ggthemes")

load("~/Projects/methyl_age/data/tissue.Rdata")
load("~/Projects/methyl_age/data/tissue.gsea.Rdata")

tissue %<>% 
  mutate(logP=-log10(pnorm(-abs(`t value`)))) %>%
  mutate(label=ifelse(logP>26, gene, ""))

volcano = ggplot(tissue, aes(Estimate, logP)) + 
  geom_point() + 
  geom_text(aes(label=label), vjust=2, size=3) +
  ylab("Log10 P-value") +
  xlab("Methylation Difference Between Liquid and Solid Tissues") + 
  xlim(-max(abs(tissue$Estimate)), max(abs(tissue$Estimate))) + 
  theme_bw() +
  ggtitle("Association between Tisse Category and Methylation for All Promoters")

grap2.data = load("~/Projects/methyl_age/data/model_data/grap2.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

grap2 = ggplot(grap2.data, aes(tissue_state, M, color=tissue_state)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-7, 5) +
  xlab("") +
  ylab("M value") +
  ggtitle("GRAP2 promoter")


rin2.data = load("~/Projects/methyl_age/data/model_data/rin2.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

rin2 = ggplot(rin2.data, aes(tissue_state, M, color=tissue_state)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-5, 5) +
  ylab("M value") +
  xlab("") +
  ggtitle("RIN2 promoter")


tiff("./figures/fig2.tiff", width=33, height=11, units="in", res=150)
grid.arrange(volcano, grap2, rin2, ncol=3)
dev.off()