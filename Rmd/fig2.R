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

tissue.gsea.reduced = tissue.gsea %>%
  arrange(-NES) %>%
  filter(FDR.q.val < 0.01) %>%
  filter(SIZE < 500) %>%
  mutate(`Log Q value`=ifelse(FDR.q.val==0, -5, log10(FDR.q.val)))


gsea = ggplot(tissue.gsea.reduced, aes(factor(NAME, levels=rev(unique(NAME))), NES, fill=`Log Q value`)) +
  geom_bar(stat="identity") + 
  coord_flip() + 
  ylab("Normalized Enrichment Score") + 
  xlab("") + 
  theme_bw() + 
  ggtitle("Gene Set Enrichment Analysis")

grap2.data = load("~/Projects/methyl_age/data/model_data/grap2.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

grap2 = ggplot(grap2.data, aes(tissue_state, M, fill=tissue_state)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-7, 5) +
  xlab("") +
  ylab("M value") +
  ggtitle("GRAP2 promoter")+
  guides(color=F)


rin2.data = load("~/Projects/methyl_age/data/model_data/rin2.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

rin2 = ggplot(rin2.data, aes(tissue_state, M, fill=tissue_state)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-5, 5) +
  ylab("M value") +
  xlab("") +
  ggtitle("RIN2 promoter") +
  guides(color=F)

tiff("./figures/fig2.tiff", width=22, height=22, units="in", res=150)
grid.arrange(volcano, gsea, grap2, rin2, ncol=2)
dev.off()