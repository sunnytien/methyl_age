library("ggplot2")
library("ggthemes")

load("~/Projects/methyl_age/data/afr.Rdata")
load("~/Projects/methyl_age/data/afr.gsea.Rdata")

afr %<>% 
  mutate(logP=-log10(pnorm(-abs(`t value`)))) %>%
  mutate(label=ifelse(logP>60 | Estimate > 0.5, gene, ""))

afr.gsea.reduced = afr.gsea %>%
  arrange(-NES) %>%
  filter(FDR.q.val == 0) %>%
  filter(SIZE < 300) %>%
  filter(SIZE > 75) %>%
  mutate(`Log Q value`=ifelse(FDR.q.val==0, -5, log10(FDR.q.val)))


gsea = ggplot(afr.gsea.reduced, aes(factor(NAME, levels=rev(unique(NAME))), NES, fill=`Log Q value`)) +
  geom_bar(stat="identity") + 
  coord_flip() + 
  ylab("Normalized Enrichment Score") + 
  xlab("") + 
  theme_bw() + 
  ggtitle("Gene Set Enrichment Analysis") + 
  theme(plot.title=element_text(size=rel(2)))+ 
  theme(axis.title.y=element_text(size=rel(1.3)))

volcano = ggplot(afr, aes(Estimate, logP)) + 
  geom_point() + 
  geom_text(aes(label=label), vjust=2, size=3) +
  ylab("Log10 P-value") +
  xlab("Methylation Difference Between Africans and Europeans") + 
  xlim(-max(abs(afr$Estimate)), max(abs(afr$Estimate))) + 
  theme_bw() +
  ggtitle("Association between Ancestry and Methylation for All Promoters") + 
  theme(plot.title=element_text(size=rel(2)))+ 
  theme(axis.title.y=element_text(size=rel(1.3)))

fmod.data = load("~/Projects/methyl_age/data/model_data/fmod.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

fmod = ggplot(fmod.data, aes(predicted.ancestry, M, fill=predicted.ancestry)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-7, 5) +
  xlab("") +
  ylab("M value") +
  ggtitle("FMOD promoter") +
  guides(color=F) +
  theme(plot.title=element_text(size=rel(2))) + 
  theme(axis.title.y=element_text(size=rel(1.3)))

or2.data = load("~/Projects/methyl_age/data/model_data/OR2L13.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

or2 = ggplot(or2.data, aes(predicted.ancestry, M, fill=predicted.ancestry)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(. ~ Probe) +
  theme_bw() +
  ylim(-7, 5) +
  xlab("") +
  ylab("M value") +
  ggtitle("OR2L13 promoter") +
  guides(color=F) + 
  theme(plot.title=element_text(size=rel(2))) + 
  theme(axis.title.y=element_text(size=rel(1.3)))

tiff("./figures/fig3.tiff", width=22, height=22, units="in", res=150)
grid.arrange(volcano, gsea, fmod, or2, ncol=2)
dev.off()