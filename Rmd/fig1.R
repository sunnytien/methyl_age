library("ggplot2")
library("ggthemes")
library("gridExtra")
library("dplyr")
library("magrittr")

load("~/Projects/methyl_age/data/age.Rdata")
load("~/Projects/methyl_age/data/age.gsea.Rdata")


age %<>% 
  mutate(logP=-log10(pnorm(-abs(`t value`)))) %>%
  mutate(label=ifelse(logP>150, gene, "")) 

age.gsea.reduced = age.gsea %>%
  arrange(-NES) %>%
  filter(FDR.q.val < 0.01) %>%
  filter(SIZE < 300) %>%
  filter(SIZE > 75) %>%
  mutate(`Log Q value`=ifelse(FDR.q.val==0, -5, log10(FDR.q.val)))
# 

gsea = ggplot(age.gsea.reduced, aes(factor(NAME, levels=rev(unique(NAME))), NES, fill=`Log Q value`)) +
  geom_bar(stat="identity") + 
  coord_flip() + 
  ylab("Normalized Enrichment Score") + 
  xlab("") + 
  theme_bw() + 
  ggtitle("Gene Set Enrichment Analysis")

# 
# grm2.data = load("~/Projects/methyl_age/data/model_data/GRM2.Rdata") %>% 
#   get %>%
#   filter(abs(M) < Inf) %>%
#   filter(series.id != "GSE56105") %>%
#   arrange(Position) %>%
#   mutate(Probe=factor(Probe, levels=unique(Probe))) 
# 
# grm2 = ggplot(grm2.data, aes(age.normed, M, color=tissue_state)) +
#   geom_point() + 
#   facet_grid(tissue_state ~ Probe) +
#   stat_smooth(method="lm",se=FALSE, color="black") + 
#   ylim(-5, 5) +
#   ylab("M value") + 
#   xlab("Normalized Age") + 
#   ggtitle("GRM2 promoter") + 
#   theme_bw()
# 
ddo.data = load("~/Projects/methyl_age/data/model_data/DDO.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

ddo = ggplot(ddo.data, aes(age.normed, M, color=tissue_state)) +
  geom_point() + 
  facet_grid(tissue_state ~ Probe) +
  stat_smooth(method="lm",se=FALSE, color="black") + 
  ylim(-5, 5) +
  xlab("Normalized Age") + 
  ylab("M value") +
  ggtitle("DDO promoter") + 
  theme_bw()

volcano = ggplot(age, aes(Estimate, logP)) + 
  geom_point(size=1) + 
  geom_text(aes(label=label), vjust=2, size=3) +
  ylab("Log10 P-value") +
  xlab("Coefficient") + 
  xlim(-max(abs(age$Estimate)), max(abs(age$Estimate))) + 
  theme_bw() +
  ggtitle("Association between Age and Methylation")

ldhd.data = load("~/Projects/methyl_age/data/model_data/ldhd.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

ldhd = ggplot(ldhd.data, aes(age.normed, M, color=tissue_state)) +
  geom_point(size=1) + 
  facet_grid(tissue_state ~ Probe) +
  stat_smooth(method="lm",se=FALSE, color="black") + 
  ylim(-5, 5) +
  xlab("Normalized Age") + 
  ylab("M value") +
  ggtitle("LDHD promoter") + 
  theme_bw() +
  scale_x_continuous(breaks=c(1,3,5)) +
  guides(color=F)

bsx.data = load("~/Projects/methyl_age/data/model_data/bsx.Rdata") %>%
  get %>%
  filter(abs(M) < Inf) %>%
  filter(series.id != "GSE56105") %>%
  arrange(Position) %>%
  mutate(Probe=factor(Probe, levels=unique(Probe)))

bsx = ggplot(bsx.data, aes(age.normed, M, color=tissue_state)) +
  geom_point(size=1) + 
  facet_grid(tissue_state ~ Probe) +
  stat_smooth(method="lm",se=FALSE, color="black") + 
  ylim(-5, 5) +
  xlab("Normalized Age") + 
  ylab("M value") +
  ggtitle("BSX promoter") + 
  theme_bw()+
  scale_x_continuous(breaks=c(1,3,5)) +
  guides(color=F)

tiff("./figures/fig1.tiff", width=22, height=22, units="in", res=150)
grid.arrange(volcano, gsea, bsx, ldhd, ncol=2)
dev.off()

tiff("./figures/ddo.tiff", width=11, height=11, units="in", res=150)
ddo
dev.off()
