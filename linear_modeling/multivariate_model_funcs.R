run.models = function(probe.info, sample.info, predicted.ancestry){

  gene = probe.info$nearestGeneSymbol[1]
  
  liquid_tissues = c("Leukocytes", "Lymphoblasts", "Lymphocytes", "Monocytes",
                     "T-cells", "Whole Blood")
  
  db = src_sqlite("./data/BMIQ.db")
  
  ### find ids present everywhere
  
  beta = tbl(db, "BMIQ") %>%
    filter(Probe %in% probe.info$Probe) %>%
    collect
  
  b = beta %>%
    select(-Probe) %>%
    as.matrix %>% 
    t
  colnames(b) = beta$Probe
  
  data = predicted.ancestry %>%
    select(gsm.id, predicted.ancestry) %>%
    inner_join(sample.info) %>%
    filter(!is.na(age)) %>%
    mutate(tissue_state=ifelse(tissue %in% liquid_tissues, "liquid", "solid"))
  rownames(data) = data$gsm.id
  
  ids = intersect(rownames(b), rownames(data)) %>%
    unique %>%
    sort
  
  b = b[ids,]
  data = data[ids, ]
  
  m.full = lm(b ~ age + predicted.ancestry + tissue_state,
              data=data)
  
  m.no.age = lm(b ~ predicted.ancestry + tissue_state,
                data=data)
  
  m.no.ancestry = lm(b ~ age + tissue_state,
              data=data)
  
  m.no.tissue = lm(b ~ age + predicted.ancestry ,
                    data=data)
  
  a = Anova(m.full, type="III")
  
  save(m.full, m.no.age, m.no.ancestry, m.no.tissue, 
       file=paste("./data/models/", gene, ".Rdata", sep=""))
  
  save(a, file=paste("./data/anovas/", gene, ".Rdata", sep=""))
}

library("nlme")
