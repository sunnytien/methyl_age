library("caret")
library("randomForest")
library("glmnet")
library("gbm")
library("kernlab")
library("dplyr")
library("doParallel")
library("purrr")
library("tidyr")
library("data.table")
library("car")

source("./model_building/get_params.R")
source("./model_building/util.R")

age.trans = function(x) sapply(x, function(y) if(y < 20) log10((y + 1)/21) + 1 else (y+1)/21)
#age.trans = function(x) yjPower(x, coef(yj))

load("./data/sample.info.Rdata")
load("./data/probe.info.Rdata")

select = dplyr::select
filter = dplyr::filter
group_by = dplyr::group_by
mutate = dplyr::mutate

### DATA PREPARATION ###

yj = powerTransform(sample.info$age, family="yjPower")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

probes = probe.info %>%
  filter(as.character(nearestGeneSymbol) %in% c(head(age$gene,5), tail(age$gene,5)))

b.tmp = beta %>%
  filter(Probe %in% probes$Probe) %>%
  collect %>%
  gather(gsm.id, beta, starts_with("GSM")) %>% # transposing
  spread(Probe, beta) # this may be more efficient with matrix stuff


data = sample.info %>%
  select(gsm.id, age) %>%
  filter(!is.na(age)) %>%
  mutate(age.normed=age.trans(age)) %>%
  select(-age) %>%
  inner_join(b.tmp)

data1 = data %>%
  sample_n(3000)

data2 = data %>%
  anti_join(data1 %>% select(gsm.id)) %>%
  sample_n(1000)

data3 = data %>%
  anti_join(data1 %>% select(gsm.id)) %>%
  anti_join(data2 %>% select(gsm.id))

#### TRAINING FIRST SET OF MODELS ###

registerDoParallel(20)

x1 = data1 %>%
  select(starts_with("cg")) %>%
  as.matrix

y1 = data1$age.normed

ctrl = trainControl(method="repeatedcv",
                    number=10,
                    repeats=1,
                    savePredictions=T,
                    allowParallel=T)

model.names = c("pcr",
                "pls",
                "ppr",
                "glmnet",
                "gaussprRadial",
                "rvmLinear",
                "rvmRadial",
                "svmLinear",
                "svmRadialCost",
                "rf",
                "cubist",
                "gbm")

models = plyr::llply(model.names,
                     train.wrapper.safe,
                     x1,       
                     y1,
                     trControl=ctrl,
                     .progress="text")

names(models) = model.names
models = models[!is.na(models)]

### TRAINING ENSEMBLE MODEL ###

x2 = data2 %>%
  select(starts_with("cg")) %>%
  as.matrix

z2 = models %>%
  lapply(predict, newdata=x2) %>%
  do.call(cbind, .)

y2 = data2$age.normed

ensemble = cv.glmnet(z2, 
                     y2,
                     alpha=0.1,
                     nfolds=10,
                     keep=T)

### TESTING ENSEMBLE MODEL ###

x3 = data3 %>%
  select(starts_with("cg")) %>%
  as.matrix

z3 = models %>%
  lapply(predict, newdata=x3) %>%
  do.call(cbind, .)

y3 = data3$age.normed

y3.pred = predict(ensemble, 
                  newx=as.matrix(z3),
                  s="lambda.min")

plot(cor(z3, y3, method="spearman"))

y3.pred = matrix(y3.pred, ncol=1)
colnames(y3.pred) = "ensemble"


cols = c(ppr="red",
         svmLinear="green",
         pls="red",
         rvmLinear="green",
         gbm="blue",
         rvmRadial="green",
         pcr="red",
         svmRadialCost="green",
         gaussprRadial="green",
         rf="blue",
         ensemble="black",
         cubist="blue")

perf = cor(cbind(z3, y3.pred), y3) %>%
  as.data.frame %>%
  mutate(method=row.names(.)) %>%
  mutate(cor=V1) %>%
  mutate(col=cols[method]) %>%
  arrange(col, cor)

par(mar=c(7.1,5.1,2.1,2.1))
barplot(perf$cor,
        col=perf$col,
        names.arg=perf$method,
        las=2,
        ylab="Correlation between Age and DNAm Age")
