library("nlme")

### adjusting for distance ###

ances = model.matrix(~predicted.ancestry + Probe - 1, data=data) %>%
  as.data.frame

data2 = data %>%
  cbind(ances) %>%
  mutate(constant.adj = 1 / (distance + 1),
         age.adj = age / (distance + 1),
         tissue_state.adj = as.numeric(tissue_state) / (distance + 1),
         predicted.ancestryAFR.adj = predicted.ancestryAFR / (distance + 1),
         predicted.ancestryASN.adj = predicted.ancestryASN / (distance + 1)) 


m = lm(M ~ age.adj + constant.adj, data=data2)



m = nls(M ~ a * age * exp(-k*distance) + distance,
        data=data2,
        start=list(a=0.3, k=0, c=0))

m = nls(M ~ a * age * exp(-k*scale(distance)), 
        data=data2,
        start=list(a=0.3, k=1))

m = nls(M ~ a * age * exp(-(Position - p)^2 / s),
        data=data2,
        start=list(a=0.3, 
                   p=mean(data2$Position), 
                   s=sd(data2$Position)))

