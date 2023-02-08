#--------------------------------------------------------------------------#
# Last update: 21-01-2023
# Author: Luca Dei Bardi (1)(2)
# Email:  luca.deibardi@uniroma1.it
# 
# The following script is part of the data analysis made for the study
# "SARS-CoV-2 spread and area deprivation in the Italian three-tier
# restrictions: a multilevel approach".
# 
# Co-authors: Anna Acampora (1), Laura Cacciani (1), Mirko Di Martino (1),
#             Nera Agabiti (1), Marina Davoli (1), Giulia Cesaroni (1)
# 
# (1) Department of Epidemiology of the Regional Health Service, ASL Roma 1, Rome, Italy.
# (2) Sapienza University of Rome, Rome, Italy.
#--------------------------------------------------------------------------#

# required packages -------------------------------------------------------

# install.packages("pacman") # <- you may need this: "PACkage MANagement tool"

rm(list = ls()) # <- this will clear your global workspace

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pacman::p_load(tidyverse,  # -- set of usefull packages
               lubridate,  # -- date handler
               janitor,    # -- data cleaner
               lmerTest,   # -- multilevel analysis with p-values
               performance
)

select <- dplyr::select    # -- not necessary but "just in case"

# > ANALYSIS --------------------------------------------------------------

# * read data -------------------------------------------------------------

full <- read_csv2("dati/full.csv") %>% 
  mutate(date=ymd(date),
         tier = ifelse(is.na(tier),"White",tier),
         nuts1 = fct_relevel(as.factor(nuts1),
                             c("North","Central","South")))

full <- full %>% 
  group_by(nuts3,tier) %>%
  mutate(foo = c(1,diff(date)),
         foo2 = rep(1:length(diff(c(1,which(foo != 1), n()+1))),
                    diff(c(1,which(foo != 1), n()+1))),
         grp = paste(foo2,tier,nuts3)) %>% 
  ungroup() %>% 
  group_by(grp) %>% 
  mutate(days = as.numeric(date-min(date)),
         rmv_short = ifelse(n()<=6,1,0)) %>% 
  ungroup() %>% 
  group_by(nuts3) %>% 
  mutate(group_id = match(grp, unique(grp)),
         group_id2 = group_id-1) %>% 
  ungroup() %>% 
  select(-foo,-foo2)

qq3 <- quantile(full$perc_low,
                probs = seq(0, 1, 1/3),na.rm=T)

full <- full %>% 
  left_join(full %>%
              select(nuts3,
                     group_id,
                     prev_tier = tier) %>%
              unique(),
            by = c("nuts3","group_id2" = "group_id")) %>% 
  mutate(prev_dummy = case_when(tier == "Yellow" & 
                                  prev_tier == "White" ~ "Lower",
                                tier == "Yellow" & 
                                  (prev_tier == "Orange"|
                                     prev_tier == "Red") ~ "Higher",
                                tier == "Orange" & 
                                  prev_tier == "Red" ~ "Higher",
                                tier == "Orange" & 
                                  (prev_tier == "White"|
                                     prev_tier == "Yellow") ~ "Lower",
                                tier == "Red" ~ "Lower"),
         ter = case_when(perc_low <  qq3[2] ~ "1 (least deprived)",
                         perc_low <  qq3[3] & perc_low >= qq3[2] ~ "2",
                         perc_low >= qq3[3] ~ "3 (most deprived)"),
         tier = fct_relevel(tier,
                            c("Yellow","Orange","Red")),
         prev_tier = fct_relevel(prev_tier,
                            c("White","Yellow","Orange","Red")),
         dens_cat = fct_relevel(as.factor(dens_cat),
                                c("low","medium-low",
                                  "medium-high","high")),
         perc_low_center = perc_low-min(full$perc_low),
         perc05_center = perc05-min(full$perc05)) %>% 
  select(-group_id,-group_id2)


# * multilevel models -----------------------------------------------------

# ** ICC estimation -------------------------------------------------------

# install.packages("performance") # <- you may need this

# -- Yellow
icc(lmer(formula=Rt ~ 1
         + (1 | grp),
         data=subset(full,
                     full$tier=="Yellow" &
                       full$date<ymd("2021-05-10"))))
## ICC = 0.418

# -- Orange
icc(lmer(formula=Rt ~ 1
         + (1 | grp),
         data=subset(full,
                     full$tier=="Orange" &
                       full$date<ymd("2021-05-10"))))

## ICC = 0.632

# -- Red
icc(lmer(formula=Rt ~ 1
         + (1 | grp),
         data=subset(full,
                     full$tier=="Red" &
                       full$date<ymd("2021-05-10"))))
## ICC = 0.553


# ** A) and B) models -----------------------------------------------------
# -- ("days models" and "SED models")

vars <- c("days","perc_low_center")
options(scipen=10) # -- scipen=0 is default

# -- Yellow days-adjusted models
y_mod <- lapply(setNames(vars, vars), function(var) {
  form = paste("Rt ~ 1 +", var, "+ (1|grp)")
  lmer(form, data=subset(full,
                         full$tier=="Yellow" &
                           full$date<ymd("2021-05-10")))
})

y_coef <- coef(summary(y_mod$days))
y_coef <- rbind(y_coef,coef(summary(y_mod$perc_low_center)))

y_coef[,c(1,2,5)]

##                    Estimate   Std. Error      Pr(>|t|)
## (Intercept)     0.989368065 0.0084491611  0.000000e+00
## days            0.004234342 0.0002416483  2.655041e-67
## (Intercept)     1.007712201 0.0133269831 7.492236e-240
## perc_low_center 0.001417095 0.0010829655  1.914454e-01

o_mod <- lapply(setNames(vars, vars), function(var) {
  form = paste("Rt ~ 1 +", var, "+ (1|grp)")
  lmer(form, data=subset(full,
                         full$tier=="Orange" &
                           full$date<ymd("2021-05-10")))
})

o_coef <- coef(summary(o_mod$days))
o_coef <- rbind(o_coef,coef(summary(o_mod$perc_low_center)))

o_coef[,c(1,2,5)]
##                  Estimate    Std. Error      Pr(>|t|)
## (Intercept)      1.032218134 0.0077266549  0.000000e+00
## days            -0.004773838 0.0002196976 1.197727e-101
## (Intercept)      0.980868578 0.0127692081  0.000000e+00
## perc_low_center  0.002288823 0.0009798750  1.978066e-02

r_mod <- lapply(setNames(vars, vars), function(var) {
  form = paste("Rt ~ 1 +", var, "+ (1|grp)")
  lmer(form, data=subset(full,
                         full$tier=="Red" &
                           full$date<ymd("2021-05-10")))
})

r_coef <- coef(summary(r_mod$days))
r_coef <- rbind(r_coef,coef(summary(r_mod$perc_low_center)))

r_coef[,c(1,2,5)]
##                  Estimate    Std. Error      Pr(>|t|)
## (Intercept)      1.049799459 0.0093298760  0.000000e+00
## days            -0.013820031 0.0003080131  0.000000e+00
## (Intercept)      0.942090738 0.0148317880 1.490219e-232
## perc_low_center  0.004440458 0.0011279075  9.489760e-05

# -- average trend for every tier
foo <- data.frame(days = seq(0,42,by=7)) %>% 
  mutate(yellow = fixef(y_mod$days)[1] + fixef(y_mod$days)[2]*days,
         orange = fixef(o_mod$days)[1] + fixef(o_mod$days)[2]*days,
         red    = fixef(r_mod$days)[1] + fixef(r_mod$days)[2]*days)

foo
## days    yellow    orange       red
##  0   0.9893681 1.0322181 1.0497995
##  7   1.0190085 0.9988013 0.9530592
## 14   1.0486489 0.9653844 0.8563190
## 21   1.0782892 0.9319675 0.7595788
## 28   1.1079296 0.8985507 0.6628386
## 35   1.1375700 0.8651338 0.5660984
## 42   1.1672104 0.8317169 0.4693582


# ** C) interaction models ------------------------------------------------

# -- Yellow interaction model
y_int <- lmer(formula=Rt ~ 1
              + days*perc_low_center
              + (1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10")))

coef(summary(y_int))[,c(1,2,5)]


##                           Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9673726774 0.01413345026 1.088264e-239
## days                  0.0057601181 0.00044603232  1.015713e-37
## perc_low_center       0.0021030122 0.00114724674  6.745654e-02
## days:perc_low_center -0.0001283412 0.00003143216  4.493599e-05


# -- Orange interaction model
o_int <- lmer(formula=Rt ~ 1
              + days*perc_low_center
              + (1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10")))

coef(summary(o_int))[,c(1,2,5)]

##                           Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9731616955 0.01291885561 0.000000e+00
## days                  0.0011272405 0.00038475559 3.402283e-03
## perc_low_center       0.0056693748 0.00099293962 1.624036e-08
## days:perc_low_center -0.0005474999 0.00002961584 1.094243e-74


# -- Red interaction model
r_int <- lmer(formula=Rt ~ 1
              + days*perc_low_center
              + (1 | grp),
              data=subset(full,
                          full$tier=="Red" &
                            full$date<ymd("2021-05-10")))

coef(summary(r_int))[,c(1,2,5)]

##                           Estimate    Std. Error      Pr(>|t|)
## (Intercept)           1.0261747518 0.01524771584 1.704836e-257
## days                 -0.0176445458 0.00050177325 2.092175e-238
## perc_low_center       0.0024428262 0.00115560383  3.500503e-02
## days:perc_low_center  0.0003394143 0.00003529912  1.135238e-21


# ** D) adjusted models ---------------------------------------------------

# -- Yellow complete model
y_adj <- lmer(formula=Rt ~ 1
              + days*perc_low_center
              + perc05_center
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10")))

coef(summary(y_adj))[,c(1,2,5)]

##                           Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9275122838 0.02775740755 8.342582e-118
## days                  0.0057864656 0.00044602614  4.784661e-38
## perc_low_center      -0.0015020932 0.00250011958  5.483070e-01
## perc05_center         0.0536214468 0.01803982764  3.135587e-03
## dens_catmedium-low    0.0173219366 0.02338342760  4.592701e-01
## dens_catmedium-high   0.0194804670 0.02355960305  4.088172e-01
## dens_cathigh         -0.0305089390 0.02490941931  2.213843e-01
## nuts1Central          0.0096985569 0.02312288834  6.751252e-01
## nuts1South            0.0550225280 0.04224406853  1.935137e-01
## days:perc_low_center -0.0001296535 0.00003142595  3.739687e-05

# -- Orange complete model
o_adj <- lmer(formula=Rt ~ 1
              + days*perc_low_center
              + perc05_center
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10")))

coef(summary(o_adj))[,c(1,2,5)]

##                          Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9849829323 0.02563665212 2.692529e-177
## days                  0.0011265534 0.00038478511  3.424338e-03
## perc_low_center       0.0046446192 0.00235473246  4.894199e-02
## perc05_center        -0.0021927940 0.01719041553  8.985339e-01
## dens_catmedium-low   -0.0009516747 0.02206686878  9.656125e-01
## dens_catmedium-high   0.0175853268 0.02208089485  4.260625e-01
## dens_cathigh         -0.0174465552 0.02361986872  4.603683e-01
## nuts1Central         -0.0184905679 0.02297213708  4.211373e-01
## nuts1South            0.0136804624 0.04017670961  7.335754e-01
## days:perc_low_center -0.0005475812 0.00002961841  1.073922e-74

# -- Red complete model
r_adj <- lmer(formula=Rt ~ 1
              + days*perc_low_center
              + perc05_center
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Red" &
                            full$date<ymd("2021-05-10")))

coef(summary(r_adj))[,c(1,2,5)]

##                          Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9932102577 0.02993962272 4.168618e-128
## days                 -0.0176547669 0.00050167623 9.597286e-239
## perc_low_center       0.0087922625 0.00284726070  2.127482e-03
## perc05_center         0.0169297847 0.02003002117  3.984011e-01
## dens_catmedium-low    0.0314064002 0.02579105921  2.239161e-01
## dens_catmedium-high   0.0419450657 0.02591093362  1.061311e-01
## dens_cathigh         -0.0161488963 0.02718480893  5.527613e-01
## nuts1Central         -0.0778890271 0.02780665805  5.294610e-03
## nuts1South           -0.1382606396 0.04906626671  5.027293e-03
## days:perc_low_center  0.0003399854 0.00003528981  9.494966e-22


# *** R2: model assumptions -----------------------------------------------

path <- "plot/"

# -- Yellow
jpeg(paste(path,"assumptions-Y.jpg",sep = ""),
     units = "in",
     width  = 7.5,
     height = 5,
     res = 1000)

par(mfrow=c(2,2))
scatter.smooth(fitted(y_adj),resid(y_adj),col="#FFE0A1") # iid
abline(h=0, lty=2)
title("Tukey-Anscombe Plot")
qqnorm(resid(y_adj),main="normality of residuals",col="#FFE0A1") # residuals~N
qqline(resid(y_adj))
scatter.smooth(fitted(y_adj),sqrt(abs(resid(y_adj))),col="#FFE0A1") # res. var vs. fitted
qqnorm(ranef(y_adj)$grp$"(Intercept)",
       main="normality of random intercepts",col="#FFE0A1") # random intercepts~N
qqline(ranef(y_adj)$grp$"(Intercept)")

dev.off()

# -- Orange, iid
jpeg(paste(path,"assumptions-O.jpg",sep = ""),
     units = "in",
     width  = 7.5,
     height = 5,
     res = 1000)

par(mfrow=c(2,2))
scatter.smooth(fitted(o_adj),resid(o_adj),col="#FFA077") # iid
abline(h=0, lty=2)
title("Tukey-Anscombe Plot")
qqnorm(resid(o_adj),main="normality of residuals",col="#FFA077") # residuals~N
qqline(resid(o_adj))
scatter.smooth(fitted(o_adj),sqrt(abs(resid(o_adj))),col="#FFA077") # res. var vs. fitted
qqnorm(ranef(o_adj)$grp$"(Intercept)",
       main="normality of random intercepts",col="#FFA077") # random intercepts~N
qqline(ranef(o_adj)$grp$"(Intercept)")

dev.off()

# -- Red, iid
jpeg(paste(path,"assumptions-R.jpg",sep = ""),
     units = "in",
     width  = 7.5,
     height = 5,
     res = 1000)

par(mfrow=c(2,2))
scatter.smooth(fitted(r_adj),resid(r_adj),col="#D0696B") # iid
abline(h=0, lty=2)
title("Tukey-Anscombe Plot")
qqnorm(resid(r_adj),main="normality of residuals",col="#D0696B") # residuals~N
qqline(resid(r_adj))
scatter.smooth(fitted(r_adj),sqrt(abs(resid(r_adj))),col="#D0696B") # res. var vs. fitted
qqnorm(ranef(r_adj)$grp$"(Intercept)",
       main="normality of random intercepts",col="#D0696B") # random intercepts~N
qqline(ranef(r_adj)$grp$"(Intercept)")

dev.off()

par(mfrow=c(1,1))


# > SENSITIVITY -----------------------------------------------------------

# * without short restrictions --------------------------------------------

sens.y.full <- lmer(formula=Rt ~ 1
                    + days*perc_low_center
                    + perc05_center
                    + dens_cat
                    + nuts1
                    + (1 | grp),
                    data=subset(full,
                                full$tier=="Yellow" &
                                  full$rmv_short==0 &
                                  full$date<ymd("2021-05-10")))

coef(summary(sens.y.full))[,c(1,2,5)]

##                           Estimate    Std. Error     Pr(>|t|)
## (Intercept)           0.8831151243 0.02805359246 7.234995e-93
## days                  0.0061260747 0.00045883206 3.995423e-40
## perc_low_center      -0.0028181955 0.00257231355 2.742253e-01
## perc05_center         0.0366509691 0.01905579263 5.552149e-02
## dens_catmedium-low    0.0532507286 0.02429665426 2.928898e-02
## dens_catmedium-high   0.0242496563 0.02484393147 3.299350e-01
## dens_cathigh          0.0111059278 0.02613737183 6.712588e-01
## nuts1Central          0.0112301627 0.02476079318 6.505313e-01
## nuts1South            0.0367584944 0.04398689953 4.041034e-01
## days:perc_low_center -0.0001210073 0.00003230268 1.812291e-04


# -- Orange
sens.o.full <- lmer(formula=Rt ~ 1
                    + days*perc_low_center
                    + perc05_center
                    + dens_cat
                    + nuts1
                    + (1 | grp),
                    data=subset(full,
                                full$tier=="Orange" &
                                  full$rmv_short==0 &
                                  full$date<ymd("2021-05-10")))

coef(summary(sens.o.full))[,c(1,2,5)]

##                           Estimate    Std. Error     Pr(>|t|)
## (Intercept)           0.9441474616 0.03506783540 4.433757e-92
## days                  0.0012947354 0.00039378445 1.014185e-03
## perc_low_center       0.0057633757 0.00313064435 6.637062e-02
## perc05_center        -0.0062100958 0.02299613227 7.872631e-01
## dens_catmedium-low    0.0076337071 0.02973789217 7.975449e-01
## dens_catmedium-high   0.0187056142 0.02968637666 5.289898e-01
## dens_cathigh          0.0146290147 0.03200799184 6.478919e-01
## nuts1Central         -0.0346722405 0.03211221249 2.809305e-01
## nuts1South            0.0324536517 0.05406503848 5.486712e-01
## days:perc_low_center -0.0005584755 0.00003032021 4.808249e-74

# -- Red
sens.r.full <- lmer(formula=Rt ~ 1
                    + days*perc_low_center
                    + perc05_center
                    + dens_cat
                    + nuts1
                    + (1 | grp),
                    data=subset(full,
                                full$tier=="Red" &
                                  full$rmv_short==0 &
                                  full$date<ymd("2021-05-10")))

coef(summary(sens.r.full))[,c(1,2,5)]

##                           Estimate     Std. Error     Pr(>|t|)
## (Intercept)           1.0907372216 0.03992514851  3.159955e-60
## days                 -0.0180152717 0.00053452983 1.286450e-214
## perc_low_center       0.0129825796 0.00406027136  1.697704e-03
## perc05_center        -0.0212038794 0.02633121324  4.219957e-01
## dens_catmedium-low    0.0171818142 0.03345163490  6.082966e-01
## dens_catmedium-high  -0.0052401601 0.03376904628  8.768982e-01
## dens_cathigh         -0.0007366732 0.03391152926  9.826986e-01
## nuts1Central         -0.1660317499 0.03625981506  9.971705e-06
## nuts1South           -0.2596153675 0.07177895093  4.114870e-04
## days:perc_low_center  0.0003540535 0.00003751578  6.914717e-21


# * weights ---------------------------------------------------------------

wgt_pop <- full %>% 
  select(nuts3,pop) %>% 
  unique() %>% 
  mutate(tot_pop=sum(pop),
         wgt = pop/tot_pop)

full <- full %>% 
  left_join(wgt_pop %>% select(nuts3,wgt),by="nuts3")


# -- Yellow complete model
w.y.int <- lmer(formula=Rt ~ 1
                + days*perc_low_center
                + perc05_center
                + dens_cat
                + nuts1
                + (1 | grp),
                data=subset(full,
                            full$tier=="Yellow" &
                              full$date<ymd("2021-05-10")),
                weights = wgt)

coef(summary(w.y.int))[,c(1,2,5)]

##                           Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9331456773 0.02744188285 1.662705e-121
## days                  0.0039034866 0.00037056294  9.458390e-26
## perc_low_center      -0.0025922199 0.00244933595  2.905304e-01
## perc05_center         0.0572004366 0.01749356221  1.172212e-03
## dens_catmedium-low    0.0219598888 0.02318135013  3.440382e-01
## dens_catmedium-high   0.0230520973 0.02320672439  3.211395e-01
## dens_cathigh         -0.0216011606 0.02412864114  3.712157e-01
## nuts1Central          0.0115838291 0.02266639281  6.095910e-01
## nuts1South            0.0649254638 0.04137004624  1.173526e-01
## days:perc_low_center -0.0000607166 0.00002896237  3.608338e-02

# -- Orange complete model
w.o.int <- lmer(formula=Rt ~ 1
                + days*perc_low_center
                + perc05_center
                + dens_cat
                + nuts1
                + (1 | grp),
                data=subset(full,
                            full$tier=="Orange" &
                              full$date<ymd("2021-05-10")),
                weights = wgt)


coef(summary(w.o.int))[,c(1,2,5)]

##                           Estimate    Std. Error      Pr(>|t|)
## (Intercept)           0.9728426336 0.02566380247 5.557407e-175
## days                  0.0023648174 0.00034082575  4.288896e-12
## perc_low_center       0.0047217953 0.00232934399  4.303208e-02
## perc05_center         0.0016857264 0.01695340947  9.208235e-01
## dens_catmedium-low   -0.0001133614 0.02207730421  9.959045e-01
## dens_catmedium-high   0.0183723579 0.02196533198  4.031987e-01
## dens_cathigh         -0.0161689459 0.02321019350  4.862704e-01
## nuts1Central         -0.0164560469 0.02273954706  4.695071e-01
## nuts1South            0.0187143388 0.03973951085  6.378423e-01
## days:perc_low_center -0.0006017817 0.00002710892 7.586876e-106

# -- Red complete model
w.r.int <- lmer(formula=Rt ~ 1
                + days*perc_low_center
                + perc05_center
                + dens_cat
                + nuts1
                + (1 | grp),
                data=subset(full,
                            full$tier=="Red" &
                              full$date<ymd("2021-05-10")),
                weights = wgt)

coef(summary(w.r.int))[,c(1,2,5)]

##                           Estimate   Std. Error      Pr(>|t|)
## (Intercept)           0.9932134271 0.0295023181 5.146010e-131
## days                 -0.0170992079 0.0004115333 5.169413e-317
## perc_low_center       0.0079177833 0.0027779487  4.552586e-03
## perc05_center         0.0154632930 0.0194521071  4.270472e-01
## dens_catmedium-low    0.0294072666 0.0254148254  2.477905e-01
## dens_catmedium-high   0.0414379736 0.0253288742  1.024917e-01
## dens_cathigh         -0.0194584212 0.0262296108  4.585595e-01
## nuts1Central         -0.0731076722 0.0270816165  7.190237e-03
## nuts1South           -0.1323703092 0.0478534847  5.886833e-03
## days:perc_low_center  0.0005036448 0.0000283491  3.630125e-68


# > REVIEWER'S COMMENTS ---------------------------------------------------

# R1: other age groups? ---------------------------------------------------


popres <- read_csv("dati/popres_prov_01-01-20.csv",skip=1) %>% 
  clean_names() %>% 
  select(nuts3 = provincia,
         age = eta,
         tot_M = totale_maschi,
         tot_F = totale_femmine) %>% 
  filter(age != "Totale") %>% 
  mutate(age = as.numeric(age),
         nuts3 = str_squish(str_replace_all(nuts3,"[[:punct:]]", " ")),
         total = tot_M+tot_F) %>% 
  select(-tot_F,-tot_M) %>% 
  group_by(nuts3) %>% 
  mutate(tot_pop = sum(total)) %>% 
  ungroup() %>% 
  mutate(age_gr = case_when(age <= 5    ~ "pop0_5",
                            age >= 6 &  age <= 13 ~ "pop6_13",
                            age >= 14 & age <= 19 ~ "pop14_19",
                            age >= 20 & age <= 29 ~ "pop20_29",
                            age >= 30 & age <= 49 ~ "pop30_49",
                            age >= 50 & age <= 64 ~ "pop50_64",
                            age >= 65   ~ "pop65_")) %>% 
  group_by(nuts3,age_gr) %>% 
  mutate(total_age_gr = sum(total)) %>% 
  ungroup() %>% 
  select(-age,-total) %>% 
  unique() %>% 
  pivot_wider(names_from = "age_gr",
              values_from = "total_age_gr") %>% 
  mutate(perc0_5   = pop0_5/tot_pop*100,
         perc6_13  = pop6_13/tot_pop*100,
         perc14_19 = pop14_19/tot_pop*100,
         perc20_29 = pop20_29/tot_pop*100,
         perc30_49 = pop30_49/tot_pop*100,
         perc50_64 = pop50_64/tot_pop*100,
         perc65    = pop65_/tot_pop*100)

cor(popres$perc0_5,popres$perc6_13)  ## =  0.9222359
cor(popres$perc0_5,popres$perc14_19) ## =  0.7654012
cor(popres$perc0_5,popres$perc20_29) ## =  0.6230641
cor(popres$perc0_5,popres$perc30_49) ## =  0.6080708
cor(popres$perc0_5,popres$perc50_64) ## = -0.8164889
cor(popres$perc0_5,popres$perc65)    ## = -0.8550843

rm(popres)


# R2: previous tier? ------------------------------------------------------

# -- Yellow complete model
y_rev <- lmer(formula=Rt ~ 1
              + days*perc_low
              + days*prev_tier
              + perc05
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10")))

coef(summary(y_rev))[,c(1,2,5)]
##                            Estimate    Std. Error      Pr(>|t|)
## (Intercept)           1.00674962138 0.07269752978  4.490278e-36
## days                  0.00213301031 0.00096300154  2.679480e-02
## perc_low             -0.00001129625 0.00187894275  9.952060e-01
## prev_tierOrange      -0.33786811606 0.01802741601  5.406924e-58
## prev_tierRed         -0.01301046886 0.02549219054  6.099051e-01
## perc05                0.04157997614 0.01353260174  2.263725e-03
## dens_catmedium-low    0.01020933541 0.01751289739  5.602431e-01
## dens_catmedium-high   0.01910119703 0.01767508895  2.804748e-01
## dens_cathigh         -0.03150574534 0.01865000581  9.192653e-02
## nuts1Central         -0.00135820211 0.01738105913  9.377526e-01
## nuts1South            0.03196200037 0.03164215344  3.130474e-01
## days:perc_low        -0.00025164014 0.00002904087  5.540630e-18
## days:prev_tierOrange  0.01592643740 0.00046966750 1.602266e-233
## days:prev_tierRed    -0.01915500277 0.02273114606  3.994402e-01


# -- Orange complete model
o_rev <- lmer(formula=Rt ~ 1
              + days*perc_low
              + days*prev_tier
              + perc05
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10")))


coef(summary(o_rev))[,c(1,2,5)]
##                           Estimate    Std. Error     Pr(>|t|)
## (Intercept)           1.2032715138 0.09183198639 4.229435e-35
## days                  0.0084374128 0.00157382755 8.515282e-08
## perc_low              0.0046058783 0.00217875600 3.488073e-02
## prev_tierYellow      -0.2739463921 0.04408586077 9.350747e-10
## prev_tierRed         -0.4331548524 0.04351652450 8.562442e-22
## perc05                0.0111149575 0.01595834126 4.863572e-01
## dens_catmedium-low   -0.0073035164 0.02045063472 7.211069e-01
## dens_catmedium-high   0.0022603626 0.02047175984 9.121142e-01
## dens_cathigh         -0.0256389095 0.02186741655 2.414199e-01
## nuts1Central         -0.0155297309 0.02124870255 4.651202e-01
## nuts1South           -0.0028076014 0.03732029046 9.400542e-01
## days:perc_low        -0.0005957501 0.00003208448 2.530889e-75
## days:prev_tierYellow  0.0032530890 0.00100275633 1.183203e-03
## days:prev_tierRed     0.0095902974 0.00100496290 1.837612e-21

r_rev <- lmer(formula=Rt ~ 1
              + days*perc_low
              + days*prev_tier
              + perc05
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Red" &
                            full$date<ymd("2021-05-10")))


coef(summary(r_rev))[,c(1,2,5)]
##                           Estimate  Std. Error      Pr(>|t|)
## (Intercept)           1.0725488050 0.104753163  2.104951e-22
## days                 -0.0355180222 0.001309953 1.921516e-149
## perc_low              0.0106053547 0.002720901  1.104362e-04
## prev_tierYellow      -0.4699758612 0.043291082  1.276284e-24
## prev_tierOrange      -0.3539695589 0.039340492  6.105481e-18
## perc05                0.0129388687 0.019035898  4.970132e-01
## dens_catmedium-low    0.0297002887 0.024513692  2.262633e-01
## dens_catmedium-high   0.0444270883 0.024622385  7.179854e-02
## dens_cathigh         -0.0137507376 0.025829664  5.947207e-01
## nuts1Central         -0.0489939150 0.026792631  6.806535e-02
## nuts1South           -0.1315845536 0.046868046  5.189287e-03
## days:perc_low         0.0001628874 0.000038528  2.409035e-05
## days:prev_tierYellow  0.0207603549 0.001098481  1.589991e-76
## days:prev_tierOrange  0.0185631972 0.000869722  5.662846e-96


# R2: relative Rt? --------------------------------------------------------

full <- full %>% 
  left_join(full %>% filter(days==0) %>% 
              select(grp,Rt_start=Rt),
            by="grp") %>% 
  mutate(Rt_rel = Rt/Rt_start)

rel.y.adj <- lmer(formula=Rt_rel ~ 1
                  + days*perc_low
                  + perc05
                  + dens_cat
                  + nuts1
                  + (1 | grp),
                  data=subset(full,
                              full$tier=="Yellow" &
                                full$date<ymd("2021-05-10")))

coef(summary(rel.y.adj))[,c(1,2,5)]
##                      Estimate      Std. Error    Pr(>|t|)
## (Intercept)          0.67544511209 0.10760296577 7.099539e-10
## days                 0.00901427831 0.00117696079 2.128498e-14
## perc_low             0.00450232415 0.00280332112 1.088381e-01
## perc05               0.04375861992 0.02022433705 3.093037e-02
## dens_catmedium-low   0.03591541569 0.02621319161 1.712228e-01
## dens_catmedium-high -0.00671056365 0.02641158695 7.995340e-01
## dens_cathigh         0.00854538159 0.02792446276 7.597102e-01
## nuts1Central        -0.03939696480 0.02592375032 1.291667e-01
## nuts1South          -0.08825436979 0.04735443327 6.291450e-02
## days:perc_low       -0.00006663375 0.00003601897 6.436128e-02

# -- Orange complete model
rel.o.adj <- lmer(formula=Rt_rel ~ 1
                  + days*perc_low
                  + perc05
                  + dens_cat
                  + nuts1
                  + (1 | grp),
                  data=subset(full,
                              full$tier=="Orange" &
                                full$date<ymd("2021-05-10")))

coef(summary(rel.o.adj))[,c(1,2,5)]
##                       Estimate    Std. Error     Pr(>|t|)
## (Intercept)          1.0563542681 0,07728198397 1.711515e-38
## days                 0.0122262867 0,00105159140 5.422910e-31
## perc_low             0.0036221587 0,00206167941 7.931140e-02
## perc05              -0.0329529106 0,01503673723 2.870192e-02
## dens_catmedium-low  -0.0014617435 0,01930264928 9.396546e-01
## dens_catmedium-high -0.0065315801 0,01931359756 7.353116e-01
## dens_cathigh         0.0083017058 0,02068195142 6.882322e-01
## nuts1Central         0.0116408895 0,02011384827 5.629183e-01
## nuts1South          -0.0022639791 0,03512049670 9.486175e-01
## days:perc_low       -0.0004853506 0,00003343134 3.802352e-47

# -- Red complete model
rel.r.adj <- lmer(formula=Rt_rel ~ 1
                  + days*perc_low
                  + perc05
                  + dens_cat
                  + nuts1
                  + (1 | grp),
                  data=subset(full,
                              full$tier=="Red" &
                                full$date<ymd("2021-05-10")))

coef(summary(rel.r.adj))[,c(1,2,5)]
##                       Estimate    Std. Error     Pr(>|t|)
## (Intercept)          0.941189094 0.05900193686  2.088089e-47
## days                -0.023198604 0.00091678045  3.19607e-132
## perc_low             0.003057633 0.00165590560  6.534067e-02
## perc05               0.004067338 0.01158800463  7.257281e-01
## dens_catmedium-low  -0.003587533 0.01491977305  8.100686e-01
## dens_catmedium-high  0.003113822 0.01499187753  8.355407e-01
## dens_cathigh        -0.006750071 0.01569822067  6.673768e-01
## nuts1Central        -0.014906695 0.01609113219  3.546537e-01
## nuts1South          -0.063834303 0.02851117969  2.555159e-02
## days:perc_low        0.000361034 0.00002838905  1.985081e-36

summary(full$Rt)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.2434  0.8181  0.9565  0.9815  1.1206  3.1783 

summary(full$Rt_rel)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.1408  0.8361  0.9785  0.9898  1.0799  4.3538 

summary(unique(full$Rt_start))
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.2780  0.8761  1.0305  1.0486  1.2019  2.4949 


# R2: continuous pop density? ---------------------------------------------

# -- Yellow complete model
y_rev <- lmer(formula=Rt ~ 1
              + days*perc_low
              + perc05
              + pop_density
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10")))

coef(summary(y_rev))[,c(1,2,5)]
##                Estimate      Std. Error     Pr(>|t|)
## (Intercept)    0.77978607494 0.09589145024 5.403875e-15
## days           0.00831005236 0.00102679155 6.807757e-16
## perc_low      -0.00121205866 0.00250167764 6.282961e-01
## perc05         0.05267933176 0.01821775296 4.044065e-03
## pop_density   -0.00002739001 0.00002462620 2.667235e-01
## nuts1Central   0.01174024284 0.02306184593 6.109812e-01
## nuts1South     0.05310096890 0.04230613242 2.101677e-01
## days:perc_low -0.00012850779 0.00003142462 4.373921e-05

# -- Orange complete model
o_rev <- lmer(formula=Rt ~ 1
              + days*perc_low
              + perc05
              + pop_density
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10")))

coef(summary(o_rev))[,c(1,2,5)]
##                Estimate      Std. Error     Pr(>|t|)
## (Intercept)    0.90345519396 0.08797332158 3.467050e-23
## days           0.01198781962 0.00093164813 1.679190e-37
## perc_low       0.00468806303 0.00233196349 4.476843e-02
## perc05        -0.00249040861 0.01721017241 8.849843e-01
## pop_density   -0.00001193956 0.00002173760 5.830002e-01
## nuts1Central  -0.01840403578 0.02293359786 4.225342e-01
## nuts1South     0.01412676942 0.04013003741 7.249241e-01
## days:perc_low -0.00054766381 0.00002961912 1.030100e-74

# -- Red complete model
r_rev <- lmer(formula=Rt ~ 1
              + days*perc_low
              + perc05
              + pop_density
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Red" &
                            full$date<ymd("2021-05-10")))

coef(summary(r_rev))[,c(1,2,5)]
##                Estimate      Std. Error     Pr(>|t|)
## (Intercept)    0.77216904122 0.10175078163 1.620329e-13
## days          -0.02439631975 0.00113991502 1.169905e-96
## perc_low       0.00897409077 0.00280596561 1.470641e-03
## perc05         0.01787669815 0.02007641925 3.736704e-01
## pop_density   -0.00002485331 0.00002403531 3.016467e-01
## nuts1Central  -0.07683382432 0.02790173988 6.110603e-03
## nuts1South    -0.13946417725 0.04908388783 4.676389e-03
## days:perc_low  0.00033999908 0.00003529345 9.551430e-22


# ------------------------------------------------------------------------- #

#  that's all folks! ------------------------------------------------------

# ------------------------------------------------------------------------- #