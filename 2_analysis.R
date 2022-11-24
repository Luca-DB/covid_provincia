#--------------------------------------------------------------------------#
# Last update: 23-11-2022
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

pacman::p_load(tidyverse,  # -- set of usefull packages
               lubridate,  # -- date handler
               lmerTest,   # -- multilevel analysis with p-values
               performance,
               openxlsx
)


# > ANALYSIS --------------------------------------------------------------

# * read data -------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
         dens_cat = fct_relevel(as.factor(dens_cat),
                                c("low","medium-low",
                                  "medium-high","high"))) %>% 
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

vars <- c("days","perc_low")
options(scipen=10) # -- scipen=0 is default

# -- Yellow days-adjusted models
y_mod <- lapply(setNames(vars, vars), function(var) {
  form = paste("Rt ~ 1 +", var, "+ (1|grp)")
  lmer(form, data=subset(full,
                         full$tier=="Yellow" &
                           full$date<ymd("2021-05-10")))
})

y_coef <- coef(summary(y_mod$days))
y_coef <- rbind(y_coef,coef(summary(y_mod$perc_low)))

y_coef[,c(1,2,5)]

##              Estimate   Std. Error      Pr(>|t|)
## (Intercept) 0.989368065 0.0084491611  0.000000e+00
## days        0.004234342 0.0002416483  2.655041e-67
## (Intercept) 0.979614824 0.0331209755 9.349648e-103
## perc_low    0.001417095 0.0010829655  1.914454e-01

o_mod <- lapply(setNames(vars, vars), function(var) {
  form = paste("Rt ~ 1 +", var, "+ (1|grp)")
  lmer(form, data=subset(full,
                         full$tier=="Orange" &
                           full$date<ymd("2021-05-10")))
})

o_coef <- coef(summary(o_mod$days))
o_coef <- rbind(o_coef,coef(summary(o_mod$perc_low)))

o_coef[,c(1,2,5)]
##              Estimate   Std. Error      Pr(>|t|)
## (Intercept)  1.032218134 0.0077266549  0.000000e+00
## days        -0.004773838 0.0002196976 1.197727e-101
## (Intercept)  0.935487056 0.0306002583 3.642874e-131
## perc_low     0.002288823 0.0009798750  1.978066e-02

r_mod <- lapply(setNames(vars, vars), function(var) {
  form = paste("Rt ~ 1 +", var, "+ (1|grp)")
  lmer(form, data=subset(full,
                         full$tier=="Red" &
                           full$date<ymd("2021-05-10")))
})

r_coef <- coef(summary(r_mod$days))
r_coef <- rbind(r_coef,coef(summary(r_mod$perc_low)))

r_coef[,c(1,2,5)]
##              Estimate   Std. Error      Pr(>|t|)
## (Intercept)  1.049799459 0.0093298760 0.000000e+00
## days        -0.013820031 0.0003080131 0.000000e+00
## (Intercept)  0.854047797 0.0353092740 1.213123e-84
## perc_low     0.004440458 0.0011279075 9.489760e-05

# ** C) interaction models ------------------------------------------------

# -- Yellow interaction model
y_int <- lmer(formula=Rt ~ 1
              + days*perc_low
              + (1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            # full$rmv_short==0 &
                            full$date<ymd("2021-05-10")))

coef(summary(y_int))[,c(1,2,5)]

##                 Estimate    Std. Error     Pr(>|t|)
## (Intercept)    0.9256753170 3.511088e-02 4.395394e-93
## days           0.0083047965 1.027031e-03 7.207190e-16
## perc_low       0.0021030122 1.147247e-03 6.745654e-02
## days:perc_low -0.0001283412 3.143216e-05 4.493599e-05


# -- Orange interaction model
o_int <- lmer(formula=Rt ~ 1
              + days*perc_low
              + (1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            # full$rmv_short==0 &
                            full$date<ymd("2021-05-10")))

coef(summary(o_int))[,c(1,2,5)]

##                 Estimate    Std. Error     Pr(>|t|)
## (Intercept)    0.8607524766 3.099760e-02 1.270732e-117
## days           0.0119827649 9.315533e-04  1.768494e-37
## perc_low       0.0056693748 9.929396e-04  1.624036e-08
## days:perc_low -0.0005474999 2.961584e-05  1.094243e-74


# -- Red interaction model
r_int <- lmer(formula=Rt ~ 1
              + days*perc_low
              + (1 | grp),
              data=subset(full,
                          full$tier=="Red" &
                            # full$rmv_short==0 &
                            full$date<ymd("2021-05-10")))

coef(summary(r_int))[,c(1,2,5)]

##                 Estimate    Std. Error     Pr(>|t|)
## (Intercept)    0.9777397484 0.03623962957  1.479416e-100
## days          -0.0243742648 0.00114006496  1.806236e-96
## perc_low       0.0024428262 0.00115560383  3.500503e-02
## days:perc_low  0.0003394143 0.00003529912  1.135238e-21


# ** D) adjusted models ---------------------------------------------------

# -- Yellow complete model
y_adj <- lmer(formula=Rt ~ 1
              + days*perc_low
              + perc05
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10")))

coef(summary(y_adj))[,c(1,2,5)]

##                       Estimate    Std. Error     Pr(>|t|)
## (Intercept)          0.7748476487 0.09597796147 8.289063e-15
## days                 0.0083571624 0.00102690785 4.717915e-16
## perc_low            -0.0015020932 0.00250011958 5.483070e-01
## perc05               0.0536214468 0.01803982763 3.135587e-03
## dens_catmedium-low   0.0173219366 0.02338342760 4.592701e-01
## dens_catmedium-high  0.0194804670 0.02355960304 4.088172e-01
## dens_cathigh        -0.0305089390 0.02490941930 2.213843e-01
## nuts1Central         0.0096985569 0.02312288834 6.751251e-01
## nuts1South           0.0550225280 0.04224406851 1.935137e-01
## days:perc_low       -0.0001296535 0.00003142595 3.739687e-05

# -- Orange complete model
o_adj <- lmer(formula=Rt ~ 1
              + days*perc_low
              + perc05
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10")))


coef(summary(o_adj))[,c(1,2,5)]

##                       Estimate    Std. Error     Pr(>|t|)
## (Intercept)          0.9003529943 0.08812791029 5.617929e-23
## days                 0.0119836894 0.00093162863 1.770163e-37
## perc_low             0.0046446192 0.00235473246 4.894199e-02
## perc05              -0.0021927940 0.01719041553 8.985339e-01
## dens_catmedium-low  -0.0009516747 0.02206686878 9.656125e-01
## dens_catmedium-high  0.0175853268 0.02208089485 4.260625e-01
## dens_cathigh        -0.0174465552 0.02361986872 4.603683e-01
## nuts1Central        -0.0184905679 0.02297213708 4.211373e-01
## nuts1South           0.0136804624 0.04017670961 7.335754e-01
## days:perc_low       -0.0005475812 0.00002961841 1.073922e-74

# -- Red complete model
r_adj <- lmer(formula=Rt ~ 1
              + days*perc_low
              + perc05
              + dens_cat
              + nuts1
              + (1 | grp),
              data=subset(full,
                          full$tier=="Red" &
                            full$date<ymd("2021-05-10")))

coef(summary(r_adj))[,c(1,2,5)]

##                       Estimate    Std. Error     Pr(>|t|)
## (Intercept)          0.7612784552 0.10173959965 3.367031e-13
## days                -0.0243958094 0.00113979261 1.125212e-96
## perc_low             0.0087922625 0.00284726070 2.127482e-03
## perc05               0.0169297847 0.02003002117 3.984011e-01
## dens_catmedium-low   0.0314064002 0.02579105921 2.239161e-01
## dens_catmedium-high  0.0419450657 0.02591093362 1.061311e-01
## dens_cathigh        -0.0161488963 0.02718480893 5.527613e-01
## nuts1Central        -0.0778890271 0.02780665805 5.294610e-03
## nuts1South          -0.1382606396 0.04906626671 5.027293e-03
## days:perc_low        0.0003399854 0.00003528981 9.494966e-22


# > SENSITIVITY -----------------------------------------------------------

# * without short restrictions --------------------------------------------

sens.y.full <- lmer(formula=Rt ~ 1
                    + days*perc_low
                    + perc05
                    + dens_cat
                    + nuts1
                    + (1 | grp),
                    data=subset(full,
                                full$tier=="Yellow" &
                                  full$rmv_short==0 &
                                  full$date<ymd("2021-05-10")))

coef(summary(sens.y.full))[,c(1,2,5)]

##                      Estimate     Std. Error     Pr(>|t|)
## (Intercept)          0.8142875910 0.09823725341 5.768087e-15
## days                 0.0085253400 0.00105590132 8.067708e-16
## perc_low            -0.0028181955 0.00257231355 2.742253e-01
## perc05               0.0366509691 0.01905579264 5.552149e-02
## dens_catmedium-low   0.0532507287 0.02429665427 2.928898e-02
## dens_catmedium-high  0.0242496563 0.02484393148 3.299350e-01
## dens_cathigh         0.0111059278 0.02613737184 6.712588e-01
## nuts1Central         0.0112301627 0.02476079319 6.505313e-01
## nuts1South           0.0367584944 0.04398689955 4.041034e-01
## days:perc_low       -0.0001210073 0.00003230268 1.812291e-04


# -- Orange
sens.o.full <- lmer(formula=Rt ~ 1
                    + days*perc_low
                    + perc05+ dens_cat
                    + nuts1
                    + (1 | grp),
                    data=subset(full,
                                full$tier=="Orange" &
                                  full$rmv_short==0 &
                                  full$date<ymd("2021-05-10")))

coef(summary(sens.o.full))[,c(1,2,5)]

##                      Estimate     Std. Error     Pr(>|t|)
## (Intercept)          0.8510043349 0.12034423324 6.892318e-12
## days                 0.0123678772 0.00095360764 4.916544e-38
## perc_low             0.0057633757 0.00313064435 6.637062e-02
## perc05              -0.0062100958 0.02299613228 7.872631e-01
## dens_catmedium-low   0.0076337071 0.02973789218 7.975449e-01
## dens_catmedium-high  0.0187056142 0.02968637668 5.289898e-01
## dens_cathigh         0.0146290147 0.03200799185 6.478919e-01
## nuts1Central        -0.0346722405 0.03211221251 2.809305e-01
## nuts1South           0.0324536517 0.05406503851 5.486712e-01
## days:perc_low       -0.0005584755 0.00003032021 4.808249e-74

# -- Red
sens.r.full <- lmer(formula=Rt ~ 1
                    + days*perc_low
                    + perc05
                    + dens_cat
                    + nuts1
                    + (1 | grp),
                    data=subset(full,
                                full$tier=="Red" &
                                  full$rmv_short==0 &
                                  full$date<ymd("2021-05-10")))

coef(summary(sens.r.full))[,c(1,2,5)]

##                      Estimate     Std. Error     Pr(>|t|)
## (Intercept)          0.9054721659 0.13608048358 5.283381e-10
## days                -0.0250352487 0.00121289638 3.350011e-89
## perc_low             0.0129825796 0.00406027136 1.697704e-03
## perc05              -0.0212038794 0.02633121325 4.219957e-01
## dens_catmedium-low   0.0171818142 0.03345163491 6.082966e-01
## dens_catmedium-high -0.0052401601 0.03376904629 8.768982e-01
## dens_cathigh        -0.0007366732 0.03391152926 9.826986e-01
## nuts1Central        -0.1660317499 0.03625981507 9.971705e-06
## nuts1South          -0.2596153675 0.07177895094 4.114870e-04
## days:perc_low        0.0003540535 0.00003751578 6.914717e-21


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
                + days*perc_low
                + perc05
                + dens_cat
                + nuts1
                + (1 | grp),
                data=subset(full,
                            full$tier=="Yellow" &
                              full$date<ymd("2021-05-10")),
                weights = wgt)

coef(summary(w.y.int))[,c(1,2,5)]

##                      Estimate     Std. Error     Pr(>|t|)
## (Intercept)          0.7899179332 0.09360357576 6.146393e-16
## days                 0.0051073416 0.00090729351 1.881584e-08
## perc_low            -0.0025922199 0.00244933599 2.905304e-01
## perc05               0.0572004367 0.01749356248 1.172212e-03
## dens_catmedium-low   0.0219598888 0.02318135047 3.440383e-01
## dens_catmedium-high  0.0230520974 0.02320672474 3.211395e-01
## dens_cathigh        -0.0216011609 0.02412864152 3.712157e-01
## nuts1Central         0.0115838289 0.02266639316 6.095910e-01
## nuts1South           0.0649254637 0.04137004687 1.173526e-01
## days:perc_low       -0.0000607166 0.00002896237 3.608338e-02

# -- Orange complete model
w.o.int <- lmer(formula=Rt ~ 1
                + days*perc_low
                # + prev_dummy
                + perc05
                # + perc_vaxI
                # + alpha
                + dens_cat
                + nuts1
                + (1 | grp),
                data=subset(full,
                            full$tier=="Orange" &
                              full$date<ymd("2021-05-10")),
                weights = wgt)


coef(summary(w.o.int))[,c(1,2,5)]

##                      Estimate     Std. Error     Pr(>|t|)
## (Intercept)          0.8734857990 0.08720118267  3.589992e-22
## days                 0.0142966103 0.00083808929  4.633088e-64
## perc_low             0.0047217953 0.00232934381  4.303206e-02
## perc05               0.0016857264 0.01695340816  9.208235e-01
## dens_catmedium-low  -0.0001133613 0.02207730254  9.959045e-01
## dens_catmedium-high  0.0183723578 0.02196533031  4.031987e-01
## dens_cathigh        -0.0161689453 0.02321019170  4.862704e-01
## nuts1Central        -0.0164560466 0.02273954532  4.695070e-01
## nuts1South           0.0187143399 0.03973950780  6.378423e-01
## days:perc_low       -0.0006017817 0.00002710892 7.586911e-106

# -- Red complete model
w.r.int <- lmer(formula=Rt ~ 1
                + days*perc_low
                + perc05
                + dens_cat
                + nuts1
                + (1 | grp),
                data=subset(full,
                            full$tier=="Red" &
                              full$date<ymd("2021-05-10")),
                weights = wgt)

coef(summary(w.r.int))[,c(1,2,5)]

##                      Estimate     Std. Error     Pr(>|t|)
## (Intercept)          0.7836100616 0.0991495064  1.858396e-14
## days                -0.0270851984 0.0009282561 9.380668e-171
## perc_low             0.0079177833 0.0027779487  4.552586e-03
## perc05               0.0154632930 0.0194521071  4.270472e-01
## dens_catmedium-low   0.0294072666 0.0254148254  2.477905e-01
## dens_catmedium-high  0.0414379736 0.0253288742  1.024917e-01
## dens_cathigh        -0.0194584212 0.0262296108  4.585595e-01
## nuts1Central        -0.0731076722 0.0270816165  7.190237e-03
## nuts1South          -0.1323703092 0.0478534847  5.886833e-03
## days:perc_low        0.0005036448 0.0000283491  3.630125e-68


# ------------------------------------------------------------------------- #

#  that's all folks! ------------------------------------------------------

# ------------------------------------------------------------------------- #