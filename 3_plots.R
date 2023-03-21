#--------------------------------------------------------------------------#
# Last update: 21-03-2023
# Author: Luca Dei Bardi (1)(2)
# Email:  luca.deibardi@uniroma1.it
# 
# The following script is part of the data analysis made for the study
# "SARS-CoV-2 spread and area economic disadvantage in the Italian three-tier
# restrictions: a multilevel approach".
# 
# Co-authors: Anna Acampora (1), Laura Cacciani (1), Mirko Di Martino (1),
#             Nera Agabiti (1), Marina Davoli (1), Giulia Cesaroni (1)
# 
# (1) Department of Epidemiology of the Regional Health Service, ASL Roma 1, Rome, Italy.
# (2) Sapienza University of Rome, Rome, Italy.
#--------------------------------------------------------------------------#

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls()) # <- this will clear your global workspace

path <- "plot/"


# required packages -------------------------------------------------------

# install.packages("pacman") # <- you may need this: "PACkage MANagement tool"

pacman::p_load(tidyverse,        # -- set of usefull packages
               viridis,          # -- beautiful color palette
               lubridate,        # -- date handler
               lmerTest,         # -- multilevel analysis
               rmapshaper,       # -- handles polygon
               sf,               # -- encodes spatial vector data
               ggthemes,         # -- extra themes for ggplot2
               janitor,          # -- data cleaner
               EpiEstim,         # -- reproduction number
               slider,           # -- sliding window functions
               eurostat,         # -- Eurostat open data
               arm,              # -- per CI Fig2 (?)
               ggpubr)           # -- arranges figures

select <- dplyr::select

# * Shapefile -------------------------------------------------------------

bord <- eurostat_geodata_60_2016 %>%        # -- regions' borders
  clean_names() %>% 
  st_transform(crs = 3035) %>% 
  filter(cntr_code == "IT", levl_code == 2) %>% 
  ms_innerlines()

gd_it <- eurostat_geodata_60_2016 %>%       # -- merging Sardinia's provinces
  clean_names() %>% 
  st_transform(crs = 3035) %>% 
  filter(cntr_code == "IT", levl_code == 3) %>% 
  ms_filter_islands(min_area = 10^8) %>%
  mutate(prov_new = ifelse(nuts_name=="Medio Campidano"|
                             nuts_name=="Carbonia-Iglesias",
                           "Sud Sardegna",
                           ifelse(nuts_name=="Ogliastra","Nuoro",
                                  ifelse(nuts_name=="Olbia-Tempio","Sassari",
                                         nuts_name)))) %>% 
  group_by(prov_new) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup() %>% 
  mutate(prov_new = str_squish(str_replace_all(prov_new,"[[:punct:]]", " ")))


# > PLOT (paper) ----------------------------------------------------------

# * read data -------------------------------------------------------------

full <- read_csv2("dati/full.csv") %>% 
  mutate(date=ymd(date),
         tier = ifelse(is.na(tier),"White",tier))

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
         ter = case_when(perc_low <  qq3[2] ~ "1 (least disadvantaged)",
                         perc_low <  qq3[3] & perc_low >= qq3[2] ~ "2",
                         perc_low >= qq3[3] ~ "3 (most disadvantaged)"),
         tier = fct_relevel(tier,
                            c("Yellow","Orange","Red")),
         dens_cat = fct_relevel(as.factor(dens_cat),
                                c("low","medium-low",
                                  "medium-high","high"))) %>% 
  select(-group_id,-group_id2)


# * Fig1: (maps) ----------------------------------------------------------

foo <- full %>% 
  select(nuts3,date,incidence,pop,
         perc_low,dens_cat,
         perc05) %>% 
  filter(date<ymd("2021-05-10"),
         date>ymd("2020-11-05")) %>% 
  group_by(nuts3) %>% 
  mutate(tot_cases = cumsum(incidence),
         cuminc = tot_cases/pop*10000,) %>%
  ungroup() %>% 
  filter(date==ymd("2021-05-09"))

p1 <-
  left_join(gd_it,foo,by=c("prov_new"="nuts3")) %>% 
  ggplot(aes(fill=perc_low))+
  geom_sf(col="grey",lwd=.333)+
  geom_sf(data = bord,col="black",
          inherit.aes = F)+
  coord_sf(datum = NA)+
  theme_map(base_size=15)+
  scale_fill_viridis_c(direction=-1,option = "D",
                       breaks=seq(10,50,by=5))+
  guides(fill=guide_legend(title="",
                           override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+
  labs(title = "% of people with income < 10,000€",
       subtitle = "2019 data")+
  theme(legend.position = c(.7,0.55),
        title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))

p1

p2 <-
  left_join(gd_it,foo,by=c("prov_new"="nuts3")) %>% 
  ggplot(aes(fill=cuminc))+
  geom_sf(col="grey",lwd=.333)+
  geom_sf(data = bord,col="black",
          inherit.aes = F)+
  coord_sf(datum = NA)+
  theme_map(base_size=15)+
  scale_fill_viridis_c(direction=-1,option="E",
                       breaks=seq(100,1300,by=200))+
  guides(fill=guide_legend(title="",
                           override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+
  labs(title = "Total cases per 10,000 residents",
       subtitle = "06 Nov 2020 - 09 May 2021")+
  theme(legend.position = c(.7,0.55),
        title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))

p2

p3 <- 
  left_join(gd_it,foo,by=c("prov_new"="nuts3")) %>% 
  ggplot(aes(fill=perc05))+
  geom_sf(col="grey",lwd=.333)+
  geom_sf(data = bord,col="black",
          inherit.aes = F)+
  coord_sf(datum = NA)+
  theme_map(base_size=15)+
  scale_fill_viridis_c(direction=1,option="B",
                       breaks=seq(0,15,by=1))+
  guides(fill=guide_legend(title="",
                           override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+
  labs(title = "% of 0-5 years old",
       subtitle = "2020 data")+
  theme(legend.position = c(.7,0.55),
        title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))

p3

p4 <- 
  left_join(gd_it, foo,
            by=c("prov_new"="nuts3")) %>% 
  ggplot(aes(fill=dens_cat))+
  geom_sf(col="grey",lwd=.333)+
  geom_sf(data = bord,col="black",
          inherit.aes = F)+
  coord_sf(datum = NA)+
  theme_map(base_size=15)+
  scale_fill_viridis_d(direction=-1,option="A")+
  guides(fill=guide_legend(title="",#pop/kmq
                           override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+
  labs(title = "Population density, quartiles",
       subtitle = "2020 data")+
  theme(legend.position = c(.7,0.55),
        title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))

p4

unico <- ggarrange(p1,p2,p3,p4,
                   ncol=2,nrow=2)

unico

ggsave(paste(path,"Fig1_maps.jpg",sep = ""),plot = unico,
       units = "in",
       width  = 12,
       height = 9)


# * Fig2 (spaghettiplot) --------------------------------------------------

# ** models ---------------------------------------------------------------

# -- multilevel linear regression by SEP and tier strata
# -- Yellow, 1st tercile
RE.y1 <- lmer(formula=Rt ~ 1
              + days
              + ( 1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10") &
                            full$ter=="1 (least disadvantaged)"))
# -- Yellow, 2nd tercile
RE.y2 <- lmer(formula=Rt ~ 1
              + days
              + ( 1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10") &
                            full$ter=="2"))

# -- Yellow 3rd tercile
RE.y3 <- lmer(formula=Rt ~ 1
              + days
              + ( 1 | grp),
              data=subset(full,
                          full$tier=="Yellow" &
                            full$date<ymd("2021-05-10") &
                            full$ter=="3 (most disadvantaged)"))

# -- Orange, 1st tercile
RE.o1 <- lmer(formula=Rt ~ 1
              + days
              + ( 1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10") &
                            full$ter=="1 (least disadvantaged)"))

# -- Orange, 2nd tercile
RE.o2 <- lmer(formula=Rt ~ 1
              + days
              + ( 1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10") &
                            full$ter=="2"))

# -- Orange, 3rd tercile
RE.o3 <- lmer(formula=Rt ~ 1
              + days
              + ( 1 | grp),
              data=subset(full,
                          full$tier=="Orange" &
                            full$date<ymd("2021-05-10") &
                            full$ter=="3 (most disadvantaged)"))

# -- Red, 1st tercile
RE.r1  <- lmer(formula=Rt ~ 1
               + days
               + ( 1 | grp),
               data=subset(full,
                           full$tier=="Red" &
                             full$date<ymd("2021-05-10") &
                             full$ter=="1 (least disadvantaged)"))

# -- Red, 2nd tercile
RE.r2  <- lmer(formula=Rt ~ 1
               + days
               + ( 1 | grp),
               data=subset(full,
                           full$tier=="Red" &
                             full$date<ymd("2021-05-10") &
                             full$ter=="2"))

# -- Red, 3rd tercile
RE.r3  <- lmer(formula=Rt ~ 1
               + days
               + ( 1 | grp),
               data=subset(full,
                           full$tier=="Red" &
                             full$date<ymd("2021-05-10") &
                             full$ter=="3 (most disadvantaged)"))
# # -- saving estimates
estimates <- data.frame(tier = c(rep("Yellow",3),
                                 rep("Orange",3),
                                 rep("Red",3)))

estimates$ter <- rep(c("1 (least disadvantaged)","2","3 (most disadvantaged)"),3)

estimates$alpha <- c(RE.y1@beta[1],RE.y2@beta[1],RE.y3@beta[1],
                     RE.o1@beta[1],RE.o2@beta[1],RE.o3@beta[1],
                     RE.r1@beta[1],RE.r2@beta[1],RE.r3@beta[1])

estimates$beta <- c(RE.y1@beta[2],RE.y2@beta[2],RE.y3@beta[2],
                    RE.o1@beta[2],RE.o2@beta[2],RE.o3@beta[2],
                    RE.r1@beta[2],RE.r2@beta[2],RE.r3@beta[2])

estimates <- estimates %>%
  mutate(tier=fct_relevel(tier,c("Yellow","Orange","Red")),
         a = round(alpha,2),
         b = round(beta,3),
         label = paste("Rt = ",a,ifelse(b<0," "," +"),b,"•days",sep=""))

models <- c("RE.y1"=RE.y1,"RE.y2"=RE.y2,"RE.y3"=RE.y3,  # -- [[1-2-3]]
            "RE.o1"=RE.o1,"RE.o2"=RE.o2,"RE.o3"=RE.o3,  # -- [[4-5-6]]
            "RE.r1"=RE.r1,"RE.r2"=RE.r2,"RE.r3"=RE.r3)  # -- [[7-8-9]]

ind <- c("RE.y1","RE.y2","RE.y3",  # -- [[1-2-3]]
         "RE.o1","RE.o2","RE.o3",  # -- [[4-5-6]]
         "RE.r1","RE.r2","RE.r3")  # -- [[7-8-9]]

newavg <- data.frame(expand_grid(days = 0:75, mod = ind))
newavg$lower <- NA
newavg$median <- NA
newavg$upper <- NA

foo <- data.frame(days = 0:75)
foo$lower <- NA
foo$median <- NA
foo$upper <- NA

for (i in ind) {
  print(i)
  # -- simulating plausible values
  sims <- sim(models[[i]],n.sim=10000) # <- this may take a while (<5 min)
  # -- saving simulations
  fs <- fixef(sims)
  # -- creating a predictor matrix
  Xmat <- model.matrix(~ 1 + days, data = foo)
  # -- matrix for fitted values
  fitmat <- matrix(ncol = nrow(fs), nrow = nrow(foo))
  # -- fitted values for each combination of the predictors,
  #    for each combination of the parameters
  for (j in 1:nrow(fs)) {
    fitmat[, j] <- Xmat %*% as.matrix(fs)[j, ]
  }
  # -- saving quantiles of the fitted values
  foo$lower <- apply(fitmat, 1, quantile, prob = 0.025)
  foo$median <- apply(fitmat, 1, quantile, prob = 0.5)
  foo$upper <- apply(fitmat, 1, quantile, prob = 0.975)
  newavg$lower[newavg$mod==i] <- foo$lower
  newavg$median[newavg$mod==i] <- foo$median
  newavg$upper[newavg$mod==i] <- foo$upper
}

newavg <- newavg %>% 
  mutate(tier = case_when(str_detect(mod,".y") ~ "Yellow",
                          str_detect(mod,".o") ~ "Orange",
                          str_detect(mod,".r") ~ "Red"),
         ter = case_when(str_detect(mod,"1") ~ "1 (least disadvantaged)",
                         str_detect(mod,"2") ~ "2",
                         str_detect(mod,"3") ~ "3 (most disadvantaged)")) %>% 
  mutate(tier = fct_relevel(tier,
                            c("Yellow","Orange","Red"))) %>% 
  left_join(full %>%
              select(ter,tier,days) %>%
              group_by(ter,tier) %>% 
              mutate(max = max(days)) %>% 
              ungroup() %>% 
              select(-days) %>% 
              unique(),
            by = c("ter","tier")) %>% 
  mutate(rmv = ifelse(days<=max,0,1)) %>% 
  filter(rmv == 0)
         
rm(fitmat,foo,fs,models,sims,Xmat,
   RE.y1,RE.y2,RE.y3,
   RE.o1,RE.o2,RE.o3,
   RE.r1,RE.r2,RE.r3)


# ** plot -----------------------------------------------------------------

colcol <- full %>% 
  filter(tier=="Yellow"|tier=="Orange"|tier=="Red") %>%
  select(tier) %>%
  unique()

full %>% 
  filter(tier=="Yellow"|tier=="Orange"|tier=="Red",
         date<ymd("2021-05-10")) %>% 
  mutate(tier = fct_relevel(tier,
                            c("Yellow","Orange","Red"))) %>% 
  ggplot(aes(x=days,y=Rt,
             group=grp))+
  geom_rect(data = colcol,
            aes(fill = tier),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,
            alpha = 0.53,inherit.aes = F,
            show.legend = F)+
  geom_line(linewidth=.1,col="grey50",show.legend = F)+
  geom_hline(yintercept = 1,col="black",linetype="dashed",linewidth=.3)+
  geom_ribbon(data=newavg,
              fill = "black", linewidth=0, col="black",
              inherit.aes = F, alpha =.5,
              aes(x=days,ymin=lower,ymax=upper,group=1))+
  geom_line(data=newavg,
            linewidth=.77,col="black",
            inherit.aes = F, linetype ="dotted",
            aes(x=days,y=median,group=1))+
  geom_text(aes(label=label, group=1,
                y = 0.25,x = 20),
            data=estimates,
            fontface = "bold",
            col="Black",size=3.3)+
  scale_fill_manual(labels=c("Yellow",
                             "Orange",
                             "Red"),
                    values=c("#FED976",
                             "#FD8D3C",   # -- Orange
                             "#BD0000"))+
  facet_grid(tier~ter)+
  coord_cartesian(ylim=c(0,2),xlim=c(0,51))+
  theme_linedraw(base_size=15)+
  scale_x_continuous(expand = c(0,0))+
  xlab("Days since tier implementation")+
  ylab("Reproduction number")

ggsave(paste(path,"Fig2_trendsCI.jpg",sep = ""),
       units = "in",
       width  = 7.5,
       height = 5)

rm(colcol,newavg,estimates)

# > PLOT (supplementary) --------------------------------------------------

# * FigS1: tier by date ---------------------------------------------------

Sys.setlocale(category = "LC_ALL", locale = "en_GB.UTF-8")

full %>% 
  select(date,nuts3,nuts2,nuts1,tier) %>% 
  mutate(new_nuts2 = ifelse(nuts3 == "Bolzano Bozen" |
                              nuts3 == "Trento",nuts3,nuts2),
         new_nuts2 = case_when(new_nuts2 == "Valle d Aosta Vallée d Aoste" ~
                                 "Valle d'Aosta",
                               T ~ new_nuts2),
         rip_geo = case_when(nuts1=="North" ~ 1,
                             nuts1=="Central" ~ 2,
                             nuts1=="South" ~ 3),
         nuts2_eng = case_when(new_nuts2 == "Valle d'Aosta" ~ "Aosta Valley",
                               new_nuts2 == "Piemonte" ~ "Piedmont",
                               new_nuts2 == "Lombardia" ~ "Lombardy",
                               new_nuts2 == "Toscana" ~ "Tuscany",
                               new_nuts2 == "Puglia" ~ "Apulia",
                               new_nuts2 == "Sardegna" ~ "Sardinia",
                               new_nuts2 == "Sicilia" ~ "Sicily",
                               new_nuts2 == "Trento" ~ "Trentino",
                               new_nuts2 == "Bolzano Bozen" ~ "South Tyrol",
                               T ~ new_nuts2),
         tier2 = as.factor(case_when(tier == "White"  ~ "Low",
                                     tier == "Yellow" ~ "Moderate",
                                     tier == "Orange" ~ "Elevated",
                                     tier == "Red"    ~ "Maximum")),
         tier2 = fct_relevel(tier2, c("Low","Moderate",
                                      "Elevated","Maximum"))) %>% 
  unique() %>% 
  filter(date>ymd("2020-11-01"),
         date<ymd("2021-05-17")) %>% 
  ggplot(aes(x=date,y=fct_reorder(nuts2_eng,rip_geo,.desc=T),col=tier2)) +
  geom_point() +
  geom_vline(xintercept = ymd("2021-01-16"),linetype="dotted",linewidth=1.1) +
  geom_text(aes(label="pandemic threat assessment update",
                y = 13.5, x = ymd("2021-01-16")),
            inherit.aes = F,show.legend=F,col=1,
            vjust=-0.7,angle=270,size=4,
            data = data.frame())+
  scale_color_manual(values=c("grey70",
                              "#FFE0A1",
                              "#FFA077",   # -- Orange
                              "#D0696B"))+
  labs(x = "", y = "") +
  scale_x_date(date_breaks = "1 month",
               date_labels =  "%b %y",
               expand = c(0.02,0.02))+
  theme_linedraw(base_size = 15)+
  guides(col=guide_legend(title="Pandemic threat"))+
  theme(legend.position="top")

Sys.setlocale(category = "LC_ALL", locale = "")

ggsave(paste(path,"FigS1_colors.jpg",sep = ""),
       units = "in",
       width  = 7.5,
       height = 5)


# * FigS2: Rt x Incidence ITA ---------------------------------------------

path_covid <- "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/"

path_it <- paste(path_covid,
                 "dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv",
                 sep = "")
shape <- 1.87
rate <- 0.28

# Serial Interval distribution (normalized)
SI <- (dgamma(0:300,shape=shape,
              rate=rate)/sum(dgamma(0:3000,shape=shape,
                                    rate=rate))) 
conf <- list(si_distr = SI, n1=10000,
             mcmc_control = make_mcmc_control(thin=1,burnin=1000000))

italia <- read_csv(path_it) %>% 
  mutate(date=as.Date(data),
         date=ymd(date)) %>% 
  select(date,nuovi = nuovi_positivi) %>% 
  mutate(inc_smt = floor(slide_index_dbl(.x = nuovi,
                                         .f = mean,
                                         .i = date,
                                         .before=3,
                                         .after=3)))

italia <- italia %>%    
  mutate(Rt_smt=c(NA,NA,NA,NA,
                  estimate_R(incid=inc_smt,
                             method="non_parametric_si",
                             config = make_config(conf))$R$"Mean(R)",
                  NA,NA,NA)) %>% 
  filter(date<ymd("2022-01-01"))

italia <- italia %>% 
  mutate(colcol = as.factor(case_when(date < ymd("2020-11-06") |
                                        date > ymd("2021-05-09") ~ 1,
                                      T ~ 2,)))


coeff <- 1/50

Sys.setlocale(category = "LC_ALL", locale = "en_GB.UTF-8")

ggplot(italia, aes(x=date, color=colcol,group=1)) +
  geom_hline(yintercept = 1,col="black",linetype="dashed")+
  geom_line(aes(y=Rt_smt), size=2) +
  geom_line(aes(y=inc_smt*coeff/1000), size=2) +
  geom_point(aes(y=Rt_smt), size=1) +
  geom_point(aes(y=inc_smt*coeff/1000), size=1) +
  geom_text(aes(label="Incidence", 
                y = 0.2, x = ymd("2020-04-01")),
            hjust=0,
            fontface = "bold",
            col="grey60",size=5,data = data.frame())+
  geom_text(aes(label="Rt", 
                y = 1.35, x = ymd("2020-04-01")),
            hjust=0,
            fontface = "bold",
            col="grey60",size=5,data = data.frame())+
  geom_text(aes(label="6 Nov 20 - 9 May 21", 
                y = 1.35, x = ymd("2021-02-27")),
            fontface = "bold",
            col="#FFA077",size=5,data = data.frame())+
  scale_y_continuous(name = "Reproduction number, Rt",
                     sec.axis = sec_axis(~./coeff,
                                         name="Incidence (x 1,000)")) +
  scale_color_manual(values = c("grey60","#FFA077"))+
  theme_linedraw(base_size = 17) +
  scale_x_date(date_breaks = "2 month",
               date_labels =  "%b %y",
               expand = c(0.02,0.02))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 33, vjust = 1, hjust=1))+
  xlab("")+
  coord_cartesian(ylim=c(0,2))

ggsave(paste(path,"FigS2_RtxInc.jpg",sep = ""),
       units = "in",
       width  = 7.5,
       height = 5)


# ------------------------------------------------------------------------- #

#  that's all folks! ------------------------------------------------------

# ------------------------------------------------------------------------- #
