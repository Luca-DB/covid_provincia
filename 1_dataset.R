#--------------------------------------------------------------------------#
# Last update: 21-01-2023
# Author: Luca Dei Bardi (1)(2)
# Email:  luca.deibardi@uniroma1.it
# 
# The following script is part of the data analysis made for the study
# "SARS-CoV-2 spread and area economic disadvantage in the italian three-tier
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls()) # <- this will clear your global workspace

pacman::p_load(tidyverse,  # -- set of useful packages
               lubridate,  # -- date handler
               janitor,    # -- data cleaner
               slider,     # -- sliding window functions (moving average)
               EpiEstim,   # -- reproduction number, among others
               readxl      # -- reads xls and xlsx files
               )

select <- dplyr::select    # -- not necessary but "just in case"


# > DATA ------------------------------------------------------------------

# * utils -----------------------------------------------------------------

path_covid <- "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/"
path_prov <- paste(path_covid,
                   "dati-province/dpc-covid19-ita-province.csv",
                   sep = "")

# -- North (East + West)
N <- c("Liguria","Lombardia","Piemonte","Valle d Aosta Vallée d Aoste",
       "Emilia Romagna","Friuli Venezia Giulia",
       "Bolzano Bozen","Trento","Veneto")
# -- Central
C <- c("Lazio","Marche","Toscana","Umbria")
# -- South (South + Insular)
S <- c("Abruzzo","Basilicata","Calabria","Campania",
       "Molise","Puglia",
       "Sicilia","Sardegna")

codes_nuts3 <- read_xls("dati/Codici-statistici-e-denominazioni-al-01_01_2020.xls",
                  sheet = "CODICI al 01012020",
                  col_types = "text") %>% 
  clean_names() %>% 
  select(nuts3 = denominazione_dell_unita_territoriale_sovracomunale_valida_a_fini_statistici,
         prov  = sigla_automobilistica) %>% 
  unique() %>% 
  mutate(nuts3 = str_squish(str_replace_all(nuts3,
                                            "[[:punct:]]", " ")),
         nuts3 = ifelse(nuts3=="Reggio Calabria",
                        "Reggio di Calabria",nuts3),
         prov = ifelse(nuts3=="Napoli", # -- Naples' code is NA,
                       "Napoli",prov))  #    better to avoid problems.

codes_mun <- read_xls("dati/Codici-statistici-e-denominazioni-al-01_01_2020.xls",
                      sheet = "CODICI al 01012020",
                      col_types = "text") %>% 
  clean_names() %>% 
  select(mun = codice_comune_formato_alfanumerico,
         municipality = denominazione_in_italiano,
         nuts3=denominazione_dell_unita_territoriale_sovracomunale_valida_a_fini_statistici) %>% 
  unique() %>% 
  mutate(nuts3 = str_squish(str_replace_all(nuts3,
                                            "[[:punct:]]", " ")),
         nuts3 = ifelse(nuts3=="Reggio Calabria",
                        "Reggio di Calabria",nuts3),
         municipality = str_squish(str_replace_all(municipality,
                                                   "[[:punct:]]", " ")))


# * COVID data ------------------------------------------------------------

# ** loading and cleaning data --------------------------------------------

# -- saving offline COVID data until May 2021
## c19_save <- read_csv(path_prov) %>%
##   mutate(date = ymd(as.Date(data))) %>%
##   filter(date <  ymd("2021-06-01"))
## 
## write_csv2(c19_save,"dati/covid_data_bck.csv")
## rm(c19_save)

# -- you can either read data from path_prov
#    or read: "dati/covid_data_bck.csv"
c19_pr <- read_csv(path_prov,
                   col_types = cols(col_character(), col_character(),
                                    col_integer(), col_character(),
                                    col_character(), col_character(),
                                    col_character(), col_skip(),
                                    col_skip(), col_integer(),
                                    .default = col_skip())) %>%
## c19_pr <- read_csv2("dati/covid_data_bck.csv") %>%
  mutate(date=as.Date(data),
         date=ymd(date),
         nuts2 = str_squish(str_replace_all(denominazione_regione,
                                            "[[:punct:]]", " ")),
         nuts2 = case_when(str_detect(nuts2,"Bolzano") ~
                             "Bolzano Bozen",
                           str_detect(nuts2,"Trento") ~
                             "Trento",
                           str_detect(nuts2,"Aosta") ~
                             "Valle d Aosta Vallée d Aoste",
                           TRUE ~ nuts2),
         nuts3 = str_squish(str_replace_all(denominazione_provincia,
                                            "[[:punct:]]", " ")),
         nuts3 = case_when(str_detect(nuts3,"Bolzano") ~ 
                             "Bolzano Bozen",
                           str_detect(nuts3,"Aosta") ~ 
                             "Valle d Aosta Vallée d Aoste",
                           TRUE ~ nuts3),
         nuts1 = case_when(nuts2%in%N ~ "North",
                           nuts2%in%C ~ "Central",
                           nuts2%in%S ~ "South"),
         nuts1 = fct_relevel(nuts1,
                             c("North","Central","South")),
         totale_casi = ifelse(nuts3=="Bolzano Bozen" &
                                date>=ymd("2021-03-22"),
# -- https://github.com/pcm-dpc/COVID-19/blob/master/note/dpc-covid19-ita-note.csv#L231
                              totale_casi-10665,
                              totale_casi)) %>% 
  filter(nuts3!="Fuori Regione Provincia Autonoma",
         nuts3!="In fase di definizione aggiornamento") %>% 
  arrange(nuts3,date) %>% 
  group_by(nuts3) %>% 
  mutate(nuovi=c(NA,diff(totale_casi)),
         nuovi = ifelse(is.na(nuovi),totale_casi,nuovi)) %>%
  ungroup() %>% 
  select(date,nuts3,nuts2,nuts1,
         totale_casi,nuovi)

# -- prob...lem solver: leveling unreasonable trends of cumulative data
prob <- c19_pr %>%                     # -- leveling drops
  mutate(yeday = date-1, 
         tmrw = date+1) %>% 
  left_join(c19_pr %>%
              select(nuts3,date,
                     tot_yeday = totale_casi),
            by=c("nuts3","yeday"="date")) %>% 
  left_join(c19_pr %>%
              select(nuts3,date,
                     tot_tmrw = totale_casi),
            by=c("nuts3","tmrw"="date")) %>%
  mutate(new_tot_casi2 = ifelse(totale_casi<tot_yeday & totale_casi<tot_tmrw&
                                  tot_yeday<=tot_tmrw, # -- IF sudden drop & increasing trend
                                floor((tot_yeday+tot_tmrw)/2), # -- DO relevel
                                totale_casi),                # -- ELSE leave as it is
         new_tot_casi2 = ifelse(is.na(new_tot_casi2), # -- to avoid NAs in first
                                totale_casi,          #    and last observation
                                new_tot_casi2))

prob <- prob %>%                       # -- leveling rises
  left_join(prob %>%
              select(nuts3,date,
                     new_tot_ieri = new_tot_casi2),
            by=c("nuts3","yeday"="date")) %>% 
  left_join(prob %>%
              select(nuts3,date,
                     new_tot_dom = new_tot_casi2),
            by=c("nuts3","tmrw"="date")) %>%
  mutate(new_tot_casi2 = ifelse(new_tot_casi2>new_tot_ieri &
                                  new_tot_casi2>new_tot_dom &
                                  new_tot_ieri<=new_tot_dom, # -- IF sudden rise & increasing trend
                                floor((new_tot_ieri+new_tot_dom)/2), # -- relevel
                                new_tot_casi2),
         new_tot_casi2 = ifelse(is.na(new_tot_casi2), # -- to avoid NAs in first
                                totale_casi,          #    and last observation
                                new_tot_casi2))     # -- ELSE leave as it is


c19_pr <- c19_pr %>% 
  left_join(prob %>% select(nuts3,date,new_tot_casi2),
            by=c("nuts3","date")) %>% 
  mutate(checky = ifelse(totale_casi==new_tot_casi2,1,0)) %>% 
  arrange(nuts3,date) %>% 
  group_by(nuts3) %>% 
  mutate(incidence=c(NA,diff(new_tot_casi2)),
         incidence = ifelse(is.na(incidence),new_tot_casi2,incidence)) %>% 
  ungroup()

# now we have new daily counts with far fewer unreasonable values
# but not every problem is solved from adjusting cumulative data.
# we will now adjust negative daily counts with a 7-days moving average.

c19_pr <- c19_pr %>%
  group_by(nuts3) %>%
  mutate(w = ifelse(incidence<0,0,1/6),
         w_incidence = w*incidence,
         avg_7day = slide_index_dbl(.x = w_incidence,
                                    .f = sum,
                                    .i = date,
                                    .before=3,
                                    .after=3),
         incidence = ifelse(incidence<0,floor(avg_7day),incidence)) %>%
  ungroup() %>% 
  select(date,nuts3,nuts2,nuts1,
         incidence)

rm(prob) # -- cleaning the environment

# no more negative values for new daily positives, we can estimate Rt:


# ** Rt estimation --------------------------------------------------------

# -- please, visit:
#    https://www.epicentro.iss.it/coronavirus/sars-cov-2-dashboard
shape <- 1.87
rate <- 0.28

# Serial Interval distribution (normalized)
SI <- (dgamma(0:300,shape=shape,
              rate=rate)/sum(dgamma(0:3000,shape=shape,
                                    rate=rate))) 
conf <- list(si_distr = SI, n1=10000,
             mcmc_control = make_mcmc_control(thin=1,burnin=1000000))

c19_pr <- c19_pr %>%    
  group_by(nuts3) %>% 
  mutate(Rt=c(NA,NA,NA,NA,
              estimate_R(incid=incidence,
                         method="non_parametric_si",
                         config = make_config(conf))$R$"Mean(R)",
              NA,NA,NA)) %>%
  ungroup()


# -- estimate_R warns us that we are using a weekly average for the Rt.
#    Also, it says we are estimating Rt too early in the epidemic.
#    It refers to the early days of the epidemic we are not analyzing.

# -- We assigned the Rt to the central day of the 7-day smoothing window.

rm(shape,rate,SI,conf) # -- cleaning the environment

# * IRPEF data ------------------------------------------------------------

# -- cleaning IRPEF data
irpef <- read_csv2(
  "dati/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2019.csv") %>% 
  clean_names() %>% 
  select(municipality = denominazione_comune,
         prov = sigla_provincia,
         nuts2 = regione,
         tot_taxpayer = numero_contribuenti,
         starts_with("reddito_complessivo")) %>% 
  pivot_longer(cols = starts_with("reddito_complessivo"),
               names_to = "income_cat",
               values_to = "value") %>% 
  mutate(tipo = ifelse(str_sub(income_cat,-9,-1)=="frequenza",
                       "taxpayers","ammontare"),
         income_cat = str_replace_all(str_replace_all(income_cat,
                                                      "[:alpha:]",""),
                                      "[:punct:]"," "),
         income_cat = str_squish(income_cat),
         income_cat = ifelse(income_cat=="","<0",income_cat),
         income_cat = str_replace_all(income_cat,"0000","0k"),
         income_cat = str_replace_all(income_cat,"000","k"),
         income_cat = str_replace_all(income_cat," ","-"),
         value = ifelse(is.na(value),0,value),    # -- this NA is missing
         prov = ifelse(is.na(prov),"Napoli",prov)) %>% # -- this is the code of Naples
  pivot_wider(names_from = tipo,
              values_from = value)

# -- grouping municipality level data at province level
#    and income classes into high-medium-low categories
irpef_prov <- irpef %>% 
  group_by(prov,income_cat) %>%
  summarise(taxpayers = sum(taxpayers)) %>% 
  unique() %>% 
  ungroup() %>%
  group_by(prov) %>% 
  mutate(tot_taxpayer = sum(taxpayers)) %>% 
  ungroup() %>% 
  left_join(codes_nuts3, by = "prov") %>% 
  filter(!is.na(nuts3),
         income_cat=="<0" | income_cat=="0-10k") %>% 
  group_by(nuts3) %>% 
  summarize(perc_low = sum(taxpayers)/tot_taxpayer*100) %>% 
  unique()

rm(irpef)


# * urbanization ----------------------------------------------------------

urb_prov <- read_xls(
  "dati/Classificazioni statistiche-e-dimensione-dei-comuni_01_01_2020.xls",
                     sheet = "Comuni al 01012020",
                     col_types = "text") %>% 
  clean_names() %>% 
  select(mun = codice_istat_del_comune_alfanumerico,
         pop = popolazione_residente_al_31_12_2019,
         kmq = superficie_territoriale_kmq_al_01_01_2020) %>% 
  left_join(codes_mun,by="mun") %>% 
  mutate(pop = as.numeric(pop),
         kmq = as.numeric(kmq)) %>% 
  group_by(nuts3) %>% 
  summarize(pop = sum(pop),
            kmq = sum(kmq),
            pop_density = pop/kmq)

qq_dens <-  quantile(urb_prov$pop_density)

urb_prov <- urb_prov %>% 
  mutate(dens_cat = case_when(pop_density <  qq_dens[2] ~ "low",
                              pop_density <  qq_dens[3] &
                                pop_density >= qq_dens[2] ~ "medium-low",
                              pop_density <  qq_dens[4] &
                                pop_density >= qq_dens[3] ~ "medium-high",
                              pop_density >= qq_dens[4] ~ "high"))


# * tiers' colors ---------------------------------------------------------

tts_day <- read_csv2("dati/tiers_nuts_complete.csv") %>% 
  mutate(nuts2 = str_squish(str_replace_all(nuts2,
                                            "[[:punct:]]", " ")),
         nuts2 = case_when(str_detect(nuts2,"Aosta") ~ 
                             "Valle d Aosta Vallée d Aoste",
                           T ~ nuts2),
         date=dmy(date)) %>% 
  filter(date >= ymd("2020-10-26") &
           date <= ymd("2021-05-23"))


# * 01/01/2020 population -------------------------------------------------
# -- https://demo.istat.it/popres/download.php?anno=2020&lingua=eng

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
  filter(age<=5) %>% 
  group_by(nuts3) %>% 
  summarize(perc05 = sum(total)/tot_pop*100) %>% 
  ungroup() %>% 
  unique()


# * join ------------------------------------------------------------------

full_day <- as.data.frame(expand_grid(c19_pr %>%
                                        select(nuts3,
                                               nuts2,
                                               nuts1) %>%
                                        unique(),
                                      c19_pr %>%
                                        filter(date >= ymd("2020-11-02"),
                                               date <= ymd("2021-05-16")) %>% 
                                        select(date) %>%
                                        unique())) %>%
  left_join(c19_pr %>% 
              select(nuts3,date,
                     incidence,Rt),
            by = c("nuts3","date")) %>% 
  left_join(irpef_prov, by="nuts3") %>% 
  left_join(urb_prov %>% 
              select(nuts3,pop,pop_density,dens_cat),
            by = "nuts3") %>% 
  left_join(tts_day, by=c("nuts2","date")) %>% 
  left_join(popres,by="nuts3")

for (i in names(full_day)) {
  cat(paste(paste("var:",i),
            paste("missing:",sum(is.na(full_day[i]))),
            "\n",
            sep = "\n"))
}

write_csv2(full_day,"dati/full.csv")

full <- read_csv2("dati/full.csv") %>% 
  mutate(date=ymd(date))


# ------------------------------------------------------------------------- #

#  that's all folks! ------------------------------------------------------

# ------------------------------------------------------------------------- #
