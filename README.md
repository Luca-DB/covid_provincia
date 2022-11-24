# covid_provincia
### This repository hosts all datasets and codes to reproduce all analyses reported in the manuscript "SARS-CoV-2 spread and area deprivation in the Italian three-tier restrictions: a multilevel approach".
It contains the following folders and files:
1. **dati** (folder): it contains the datasets needed for the analyses:
   * **Codici-statistici-e-denominazioni-al-01_01_2020.xls**: data downloaded from the [Italian National Institute of Statistics][1] (ISTAT). The file contains the territorial codes, names, provinces (NUTS-3), regions (NUTS-2), and geographical repartition (NUTS-1) of all Italian municipalities as of 01/01/2020.
   * **Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2019.csv**: data downloaded from the [Italian Ministry of Economy and Finance][2] (MEF). The file contains information grouped at the municipality level about 2019 Italian taxpayers.
   * **Classificazioni statistiche-e-dimensione-dei-comuni_01_01_2020.xls**: data downloaded from [ISTAT][3]. The file contains land areas and populations of all Italian Municipalities as of 01/01/2020.
   * **tiers_nuts_complete.csv**: data gathered from [SKY TG24 news archive][4]. The file contains daily data about the Italian four-tier restriction system implemented from 06/11/2020 to 31/03/2022.
   * **popres_prov_01-01-20.csv**: data downloaded from [ISTAT][5]. The file contains the Italian resident population by province and age.
   * **covid_data_bck.csv**: data from the GitHub of the [Italian Civil Protection Department][6]. The file is generated from the R code "1_dataset.R" to provide an offline version of the data until 31/05/2021.
   * **full.csv**: dataset generated from the R code "1_dataset.R". It summarizes the previous datasets, containing only the data relevant to the analyses.
   
2. **plot** (folder): it contains the outputs of the R code "3_plots.R".

3. **1_dataset.R**: This R code loads, cleans, and merges all the data needed for the analyses. Its main output is the file "full.csv" that is loaded and used in the other two scripts.

4. **2_analysis.R**: This R code contains all the analyses of the study.

5. **3_plots.R**: This R code generates all the figures of the study and saves the output in the "plot" folder.

The three codes were made according to their numerical order (1_dataset.R > 2_analysis.R > 3_plots.R). However, they are independent of the others and work in any order.




[1]: https://www.istat.it/it/archivio/6789
[2]: https://www1.finanze.gov.it/finanze/analisi_stat/public/index.php?search_class[0]=cCOMUNE&opendata=yes
[3]: https://www.istat.it/it/archivio/156224
[4]: https://tg24.sky.it/archivio
[5]: https://demo.istat.it/popres/download.php?anno=2020&lingua=eng
[6]: https://github.com/pcm-dpc/COVID-19
