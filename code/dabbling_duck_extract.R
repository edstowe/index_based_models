library(readxl)
library(tidyverse)

nfs_sheets <- readxl::excel_sheets("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx") %>%
  as.data.frame()

# Dabbling Duck - North Ferry Slough
nfs_dd_sheets <- data.frame(
  sheet_name = readxl::excel_sheets("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx")) %>%
  filter(str_detect(sheet_name, "DD NFS"))

all_dd_nfs <- data.frame()

for(i in 1:nrow(nfs_dd_sheets)){
  sheet_i <- read_excel("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx", sheet = nfs_dd_sheets$sheet_name[i]) %>%
    rename(var_num = `...2`, var= `...3`, score= `...10`, comments= `...11`)
  vars <- select(sheet_i, var_num, var) %>% filter(!is.na(var_num))
  vals <- select(sheet_i, score) %>% filter(!is.na(score))
  alt <- sheet_i$comments[8]
  joined <- bind_cols(vars, vals[1:11,]) %>% 
    mutate(alternative = alt,
           sheet_name = nfs_dd_sheets$sheet_name[i])
  all_dd_nfs <- bind_rows(all_dd_nfs, joined)
}


dd_nfs_wide <-all_dd_nfs %>% select(-var_num) %>% pivot_wider(names_from = var, values_from = score)

all_dd_nfs %>%
  mutate(fwop = case_when(str_detect(alternative, "FWOP") ~ "FWOP", TRUE ~ "Not FWOP")) %>%
  ggplot()+
  geom_histogram(aes(score, fill = fwop))+
  facet_wrap(~var)

dd_times <- all_dd_nfs  %>% 
  mutate(alt = str_remove(sheet_name, "DD NFS "),
         year = case_when(str_detect(alt, "FWOP")~50,
                          str_detect(alt, "Existing")~0,
                          str_detect(alt, "TY1")~1,
                          str_detect(alt, "TY3")~3,
                          str_detect(alt, "TY50")~50),
         proj = case_when(alt %in% c("TY1", "TY3", "TY50") ~ "Main",
                          str_detect(alt, "(Base)")~"Base",
                          str_detect(alt, "(Wind)")~"Wind",
                          str_detect(alt, "(Islands)")~"Islands",
                          str_detect(alt, "(0.5 IsH9)")~"IsH9",
                          TRUE ~ alt)) %>%
  relocate(c(proj, year), .before = var_num) %>%
  select(-alt)

dd_proj <- filter(dd_times, !proj %in% c("FWOP", "Existing"))

pd <- position_dodge(4)

nfs <- ggplot(dd_proj, aes(year, score, color = proj))+
  geom_line(position = pd)+
  geom_point(position = pd)+
  facet_wrap(~var, scales = "free")+
  ggtitle("NFS Dabbling Duck")

# Dabbling Duck - McMillan
mcm_dd_sheets <- data.frame(
  sheet_name = readxl::excel_sheets("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx")) %>%
  filter(str_detect(sheet_name, "DD McMillan"))

all_dd_mcm <- data.frame()

for(i in 1:nrow(mcm_dd_sheets)){
  sheet_i <- read_excel("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx", sheet = mcm_dd_sheets$sheet_name[i]) %>%
    rename(var_num = `...2`, var= `...3`, score= `...10`, comments= `...11`)
  vars <- select(sheet_i, var_num, var) %>% filter(!is.na(var_num))
  vals <- select(sheet_i, score) %>% filter(!is.na(score))
  alt <- sheet_i$comments[8]
  joined <- bind_cols(vars, vals[1:11,]) %>% 
    mutate(alternative = alt,
           sheet_name = mcm_dd_sheets$sheet_name[i])
  all_dd_mcm <- bind_rows(all_dd_mcm, joined)
}

dd_times_mcm <- all_dd_mcm  %>% 
  filter(!sheet_name %in% c("DD McMillan FWOP (Base)","DD McMillan FWOP (Top)")) %>%
  mutate(alt = str_remove(sheet_name, "DD McMillan "),
         year = case_when(str_detect(alt, "FWOP")~50,
                          str_detect(alt, "Existing")~0,
                          str_detect(alt, "TY1")~1,
                          str_detect(alt, "TY3")~3,
                          str_detect(alt, "TY50")~50),
         proj = case_when(alt %in% c("TY1", "TY3", "TY50") ~ "Main",
                          str_detect(alt, "(Base)")~"Base",
                          str_detect(alt, "(SP4 Ext)")~"SP4 Ext",
                          str_detect(alt, "(Top)")~"Top",
                          TRUE ~ alt)) %>%
  relocate(c(proj, year), .before = var_num) %>%
  select(-alt)

dd_proj_mcm <- filter(dd_times_mcm, !proj %in% c("FWOP", "Existing"))

pd <- position_dodge(4)

mcm <- ggplot(dd_proj_mcm, aes(year, score, color = proj))+
  geom_line(position = pd)+
  geom_point(position = pd)+
  facet_wrap(~var, scales = "free")+
  ggtitle("McMillan Dabbling Duck")

# Dabbling Duck - South Ferry Slough
sfs_dd_sheets <- data.frame(
  sheet_name = readxl::excel_sheets("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx")) %>%
  filter(str_detect(sheet_name, "DD SFS"))

all_dd_sfs <- data.frame()

for(i in 1:nrow(sfs_dd_sheets)){
  sheet_i <- read_excel("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx", sheet = sfs_dd_sheets$sheet_name[i]) %>%
    rename(var_num = `...2`, var= `...3`, score= `...10`, comments= `...11`)
  vars <- select(sheet_i, var_num, var) %>% filter(!is.na(var_num))
  vals <- select(sheet_i, score) %>% filter(!is.na(score))
  alt <- sheet_i$comments[8]
  joined <- bind_cols(vars, vals[1:11,]) %>% 
    mutate(alternative = alt,
           sheet_name = sfs_dd_sheets$sheet_name[i])
  all_dd_sfs <- bind_rows(all_dd_sfs, joined)
}

dd_times_sfs <- all_dd_sfs  %>% 
  mutate(alt = str_remove(sheet_name, "DD SFS "),
         year = case_when(str_detect(alt, "FWOP")~50,
                          str_detect(alt, "Existing")~0,
                          str_detect(alt, "TY1")~1,
                          str_detect(alt, "TY3")~3,
                          str_detect(alt, "TY50")~50),
         proj = case_when(str_detect(alt, "AG3")~"AG3",
                          str_detect(alt, "AG4")~"AG4",
                          str_detect(alt, "AG5")~"AG5",
                          str_detect(alt, "AG6")~"AG6",
                          str_detect(alt, "AG7")~"AG7")) %>%
  relocate(c(proj, year), .before = var_num)

dd_proj_sfs <-dd_times_sfs %>% filter(!str_detect(alt, "FWOP"))

pd <- position_dodge(4)

sfs <- ggplot(dd_proj_sfs, aes(year, score, color = proj))+
  geom_line(position = pd)+
  geom_point(position = pd)+
  facet_wrap(~var, scales = "free")+
  ggtitle("SFS Dabbling Duck")



nfs
mcm
sfs


# Veery

nfs_veery_sheets <- data.frame(
  sheet_name = readxl::excel_sheets("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx")) %>%
  filter(str_detect(sheet_name, "Veery NFS"))

all_veery_nfs <- data.frame()

for(i in 1:nrow(nfs_veery_sheets)){
  
  sheet_i <- read_excel("LP10_Iteration_4_HEP_Analysis_May2021_Final.xlsx", 
                        sheet = nfs_veery_sheets$sheet_name[i]) %>%
    select(-c(`...3`,`...7`,`...8`))

  col_names <- sheet_i[7,] %>% unlist()
    
  sheet_i <- set_names(sheet_i, col_names)
  
  alt <- filter(sheet_i, Variable == "Enter Condition:") %>% pull(Description)
  
  vars <- filter(sheet_i, Variable %in% paste0("V", 1:6)) %>%
    mutate(alt = rep(alt, each = 6)) %>% relocate(alt, .before = Variable)
  
  all_veery_nfs <- bind_rows(all_veery_nfs, vars)
}

veery_times <- all_veery_nfs %>% 
  mutate(alt = str_remove(alt, "LP 10 - |LP 10 -"),
         alt = str_replace(alt, "North of Ferry Slough", "NFS"),
         year = as.numeric(str_extract(alt, "\\d+")),
         Alternative = str_replace_all(alt, "\\s*\\([^\\)]+\\)", "")) %>%
  relocate(c(Alternative, year), .before = Variable) %>%
  select(-alt) %>%
  filter(Variable != "V2") %>%
  mutate(Data = as.numeric(Data))

filter(veery_times, Variable == "V3") %>%
  ggplot()+
  geom_point(aes(year, as.numeric(Data)))+
  facet_wrap(~Alternative)

table(veery_times$Alternative, veery_times$year)

pd <- position_dodge(2)

veery_times %>%
  filter(Alternative %in% c("NFS Base FWP", "NFS FWP")) %>%
ggplot(aes(year, Data, color = Alternative))+
  geom_line(position = pd)+
  geom_point(position = pd)+
  facet_wrap(~Description, 
             labeller = label_wrap_gen(width = 20),
             scales = "free")







