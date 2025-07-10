library(ggrepel)
library(readr)
library(dplyr)
library(tidyr)
library(knitr)
library(purrr)
library(fs)
library(tidyverse)
library(janitor)
library(lsr)
library(readxl)
library(lme4)
library(ggpubr)
library(viridis)
library(car)
library(cowplot)
library(ggh4x)
library(broom)
library(readr)
library(reshape2)
library(pheatmap)
library(httk)
library(tcpl)
library(data.table)
library(xlsx)
library(bio3d)

####Statistics----
quant_concordance_stats <- function(x,y,dataset){
  reg<-lm(log10(y)~log10(x))
  r2 = summary(reg)$r.squared
  df = summary(reg)$df[2]
  rmse = sqrt(mean(reg$residuals^2))
  r = cor(log10(y), log10(x), method='pearson',use="complete.obs")
  rho= cor(log10(y), log10(x), method='spearman',use="complete.obs")
  d<-(log10(y)-log10(x)) %>% na.omit()
  rmsd = sqrt(mean(d^2))
  absdif = mean(abs(log10(y)-log10(x)),na.rm=T)
  mad = median(abs(log10(y)-log10(x)),na.rm=T)
  signdif = mean(log10(x)-log10(y),na.rm=T)
  cat(paste("\nn = ", df+2, "\absdif = ", absdif, "\nsigndif = ", signdif, "\nrmsd = ", rmsd, "\nrho = ", rho, "\nr = ", r, "\nrsq = ", r2, "\nrmse = ", rmse))
}


###Read in data files----
####Nonclinical data (dataset #1)----
#process IUCLID file
df<-read_xlsx("IUCLID_RepeatedDoseToxicityOral.xlsx")%>%
  clean_names()%>%
  rename(species=materials_and_methods_test_animals_species_value)%>%
  filter(species %in% c("rat","mouse"))%>%
  rename(sample_size=materials_and_methods_administration_exposure_no_of_animals_per_sex_per_dose)%>%
  rename(critical_effect=results_and_discussion_effect_levels_efflevel_basis_remarks)%>%
  rename(route=materials_and_methods_administration_exposure_route_of_administration_value)%>%
  rename(studyID=document_key)%>%
  rename(design=materials_and_methods_administration_exposure_details_on_study_design)%>%
  rename(sex=results_and_discussion_effect_levels_efflevel_sex_value)%>%
  separate(route,c("route",NA))%>%
  rename(POD_type=results_and_discussion_effect_levels_efflevel_endpoint_value)%>%
  rename(lower_dose=results_and_discussion_effect_levels_efflevel_effect_level_lower_value)%>%
  mutate_at('lower_dose',as.numeric)%>%
  rename(upper_dose=results_and_discussion_effect_levels_efflevel_effect_level_upper_value)%>%
  filter(grepl("NDA",chemical_name_from_substance_names))%>%
  separate(chemical_name_from_substance_names,c(NA,"drug"),sep="-")%>%
  separate(materials_and_methods_administration_exposure_duration_of_treatment_exposure,c("duration_value","duration_units"))%>%
  rename(flag=administrative_data_purpose_flag_value)%>%
  rename(NOAELnotes= results_and_discussion_effect_levels_efflevel_basis_other)%>%
  rename(result_remarks= results_and_discussion_effect_levels_efflevel_remarks_on_results_other)%>%
  rename(dose_units=results_and_discussion_effect_levels_efflevel_effect_level_unit_code)%>%
  rename(hazard_category=results_and_discussion_effect_levels_efflevel_basis_value)%>%
  filter(!duration_value=="PND")%>%
  filter(!duration_units=="postnatal")%>%
  filter(!hazard_category=="histopathology: neoplastic")%>%
  mutate(duration_clean=case_when(duration_value=="twelve"~ "12",
                                  duration_value=="one" ~ "1",
                                  duration_value=="not" ~ "NR"))%>%
  mutate(duration_clean = case_when(!is.na(duration_clean) ~ duration_clean,
                                    TRUE ~ duration_value))%>%
  mutate_at('duration_clean',as.numeric) %>%
  mutate(days = case_when(duration_units=="day" ~ duration_clean,
                          duration_units=="days" ~ duration_clean,
                          duration_units=="weeks" ~ duration_clean*7,
                          duration_units=="week" ~ duration_clean*7,
                          duration_units=="Week" ~ duration_clean*7,
                          duration_units=="months" ~ duration_clean*30.42,
                          duration_units=="month" ~ duration_clean*30.42,
                          duration_units=="Month" ~ duration_clean*30.42,
                          duration_units=="Months" ~ duration_clean*30.42,
                          duration_units=="years" ~ duration_clean*365,
                          duration_units=="year" ~ duration_clean*365)) %>%
  mutate(duration_type = case_when(days >= 1 & days <= 30 ~ "short-term",
                                   days > 30 & days <= 90 ~ "subchronic",
                                   days > 90 ~ "chronic"))%>%
  mutate(HED_fda = case_when(species== "rat" ~ lower_dose/6.2,
                             species=="mouse" ~ lower_dose/12.3))%>%
  mutate(HED = case_when(species== "rat" ~ lower_dose*0.24,
                         species=="mouse" ~ lower_dose*0.13))%>%
  mutate(drug=tolower(drug))%>%
  mutate(drug=case_when(drug=="abiraterone"~"abiraterone acetate",
                        drug=="afatinib"~"afatinib dimaleate",
                        drug=="alogliptin"~"alogliptin benzoate",
                        drug=="bosentan hydrate"~"bosentan",
                        drug=="cariprazine"~"cariprazine hydrochloride",
                        drug=="citalopram"~"citalopram hydrobromide",
                        drug=="eliglustat"~"eliglustat tartrate",
                        drug=="lenvatinib"~"lenvatinib mesylate",
                        drug=="nintedanib esilate"~"nintedanib esylate",
                        drug=="osimertinib"~"osimertinib mesylate",
                        drug=="patiromer"~"patiromer sorbitex calcium",
                        drug=="simeprevir"~"simeprevir sodium",
                        drug=="topotecan"~"topotecan hydrochloride",
                        drug=="valbenazine"~"valbenazine tosylate",
                        drug=="vorapaxar"~"vorapaxar sulfate",
                        drug=="vortioxetine"~"vortioxetine hydrobromide",
                        drug=="cevimeline hcl"~"cevimeline hydrochloride",
                        drug=="ivabradine"~"ivabradine hydrochloride",
                        drug=="ruxolitinib"~"ruxolitinib phosphate",
                        drug=="tolterodine"~"tolterodine tartrate",
                        TRUE~drug))%>%
  filter(!POD_type =="NOEL")%>%
  filter(!POD_type=="NOAEL")%>%
  filter(species %in% c("rat","mouse"))%>%
  filter(dose_units=="mg/kg bw/day (actual dose received)")%>%
  rename(remarks=materials_and_methods_administration_exposure_doses_concentrations_entry_0_remarks)%>%
  filter(!studyID %in% c("461bafeb-6436-41aa-951a-788a851f09e9","84a20c05-26fc-45aa-a7fd-7ce399aa88f5","89ee60d7-1129-4f1d-8180-77950b50f64e","bd79e3a1-6b91-44ae-a215-840d4b94b5ae","c514e9fa-6189-4860-ae7c-c1b99db15b11","736c579d-a47f-410f-83ed-152f93b1fb81","59fd8801-e8fc-40c7-babd-bfb64597690e","bd79e3a1-6b91-44ae-a215-840d4b94b5ae","461bafeb-6436-41aa-951a-788a851f09e9","7b7edcaa-921b-42af-a02e-c9fd466c7229","c514e9fa-6189-4860-ae7c-c1b99db15b11","9a0f3c47-4f39-403d-b12a-0976104e4bdc","4a807546-a7ea-467a-8a95-47a3484ec84c"))%>% #filter for studies that include coexposures
  select(drug,species,route,POD_type,NOAELnotes,lower_dose,upper_dose,dose_units,HED,sex,duration_clean,duration_units,days,duration_type,critical_effect,hazard_category,sample_size,studyID,remarks,result_remarks)%>%
  rename(drugname=drug)%>%
  separate(drugname, sep = " ", into = c("drug", "drug_form"),remove=F)%>%
  mutate(drug=case_when(drugname=="sodium zirconium cyclosilicate"~"sodium zirconium",
                        drugname=="sodium gamma hydroxybutyrate"~"sodium hydroxybutyrate",
                        T~drug))%>%
  rename(effect=critical_effect)

#####Table S2 of clinical safety data----
#write.csv(df,"Supplemental Table S2 nonclinical safety dataset1.csv", row.names = F)

####Clinical dataset (dataset #2)----
##read in clinical safety data file
h<-read_xlsx("human_loael_clinical_dataset2.xlsx")%>%
  clean_names()%>%
  filter(monotherapy=="yes")%>%
  filter(placebo_controlled=="yes")%>%
  filter(repeat_dose=="yes")%>%
  mutate_at('adverse_effect_dose_mg_d', as.numeric)%>%
  mutate_at('median_weight_kg_of_treatment_group',as.numeric)%>%
  mutate_at('highest_dose_safety_tested_mg_d',as.numeric)%>%
  mutate(mg_kg_d_study=adverse_effect_dose_mg_d/median_weight_kg_of_treatment_group)%>%
  mutate(mg_kg_d = coalesce(mg_kg_d_study,adverse_effect_dose_mg_d/80))%>%
  #mutate(mg_kg_d_withNE = coalesce(mg_kg_d,highest_dose_safety_tested_mg_d/80))%>%
  rename(drugname=drug)%>%
  separate(drugname, sep = " ", into = c("drug", "drug_form"),remove=F)%>%
  mutate(drug=case_when(drugname=="sodium zirconium cyclosilicate"~"sodium zirconium",
                        drugname=="sodium oxybate"~"sodium oxybate",
                        T~drug))

longh<-h%>%
  mutate(effect=str_replace_all(adverse_effects," \\s*\\([^\\)]+\\)",""))%>%
  separate_longer_delim(effect, delim = ", ")%>%
  mutate(effect=tolower(effect))%>%
  mutate(effect = case_when(grepl("alt increase|increased alt|increase in alt|alt elevation|alt >1-3x uln|alt increase,|increase ALT|\\balt\\b|alanine aminotransferase|sgpt",ignore.case=T,effect) ~"ALT increase",
                            grepl("ast increase|aspartate aminotransferase|\\bast\\b|increase ast|increase in ast|increased ast|sgot",ignore.case=T,effect) ~ "AST increase",
                            grepl("weight decrease|decreased weight|weight decreased|weight loss",ignore.case = T,effect)~"weight decrease",
                            grepl("weight increase|weight gain|weight increased|increase in body weight",ignore.case = T,effect)~"weight increase",
                            grepl("bun|blood urea nitrogen",ignore.case = T,effect)~ "blood urea nitrogen change",
                            grepl("hemoglobin decrease|low hemoglobin|hemoglobin below|decrease hemoglobin|haemaglobin decrease|haemoglobin decrease|low haemoglobin",ignore.case = T,effect)~"hemoglobin decrease",
                            effect=="ne"~"no effects meet criteria",
                            grepl("bilirubin|hyperbili",ignore.case = T,effect)~"bilirubin increase",
                            TRUE ~ effect)) %>%
   select(drug,effect,mg_kg_d)

longh_multi<-h %>%
  filter(multiple_dose=="Y")%>% ###filter to select only multiple dose studies
  mutate(effect=str_replace_all(adverse_effects," \\s*\\([^\\)]+\\)",""))%>%
  separate_longer_delim(effect, delim = ", ")%>%
  mutate(effect=tolower(effect))%>%
  mutate(effect = case_when(grepl("alt increase|increased alt|increase in alt|alt elevation|alt >1-3x uln|alt increase,|increase ALT|\\balt\\b|alanine aminotransferase|sgpt",ignore.case=T,effect) ~"ALT increase",
                            grepl("ast increase|aspartate aminotransferase|\\bast\\b|increase ast|increase in ast|increased ast|sgot",ignore.case=T,effect) ~ "AST increase",
                            grepl("weight decrease|decreased weight|weight decreased|weight loss",ignore.case = T,effect)~"weight decrease",
                            grepl("weight increase|weight gain|weight increased|increase in body weight",ignore.case = T,effect)~"weight increase",
                            grepl("bun|blood urea nitrogen",ignore.case = T,effect)~ "blood urea nitrogen change",
                            grepl("hemoglobin decrease|low hemoglobin|hemoglobin below|decrease hemoglobin|haemaglobin decrease|haemoglobin decrease|low haemoglobin",ignore.case = T,effect)~"hemoglobin decrease",
                            effect=="ne"~"no effects meet criteria",
                            grepl("bilirubin|hyperbili",ignore.case = T,effect)~"bilirubin increase",
                            TRUE ~ effect)) %>%
   select(drug,effect,mg_kg_d)
hdrugs<-h%>%select(drug)%>%distinct()

#####Table S3 of clinical safety data----
#write.csv(longh,"Supplemental Table S3 clinical safety dataset1.csv", row.names = F)

####Clinical safety (dataset #3)----
boxed_all<-read.csv("chEMBL_blackbox_warnings_20240819.csv")%>%
  clean_names()%>%
  rename(drug=parent_molecule_name)%>%
  mutate(drug=tolower(drug))%>%
  rename(drugname=drug)%>%
  separate(drugname, sep = " ", into = c("drug", "drug_form"),remove=F)%>%
  select(drug,warning_class)%>%
  rename(effect=warning_class)%>%
  filter(!effect %in% c("None","misuse","teratogenicity","infectious disease","carcinogenicity"))

sider<-read_excel("SIDER_dataextract.xlsx") %>%
  rename(effect_original=effect)%>%
  mutate(effect=str_replace_all(effect_original," \\s*\\([^\\)]+\\)",""))%>%
  rename(drugname=drug)%>%
  separate(drugname, sep = " ", into = c("drug", "drug_form"),remove=F)%>%
  mutate(drug=case_when(drugname=="sodium zirconium cyclosilicate"~"sodium zirconium",
                        T~drug))%>%
  select(drug,effect,RR,dose_mg,p_value)%>%
  na.omit(RR)%>%
  filter(RR>1.0)%>%
  filter(p_value<0.05)%>%
  mutate(effect = case_when(grepl("alt increase|increased alt|increase in alt|alt elevation|alt >1-3x uln|alt increase,|increase ALT|\\balt\\b|alanine aminotransferase|sgpt",ignore.case=T,effect) ~"ALT increase",
                            grepl("ast increase|aspartate aminotransferase|\\bast\\b|increase ast|increase in ast|increased ast|sgot",ignore.case=T,effect) ~ "AST increase",
                            grepl("weight decrease|decreased weight|weight decreased|weight loss",ignore.case = T,effect)~"weight decrease",
                            grepl("weight increase|weight gain|weight increased|increase in body weight",ignore.case = T,effect)~"weight increase",
                            grepl("bun|blood urea nitrogen",ignore.case = T,effect)~ "blood urea nitrogen change",
                            grepl("hemoglobin decrease|low hemoglobin|hemoglobin below|decrease hemoglobin|haemaglobin decrease|haemoglobin decrease|low haemoglobin",ignore.case = T,effect)~"hemoglobin decrease",
                            grepl("bilirubin|hyperbili",ignore.case = T,effect)~"bilirubin increase",
                            TRUE ~ effect)) %>%
  mutate_at('dose_mg',as.numeric)%>%
  mutate(mg_kg_d = dose_mg/80)
sider_multidose<-read_excel("SIDER_dataextract.xlsx") %>%
  rename(effect_original=effect)%>%
  mutate(effect=str_replace_all(effect_original," \\s*\\([^\\)]+\\)",""))%>%
  rename(drugname=drug)%>%
  separate(drugname, sep = " ", into = c("drug", "drug_form"),remove=F)%>%
  mutate(drug=case_when(drugname=="sodium zirconium cyclosilicate"~"sodium zirconium",
                        T~drug))%>%
  group_by(drug) %>%
  mutate_at('dose_mg',as.numeric) %>%
  summarise(max=max(dose_mg),
            min=min(dose_mg),
            range=max-min)%>%
  mutate(multiple_dose= case_when(range>0 ~ "Y",
                                  T ~ "N"))%>%
  select(drug,multiple_dose)
sider_multi<-merge(sider,sider_multidose,by="drug")%>%
  filter(multiple_dose=="Y") ###select only drugs with multiple doses

clinset2_comb1<- lst(boxed_all, sider) %>%
  bind_rows(.id = "source")%>%
  distinct()%>%
  mutate(effect = case_when(effect=="psychiatric toxicity" ~ "neurotoxicity",
                            effect=="vascular toxicity" ~ "cardiovascular toxicity",
                            effect=="cardiotoxicity" ~ "cardiovascular toxicity",
                            TRUE ~ effect))%>%
  mutate(effect=tolower(effect))
sider_drugs<-sider%>%
  distinct(drug)
clinset2_comb<-merge(sider_drugs,clinset2_comb1,by="drug") #filter to drugs with dose data

clinset2_comb1_multi<- lst(boxed_all, sider_multi) %>%
  bind_rows(.id = "source")%>%
  distinct()%>%
  mutate(effect = case_when(effect=="psychiatric toxicity" ~ "neurotoxicity",
                            effect=="vascular toxicity" ~ "cardiovascular toxicity",
                            effect=="cardiotoxicity" ~ "cardiovascular toxicity",
                            TRUE ~ effect))%>%
  mutate(effect=tolower(effect))
clinset2_comb_multi<-merge(sider_drugs,clinset2_comb1_multi,by="drug") #filter to drugs with dose data

#####Table S4 of clinical safety data----
#write.csv(clinset2_comb,"Supplemental Table S4 clinical safety dataset2.csv", row.names=F)

###Hazard mapping----

######Table S5 all categorized hazards----
all_effects<-read.csv("hazard_assignments.csv")
all_effects_clin<- all_effects %>%
  filter(type=="clinical")%>%
  select(-c(n_clinset1,n_clinset2,n_clinsetnon))
all_effects_nonclin<- all_effects %>%
  filter(type=="nonclinical")%>%
  select(-c(n_clinset1,n_clinset2,n_clinsetnon))
#clinical
#dataset2
clinset1<-merge(longh,all_effects_clin,by="effect",all=T)%>%
  filter(hazard_sys_norm!="NA")
clinset1_multi<-merge(longh_multi,all_effects_clin,by="effect")
clinset1_drop_multi<- clinset1_multi %>%
  filter(is.na(category)|category!="incidental")
#dataset3
clinset2<-merge(clinset2_comb,all_effects_clin,by="effect")%>%
  filter(hazard_sys_norm!="NA")
clinset2_multi<-merge(clinset2_comb_multi,all_effects_clin,by="effect")
clinset2_drop_multi<-clinset2_multi %>%
  filter(is.na(category)|category!="incidental")
#nonclinical
nonclin<- merge(df,all_effects_nonclin,by="effect")

####Figure 2----
ceh<-clinset1%>%
  group_by(hazard_sys_norm)%>%
  summarise(drugs=n_distinct(drug))%>%
  na.omit()%>%
  filter(hazard_sys_norm %in% c("gastrointestinal","metabolic","cardiovascular","nervous","hematologic","hepatic","musculoskeletal","dermal","renal"))
ceh_drop_multi<-clinset1_drop_multi%>%
  group_by(hazard_sys_norm)%>%
  summarise(drugs=n_distinct(drug))%>%
  na.omit()%>%
  filter(hazard_sys_norm %in% c("gastrointestinal","metabolic","cardiovascular","nervous","hematologic","hepatic","musculoskeletal","dermal","renal"))
cehset2<-clinset2%>%
  group_by(hazard_sys_norm)%>%
  summarise(drugs=n_distinct(drug))%>%
  na.omit()%>%
  filter(hazard_sys_norm %in% c("gastrointestinal","metabolic","cardiovascular","nervous","hematologic","hepatic","musculoskeletal","dermal","renal"))
ceh_set2_drop_multi<-clinset2_drop_multi%>%
  group_by(hazard_sys_norm)%>%
  summarise(drugs=n_distinct(drug))%>%
  na.omit()%>%
  filter(hazard_sys_norm %in% c("gastrointestinal","metabolic","cardiovascular","nervous","hematologic","hepatic","musculoskeletal","dermal","renal"))
nonclin_m<-merge(hdrugs,nonclin,by="drug")
ceh_nonclin_rat<-nonclin_m%>%
  filter(species=="rat")%>%
  group_by(hazard_sys_norm)%>%
  summarise(drugs=n_distinct(drug))%>%
  na.omit()%>%
  filter(hazard_sys_norm %in% c("gastrointestinal","metabolic","cardiovascular","nervous","hematologic","hepatic","musculoskeletal","dermal","renal"))
ceh_nonclin_mouse<-nonclin_m%>%
  filter(species=="mouse")%>%
  group_by(hazard_sys_norm)%>%
  summarise(drugs=n_distinct(drug))%>%
  na.omit()%>%
  filter(hazard_sys_norm %in% c("gastrointestinal","metabolic","cardiovascular","nervous","hematologic","hepatic","musculoskeletal","dermal","renal"))

fig2data<-list(ceh, ceh_drop_multi,cehset2,ceh_set2_drop_multi,ceh_nonclin_rat,ceh_nonclin_mouse)
fig2data<-fig2data %>% reduce(inner_join,by='hazard_sys_norm')%>%
  rename(clinset1=drugs.x)%>%
  rename(clinset1_filt=drugs.y)%>%
  rename(clinset2=drugs.x.x)%>%
  rename(clinset2_filt=drugs.y.y)%>%
  rename(rat=drugs.x.x.x)%>%
  rename(mouse=drugs.y.y.y)%>%
  pivot_longer(cols=!hazard_sys_norm,names_to="dataset",values_to="drugs")

fig2data$hazard_sys_norm = factor(fig2data$hazard_sys_norm, levels = c("musculoskeletal","renal","hepatic","dermal","hematologic","cardiovascular","metabolic","nervous","gastrointestinal"), ordered = TRUE)

labs<-c("Nonclinical, Rat (#1)","Clinical (#2)", "Clinical (#3)","Nonclinical, Mouse (#1)","Clinical, filtered (#2.1)", "Clinical, filtered (#3.1)")
names(labs)<-c("rat","clinset1","clinset2","mouse","clinset1_filt","clinset2_filt")

fig2data %>%
  mutate(dataset=factor(dataset, levels=c("rat","clinset1","clinset2","mouse","clinset1_filt","clinset2_filt")))%>%
  ggplot(aes(x=hazard_sys_norm,y=drugs,fill=dataset))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+
  theme_minimal()+
  geom_text(aes(label=drugs), hjust=-0.1)+
  coord_flip()+
  theme(axis.title.y = element_blank())+
  facet_wrap(~dataset,
             labeller = labeller(dataset=labs))+
  ylim(0,165)+
  ylab("Count of drugs")+
  theme(text=element_text(size=15))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("lightgrey","#1fa187","#1fa187","lightgrey","#1fa187","#1fa187"))


###Pharmacokinetics dataset----
pk <- read_xlsx("cmax_auc_dataextraction.xlsx") %>%
  clean_names() %>%
  mutate_at('cmax', as.numeric) %>%
  mutate(cmax_ugml = case_when(cmax_units=="ug/mL" ~ cmax,
                               cmax_units=="ng/mL" ~ cmax*0.001,
                               cmax_units=="mg/mL" ~ cmax*1000,
                               cmax_units=="mg/L" ~ cmax,
                               cmax_units=="pg/mL" ~ cmax*0.000001,
                               cmax_units=="ug/L" ~ cmax*0.001)) %>%
  mutate_at('drug_dose', as.numeric) %>%
  mutate_at('auc_inf', as.numeric) %>%
  mutate(auc_mghml = case_when(auc_units=="mg*h/mL" ~ auc_inf,
                               auc_units=="ng*h/mL" ~ auc_inf*0.000001,
                               auc_units=="mg*h/L" ~ auc_inf*1000,
                               auc_units=="ug*h/mL" ~ auc_inf*0.001,
                               auc_units=="ug*h/L" ~ auc_inf*0.000001)) %>%
  mutate(drug_mgkg = case_when(species=="human" & drug_dose_units=="mg" ~ drug_dose/80,
                               species=="human" & drug_dose_units=="g" ~ (drug_dose*1000)/80,
                               species=="rat" & drug_dose_units=="mg/kg" ~ drug_dose,
                               species=="human" & drug_dose_units=="mg/kg" ~ drug_dose,
                               species=="mouse" & drug_dose_units=="mg/kg" ~ drug_dose))%>%
  mutate(cmax_norm = cmax_ugml/drug_mgkg)%>%
  mutate(auc_norm = auc_mghml/drug_mgkg)%>%
  group_by(drug,species)%>%
  summarise(cmax_norm_avg=mean(cmax_norm,na.rm=T),
            auc_norm_avg=mean(auc_norm,na.rm=T))%>%
  separate(drug, sep = " ", into = c("drug"))

###Bioactivity dataset----
invitrodb<-read_excel("concordance_invitrodb_v4_1.xlsx",sheet=1)
k<-read_excel("concordance_invitrodb_v4_1.xlsx",sheet=2)%>%
  clean_names()%>%
  filter(in_ad>0)%>%
  select(casrn)%>%
  distinct()%>%
  rename(casn=casrn)
k1<-read_excel("concordance_invitrodb_v4_1.xlsx",sheet=1)%>%
  group_by(casn)%>%
  tally()%>%
  filter(n>2)%>%
  select(casn)
k_all<- read_excel("concordance_invitrodb_v4_1.xlsx") %>%
  select(casn,ac50)%>%
  mutate_at('ac50',as.numeric)
k2<-merge(k_all,k,by="casn")
kfinal<-merge(k2,k1,by="casn")

kfinal %>%
  group_by(casn)%>%
  tally()%>%
  summarise(min=min(n),
            max=max(n),
            median=median(n))

p = c(0.05, 0.10, 0.25)

tcdf<- kfinal %>%
  group_by(casn) %>%
  summarise(q5 = quantile(ac50, probs = p[1], na.rm = T,type=7)) %>%
  na.omit()

#load predictions
load_sipes2017()
load_pradeep2020()
load_dawson2021()

####Estimate AED with httk----

#subset for drugs with httk data
invitro.dist.subset<-as.data.frame(subset(tcdf, casn %in% get_cheminfo(species='Human')))
length(unique(invitro.dist.subset$casn))
colnames(invitro.dist.subset)
colnames(invitro.dist.subset) <- c('chemcas','AC50p5')

#non-restrictive clearance
# aed.df=data.frame()
# 
# for (casn in unique(invitro.dist.subset[,'chemcas']))
# {
#   aed_50<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn,which.quantile=c(0.5),restrictive.clearance = F,output.units='mgpkgpday',species="Human",model='pbtk',Caco2.options = list(Caco2.Pab.default = 1.6, Caco2.Fabs = TRUE, Caco2.Fgut = TRUE, overwrite.invivo = FALSE, keepit100 = FALSE)))
#   aed_95<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn,which.quantile=c(0.95),restrictive.clearance = F,output.units='mgpkgpday',species="Human",model='pbtk',Caco2.options = list(Caco2.Pab.default = 1.6, Caco2.Fabs = TRUE, Caco2.Fgut = TRUE, overwrite.invivo = FALSE, keepit100 = FALSE)))
#   aed.df<-rbind(aed.df,cbind(casn, aed_50,aed_95))
# }

# colnames(aed.df)<-c('casn','AED50','AED95')
# new<- merge(aed.df,tcdf,by=("casn"))
# names<-read.csv("batchdownload.csv") %>%
#   separate(TOXCAST_NUMBER_OF_ASSAYS.TOTAL, sep = "/", c("positive","tested"))
# colnames(names)<-c('drug','casn','positive','tested')
# btemp<- merge(new,names,by="casn") %>%
#   distinct(casn,.keep_all=T) %>%
#   mutate(drug=tolower(drug))%>%
#   separate(drug, sep = " ", into = c("drug"))
# btemp$positive=as.numeric(btemp$positive)
# btemp$tested=as.numeric(btemp$tested)
# btemp2<- btemp %>%
#   mutate(peract = (positive/tested)*100)%>%
#   mutate(drug=tolower(drug))%>%
#   mutate_at('AED50', as.numeric) %>%
#   mutate_at('AED95', as.numeric)%>%
#   mutate_at('tested',as.numeric)%>%
#   filter(tested>300)
#write.csv(btemp2,"btemp2_20250103.csv")
btemp2<-read.csv("btemp2_20250103.csv")

###Create merged files for analyses----
q = 0.05

#5th%ile HED LOAEL rat and mouse all records
rodent_p5_study_smin<- df %>%
  group_by(studyID,drug,species)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug, species) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))%>%
  spread(species,q5)

#clinical dataset2 min human LOAEL all records
h5_set1<- longh %>%
  group_by(drug)%>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
p5_hed_s_set1 <-inner_join(rodent_p5_study_smin,h5_set1,by="drug")

#clinical dataset3 min human LOAEL all records
h5_set2<-clinset2_comb %>%
  group_by(drug)%>%
  summarise(humanmin = min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
p5_hed_s_set2 <-inner_join(rodent_p5_study_smin,h5_set2,by="drug")

#drop single dose drugs & incidental effects
h5_set1_drop_multi<- clinset1_drop_multi %>%
  group_by(drug)%>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
p5_hed_s_set1_drop_multi <-inner_join(rodent_p5_study_smin,h5_set1_drop_multi,by="drug")
h5_set2_drop_multi<- clinset2_drop_multi %>%
  group_by(drug)%>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
p5_hed_s_set2_drop_multi <-inner_join(rodent_p5_study_smin,h5_set2_drop_multi,by="drug")

#bioactivity
aed_set1<-merge(btemp2,h5_set1,by="drug")
aed_set2<-merge(btemp2,h5_set2,by="drug")
aed_set1_drop_multi<-merge(btemp2,h5_set1_drop_multi,by="drug")
aed_set2_drop_multi<-merge(btemp2,h5_set2_drop_multi,by="drug")
aed_nonclin<-merge(btemp2,rodent_p5_study_smin,by="drug")

#resulting data files
#p5_hed_s_set1 = 5th%ile rodent LOAEL & human LOAEL clinical dataset1
#p5_hed_s_set2 = 5th%ile rodent LOAEL & human LOAEL clinical dataset2
#p5_hed_s_set1_drop  = excludes clinical signs from 5th%ile rodent LOAEL & human LOAEL clinical dataset1
#p5_hed_s_set2_drop  = excludes clinical signs from 5th%ile rodent LOAEL & human LOAEL clinical dataset2
#aed_set1 = AED & human LOAEL clinical dataset1
#aed_set2 = AED & human LOAEL clinical dataset2
#aed_set1_drop = excludes clinical signs; AED & human LOAEL clinical dataset1
#aed_set2_drop = excludes clinical signs; AED & human LOAEL clinical dataset2

#PK adjusted
rodent_p5_study_smin_adm<- df %>%
  group_by(studyID,drug,species)%>%
  summarise(smin = min(lower_dose)) %>%
  group_by(drug, species) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))%>%
  spread(species,q5)
#clin dataset #2
p5_adm_s_nda <-inner_join(rodent_p5_study_smin_adm,h5_set1,by="drug")
df_g<-p5_adm_s_nda %>% rename(human=humanmin)%>%
  pivot_longer(cols=!drug,names_to="species",values_to="p5")
df_pk<- merge(df_g, pk, by=c("drug","species")) %>%
  mutate(cmax_corrected = cmax_norm_avg*p5)%>%
  select(drug,species,cmax_corrected)%>%
  spread(species,cmax_corrected)
df_pk_auc<- merge(df_g, pk, by=c("drug","species")) %>%
  mutate(auc_corrected = auc_norm_avg*p5)%>%
  select(drug,species,auc_corrected)%>%
  spread(species,auc_corrected)
#clin dataset #2.1
p5_adm_s_nda_2.2 <-inner_join(rodent_p5_study_smin_adm,h5_set1_drop_multi,by="drug")
df_g_2.2<-p5_adm_s_nda_2.2 %>% rename(human=humanmin)%>%
  pivot_longer(cols=!drug,names_to="species",values_to="p5")
df_pk_2.2<- merge(df_g_2.2, pk, by=c("drug","species")) %>%
  mutate(cmax_corrected = cmax_norm_avg*p5)%>%
  select(drug,species,cmax_corrected)%>%
  spread(species,cmax_corrected)
df_pk_auc_2.2<- merge(df_g_2.2, pk, by=c("drug","species")) %>%
  mutate(auc_corrected = auc_norm_avg*p5)%>%
  select(drug,species,auc_corrected)%>%
  spread(species,auc_corrected)
#clin dataset #3
p5_adm_s_nda_set2 <-inner_join(rodent_p5_study_smin_adm,h5_set2,by="drug")
df_g_set2<-p5_adm_s_nda_set2 %>% rename(human=humanmin)%>%
  pivot_longer(cols=!drug,names_to="species",values_to="p5")
df_pk_set2<- merge(df_g_set2, pk, by=c("drug","species")) %>%
  mutate(cmax_corrected = cmax_norm_avg*p5)%>%
  select(drug,species,cmax_corrected)%>%
  spread(species,cmax_corrected)
df_pk_auc_set2<- merge(df_g_set2, pk, by=c("drug","species")) %>%
  mutate(auc_corrected = auc_norm_avg*p5)%>%
  select(drug,species,auc_corrected)%>%
  spread(species,auc_corrected)
#rat to mouse
df_gr<-rodent_p5_study_smin_adm %>%
  pivot_longer(cols=!drug,names_to="species",values_to="p5")
df_pkr<- merge(df_gr, pk, by=c("drug","species")) %>%
  mutate(cmax_corrected = cmax_norm_avg*p5)%>%
  select(drug,species,cmax_corrected)%>%
  spread(species,cmax_corrected)

###Dataset summary values----
####Clinical dataset #2----
h %>% select(drug) %>% distinct() #214 drugs meet inclusion criteria
h %>% select(drug,adverse_effects) %>%
  filter(adverse_effects=="NE")%>% distinct() #38 drugs with "no effects" that meet criteria
clinset1%>%select(drug)%>%distinct()
longh%>%select(drug)%>%distinct()
clinset1_drop_multi%>%select(drug)%>%distinct()
#duration
dur<-h %>%
  separate(duration, sep = " ", into = c("duration_value", "duration_units"))%>%
  mutate_at('duration_value',as.numeric) %>%
  mutate(days = case_when(duration_units=="day" ~ duration_value,
                          duration_units=="days" ~ duration_value,
                          duration_units=="d" ~ duration_value,
                          duration_units=="w" ~ duration_value*7,
                          duration_units=="weeks" ~ duration_value*7,
                          duration_units=="week" ~ duration_value*7,
                          duration_units=="m" ~ duration_value*30.42,
                          duration_units=="months" ~ duration_value*30.42,
                          duration_units=="month" ~ duration_value*30.42,
                          duration_units=="years" ~ duration_value*365,
                          duration_units=="y" ~ duration_value*365,
                          duration_units=="year" ~ duration_value*365))
dur %>%
  summarise(median_duration=median(days,na.rm=T),
            min_dur=min(days,na.rm=T),
            max_dur=max(days,na.rm=T))
#treatment group size
sz<- h %>% mutate_at('treatment_group_size',as.numeric)
sz%>%
  summarise(median_sz=median(treatment_group_size,na.rm=T),
            min_sz=min(treatment_group_size,na.rm=T),
            max_sz=max(treatment_group_size,na.rm=T))
#sex
h %>% select(drug,sex)%>% distinct()%>%count(sex)
#multiple doses tested
h %>%
  filter(multiple_dose=="Y")%>%
  distinct(drug) #count
h %>%
  filter(multiple_dose=="Y")%>%
  mutate_at('highest_dose_safety_tested_mg_d', as.numeric)%>%
  mutate_at('lowest_dose_safety_tested_mg_d', as.numeric)%>%
  mutate(test_range= (highest_dose_safety_tested_mg_d)-(lowest_dose_safety_tested_mg_d))%>%
  distinct(drug,test_range)%>%
  summarise(median=median(test_range),
            min=min(test_range),
            max=max(test_range))

#multiple doses with varied adverse effects - dose range
range_c<- h %>%
  group_by(drug) %>%
  summarise(max=max(adverse_effect_dose_mg_d),
            min=min(adverse_effect_dose_mg_d),
            range=max-min)%>%
  #summarise(range = diff(range(mg_kg_d, na.rm = TRUE)))%>%
  filter(!range==0)%>%
  filter(!range=="-Inf")
range_c %>%
  summarise(median=median(range),
            min=min(range),
            max=max(range))
####Clinical dataset #3----
sider %>% select(drug) %>% distinct() #117 drugs meet inclusion criteria
#multiple doses - dose range
range_c<- sider %>%
  group_by(drug) %>%
  summarise(max=max(dose_mg),
            min=min(dose_mg),
            range=max-min)%>%
  #summarise(range = diff(range(mg_kg_d, na.rm = TRUE)))%>%
  filter(!range==0)%>%
  filter(!range=="-Inf")
range_c %>%
  summarise(median=median(range),
            min=min(range),
            max=max(range))
####Correlation between clinical datasets----
clincom<-merge(h5_set1,h5_set2,by="drug")
quant_concordance_stats(x = clincom$humanmin.x, y =  clincom$humanmin.y, climcom)
####Supplemental Figure S1----
clincom %>%
  ggplot(aes(x=log10(humanmin.x),y=log10(humanmin.y)))+
  geom_point()+
  theme_minimal()+
  geom_text_repel(aes(label=drug))+
  xlab(expression(Clinical ~ dataset ~ '#2' ~ min. ~ LOAEL  ~ log[10] ~ 'mg/kg-d'))+
  ylab(expression(Clinical ~ dataset ~ '#3' ~ min. ~ LOAEL  ~ log[10] ~ 'mg/kg-d'))+
  theme(text=element_text(size=15))+
  geom_smooth(method=lm, se=FALSE)+
  draw_label('r = 0.92, RMSE = 0.36 log10-mg/kg-day, MAD = 0.07', -2.5, 2, hjust = 0, vjust = 0)

h5_set1_multi<- clinset1_multi %>%
  group_by(drug)%>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
h5_set2_multi<- clinset2_multi %>%
  group_by(drug)%>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
clincom_multi<-merge(h5_set1_multi,h5_set2_multi,by="drug")
quant_concordance_stats(x = clincom_multi$humanmin.x, y =  clincom_multi$humanmin.y, climcom)


###Protective quantitative concordance analysis----

#clin set2 vs nonclin - all effects
quant_concordance_stats(x = p5_hed_s_set1$rat, y =  p5_hed_s_set1$humanmin)
quant_concordance_stats(x = p5_hed_s_set1$mouse, y =  p5_hed_s_set1$humanmin)
#clin set2.2 vs nonclin - multi dose only AND drop incidental
quant_concordance_stats(x = p5_hed_s_set1_drop_multi$rat, y =  p5_hed_s_set1_drop_multi$humanmin)
quant_concordance_stats(x = p5_hed_s_set1_drop_multi$mouse, y =  p5_hed_s_set1_drop_multi$humanmin)
#clin set3 vs nonclin - all effects
quant_concordance_stats(x = p5_hed_s_set2$rat, y =  p5_hed_s_set2$humanmin)
quant_concordance_stats(x = p5_hed_s_set2$mouse, y =  p5_hed_s_set2$humanmin)
#clin set3.1 vs nonclin - multi dose only AND drop incidental
quant_concordance_stats(x = p5_hed_s_set2_drop_multi$rat, y =  p5_hed_s_set2_drop_multi$humanmin)
quant_concordance_stats(x = p5_hed_s_set2_drop_multi$mouse, y =  p5_hed_s_set2_drop_multi$humanmin)
#clin set2 vs bioactivity - all effects
quant_concordance_stats(x = aed_set1$AED50, y =  aed_set1$humanmin)
#clin set2.1 vs nonclin - multi dose only AND drop incidental
quant_concordance_stats(x = aed_set1_drop_multi$AED50, y =  aed_set1_drop_multi$humanmin)
#clin set3 vs bioactivity - all effects
quant_concordance_stats(x = aed_set2$AED50, y =  aed_set2$humanmin)
#clin set3.1 vs bioactivity - multi dose only AND drop incidental
quant_concordance_stats(x = aed_set2_drop_multi$AED50, y =  aed_set2_drop_multi$humanmin)
#rat vs mouse
quant_concordance_stats(x = rodent_p5_study_smin$mouse, y =  rodent_p5_study_smin$rat)
#nonclin vs bioactivity
quant_concordance_stats(x = aed_nonclin$AED50, y =  aed_nonclin$rat)
quant_concordance_stats(x = aed_nonclin$AED50, y =  aed_nonclin$mouse)
#Cmax adjusted - clin set2
quant_concordance_stats(x = df_pk$rat, y =  df_pk$human, df_pk)
quant_concordance_stats(x = df_pk$mouse, y =  df_pk$human, df_pk)
#AUC adjusted - clin set2
quant_concordance_stats(x = df_pk_auc$rat, y =  df_pk_auc$human)
quant_concordance_stats(x = df_pk_auc$mouse, y =  df_pk_auc$human)
#Cmax adjusted - clin set2.1
quant_concordance_stats(x = df_pk_2.2$rat, y =  df_pk_2.2$human)
quant_concordance_stats(x = df_pk_2.2$mouse, y =  df_pk_2.2$human)
#AUC adjusted - clin set2.1
quant_concordance_stats(x = df_pk_auc_2.2$rat, y =  df_pk_auc_2.2$human)
quant_concordance_stats(x = df_pk_auc_2.2$mouse, y =  df_pk_auc_2.2$human)
#Cmax adjusted - clin set3
quant_concordance_stats(x = df_pk_set2$rat, y =  df_pk_set2$human)
quant_concordance_stats(x = df_pk_set2$mouse, y =  df_pk_set2$human)
#AUC adjusted - clin set3.1
quant_concordance_stats(x = df_pk_auc_set2$rat, y =  df_pk_auc_set2$human)
quant_concordance_stats(x = df_pk_auc_set2$mouse, y =  df_pk_auc_set2$human)
#Cmax rat to mouse
quant_concordance_stats(x = df_pkr$mouse, y =  df_pkr$rat, df_pk)

####Figures----
#####Human to rodent scatterplots----
p1<-ggplot() +
  geom_point(data = p5_hed_s_set1, aes(y=log10(humanmin), x=log10(rat),fill="rat"),color="#440154", size=3,alpha=0.6)+
  theme_minimal()+
  geom_smooth(data = p5_hed_s_set1, method=lm, aes(y=log10(humanmin), x=log10(rat)),se=F,color="#440154")+
  theme(text=element_text(size=12))+
  geom_smooth(data = p5_hed_s_set1, method=lm, aes(y=log10(humanmin), x=log10(mouse)),se=F,color="#26828e")+
  xlab(expression(Rodent ~ p5 ~ LOAEL[HED] ~ log[10] ~ `mg/kg-d`))+
  ylab(expression(Human ~ LOAEL ~ log[10] ~ 'mg/kg-d'))+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  geom_point(data = p5_hed_s_set1, aes(y=log10(humanmin), x=log10(mouse),fill="mouse"),color="#26828e", size=3,shape=17,alpha=0.6)+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(-2.3,2.8),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_y_continuous(limits=c(-2.8,2.3),breaks=c(-3,-2,-1,0,1,2,3))
p2<-ggplot() +
  geom_point(data = p5_hed_s_set1_drop_multi, aes(y=log10(humanmin), x=log10(rat),fill="rat"),color="#440154", size=3,alpha=0.6)+
  theme_minimal()+
  geom_smooth(data = p5_hed_s_set1_drop_multi, method=lm, aes(y=log10(humanmin), x=log10(rat)),se=F,color="#440154")+
  theme(text=element_text(size=12))+
  geom_smooth(data = p5_hed_s_set1_drop_multi, method=lm, aes(y=log10(humanmin), x=log10(mouse)),se=F,color="#26828e")+
  xlab(expression(Rodent ~ p5 ~ LOAEL[HED] ~ log[10] ~ `mg/kg-d`))+
  ylab(expression(Human ~ LOAEL ~ log[10] ~ 'mg/kg-d'))+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  geom_point(data = p5_hed_s_set1_drop_multi, aes(y=log10(humanmin), x=log10(mouse),fill="mouse"),color="#26828e", size=3,shape=17,alpha=0.6)+
  scale_fill_discrete(name = "Species")+
  scale_x_continuous(limits=c(-2.3,2.8),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_y_continuous(limits=c(-2.8,2.3),breaks=c(-3,-2,-1,0,1,2,3))
p3<-ggplot() +
  geom_point(data = p5_hed_s_set2, aes(y=log10(humanmin), x=log10(rat),fill="rat"),color="#440154", size=3,alpha=0.6)+
  theme_minimal()+
  geom_smooth(data = p5_hed_s_set2, method=lm, aes(y=log10(humanmin), x=log10(rat)),se=F,color="#440154")+
  theme(text=element_text(size=12))+
  geom_smooth(data = p5_hed_s_set2, method=lm, aes(y=log10(humanmin), x=log10(mouse)),se=F,color="#26828e")+
  xlab(expression(Rodent ~ p5 ~ LOAEL[HED] ~ log[10] ~ `mg/kg-d`))+
  ylab(expression(Human ~ LOAEL ~ log[10] ~ 'mg/kg-d'))+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  geom_point(data = p5_hed_s_set2, aes(y=log10(humanmin), x=log10(mouse),fill="mouse"),color="#26828e", size=3,shape=17,alpha=0.6)+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(-2.3,2.8),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_y_continuous(limits=c(-2.8,2.3),breaks=c(-3,-2,-1,0,1,2,3))
p4<-ggplot() +
  geom_point(data = p5_hed_s_set2_drop_multi, aes(y=log10(humanmin), x=log10(rat),fill="rat"),color="#440154", size=3,alpha=0.6)+
  theme_minimal()+
  geom_smooth(data = p5_hed_s_set2_drop_multi, method=lm, aes(y=log10(humanmin), x=log10(rat)),se=F,color="#440154")+
  theme(text=element_text(size=12))+
  geom_smooth(data = p5_hed_s_set2_drop_multi, method=lm, aes(y=log10(humanmin), x=log10(mouse)),se=F,color="#26828e")+
  xlab(expression(Rodent ~ p5 ~ LOAEL[HED] ~ log[10] ~ `mg/kg-d`))+
  ylab(expression(Human ~ LOAEL ~ log[10] ~ 'mg/kg-d'))+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  geom_point(data = p5_hed_s_set2_drop_multi, aes(y=log10(humanmin), x=log10(mouse),fill="mouse"),color="#26828e", size=3,shape=17,alpha=0.6)+
  scale_fill_discrete(name = "Species",guide = "none")+
  scale_x_continuous(limits=c(-2.3,2.8),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_y_continuous(limits=c(-2.8,2.3),breaks=c(-3,-2,-1,0,1,2,3))
######Figure 3----
prow<-plot_grid(p1+ theme(legend.position="none"),p2+ theme(legend.position="none"),p3+ theme(legend.position="none"),p4 + theme(legend.position="none"),align = "vh",nrow = 2,labels="AUTO")
legend <- get_legend(
  p2 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave2("scatterplots.png",width=9,height=7)

#####Bioactivity scatterplots----
p1<-aed_set1%>%
  ggplot(aes(y=log10(humanmin), x=log10(AED50))) +
  geom_point(size=3,color="#3b528b")+
  theme_minimal()+
  geom_smooth(method=lm, se=FALSE,color="#3b528b")+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  ylab(expression(Human ~ LOAEL ~ log[10] ~ 'mg/kg-d'))+
  xlab(expression(italic('in') ~ italic(vitro) ~ p5 ~ AED50 ~ log[10] ~ 'mg/kg-d'))+
  xlim(-3.5,3.5)+
  ylim(-3,3)+
  theme(text=element_text(size=15))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  draw_label('n=38', -2.5, 2, hjust = 0, vjust = 0)+
  scale_y_continuous(limits=c(-2.6,2.1),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_x_continuous(limits=c(-3,2.8),breaks=c(-3,-2,-1,0,1,2,3))

p2<-aed_set1_drop_multi%>%
  ggplot(aes(y=log10(humanmin), x=log10(AED50))) +
  geom_point(size=3,color="#3b528b")+
  theme_minimal()+
  geom_smooth(method=lm, se=FALSE,color="#3b528b")+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  ylab(expression(Human ~ LOAEL ~ log[10] ~ 'mg/kg-d'))+
  xlab(expression(italic('in') ~ italic(vitro) ~ p5 ~ AED50 ~ log[10] ~ 'mg/kg-d'))+
  xlim(-3.5,3.5)+
  ylim(-3,3)+
  theme(text=element_text(size=15))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  draw_label('n=21', -2.5, 2, hjust = 0, vjust = 0)+
  scale_y_continuous(limits=c(-2.6,2.1),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_x_continuous(limits=c(-3,2.8),breaks=c(-3,-2,-1,0,1,2,3))

p3<-ggplot() +
  geom_point(data = aed_nonclin, aes(y=log10(rat), x=log10(AED50),fill="rat"),color="#440154", size=3,alpha=0.6)+
  theme_minimal()+
  geom_smooth(data = aed_nonclin, method=lm, aes(y=log10(rat), x=log10(AED50)),se=F,color="#440154")+
  theme(text=element_text(size=15))+
  geom_smooth(data = aed_nonclin, method=lm, aes(y=log10(mouse), x=log10(AED50)),se=F,color="#26828e")+
  xlab(expression(italic('in') ~ italic(vitro) ~ p5 ~ AED50 ~ log[10] ~ 'mg/kg-d'))+
  ylab(expression(p5 ~ Rodent ~ LOAEL[HED] ~ log[10] ~ `mg/kg-d`))+
  geom_abline(intercept = 0, slope = 1,col="red",linetype="dashed")+
  geom_point(data = aed_nonclin, aes(y=log10(mouse), x=log10(AED50),fill="mouse"),color="#26828e", size=3,shape=17,alpha=0.6)+
  scale_fill_discrete(name = "Species")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  draw_label('n=14 (mouse)', -2.5, 2, hjust = 0, vjust = 0)+
  draw_label('n=26 (rat)', -2.5, 1.7, hjust = 0, vjust = 0)+
  scale_x_continuous(limits=c(-3,2.8),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_y_continuous(limits=c(-2.6,2.1),breaks=c(-3,-2,-1,0,1,2,3))
######Figure 4----
prow<-plot_grid(p1,p2,p3+ theme(legend.position="none"),align = "vh",nrow = 1,labels="AUTO")
legend <- get_legend(
  p3 + theme(legend.box.margin = margin(0, 0, 0, 8))
)
plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave2("bioactivity_scatterplots.png",width=17,height=5)

###Bland-Altman plot----
ba<- p5_hed_s_set1 %>% select(drug,rat,humanmin)%>% na.omit()
#log10
ba$diff <- log10(ba$humanmin) - log10(ba$rat)
mean_diff <- mean(ba$diff)
mean_diff
ba$avg <- (log10(ba$humanmin) + log10(ba$rat)) / 2
lower <- mean_diff - 1.96*sd(ba$diff)
upper <- mean_diff + 1.96*sd(ba$diff)
p1<-ggplot(ba, aes(x = avg, y = diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ylab("Difference Between Measurements") +
  xlab("Average Measurement")+
  theme_minimal()
ba<- p5_hed_s_set1 %>% select(drug,mouse,humanmin)%>% na.omit()
#log10
ba$diff <- log10(ba$humanmin) - log10(ba$mouse)
mean_diff <- mean(ba$diff)
mean_diff
ba$avg <- (log10(ba$humanmin) + log10(ba$mouse)) / 2
lower <- mean_diff - 1.96*sd(ba$diff)
upper <- mean_diff + 1.96*sd(ba$diff)
p2<-ggplot(ba, aes(x = avg, y = diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ylab("Difference Between Measurements") +
  xlab("Average Measurement")+
  theme_minimal()
ba<- p5_hed_s_set2 %>% select(drug,rat,humanmin)%>% na.omit()
#log10
ba$diff <- log10(ba$humanmin) - log10(ba$rat)
mean_diff <- mean(ba$diff)
mean_diff
ba$avg <- (log10(ba$humanmin) + log10(ba$rat)) / 2
lower <- mean_diff - 1.96*sd(ba$diff)
upper <- mean_diff + 1.96*sd(ba$diff)
p3<-ggplot(ba, aes(x = avg, y = diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ylab("Difference Between Measurements") +
  xlab("Average Measurement")+
  theme_minimal()
ba<- p5_hed_s_set2 %>% select(drug,mouse,humanmin)%>% na.omit()
#log10
ba$diff <- log10(ba$humanmin) - log10(ba$mouse)
mean_diff <- mean(ba$diff)
mean_diff
ba$avg <- (log10(ba$humanmin) + log10(ba$mouse)) / 2
lower <- mean_diff - 1.96*sd(ba$diff)
upper <- mean_diff + 1.96*sd(ba$diff)
p4<-ggplot(ba, aes(x = avg, y = diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ylab("Difference Between Measurements") +
  xlab("Average Measurement")+
  theme_minimal()
#bioactivity
ba_aed<- aed_set1 %>% select(drug,humanmin,AED50)%>%na.omit()
ba_aed$diff <- log10(ba_aed$humanmin) - log10(ba_aed$AED50)
mean_diff <- mean(ba_aed$diff)
mean_diff
ba_aed$avg <- (log10(ba_aed$humanmin) + log10(ba_aed$AED50)) / 2
lower <- mean_diff - 1.96*sd(ba_aed$diff)
upper <- mean_diff + 1.96*sd(ba_aed$diff)
p5<-ggplot(ba_aed, aes(x = avg, y = diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean_diff) +
  geom_hline(yintercept = lower, color = "red", linetype="dashed") +
  geom_hline(yintercept = upper, color = "red", linetype="dashed") +
  ylab("Difference Between Measurements") +
  xlab("Average Measurement")+
  theme_minimal()
######Supplemental Figure S2----
plot_grid(p1,p2,p3,p4,p5,labels="AUTO",nrow=1)

###Predictive quantitative concordance analysis----
####Clinical dataset #2 All effects and drugs----
#####Metabolic---
#matching human~rat
rat_met <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="metabolic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_met <- clinset1 %>%
  filter(hazard_sys_norm=="metabolic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmet<-merge(rat_met,h_met,by="drug")
quant_concordance_stats(x = hmet$rp5, y =  hmet$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmet_h<-merge(r_all,h_met,by="drug")
quant_concordance_stats(x = hmet_h$q5, y =  hmet_h$humanmin)

#####Gastrointestinal---
#matching human~rat
rat_gas <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="gastrointestinal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_gas <- clinset1 %>%
  filter(hazard_sys_norm=="gastrointestinal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hgast<-merge(rat_gas,h_gas,by="drug")
quant_concordance_stats(x = hgast$rp5, y =  hgast$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hgast_h<-merge(r_all,h_gas,by="drug")
quant_concordance_stats(x = hgast_h$q5, y =  hgast_h$humanmin)

#####Cardiovascular---
#matching human~rat
rat_car <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="cardiovascular")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_car <- clinset1 %>%
  filter(hazard_sys_norm=="cardiovascular") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hcar<-merge(rat_car,h_car,by="drug")
quant_concordance_stats(x = hcar$rp5, y =  hcar$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hcar_h<-merge(r_all,h_car,by="drug")
quant_concordance_stats(x = hcar_h$q5, y =  hcar_h$humanmin)

#####Nervous---
#matching human~rat
rat_ner <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="nervous")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_ner <- clinset1 %>%
  filter(hazard_sys_norm=="nervous") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hner<-merge(rat_ner,h_ner,by="drug")
quant_concordance_stats(x = hner$rp5, y =  hner$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hner_h<-merge(r_all,h_ner,by="drug")
quant_concordance_stats(x = hner_h$q5, y =  hner_h$humanmin)

#####Hematologic---
#matching human~rat
rat_hem <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="hematologic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hem <- clinset1 %>%
  filter(hazard_sys_norm=="hematologic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhem<-merge(rat_hem,h_hem,by="drug")
quant_concordance_stats(x = hhem$rp5, y =  hhem$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhem_h<-merge(r_all,h_hem,by="drug")
quant_concordance_stats(x = hhem_h$q5, y =  hhem_h$humanmin)

#####Hepatic---
#matching human~rat
rat_hep <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="hepatic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hep <- clinset1 %>%
  filter(hazard_sys_norm=="hepatic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhr<-merge(rat_hep,h_hep,by="drug")
quant_concordance_stats(x = hhr$rp5, y =  hhr$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhr_h<-merge(r_all,h_hep,by="drug")
quant_concordance_stats(x = hhr_h$q5, y =  hhr_h$humanmin)

#####Dermal---
#matching human~rat
rat_derm <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="dermal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_derm <- clinset1 %>%
  filter(hazard_sys_norm=="dermal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hderm<-merge(rat_derm,h_derm,by="drug")
quant_concordance_stats(x = hderm$rp5, y =  hderm$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hderm_h<-merge(r_all,h_derm,by="drug")
quant_concordance_stats(x = hderm_h$q5, y =  hderm_h$humanmin)

#####Musculoskeletal---
#matching human~rat
rat_musc <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="musculoskeletal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_musc <- clinset1 %>%
  filter(hazard_sys_norm=="musculoskeletal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmusc<-merge(rat_musc,h_musc,by="drug")
quant_concordance_stats(x = hmusc$rp5, y =  hmusc$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmusc_h<-merge(r_all,h_musc,by="drug")
quant_concordance_stats(x = hmusc_h$q5, y =  hmusc_h$humanmin)

#####Renal---
#matching human~rat
rat_renal <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="renal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_renal <- clinset1 %>%
  filter(hazard_sys_norm=="renal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hrenal<-merge(rat_renal,h_renal,by="drug")
quant_concordance_stats(x = hrenal$rp5, y =  hrenal$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hrenal_h<-merge(r_all,h_renal,by="drug")
quant_concordance_stats(x = hrenal_h$q5, y =  hrenal_h$humanmin)

####Clinical dataset #2.1 Filtered----
#####Metabolic---
#matching human~rat
rat_met <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="metabolic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_met <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="metabolic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmet<-merge(rat_met,h_met,by="drug")
quant_concordance_stats(x = hmet$rp5, y =  hmet$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmet_h<-merge(r_all,h_met,by="drug")
quant_concordance_stats(x = hmet_h$q5, y =  hmet_h$humanmin)

#####Gastrointestinal---
#matching human~rat
rat_gas <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="gastrointestinal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_gas <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="gastrointestinal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hgast<-merge(rat_gas,h_gas,by="drug")
quant_concordance_stats(x = hgast$rp5, y =  hgast$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hgast_h<-merge(r_all,h_gas,by="drug")
quant_concordance_stats(x = hgast_h$q5, y =  hgast_h$humanmin)

#####Cardiovascular---
#matching human~rat
rat_car <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="cardiovascular")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_car <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="cardiovascular") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hcar<-merge(rat_car,h_car,by="drug")
quant_concordance_stats(x = hcar$rp5, y =  hcar$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hcar_h<-merge(r_all,h_car,by="drug")
quant_concordance_stats(x = hcar_h$q5, y =  hcar_h$humanmin)

#####Nervous---
#matching human~rat
rat_ner <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="nervous")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_ner <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="nervous") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hner<-merge(rat_ner,h_ner,by="drug")
quant_concordance_stats(x = hner$rp5, y =  hner$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hner_h<-merge(r_all,h_ner,by="drug")
quant_concordance_stats(x = hner_h$q5, y =  hner_h$humanmin)

#####Hematologic---
#matching human~rat
rat_hem <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="hematologic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hem <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="hematologic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhem<-merge(rat_hem,h_hem,by="drug")
quant_concordance_stats(x = hhem$rp5, y =  hhem$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhem_h<-merge(r_all,h_hem,by="drug")
quant_concordance_stats(x = hhem_h$q5, y =  hhem_h$humanmin)

#####Hepatic---
#matching human~rat
rat_hep <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="hepatic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hep <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="hepatic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhr<-merge(rat_hep,h_hep,by="drug")
quant_concordance_stats(x = hhr$rp5, y =  hhr$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhr_h<-merge(r_all,h_hep,by="drug")
quant_concordance_stats(x = hhr_h$q5, y =  hhr_h$humanmin)

#####Dermal---
#matching human~rat
rat_derm <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="dermal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_derm <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="dermal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hderm<-merge(rat_derm,h_derm,by="drug")
quant_concordance_stats(x = hderm$rp5, y =  hderm$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hderm_h<-merge(r_all,h_derm,by="drug")
quant_concordance_stats(x = hderm_h$q5, y =  hderm_h$humanmin)

#####Musculoskeletal---
#matching human~rat
rat_musc <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="musculoskeletal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_musc <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="musculoskeletal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmusc<-merge(rat_musc,h_musc,by="drug")
quant_concordance_stats(x = hmusc$rp5, y =  hmusc$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmusc_h<-merge(r_all,h_musc,by="drug")
quant_concordance_stats(x = hmusc_h$q5, y =  hmusc_h$humanmin)

#####Renal---
#matching human~rat
rat_renal <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="renal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_renal <- clinset1_drop_multi %>%
  filter(hazard_sys_norm=="renal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hrenal<-merge(rat_renal,h_renal,by="drug")
quant_concordance_stats(x = hrenal$rp5, y =  hrenal$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hrenal_h<-merge(r_all,h_renal,by="drug")
quant_concordance_stats(x = hrenal_h$q5, y =  hrenal_h$humanmin)

####Predictive concordance bar plots----
dat_r <- data.frame(
  Not_matched = c(0.43,0.39,0.55,0.36,0.46,0.37,0.53,0.61,0.54),
  Matched = c(0.39,0.5,0.58,0.38,0.33,0.39,0.51,0.5,0.77),
  sys = as.factor(c("Gastrointestinal","Nervous","Metabolic","Cardiovascular","Hematologic","Dermal","Hepatic","Renal","Musculoskeletal"))
)
dat_long_r <- dat_r %>%
  gather("Data", "value", -sys)
dat_long_r$n <- c(78,77,55,41,36,34,36,26,25,63,49,55,18,30,20,33,24,12)
p1<-dat_long_r%>%
  ggplot(aes(x = factor(sys, levels=c("Gastrointestinal","Nervous","Metabolic","Cardiovascular","Hematologic","Dermal","Hepatic","Renal","Musculoskeletal")), y = value, fill = Data)) +
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("lightgrey", "#1fa187"),labels = c("Matched","Not matched"))+
  theme_minimal()+
  theme(axis.title.x=element_blank())+
  ylab("r")+
  theme(text=element_text(size=15))+
  ylim(0,0.9)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.title=element_blank())
#RMSD
dat <- data.frame(
  Not_matched = c(1,1,0.97,1,00.98,1,0.99,0.83,0.94),
  Matched = c(1.4,1.4,0.8,1.3,1.1,0.98,1.2,0.96,1.1),
  sys = as.factor(c("Gastrointestinal","Nervous","Metabolic","Cardiovascular","Hematologic","Dermal","Hepatic","Renal","Musculoskeletal"))
)
dat_long <- dat %>%
  gather("Data", "value", -sys)
dat_long$n <- c(78,77,55,41,36,34,36,26,25,63,49,55,18,30,20,33,24,12)
p2<-ggplot(dat_long, aes(x = factor(sys, levels=c("Gastrointestinal","Nervous","Metabolic","Cardiovascular","Hematologic","Dermal","Hepatic","Renal","Musculoskeletal")), y = value, fill = Data)) +
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("lightgrey", "#1fa187"),labels = c("Matched","Not matched"))+
  theme_minimal()+
  theme(axis.title.x=element_blank())+
  ylab("RMSD")+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.title=element_blank())
#####Figure 5----
plot_grid(p1,p2,align = "v",nrow = 2,labels="AUTO")

###Qualitative concordance----
####Clinical dataset #2 All effects and drugs----
#####Hepatic---
hr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhep=any(r_hepatic=="Y"))%>%
  select(drug,rhep)%>%
  distinct()
hh<- clinset1 %>%
  mutate(h_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhep=any(h_hepatic=="Y"))%>%
  select(drug,humanhep)%>%
  distinct()
hepmerge<-merge(hh,hr,by="drug")
hepsum<- hepmerge %>%
  mutate(hepatic = case_when(rhep=="TRUE" & humanhep=="TRUE"~ "TP",
                             rhep=="TRUE" & humanhep=="FALSE" ~ "FP",
                             rhep=="FALSE" & humanhep=="TRUE" ~ "FN",
                             rhep=="FALSE" & humanhep=="FALSE" ~ "TN"))
hepsum %>% count(hepatic)
Hepatic<-hepsum %>%
  count(hepatic)%>%
  spread(hepatic,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Cardio---
cr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rcar=any(r_cardio=="Y"))%>%
  select(drug,rcar)%>%
  distinct()
hc<- clinset1 %>%
  mutate(h_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humancar=any(h_cardio=="Y"))%>%
  select(drug,humancar)%>%
  distinct()
carmerge<-merge(hc,cr,by="drug")
carsum<- carmerge %>%
  mutate(cardio = case_when(rcar=="TRUE" & humancar=="TRUE"~ "TP",
                            rcar=="TRUE" & humancar=="FALSE" ~ "FP",
                            rcar=="FALSE" & humancar=="TRUE" ~ "FN",
                            rcar=="FALSE" & humancar=="FALSE" ~ "TN"))
Cardiovascular<-carsum %>%
  count(cardio)%>%
  spread(cardio,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Renal---
rer <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rren=any(r_renal=="Y"))%>%
  select(drug,rren)%>%
  distinct()
hren<- clinset1 %>%
  mutate(h_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanren=any(h_renal=="Y"))%>%
  select(drug,humanren)%>%
  distinct()
renmerge<-merge(hren,rer,by="drug")
rensum<- renmerge %>%
  mutate(renal = case_when(rren=="TRUE" & humanren=="TRUE"~ "TP",
                           rren=="TRUE" & humanren=="FALSE" ~ "FP",
                           rren=="FALSE" & humanren=="TRUE" ~ "FN",
                           rren=="FALSE" & humanren=="FALSE" ~ "TN"))
rensum %>% count(renal)
Renal<- rensum %>%
  count(renal)%>%
  spread(renal,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Hemato---
hemr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhem=any(r_hem=="Y"))%>%
  select(drug,rhem)%>%
  distinct()
hhem<- clinset1 %>%
  mutate(h_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhem=any(h_hem=="Y"))%>%
  select(drug,humanhem)%>%
  distinct()
hemmerge<-merge(hemr,hhem,by="drug")
hemsum<- hemmerge %>%
  mutate(hemat = case_when(rhem=="TRUE" & humanhem=="TRUE"~ "TP",
                           rhem=="TRUE" & humanhem=="FALSE" ~ "FP",
                           rhem=="FALSE" & humanhem=="TRUE" ~ "FN",
                           rhem=="FALSE" & humanhem=="FALSE" ~ "TN"))
hemsum %>% count(hemat)
Hematologic<-hemsum %>%
  count(hemat)%>%
  spread(hemat,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Metabolic---
metr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmet=any(r_met=="Y"))%>%
  select(drug,rmet)%>%
  distinct()
hmet<- clinset1 %>%
  mutate(h_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmet=any(h_met=="Y"))%>%
  select(drug,humanmet)%>%
  distinct()
metmerge<-merge(metr,hmet,by="drug")
metsum<- metmerge %>%
  mutate(met = case_when(rmet=="TRUE" & humanmet=="TRUE"~ "TP",
                         rmet=="TRUE" & humanmet=="FALSE" ~ "FP",
                         rmet=="FALSE" & humanmet=="TRUE" ~ "FN",
                         rmet=="FALSE" & humanmet=="FALSE" ~ "TN"))
metsum %>% count(met)
Metabolic<- metsum %>%
  count(met)%>%
  spread(met,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Psychiatric---
psyr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rpsy=any(r_psy=="Y"))%>%
  select(drug,rpsy)%>%
  distinct()
hpsy<- clinset1 %>%
  mutate(h_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanpsy=any(h_psy=="Y"))%>%
  select(drug,humanpsy)%>%
  distinct()
psymerge<-merge(psyr,hpsy,by="drug")
psysum<- psymerge %>%
  mutate(neuro = case_when(rpsy=="TRUE" & humanpsy=="TRUE"~ "TP",
                           rpsy=="TRUE" & humanpsy=="FALSE" ~ "FP",
                           rpsy=="FALSE" & humanpsy=="TRUE" ~ "FN",
                           rpsy=="FALSE" & humanpsy=="FALSE" ~ "TN"))
psysum %>% count(neuro)
Nervous<- psysum %>%
  count(neuro)%>%
  spread(neuro,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Gastrointestinal---
gasr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rgas=any(r_gas=="Y"))%>%
  select(drug,rgas)%>%
  distinct()
hgas<- clinset1 %>%
  mutate(h_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humangas=any(h_gas=="Y"))%>%
  select(drug,humangas)%>%
  distinct()
gasmerge<-merge(gasr,hgas,by="drug")
gassum<- gasmerge %>%
  mutate(gas = case_when(rgas=="TRUE" & humangas=="TRUE"~ "TP",
                         rgas=="TRUE" & humangas=="FALSE" ~ "FP",
                         rgas=="FALSE" & humangas=="TRUE" ~ "FN",
                         rgas=="FALSE" & humangas=="FALSE" ~ "TN"))
gassum %>% count(gas)
Gastrointestinal<- gassum %>%
  count(gas)%>%
  spread(gas,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Musculo---
musr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmus=any(r_mus=="Y"))%>%
  select(drug,rmus)%>%
  distinct()
hmus<- clinset1 %>%
  mutate(h_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmus=any(h_mus=="Y"))%>%
  select(drug,humanmus)%>%
  distinct()
musmerge<-merge(musr,hmus,by="drug")
mussum<- musmerge %>%
  mutate(mus = case_when(rmus=="TRUE" & humanmus=="TRUE"~ "TP",
                         rmus=="TRUE" & humanmus=="FALSE" ~ "FP",
                         rmus=="FALSE" & humanmus=="TRUE" ~ "FN",
                         rmus=="FALSE" & humanmus=="FALSE" ~ "TN"))
mussum %>% count(mus)
Musculoskeletal<- mussum %>%
  count(mus)%>%
  spread(mus,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Derm---
dermr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rderm=any(r_derm=="Y"))%>%
  select(drug,rderm)%>%
  distinct()
hderm<- clinset1 %>%
  mutate(h_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanderm=any(h_derm=="Y"))%>%
  select(drug,humanderm)%>%
  distinct()
dermmerge<-merge(dermr,hderm,by="drug")
dermsum<- dermmerge %>%
  mutate(derm = case_when(rderm=="TRUE" & humanderm=="TRUE"~ "TP",
                          rderm=="TRUE" & humanderm=="FALSE" ~ "FP",
                          rderm=="FALSE" & humanderm=="TRUE" ~ "FN",
                          rderm=="FALSE" & humanderm=="FALSE" ~ "TN"))
dermsum %>% count(derm)
Dermal<- dermsum %>%
  count(derm)%>%
  spread(derm,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

######Qualitative summary table----
qualtable_clinset1<- lst(Gastrointestinal,Nervous,Metabolic,Cardiovascular,Hematologic,Dermal,Hepatic,Renal,Musculoskeletal) %>%
  bind_rows(.id = "hazard_sys_norm")

#####Qualitative bar plots----
qualtable_clinset1_pred<- lst(Gastrointestinal,Nervous,Metabolic,Cardiovascular,Hematologic,Dermal,Hepatic,Renal,Musculoskeletal) %>%
  bind_rows(.id = "hazard_sys_norm")%>%
  select(hazard_sys_norm,PPV,NPV)%>%
  gather("data", "value", -hazard_sys_norm)
p1<-ggplot(qualtable_clinset1_pred, aes(x =factor(hazard_sys_norm, levels=c("Gastrointestinal","Nervous","Metabolic","Cardiovascular","Hematologic","Dermal","Hepatic","Renal","Musculoskeletal")), y = value, fill = data)) +
  geom_col(position = "dodge")+
  theme_minimal()+
  theme(text=element_text(size=18))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  #geom_text(aes(label=n),position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Percent")+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("lightgrey", "#1fa187"))+
  ylim(0,100)
qualtable_clinset1_lr<- lst(Gastrointestinal,Nervous,Metabolic,Cardiovascular,Hematologic,Dermal,Hepatic,Renal,Musculoskeletal) %>%
  bind_rows(.id = "hazard_sys_norm")%>%
  select(hazard_sys_norm,LRneg,LRpos)%>%
  gather("data", "value", -hazard_sys_norm)
p2<-ggplot(qualtable_clinset1_lr, aes(x =factor(hazard_sys_norm, levels=c("Gastrointestinal","Nervous","Metabolic","Cardiovascular","Hematologic","Dermal","Hepatic","Renal","Musculoskeletal")), y = value, fill = data)) +
  geom_col(position = "dodge")+
  theme_minimal()+
  theme(text=element_text(size=18))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  #geom_text(aes(label=n),position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Ratio")+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("lightgrey", "#1fa187"),labels=c('iLR-','LR+'))+
  ylim(0,3.5)
#####Figure 6----
plot_grid(p1,p2,align = "v",nrow = 2,labels="AUTO")

####Clinical dataset #2 Filtered----
#####Hepatic---
hr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhep=any(r_hepatic=="Y"))%>%
  select(drug,rhep)%>%
  distinct()
hh<- clinset1_drop_multi %>%
  mutate(h_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhep=any(h_hepatic=="Y"))%>%
  select(drug,humanhep)%>%
  distinct()
hepmerge<-merge(hh,hr,by="drug")
hepsum<- hepmerge %>%
  mutate(hepatic = case_when(rhep=="TRUE" & humanhep=="TRUE"~ "TP",
                             rhep=="TRUE" & humanhep=="FALSE" ~ "FP",
                             rhep=="FALSE" & humanhep=="TRUE" ~ "FN",
                             rhep=="FALSE" & humanhep=="FALSE" ~ "TN"))
hepsum %>% count(hepatic)
hepatic<-hepsum %>%
  count(hepatic)%>%
  spread(hepatic,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Cardio---
cr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rcar=any(r_cardio=="Y"))%>%
  select(drug,rcar)%>%
  distinct()
hc<- clinset1_drop_multi %>%
  mutate(h_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humancar=any(h_cardio=="Y"))%>%
  select(drug,humancar)%>%
  distinct()
carmerge<-merge(hc,cr,by="drug")
carsum<- carmerge %>%
  mutate(cardio = case_when(rcar=="TRUE" & humancar=="TRUE"~ "TP",
                            rcar=="TRUE" & humancar=="FALSE" ~ "FP",
                            rcar=="FALSE" & humancar=="TRUE" ~ "FN",
                            rcar=="FALSE" & humancar=="FALSE" ~ "TN"))
cardio<-carsum %>%
  count(cardio)%>%
  spread(cardio,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Renal---
rer <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rren=any(r_renal=="Y"))%>%
  select(drug,rren)%>%
  distinct()
hren<- clinset1_drop_multi %>%
  mutate(h_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanren=any(h_renal=="Y"))%>%
  select(drug,humanren)%>%
  distinct()
renmerge<-merge(hren,rer,by="drug")
rensum<- renmerge %>%
  mutate(renal = case_when(rren=="TRUE" & humanren=="TRUE"~ "TP",
                           rren=="TRUE" & humanren=="FALSE" ~ "FP",
                           rren=="FALSE" & humanren=="TRUE" ~ "FN",
                           rren=="FALSE" & humanren=="FALSE" ~ "TN"))
rensum %>% count(renal)
nephro<- rensum %>%
  count(renal)%>%
  spread(renal,n)%>%
  mutate(prev= (TP+0)/(0+FP+TN+TP),
         sens= (TP/(TP+0)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+0))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Hemato---
hemr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhem=any(r_hem=="Y"))%>%
  select(drug,rhem)%>%
  distinct()
hhem<- clinset1_drop_multi %>%
  mutate(h_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhem=any(h_hem=="Y"))%>%
  select(drug,humanhem)%>%
  distinct()
hemmerge<-merge(hemr,hhem,by="drug")
hemsum<- hemmerge %>%
  mutate(hemat = case_when(rhem=="TRUE" & humanhem=="TRUE"~ "TP",
                           rhem=="TRUE" & humanhem=="FALSE" ~ "FP",
                           rhem=="FALSE" & humanhem=="TRUE" ~ "FN",
                           rhem=="FALSE" & humanhem=="FALSE" ~ "TN"))
hemsum %>% count(hemat)
hemat<-hemsum %>%
  count(hemat)%>%
  spread(hemat,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Metabolic---
metr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmet=any(r_met=="Y"))%>%
  select(drug,rmet)%>%
  distinct()
hmet<- clinset1_drop_multi %>%
  mutate(h_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmet=any(h_met=="Y"))%>%
  select(drug,humanmet)%>%
  distinct()
metmerge<-merge(metr,hmet,by="drug")
metsum<- metmerge %>%
  mutate(met = case_when(rmet=="TRUE" & humanmet=="TRUE"~ "TP",
                         rmet=="TRUE" & humanmet=="FALSE" ~ "FP",
                         rmet=="FALSE" & humanmet=="TRUE" ~ "FN",
                         rmet=="FALSE" & humanmet=="FALSE" ~ "TN"))
metsum %>% count(met)
metab<- metsum %>%
  count(met)%>%
  spread(met,n)%>%
  mutate(prev= (TP+0)/(0+FP+TN+TP),
         sens= (TP/(TP+0)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+0))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Psychiatric---
psyr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rpsy=any(r_psy=="Y"))%>%
  select(drug,rpsy)%>%
  distinct()
hpsy<- clinset1_drop_multi %>%
  mutate(h_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanpsy=any(h_psy=="Y"))%>%
  select(drug,humanpsy)%>%
  distinct()
psymerge<-merge(psyr,hpsy,by="drug")
psysum<- psymerge %>%
  mutate(neuro = case_when(rpsy=="TRUE" & humanpsy=="TRUE"~ "TP",
                           rpsy=="TRUE" & humanpsy=="FALSE" ~ "FP",
                           rpsy=="FALSE" & humanpsy=="TRUE" ~ "FN",
                           rpsy=="FALSE" & humanpsy=="FALSE" ~ "TN"))
psysum %>% count(neuro)
psych<- psysum %>%
  count(neuro)%>%
  spread(neuro,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Gastrointestinal---
gasr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rgas=any(r_gas=="Y"))%>%
  select(drug,rgas)%>%
  distinct()
hgas<- clinset1_drop_multi %>%
  mutate(h_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humangas=any(h_gas=="Y"))%>%
  select(drug,humangas)%>%
  distinct()
gasmerge<-merge(gasr,hgas,by="drug")
gassum<- gasmerge %>%
  mutate(gas = case_when(rgas=="TRUE" & humangas=="TRUE"~ "TP",
                         rgas=="TRUE" & humangas=="FALSE" ~ "FP",
                         rgas=="FALSE" & humangas=="TRUE" ~ "FN",
                         rgas=="FALSE" & humangas=="FALSE" ~ "TN"))
gassum %>% count(gas)
gastro<- gassum %>%
  count(gas)%>%
  spread(gas,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Musculo---
musr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmus=any(r_mus=="Y"))%>%
  select(drug,rmus)%>%
  distinct()
hmus<- clinset1_drop_multi %>%
  mutate(h_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmus=any(h_mus=="Y"))%>%
  select(drug,humanmus)%>%
  distinct()
musmerge<-merge(musr,hmus,by="drug")
mussum<- musmerge %>%
  mutate(mus = case_when(rmus=="TRUE" & humanmus=="TRUE"~ "TP",
                         rmus=="TRUE" & humanmus=="FALSE" ~ "FP",
                         rmus=="FALSE" & humanmus=="TRUE" ~ "FN",
                         rmus=="FALSE" & humanmus=="FALSE" ~ "TN"))
mussum %>% count(mus)
musculo<- mussum %>%
  count(mus)%>%
  spread(mus,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Derm---
dermr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rderm=any(r_derm=="Y"))%>%
  select(drug,rderm)%>%
  distinct()
hderm<- clinset1_drop_multi %>%
  mutate(h_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanderm=any(h_derm=="Y"))%>%
  select(drug,humanderm)%>%
  distinct()
dermmerge<-merge(dermr,hderm,by="drug")
dermsum<- dermmerge %>%
  mutate(derm = case_when(rderm=="TRUE" & humanderm=="TRUE"~ "TP",
                          rderm=="TRUE" & humanderm=="FALSE" ~ "FP",
                          rderm=="FALSE" & humanderm=="TRUE" ~ "FN",
                          rderm=="FALSE" & humanderm=="FALSE" ~ "TN"))
dermsum %>% count(derm)
derm<- dermsum %>%
  count(derm)%>%
  spread(derm,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Qualitative summary table----
qualtable_clinset1_drop_multi <- lst(gastro,metab,cardio,psych,hemat,hepatic,derm,musculo,nephro) %>%
  bind_rows(.id = "hazard_sys_norm")

###CDF plots----
rodent_p5_admin<- df %>%
  group_by(studyID,drug,species)%>%
  summarise(smin = min(lower_dose)) %>%
  group_by(drug, species) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))%>%
  spread(species,q5)
p5_hed_admin <-inner_join(rodent_p5_admin,h5_set1,by="drug")

cdf<-p5_hed_s_set1 %>%
  mutate(HEDrat= log10(rat)-log10(humanmin))%>%
  mutate(HEDmouse=log10(mouse)-log10(humanmin))%>%
  #mutate(HEDmouse_1000x = log10(mouse/1000)-log10(humanmin))%>%
  select(drug, HEDrat,HEDmouse)
cdf_ad<-p5_hed_admin %>%
  mutate(admin_doserat= log10(rat)-log10(humanmin)) %>%
  mutate(admin_dosemouse=log10(mouse)-log10(humanmin))%>%
  select(admin_doserat,admin_dosemouse,drug)
cdf_rh<-merge(cdf_ad,cdf,by="drug")%>%
  pivot_longer(cols=!drug,names_to="type",values_to="value")

p1<-ggplot(cdf_rh, aes(x=value,color=type)) +
  stat_ecdf(geom = "point",alpha=0.6)+
  theme_minimal()+
  xlab(expression(log[10] ~ Ratio ~ p5 ~ 'Rodent/Human' ~ LOAEL))+
  geom_vline(xintercept = 0,linetype="dashed",color="red")+
  theme(text=element_text(size=11))+
  theme(legend.title=element_blank())+
  scale_x_continuous(limits=c(-3,3),breaks=c(-3,-2,-1,0,1,2,3))+
  ylab("Cumulative Frequency")+
  scale_color_viridis_d(labels=c('Mouse p5 LOAEL','Rat p5 LOAEL',expression(Mouse ~ p5 ~ LOAEL[HED]),expression(Rat ~ p5 ~ LOAEL[HED])))+
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow = 2)) +
  geom_vline(xintercept=2.6,color="#440154",linetype="dashed")+
  geom_vline(xintercept=2.38,color="#31688e",linetype="dashed")+
  geom_vline(xintercept=2.24,color="#35b779",linetype="dashed")+
  geom_vline(xintercept=1.89,color="#fde725",linetype="dashed")

cdf_iv<-aed_set1 %>% select(drug,humanmin, AED50, AED95)%>%
  mutate(AED_50= log10(AED50)-log10(humanmin))%>%
  mutate(AED_95=log10(AED95)-log10(humanmin))%>%
  #mutate(AED_95X=log10(AED95/10)-log10(humanmin))%>%
  select(drug,AED_50,AED_95)%>%
  pivot_longer(cols=!drug,names_to="type",values_to="value")

p2<-ggplot(cdf_iv, aes(x=value,color=type)) +
  stat_ecdf(geom = "point")+
  theme_minimal()+
  theme(text=element_text(size=11))+
  xlab(expression(log[10] ~ Ratio ~ p5 ~ 'AED/Human' ~ LOAEL))+
  geom_vline(xintercept = 0,linetype="dashed",color="red")+
  theme(legend.title=element_blank())+
  scale_color_manual(values=c("#35b779","#440154"),labels=c('p5 AED50', 'p5 AED95'))+
  scale_x_continuous(limits=c(-3,3),breaks=c(-3,-2,-1,0,1,2,3))+
  ylab("Cumulative Frequency")+
  theme(legend.position="bottom")+
  geom_vline(xintercept=1.2,color="#440154",linetype="dashed")+
  geom_vline(xintercept=1.92,color="#35b779",linetype="dashed")

cdf_iv_r<-aed_nonclin %>% select(drug,rat,mouse, AED50)%>%
  mutate(AED50_ratHED= log10(AED50)-log10(rat))%>%
  mutate(AED50_mouseHED=log10(AED50)-log10(mouse))%>%
  #mutate(AED_95X=log10(AED95/10)-log10(humanmin))%>%
  select(drug,AED50_ratHED,AED50_mouseHED)%>%
  pivot_longer(cols=!drug,names_to="type",values_to="value")

p3<-ggplot(cdf_iv_r, aes(x=value,color=type)) +
  stat_ecdf(geom = "point")+
  theme_minimal()+
  theme(text=element_text(size=11))+
  xlab(expression(log[10] ~ Ratio ~ p5 ~ 'AED50/p5' ~ Rodent ~ LOAEL))+
  geom_vline(xintercept = 0,linetype="dashed",color="red")+
  theme(legend.title=element_blank())+
  scale_color_manual(values=c("#35b779","#440154"),labels=c(expression(Mouse ~ p5 ~ LOAEL[HED]),expression(Rat ~ p5 ~ LOAEL[HED])))+
  scale_x_continuous(limits=c(-3,3),breaks=c(-3,-2,-1,0,1,2,3))+
  ylab("Cumulative Frequency")+
  theme(legend.position="bottom")+
  geom_vline(xintercept=2.1,color="#440154",linetype="dashed")+
  geom_vline(xintercept=1.58,color="#35b779",linetype="dashed")

#####Figure 7----
plot_grid(p1,p2,p3,align = "vh",labels="AUTO",nrow=1)
ggsave2("CDF_plots.png",width=9,height=3)

###Supplemental analyses----

####Class----
class<-read_xlsx("Evangelisti supp1.xlsx")%>%
  clean_names()%>%
  select(active_ingredient, atc_anatomical_class,pharmacological_class)%>%
  rename(drugname=active_ingredient)%>%
  mutate(drugname=tolower(drugname))%>%
  separate(drugname, sep = " ", into = c("drug", "drug_form"),remove=F)%>%
  mutate(drug=case_when(drugname=="sodium zirconium cyclosilicate"~"sodium zirconium",
                        drugname=="sodium oxybate"~"sodium oxybate",
                        T~drug))
p5_hed_s_set1_anti<-merge(p5_hed_s_set1,class,by="drug")%>%
  filter(atc_anatomical_class=="ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS")
p5_hed_s_set1_nonanti<-merge(p5_hed_s_set1,class,by="drug")%>%
  filter(atc_anatomical_class!="ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS")
p5_hed_s_set2_anti<-merge(p5_hed_s_set2,class,by="drug")%>%
  filter(atc_anatomical_class=="ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS")
p5_hed_s_set2_nonanti<-merge(p5_hed_s_set2,class,by="drug")%>%
  filter(atc_anatomical_class!="ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS")
#clin set2 vs nonclin - antineoplastics
quant_concordance_stats(x = p5_hed_s_set1_anti$rat, y =  p5_hed_s_set1_anti$humanmin)
quant_concordance_stats(x = p5_hed_s_set1_anti$mouse, y =  p5_hed_s_set1_anti$humanmin)
#clin set2 vs nonclin - non-antineoplastics
quant_concordance_stats(x = p5_hed_s_set1_nonanti$rat, y =  p5_hed_s_set1_nonanti$humanmin)
quant_concordance_stats(x = p5_hed_s_set1_nonanti$mouse, y =  p5_hed_s_set1_nonanti$humanmin)
#clin set3 vs nonclin - antineoplastics
quant_concordance_stats(x = p5_hed_s_set2_anti$rat, y =  p5_hed_s_set2_anti$humanmin)
quant_concordance_stats(x = p5_hed_s_set2_anti$mouse, y =  p5_hed_s_set2_anti$humanmin)
#clin set3 vs nonclin - non-antineoplastics
quant_concordance_stats(x = p5_hed_s_set2_nonanti$rat, y =  p5_hed_s_set2_nonanti$humanmin)
quant_concordance_stats(x = p5_hed_s_set2_nonanti$mouse, y =  p5_hed_s_set2_nonanti$humanmin)
#bioactivity
aed_set1_anti<-merge(aed_set1,class,by="drug")%>%
  filter(atc_anatomical_class=="ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS")
quant_concordance_stats(x = aed_set1_anti$AED50, y =  aed_set1_anti$humanmin)

####Duration----
rodent_p5_ch<- df %>%
  filter(duration_type=="chronic")%>%
  group_by(drug, species) %>%
  summarise(quant05 = quantile(HED, probs = q[1], na.rm = T,type=7)) %>%
  spread(species, quant05)
p5_hed_s_set1_chronic <-inner_join(rodent_p5_ch,h5_set1,by="drug")
#clin set1 vs nonclin - chronic
quant_concordance_stats(x = p5_hed_s_set1_chronic$rat, y =  p5_hed_s_set1_chronic$humanmin)
quant_concordance_stats(x = p5_hed_s_set1_chronic$mouse, y =  p5_hed_s_set1_chronic$humanmin)
p5_hed_s_set2_chronic <-inner_join(rodent_p5_ch,h5_set2,by="drug")
#clin set1 vs nonclin - chronic
quant_concordance_stats(x = p5_hed_s_set2_chronic$rat, y =  p5_hed_s_set2_chronic$humanmin)
quant_concordance_stats(x = p5_hed_s_set2_chronic$mouse, y =  p5_hed_s_set2_chronic$humanmin)

####Percentile - median----
q = c(0.05,0.5)
rodent_p5_study_med<- df %>%
  group_by(studyID,drug,species)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug, species) %>%
  summarise(q5 = quantile(smin, probs = q[2], na.rm = T,type=7))%>%
  spread(species,q5)
#clinical dataset1 human LOAEL all records
h5_set1_med<- clinset1 %>% #clinset1_drop_multi
  group_by(drug) %>%
  summarise(q5 = quantile(mg_kg_d, probs = q[2], na.rm = T,type=7))
p5_hed_s_set1_med <-inner_join(rodent_p5_study_med,h5_set1_med,by="drug")
#clin set1 vs nonclin - median
quant_concordance_stats(x = p5_hed_s_set1_med$rat, y =  p5_hed_s_set1_med$q5)
quant_concordance_stats(x = p5_hed_s_set1_med$mouse, y =  p5_hed_s_set1_med$q5)
#clinical dataset2 human LOAEL all records
h5_set2_med<- clinset2 %>% #clinset2_drop_multi
  group_by(drug) %>%
  summarise(q5 = quantile(mg_kg_d, probs = q[2], na.rm = T,type=7))
p5_hed_s_set2_med <-inner_join(rodent_p5_study_med,h5_set2_med,by="drug")
#clin set2 vs nonclin - median
quant_concordance_stats(x = p5_hed_s_set2_med$rat, y =  p5_hed_s_set2_med$q5)
quant_concordance_stats(x = p5_hed_s_set2_med$mouse, y =  p5_hed_s_set2_med$q5)

###Predictive quantitative concordance analysis----
####Clinical dataset #2 All effects and drugs for mouse----
#####Metabolic---
#matching human~mouse
rat_met <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="metabolic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_met <- clinset1 %>%
  filter(hazard_sys_norm=="metabolic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmet<-merge(rat_met,h_met,by="drug")
quant_concordance_stats(x = hmet$rp5, y =  hmet$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmet_h<-merge(r_all,h_met,by="drug")
quant_concordance_stats(x = hmet_h$q5, y =  hmet_h$humanmin)

#####Gastrointestinal---
#matching human~mouse
rat_gas <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="gastrointestinal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_gas <- clinset1 %>%
  filter(hazard_sys_norm=="gastrointestinal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hgast<-merge(rat_gas,h_gas,by="drug")
quant_concordance_stats(x = hgast$rp5, y =  hgast$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hgast_h<-merge(r_all,h_gas,by="drug")
quant_concordance_stats(x = hgast_h$q5, y =  hgast_h$humanmin)

#####Cardiovascular---
#matching human~mouse
rat_car <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="cardiovascular")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_car <- clinset1 %>%
  filter(hazard_sys_norm=="cardiovascular") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hcar<-merge(rat_car,h_car,by="drug")
quant_concordance_stats(x = hcar$rp5, y =  hcar$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hcar_h<-merge(r_all,h_car,by="drug")
quant_concordance_stats(x = hcar_h$q5, y =  hcar_h$humanmin)

#####Nervous---
#matching human~mouse
rat_ner <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="nervous")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_ner <- clinset1 %>%
  filter(hazard_sys_norm=="nervous") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hner<-merge(rat_ner,h_ner,by="drug")
quant_concordance_stats(x = hner$rp5, y =  hner$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hner_h<-merge(r_all,h_ner,by="drug")
quant_concordance_stats(x = hner_h$q5, y =  hner_h$humanmin)

#####Hematologic---
#matching human~mouse
rat_hem <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="hematologic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hem <- clinset1 %>%
  filter(hazard_sys_norm=="hematologic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhem<-merge(rat_hem,h_hem,by="drug")
quant_concordance_stats(x = hhem$rp5, y =  hhem$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhem_h<-merge(r_all,h_hem,by="drug")
quant_concordance_stats(x = hhem_h$q5, y =  hhem_h$humanmin)

#####Hepatic---
#matching human~mouse
rat_hep <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="hepatic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hep <- clinset1 %>%
  filter(hazard_sys_norm=="hepatic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhr<-merge(rat_hep,h_hep,by="drug")
quant_concordance_stats(x = hhr$rp5, y =  hhr$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhr_h<-merge(r_all,h_hep,by="drug")
quant_concordance_stats(x = hhr_h$q5, y =  hhr_h$humanmin)

#####Dermal---
#matching human~mouse
rat_derm <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="dermal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_derm <- clinset1 %>%
  filter(hazard_sys_norm=="dermal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hderm<-merge(rat_derm,h_derm,by="drug")
quant_concordance_stats(x = hderm$rp5, y =  hderm$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hderm_h<-merge(r_all,h_derm,by="drug")
quant_concordance_stats(x = hderm_h$q5, y =  hderm_h$humanmin)

#####Musculoskeletal---
#matching human~mouse
rat_musc <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="musculoskeletal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_musc <- clinset1 %>%
  filter(hazard_sys_norm=="musculoskeletal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmusc<-merge(rat_musc,h_musc,by="drug")
quant_concordance_stats(x = hmusc$rp5, y =  hmusc$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmusc_h<-merge(r_all,h_musc,by="drug")
quant_concordance_stats(x = hmusc_h$q5, y =  hmusc_h$humanmin)

#####Renal---
#matching human~mouse
rat_renal <- nonclin %>%
  filter(species=="mouse")%>%
  filter(hazard_sys_norm=="renal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_renal <- clinset1 %>%
  filter(hazard_sys_norm=="renal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hrenal<-merge(rat_renal,h_renal,by="drug")
quant_concordance_stats(x = hrenal$rp5, y =  hrenal$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="mouse")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hrenal_h<-merge(r_all,h_renal,by="drug")
quant_concordance_stats(x = hrenal_h$q5, y =  hrenal_h$humanmin)

####Clinical dataset #3 All effects and drugs----
#####Metabolic---
#matching human~rat
rat_met <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="metabolic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_met <- clinset2 %>%
  filter(hazard_sys_norm=="metabolic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmet<-merge(rat_met,h_met,by="drug")
quant_concordance_stats(x = hmet$rp5, y =  hmet$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmet_h<-merge(r_all,h_met,by="drug")
quant_concordance_stats(x = hmet_h$q5, y =  hmet_h$humanmin)

#####Gastrointestinal---
#matching human~rat
rat_gas <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="gastrointestinal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_gas <- clinset2 %>%
  filter(hazard_sys_norm=="gastrointestinal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hgast<-merge(rat_gas,h_gas,by="drug")
quant_concordance_stats(x = hgast$rp5, y =  hgast$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hgast_h<-merge(r_all,h_gas,by="drug")
quant_concordance_stats(x = hgast_h$q5, y =  hgast_h$humanmin)

#####Cardiovascular---
#matching human~rat
rat_car <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="cardiovascular")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_car <- clinset2 %>%
  filter(hazard_sys_norm=="cardiovascular") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hcar<-merge(rat_car,h_car,by="drug")
quant_concordance_stats(x = hcar$rp5, y =  hcar$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hcar_h<-merge(r_all,h_car,by="drug")
quant_concordance_stats(x = hcar_h$q5, y =  hcar_h$humanmin)

#####Nervous---
#matching human~rat
rat_ner <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="nervous")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_ner <- clinset2 %>%
  filter(hazard_sys_norm=="nervous") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hner<-merge(rat_ner,h_ner,by="drug")
quant_concordance_stats(x = hner$rp5, y =  hner$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hner_h<-merge(r_all,h_ner,by="drug")
quant_concordance_stats(x = hner_h$q5, y =  hner_h$humanmin)

#####Hematologic---
#matching human~rat
rat_hem <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="hematologic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hem <- clinset2 %>%
  filter(hazard_sys_norm=="hematologic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhem<-merge(rat_hem,h_hem,by="drug")
quant_concordance_stats(x = hhem$rp5, y =  hhem$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhem_h<-merge(r_all,h_hem,by="drug")
quant_concordance_stats(x = hhem_h$q5, y =  hhem_h$humanmin)

#####Hepatic---
#matching human~rat
rat_hep <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="hepatic")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_hep <- clinset2 %>%
  filter(hazard_sys_norm=="hepatic") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hhr<-merge(rat_hep,h_hep,by="drug")
quant_concordance_stats(x = hhr$rp5, y =  hhr$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hhr_h<-merge(r_all,h_hep,by="drug")
quant_concordance_stats(x = hhr_h$q5, y =  hhr_h$humanmin)

#####Dermal---
#matching human~rat
rat_derm <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="dermal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_derm <- clinset2 %>%
  filter(hazard_sys_norm=="dermal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hderm<-merge(rat_derm,h_derm,by="drug")
quant_concordance_stats(x = hderm$rp5, y =  hderm$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hderm_h<-merge(r_all,h_derm,by="drug")
quant_concordance_stats(x = hderm_h$q5, y =  hderm_h$humanmin)

#####Musculoskeletal---
#matching human~rat
rat_musc <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="musculoskeletal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_musc <- clinset2 %>%
  filter(hazard_sys_norm=="musculoskeletal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hmusc<-merge(rat_musc,h_musc,by="drug")
quant_concordance_stats(x = hmusc$rp5, y =  hmusc$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hmusc_h<-merge(r_all,h_musc,by="drug")
quant_concordance_stats(x = hmusc_h$q5, y =  hmusc_h$humanmin)

#####Renal---
#matching human~rat
rat_renal <- nonclin %>%
  filter(species=="rat")%>%
  filter(hazard_sys_norm=="renal")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(rp5 = quantile(smin, probs = q[1], na.rm = T,type=7))
h_renal <- clinset2 %>%
  filter(hazard_sys_norm=="renal") %>%
  group_by(drug) %>%
  summarise(humanmin= min(mg_kg_d, na.rm=T))%>%
  filter(!humanmin==Inf)
hrenal<-merge(rat_renal,h_renal,by="drug")
quant_concordance_stats(x = hrenal$rp5, y =  hrenal$humanmin)
#human only
r_all <- nonclin %>%
  filter(species=="rat")%>%
  group_by(studyID,drug)%>%
  summarise(smin = min(HED)) %>%
  group_by(drug) %>%
  summarise(q5 = quantile(smin, probs = q[1], na.rm = T,type=7))
hrenal_h<-merge(r_all,h_renal,by="drug")
quant_concordance_stats(x = hrenal_h$q5, y =  hrenal_h$humanmin)

###Qualitative concordance----
####Clinical dataset #3 All effects and drugs----
#####Hepatic---
hr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhep=any(r_hepatic=="Y"))%>%
  select(drug,rhep)%>%
  distinct()
hh<- clinset2 %>%
  mutate(h_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhep=any(h_hepatic=="Y"))%>%
  select(drug,humanhep)%>%
  distinct()
hepmerge<-merge(hh,hr,by="drug")
hepsum<- hepmerge %>%
  mutate(hepatic = case_when(rhep=="TRUE" & humanhep=="TRUE"~ "TP",
                             rhep=="TRUE" & humanhep=="FALSE" ~ "FP",
                             rhep=="FALSE" & humanhep=="TRUE" ~ "FN",
                             rhep=="FALSE" & humanhep=="FALSE" ~ "TN"))
hepsum %>% count(hepatic)
Hepatic<-hepsum %>%
  count(hepatic)%>%
  spread(hepatic,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Cardio---
cr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rcar=any(r_cardio=="Y"))%>%
  select(drug,rcar)%>%
  distinct()
hc<- clinset2 %>%
  mutate(h_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humancar=any(h_cardio=="Y"))%>%
  select(drug,humancar)%>%
  distinct()
carmerge<-merge(hc,cr,by="drug")
carsum<- carmerge %>%
  mutate(cardio = case_when(rcar=="TRUE" & humancar=="TRUE"~ "TP",
                            rcar=="TRUE" & humancar=="FALSE" ~ "FP",
                            rcar=="FALSE" & humancar=="TRUE" ~ "FN",
                            rcar=="FALSE" & humancar=="FALSE" ~ "TN"))
Cardiovascular<-carsum %>%
  count(cardio)%>%
  spread(cardio,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Renal---
rer <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rren=any(r_renal=="Y"))%>%
  select(drug,rren)%>%
  distinct()
hren<- clinset2 %>%
  mutate(h_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanren=any(h_renal=="Y"))%>%
  select(drug,humanren)%>%
  distinct()
renmerge<-merge(hren,rer,by="drug")
rensum<- renmerge %>%
  mutate(renal = case_when(rren=="TRUE" & humanren=="TRUE"~ "TP",
                           rren=="TRUE" & humanren=="FALSE" ~ "FP",
                           rren=="FALSE" & humanren=="TRUE" ~ "FN",
                           rren=="FALSE" & humanren=="FALSE" ~ "TN"))
rensum %>% count(renal)
Renal<- rensum %>%
  count(renal)%>%
  spread(renal,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Hemato---
hemr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhem=any(r_hem=="Y"))%>%
  select(drug,rhem)%>%
  distinct()
hhem<- clinset2 %>%
  mutate(h_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhem=any(h_hem=="Y"))%>%
  select(drug,humanhem)%>%
  distinct()
hemmerge<-merge(hemr,hhem,by="drug")
hemsum<- hemmerge %>%
  mutate(hemat = case_when(rhem=="TRUE" & humanhem=="TRUE"~ "TP",
                           rhem=="TRUE" & humanhem=="FALSE" ~ "FP",
                           rhem=="FALSE" & humanhem=="TRUE" ~ "FN",
                           rhem=="FALSE" & humanhem=="FALSE" ~ "TN"))
hemsum %>% count(hemat)
Hematologic<-hemsum %>%
  count(hemat)%>%
  spread(hemat,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Metabolic---
metr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmet=any(r_met=="Y"))%>%
  select(drug,rmet)%>%
  distinct()
hmet<- clinset2 %>%
  mutate(h_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmet=any(h_met=="Y"))%>%
  select(drug,humanmet)%>%
  distinct()
metmerge<-merge(metr,hmet,by="drug")
metsum<- metmerge %>%
  mutate(met = case_when(rmet=="TRUE" & humanmet=="TRUE"~ "TP",
                         rmet=="TRUE" & humanmet=="FALSE" ~ "FP",
                         rmet=="FALSE" & humanmet=="TRUE" ~ "FN",
                         rmet=="FALSE" & humanmet=="FALSE" ~ "TN"))
metsum %>% count(met)
Metabolic<- metsum %>%
  count(met)%>%
  spread(met,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Psychiatric---
psyr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rpsy=any(r_psy=="Y"))%>%
  select(drug,rpsy)%>%
  distinct()
hpsy<- clinset2 %>%
  mutate(h_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanpsy=any(h_psy=="Y"))%>%
  select(drug,humanpsy)%>%
  distinct()
psymerge<-merge(psyr,hpsy,by="drug")
psysum<- psymerge %>%
  mutate(neuro = case_when(rpsy=="TRUE" & humanpsy=="TRUE"~ "TP",
                           rpsy=="TRUE" & humanpsy=="FALSE" ~ "FP",
                           rpsy=="FALSE" & humanpsy=="TRUE" ~ "FN",
                           rpsy=="FALSE" & humanpsy=="FALSE" ~ "TN"))
psysum %>% count(neuro)
Nervous<- psysum %>%
  count(neuro)%>%
  spread(neuro,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Gastrointestinal---
gasr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rgas=any(r_gas=="Y"))%>%
  select(drug,rgas)%>%
  distinct()
hgas<- clinset2 %>%
  mutate(h_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humangas=any(h_gas=="Y"))%>%
  select(drug,humangas)%>%
  distinct()
gasmerge<-merge(gasr,hgas,by="drug")
gassum<- gasmerge %>%
  mutate(gas = case_when(rgas=="TRUE" & humangas=="TRUE"~ "TP",
                         rgas=="TRUE" & humangas=="FALSE" ~ "FP",
                         rgas=="FALSE" & humangas=="TRUE" ~ "FN",
                         rgas=="FALSE" & humangas=="FALSE" ~ "TN"))
gassum %>% count(gas)
Gastrointestinal<- gassum %>%
  count(gas)%>%
  spread(gas,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Musculo---
musr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmus=any(r_mus=="Y"))%>%
  select(drug,rmus)%>%
  distinct()
hmus<- clinset2 %>%
  mutate(h_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmus=any(h_mus=="Y"))%>%
  select(drug,humanmus)%>%
  distinct()
musmerge<-merge(musr,hmus,by="drug")
mussum<- musmerge %>%
  mutate(mus = case_when(rmus=="TRUE" & humanmus=="TRUE"~ "TP",
                         rmus=="TRUE" & humanmus=="FALSE" ~ "FP",
                         rmus=="FALSE" & humanmus=="TRUE" ~ "FN",
                         rmus=="FALSE" & humanmus=="FALSE" ~ "TN"))
mussum %>% count(mus)
Musculoskeletal<- mussum %>%
  count(mus)%>%
  spread(mus,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Derm---
dermr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rderm=any(r_derm=="Y"))%>%
  select(drug,rderm)%>%
  distinct()
hderm<- clinset2 %>%
  mutate(h_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanderm=any(h_derm=="Y"))%>%
  select(drug,humanderm)%>%
  distinct()
dermmerge<-merge(dermr,hderm,by="drug")
dermsum<- dermmerge %>%
  mutate(derm = case_when(rderm=="TRUE" & humanderm=="TRUE"~ "TP",
                          rderm=="TRUE" & humanderm=="FALSE" ~ "FP",
                          rderm=="FALSE" & humanderm=="TRUE" ~ "FN",
                          rderm=="FALSE" & humanderm=="FALSE" ~ "TN"))
dermsum %>% count(derm)
Dermal<- dermsum %>%
  count(derm)%>%
  spread(derm,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

######Qualitative summary table----
qualtable_clinset2<- lst(Gastrointestinal,Metabolic,Cardiovascular,Nervous,Hematologic,Hepatic,Musculoskeletal,Dermal,Renal) %>%
  bind_rows(.id = "hazard_sys_norm")

####Clinical dataset #2 All effects and drugs vs mouse----
#####Hepatic---
hr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhep=any(r_hepatic=="Y"))%>%
  select(drug,rhep)%>%
  distinct()
hh<- clinset1 %>%
  mutate(h_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhep=any(h_hepatic=="Y"))%>%
  select(drug,humanhep)%>%
  distinct()
hepmerge<-merge(hh,hr,by="drug")
hepsum<- hepmerge %>%
  mutate(hepatic = case_when(rhep=="TRUE" & humanhep=="TRUE"~ "TP",
                             rhep=="TRUE" & humanhep=="FALSE" ~ "FP",
                             rhep=="FALSE" & humanhep=="TRUE" ~ "FN",
                             rhep=="FALSE" & humanhep=="FALSE" ~ "TN"))
hepsum %>% count(hepatic)
Hepatic<-hepsum %>%
  count(hepatic)%>%
  spread(hepatic,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Cardio---
cr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rcar=any(r_cardio=="Y"))%>%
  select(drug,rcar)%>%
  distinct()
hc<- clinset1 %>%
  mutate(h_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humancar=any(h_cardio=="Y"))%>%
  select(drug,humancar)%>%
  distinct()
carmerge<-merge(hc,cr,by="drug")
carsum<- carmerge %>%
  mutate(cardio = case_when(rcar=="TRUE" & humancar=="TRUE"~ "TP",
                            rcar=="TRUE" & humancar=="FALSE" ~ "FP",
                            rcar=="FALSE" & humancar=="TRUE" ~ "FN",
                            rcar=="FALSE" & humancar=="FALSE" ~ "TN"))
Cardiovascular<-carsum %>%
  count(cardio)%>%
  spread(cardio,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Renal---
rer <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rren=any(r_renal=="Y"))%>%
  select(drug,rren)%>%
  distinct()
hren<- clinset1 %>%
  mutate(h_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanren=any(h_renal=="Y"))%>%
  select(drug,humanren)%>%
  distinct()
renmerge<-merge(hren,rer,by="drug")
rensum<- renmerge %>%
  mutate(renal = case_when(rren=="TRUE" & humanren=="TRUE"~ "TP",
                           rren=="TRUE" & humanren=="FALSE" ~ "FP",
                           rren=="FALSE" & humanren=="TRUE" ~ "FN",
                           rren=="FALSE" & humanren=="FALSE" ~ "TN"))
rensum %>% count(renal)
Renal<- rensum %>%
  count(renal)%>%
  spread(renal,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Hemato---
hemr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhem=any(r_hem=="Y"))%>%
  select(drug,rhem)%>%
  distinct()
hhem<- clinset1 %>%
  mutate(h_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhem=any(h_hem=="Y"))%>%
  select(drug,humanhem)%>%
  distinct()
hemmerge<-merge(hemr,hhem,by="drug")
hemsum<- hemmerge %>%
  mutate(hemat = case_when(rhem=="TRUE" & humanhem=="TRUE"~ "TP",
                           rhem=="TRUE" & humanhem=="FALSE" ~ "FP",
                           rhem=="FALSE" & humanhem=="TRUE" ~ "FN",
                           rhem=="FALSE" & humanhem=="FALSE" ~ "TN"))
hemsum %>% count(hemat)
Hematologic<-hemsum %>%
  count(hemat)%>%
  spread(hemat,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Metabolic---
metr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmet=any(r_met=="Y"))%>%
  select(drug,rmet)%>%
  distinct()
hmet<- clinset1 %>%
  mutate(h_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmet=any(h_met=="Y"))%>%
  select(drug,humanmet)%>%
  distinct()
metmerge<-merge(metr,hmet,by="drug")
metsum<- metmerge %>%
  mutate(met = case_when(rmet=="TRUE" & humanmet=="TRUE"~ "TP",
                         rmet=="TRUE" & humanmet=="FALSE" ~ "FP",
                         rmet=="FALSE" & humanmet=="TRUE" ~ "FN",
                         rmet=="FALSE" & humanmet=="FALSE" ~ "TN"))
metsum %>% count(met)
Metabolic<- metsum %>%
  count(met)%>%
  spread(met,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Psychiatric---
psyr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rpsy=any(r_psy=="Y"))%>%
  select(drug,rpsy)%>%
  distinct()
hpsy<- clinset1 %>%
  mutate(h_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanpsy=any(h_psy=="Y"))%>%
  select(drug,humanpsy)%>%
  distinct()
psymerge<-merge(psyr,hpsy,by="drug")
psysum<- psymerge %>%
  mutate(neuro = case_when(rpsy=="TRUE" & humanpsy=="TRUE"~ "TP",
                           rpsy=="TRUE" & humanpsy=="FALSE" ~ "FP",
                           rpsy=="FALSE" & humanpsy=="TRUE" ~ "FN",
                           rpsy=="FALSE" & humanpsy=="FALSE" ~ "TN"))
psysum %>% count(neuro)
Nervous<- psysum %>%
  count(neuro)%>%
  spread(neuro,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Gastrointestinal---
gasr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rgas=any(r_gas=="Y"))%>%
  select(drug,rgas)%>%
  distinct()
hgas<- clinset1 %>%
  mutate(h_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humangas=any(h_gas=="Y"))%>%
  select(drug,humangas)%>%
  distinct()
gasmerge<-merge(gasr,hgas,by="drug")
gassum<- gasmerge %>%
  mutate(gas = case_when(rgas=="TRUE" & humangas=="TRUE"~ "TP",
                         rgas=="TRUE" & humangas=="FALSE" ~ "FP",
                         rgas=="FALSE" & humangas=="TRUE" ~ "FN",
                         rgas=="FALSE" & humangas=="FALSE" ~ "TN"))
gassum %>% count(gas)
Gastrointestinal<- gassum %>%
  count(gas)%>%
  spread(gas,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Musculo---
musr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmus=any(r_mus=="Y"))%>%
  select(drug,rmus)%>%
  distinct()
hmus<- clinset1 %>%
  mutate(h_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmus=any(h_mus=="Y"))%>%
  select(drug,humanmus)%>%
  distinct()
musmerge<-merge(musr,hmus,by="drug")
mussum<- musmerge %>%
  mutate(mus = case_when(rmus=="TRUE" & humanmus=="TRUE"~ "TP",
                         rmus=="TRUE" & humanmus=="FALSE" ~ "FP",
                         rmus=="FALSE" & humanmus=="TRUE" ~ "FN",
                         rmus=="FALSE" & humanmus=="FALSE" ~ "TN"))
mussum %>% count(mus)
Musculoskeletal<- mussum %>%
  count(mus)%>%
  spread(mus,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Derm---
dermr <- nonclin %>%
  filter(species=="mouse")%>%
  mutate(r_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rderm=any(r_derm=="Y"))%>%
  select(drug,rderm)%>%
  distinct()
hderm<- clinset1 %>%
  mutate(h_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanderm=any(h_derm=="Y"))%>%
  select(drug,humanderm)%>%
  distinct()
dermmerge<-merge(dermr,hderm,by="drug")
dermsum<- dermmerge %>%
  mutate(derm = case_when(rderm=="TRUE" & humanderm=="TRUE"~ "TP",
                          rderm=="TRUE" & humanderm=="FALSE" ~ "FP",
                          rderm=="FALSE" & humanderm=="TRUE" ~ "FN",
                          rderm=="FALSE" & humanderm=="FALSE" ~ "TN"))
dermsum %>% count(derm)
Dermal<- dermsum %>%
  count(derm)%>%
  spread(derm,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

######Qualitative summary table----
qualtable_clinset1_mouse<- lst(Gastrointestinal,Metabolic,Cardiovascular,Nervous,Hematologic,Hepatic,Musculoskeletal,Dermal,Renal) %>%
  bind_rows(.id = "hazard_sys_norm")

####Rat vs mouse----
#####Hepatic---
hr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhep=any(r_hepatic=="Y"))%>%
  select(drug,rhep)%>%
  distinct()
hh<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_hepatic= ifelse(grepl("hepat", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hepatic)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhep=any(h_hepatic=="Y"))%>%
  select(drug,humanhep)%>%
  distinct()
hepmerge<-merge(hh,hr,by="drug")
hepsum<- hepmerge %>%
  mutate(hepatic = case_when(rhep=="TRUE" & humanhep=="TRUE"~ "TP",
                             rhep=="TRUE" & humanhep=="FALSE" ~ "FP",
                             rhep=="FALSE" & humanhep=="TRUE" ~ "FN",
                             rhep=="FALSE" & humanhep=="FALSE" ~ "TN"))
hepsum %>% count(hepatic)
Hepatic<-hepsum %>%
  count(hepatic)%>%
  spread(hepatic,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Cardio---
cr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rcar=any(r_cardio=="Y"))%>%
  select(drug,rcar)%>%
  distinct()
hc<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_cardio= ifelse(grepl("cardio", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_cardio)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humancar=any(h_cardio=="Y"))%>%
  select(drug,humancar)%>%
  distinct()
carmerge<-merge(hc,cr,by="drug")
carsum<- carmerge %>%
  mutate(cardio = case_when(rcar=="TRUE" & humancar=="TRUE"~ "TP",
                            rcar=="TRUE" & humancar=="FALSE" ~ "FP",
                            rcar=="FALSE" & humancar=="TRUE" ~ "FN",
                            rcar=="FALSE" & humancar=="FALSE" ~ "TN"))
Cardiovascular<-carsum %>%
  count(cardio)%>%
  spread(cardio,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Renal---
rer <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rren=any(r_renal=="Y"))%>%
  select(drug,rren)%>%
  distinct()
hren<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_renal= ifelse(grepl("renal", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_renal)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanren=any(h_renal=="Y"))%>%
  select(drug,humanren)%>%
  distinct()
renmerge<-merge(hren,rer,by="drug")
rensum<- renmerge %>%
  mutate(renal = case_when(rren=="TRUE" & humanren=="TRUE"~ "TP",
                           rren=="TRUE" & humanren=="FALSE" ~ "FP",
                           rren=="FALSE" & humanren=="TRUE" ~ "FN",
                           rren=="FALSE" & humanren=="FALSE" ~ "TN"))
rensum %>% count(renal)
Renal<- rensum %>%
  count(renal)%>%
  spread(renal,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Hemato---
hemr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rhem=any(r_hem=="Y"))%>%
  select(drug,rhem)%>%
  distinct()
hhem<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_hem= ifelse(grepl("hematologic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_hem)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanhem=any(h_hem=="Y"))%>%
  select(drug,humanhem)%>%
  distinct()
hemmerge<-merge(hemr,hhem,by="drug")
hemsum<- hemmerge %>%
  mutate(hemat = case_when(rhem=="TRUE" & humanhem=="TRUE"~ "TP",
                           rhem=="TRUE" & humanhem=="FALSE" ~ "FP",
                           rhem=="FALSE" & humanhem=="TRUE" ~ "FN",
                           rhem=="FALSE" & humanhem=="FALSE" ~ "TN"))
hemsum %>% count(hemat)
Hematologic<-hemsum %>%
  count(hemat)%>%
  spread(hemat,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Metabolic---
metr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmet=any(r_met=="Y"))%>%
  select(drug,rmet)%>%
  distinct()
hmet<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_met= ifelse(grepl("metabolic", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_met)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmet=any(h_met=="Y"))%>%
  select(drug,humanmet)%>%
  distinct()
metmerge<-merge(metr,hmet,by="drug")
metsum<- metmerge %>%
  mutate(met = case_when(rmet=="TRUE" & humanmet=="TRUE"~ "TP",
                         rmet=="TRUE" & humanmet=="FALSE" ~ "FP",
                         rmet=="FALSE" & humanmet=="TRUE" ~ "FN",
                         rmet=="FALSE" & humanmet=="FALSE" ~ "TN"))
metsum %>% count(met)
Metabolic<- metsum %>%
  count(met)%>%
  spread(met,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Psychiatric---
psyr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rpsy=any(r_psy=="Y"))%>%
  select(drug,rpsy)%>%
  distinct()
hpsy<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_psy= ifelse(grepl("nervous", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_psy)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanpsy=any(h_psy=="Y"))%>%
  select(drug,humanpsy)%>%
  distinct()
psymerge<-merge(psyr,hpsy,by="drug")
psysum<- psymerge %>%
  mutate(neuro = case_when(rpsy=="TRUE" & humanpsy=="TRUE"~ "TP",
                           rpsy=="TRUE" & humanpsy=="FALSE" ~ "FP",
                           rpsy=="FALSE" & humanpsy=="TRUE" ~ "FN",
                           rpsy=="FALSE" & humanpsy=="FALSE" ~ "TN"))
psysum %>% count(neuro)
Nervous<- psysum %>%
  count(neuro)%>%
  spread(neuro,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Gastrointestinal---
gasr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rgas=any(r_gas=="Y"))%>%
  select(drug,rgas)%>%
  distinct()
hgas<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_gas= ifelse(grepl("gastro", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_gas)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humangas=any(h_gas=="Y"))%>%
  select(drug,humangas)%>%
  distinct()
gasmerge<-merge(gasr,hgas,by="drug")
gassum<- gasmerge %>%
  mutate(gas = case_when(rgas=="TRUE" & humangas=="TRUE"~ "TP",
                         rgas=="TRUE" & humangas=="FALSE" ~ "FP",
                         rgas=="FALSE" & humangas=="TRUE" ~ "FN",
                         rgas=="FALSE" & humangas=="FALSE" ~ "TN"))
gassum %>% count(gas)
Gastrointestinal<- gassum %>%
  count(gas)%>%
  spread(gas,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Musculo---
musr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rmus=any(r_mus=="Y"))%>%
  select(drug,rmus)%>%
  distinct()
hmus<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_mus= ifelse(grepl("muscul", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_mus)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanmus=any(h_mus=="Y"))%>%
  select(drug,humanmus)%>%
  distinct()
musmerge<-merge(musr,hmus,by="drug")
mussum<- musmerge %>%
  mutate(mus = case_when(rmus=="TRUE" & humanmus=="TRUE"~ "TP",
                         rmus=="TRUE" & humanmus=="FALSE" ~ "FP",
                         rmus=="FALSE" & humanmus=="TRUE" ~ "FN",
                         rmus=="FALSE" & humanmus=="FALSE" ~ "TN"))
mussum %>% count(mus)
Musculoskeletal<- mussum %>%
  count(mus)%>%
  spread(mus,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

####Derm---
dermr <- nonclin %>%
  filter(species=="rat")%>%
  mutate(r_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, r_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(rderm=any(r_derm=="Y"))%>%
  select(drug,rderm)%>%
  distinct()
hderm<- nonclin %>%
  filter(species=="mouse")%>%
  mutate(h_derm= ifelse(grepl("derm", hazard_sys_norm), "Y","N")) %>%
  select(drug, h_derm)%>%
  distinct() %>%
  group_by(drug)%>%
  mutate(humanderm=any(h_derm=="Y"))%>%
  select(drug,humanderm)%>%
  distinct()
dermmerge<-merge(dermr,hderm,by="drug")
dermsum<- dermmerge %>%
  mutate(derm = case_when(rderm=="TRUE" & humanderm=="TRUE"~ "TP",
                          rderm=="TRUE" & humanderm=="FALSE" ~ "FP",
                          rderm=="FALSE" & humanderm=="TRUE" ~ "FN",
                          rderm=="FALSE" & humanderm=="FALSE" ~ "TN"))
dermsum %>% count(derm)
Dermal<- dermsum %>%
  count(derm)%>%
  spread(derm,n)%>%
  mutate(prev= (TP+FN)/(FN+FP+TN+TP),
         sens= (TP/(TP+FN)),
         spec= (TN/(TN+FP)),
         PPV= (TP/(TP+FP))*100,
         NPV= (TN/(TN+FN))*100,
         LRneg= spec/(1-sens),
         LRpos= sens/(1-spec),
         BA= (sens+spec)/2)

######Qualitative summary table----
qualtable_rat_mouse<- lst(Gastrointestinal,Metabolic,Cardiovascular,Nervous,Hematologic,Hepatic,Musculoskeletal,Dermal,Renal) %>%
  bind_rows(.id = "hazard_sys_norm")

