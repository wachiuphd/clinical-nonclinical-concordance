library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(janitor)
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

#####Table S2 of nonclinical safety data----
write.csv(df,"Supplemental Table S2 nonclinical safety dataset1.csv", row.names = F)

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
write.csv(longh,"Supplemental Table S3 clinical safety dataset1.csv", row.names = F)

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
write.csv(clinset2_comb,"Supplemental Table S4 clinical safety dataset2.csv", row.names=F)
