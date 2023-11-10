# .............................
# .............................
# Packages ####
library("readxl")
library("ClimMobTools")
library("gosset")
library("janitor")

list.files("data")

dat1 = read.csv("data/raw/tricot-consumer-cassava-nigeria.csv")

dat2 = read.csv("data/raw/tricot-consumer-cassava-cameroon.csv")

dat3 = read.csv("data/raw/tricot-consumer-sweetpotato-uganda.csv")

# deal with dat1 
names(dat1)

dat1$technology = "cassava-eba"
dat1$country = "Nigeria"

names(dat1)[grepl("_overall_best_", names(dat1))] = "best_overall"
names(dat1)[grepl("_overall_worst_", names(dat1))] = "worst_overall"
names(dat1)[grepl("describe_the_taste_of_eba_you_prefer_like", names(dat1))] = "reason_like"
names(dat1)[grepl("describe_the_taste_of_eba_you_detest_dislike", names(dat1))] = "reason_dislike"
names(dat1)[grepl("gps_coordinate_n", names(dat1))] = "latitude"
names(dat1)[grepl("gps_coordinate_e", names(dat1))] = "longitude"
names(dat1)[grepl("consumers_gender", names(dat1))] = "gender"

names(dat1) = gsub("strechability", "stretchability", names(dat1))

names(dat1)

# remove farmer name and phone
rmv = !grepl("farmer_name|consumer_name|phone|gps_point_|which_attribute_was|are_there_other_attribute|
              unique_row_|do_you_consent_to_provide|submission_id|x_form_id|enumerator_id|package_code", names(dat1))

dat1 = dat1[rmv]

names(dat1)

sel = c("best_", "worst_", "id", "option_", "package_", 
         "reason_", "long", "lat", "country", "technology", 
         "ethinic", "gender", "age$", "marital_", "occupation",
         "how_often", "what_accasion_")

sel = paste(sel, collapse = "|")         

sel = grepl(sel, names(dat1))

dat1 = dat1[sel]

names(dat1)

# fix names with tricot traits
sel = grepl("best_", names(dat1))
newname = gsub("best_", "", names(dat1)[sel])
newname = paste0(newname, "_pos")
names(dat1)[sel] = newname

sel = grepl("worst_", names(dat1))
newname = gsub("worst_", "", names(dat1)[sel])
newname = paste0(newname, "_neg")
names(dat1)[sel] = newname

names(dat1) = gsub("_brightness_", "_", names(dat1))
names(dat1) = gsub("how_often_do_you_consume_eba", "consumption", names(dat1))

# now dat2
dat2$technology = "cassava-eba"
dat2$country = "Cameroon"

names(dat2) = make_clean_names(names(dat2))

names(dat2) = gsub("registration_|qst_", "", names(dat2))

names(dat2)

names(dat2)[grepl("gender1", names(dat2))] = "gender"
names(dat2)[grepl("ethnicgroup", names(dat2))] = "ethnic_group"

names(dat2) = gsub("package_item_", "option_", names(dat2))
names(dat2) = gsub("overallappreciation_", "overall_", names(dat2))
names(dat2) = gsub("notpreferredtaste", "reason_dislike", names(dat2))
names(dat2) = gsub("preferredtaste", "reason_like", names(dat2))
names(dat2) = gsub("strechability", "stretchability", names(dat2))


sel = c("_pos", "_neg", "id", "option_", "package_", 
         "reason_", "long", "lat", "country", "technology", 
         "ethinic", "gender", "age$", "marital_", "occupation",
         "consumption", "what_accasion_")

sel = paste(sel, collapse = "|")         

sel = grepl(sel, names(dat2))

names(dat2)[sel]

dat2 = dat2[sel]

names(dat2)

rmv = !grepl("participant_name|phone|gps_point_|unique_row_identifier|
              xform_id_|enumeratorid|surveyid|xform_id_string", names(dat2))

dat2 = dat2[rmv]

dat = rowbind(dat1, dat2)

names(dat)

# now organise the technology names
items = sort(unique(unlist(dat[,paste0("option_", letters[1:3])])))

items

# the local variates have switched names when they were tested in the field
dat[,paste0("option_", letters[1:3])] = 
  lapply(dat[,paste0("option_", letters[1:3])], function(x){
  x = gsub("Local 1_SAPE", "Sape", x)
  x = gsub("Local 2 _Madame", "Madame", x)
  x = gsub(" ", "", x)
  x
})

items = sort(unique(unlist(dat[,paste0("option_", letters[1:3])])))

items

table(unlist(dat[,paste0("option_", letters[1:3])]))

# now do the same for the genotype features 
excel_sheets("data/raw/RTBFoods-WP5_Lab-gari.xlsx")

# first the biophysical
feat1  = read_excel("data/raw/RTBFoods-WP5_Lab-gari.xlsx", 
                  sheet = 1, 
                  na = c(".", " "))
names(feat1) = make_clean_names(names(feat1))

# the other region
feat2  = read_excel("data/raw/RTBFoods-WP5_Lab-gari.xlsx", 
                     sheet = 3, 
                     na = c(".", " "))
names(feat2) = make_clean_names(names(feat2))

feat = rowbind(feat1, feat2)

names(feat)

rmv = c("processor", "sample", "batch", "trough_1", "peak_1")
rmv = !names(feat) %in% rmv

feat = feat[rmv]

names(feat)

# now from Cameroon 
list.files("data/raw", full.names = TRUE)

excel_sheets("data/raw/ITPA STPA Gari Eba_Cameroon Samples_2022_1.xlsx")

feat3 = read_excel("data/raw/ITPA STPA Gari Eba_Cameroon Samples_2022_1.xlsx", 
                    sheet = 2, 
                    na = c(".", " "))
names(feat3) = make_clean_names(names(feat3))

feat4 = read_excel("data/raw/ITPA STPA Gari Eba_Cameroon Samples_2022_1.xlsx", 
                    sheet = 4, 
                    na = c(".", " "))
names(feat4) = make_clean_names(names(feat4))

feat5 = read_excel("data/raw/ITPA STPA Gari Eba_Cameroon Samples_2022_1.xlsx", 
                    sheet = 5, 
                    na = c(".", " "))
names(feat5) = make_clean_names(names(feat5))

names(feat)

names(feat3)

sel = c("variety_name", "color_l_gari", "color_a_gari", "color_b_gari")

feat3 = feat3[sel]

names(feat3) = gsub("_gari", "", names(feat3))

feat = rowbind(feat, feat3)

names(feat)

names(feat4)

rmv = c("sn",  "sample_code" , "date_of_analysis")
rmv = !names(feat4) %in% rmv

feat4 = feat4[rmv]

names(feat)

rmv = c("s_n",  "sample_codes" , "date_of_analysis")
rmv = !names(feat5) %in% rmv

feat5 = feat5[rmv]

feat = rowbind(feat, feat4)

feat = rowbind(feat, feat5)

# now for the othe traits
excel_sheets("data/raw/RTBFoods-WP5_Lab-gari.xlsx")

ebatex1  = read_excel("data/raw/RTBFoods-WP5_Lab-gari.xlsx", 
                     sheet = 2, 
                     na = c(".", " "))
names(ebatex1) = make_clean_names(names(ebatex1))


ebatex2  = read_excel("data/raw/RTBFoods-WP5_Lab-gari.xlsx", 
                       sheet = 4, 
                       na = c(".", " "))
names(ebatex2) = make_clean_names(names(ebatex2))


head(ebatex1)

ebatex = rowbind(ebatex1, ebatex2)

rmv = c("processor", "sample", "batch")
rmv = !names(ebatex) %in% rmv

ebatex = ebatex[rmv]

names(ebatex)

sel = names(feat4)[names(feat4) %in% names(ebatex)]

#feat4 = feat4[sel]

#ebatex = rowbind(ebatex, feat)

names(feat)[names(feat) == "variety_name"] = "genotype"

names(ebatex)[names(ebatex) == "variety_name"] = "genotype"

feat = rowbind(feat, ebatex)

# check whether the varieties names have the same names as the lab data
pack = paste0("option_", letters[1:3])

sort(unique(unlist(dat[pack])))
sort(unique(feat$genotype))

feat$genotype = gsub("TMS92/0326", "IITATMSIBA920326", feat$genotype)
feat$genotype = gsub("TMS92/0326", "IITATMSIBA920326", feat$genotype)
feat$genotype = gsub("LOCAL 1", "Sape", feat$genotype)
feat$genotype = gsub("LOCAL 2", "Madame", feat$genotype)
feat$genotype = gsub("SAPE", "Sape", feat$genotype)

all(unique(unlist(dat[pack])) %in% as.character(feat$genotype))
all(as.character(feat$genotype) %in% unique(unlist(dat[pack])))

feat = feat[as.character(feat$genotype) %in% unique(unlist(dat[pack])), ]

all(unique(unlist(dat[pack])) %in% as.character(feat$genotype))
all(as.character(feat$genotype) %in% unique(unlist(dat[pack])))

feat = feat[!is.na(feat$genotype), ]

# write files
write.csv(feat, "data/cassava-biophysical-features.csv", row.names = FALSE)

write.csv(dat, "data/cassava-tricot-data.csv", row.names = FALSE)

# .........................................
# .........................................
# now work with the sweetpotato data

names(dat3)

dat3$technology = "sweetpotato"
dat3$country = "Uganda"

names(dat3) = make_clean_names(names(dat3))
names(dat3) = gsub("registration_|consumer_testing_4map_", "", names(dat3))
names(dat3) = gsub("sexparticipant", "gender", names(dat3))
names(dat3) = gsub("rank", "", names(dat3))
names(dat3) = gsub("descriptndislike", "reason_dislike", names(dat3))
names(dat3) = gsub("descriptnlike", "reason_like", names(dat3))
names(dat3) = gsub("mostpreferedsample_", "overall_", names(dat3))
names(dat3) = gsub("mostpreferedsample_", "overall_", names(dat3))
names(dat3) = gsub("package_item_", "option_", names(dat3))
names(dat3) = gsub("districtparticipant", "district", names(dat3))


names(dat3)


sel = c("_pos", "_neg", "id", "option_", "district", 
         "reason_", "country", "technology", 
         "ethinic", "gender")

sel = paste(sel, collapse = "|")         

sel = grepl(sel, names(dat3))

dat3 = dat3[sel]

names(dat3)

rmv = !grepl("package_|surveyid|xform_id_string|[[:digit:]]+", names(dat3))
# 
dat3 = dat3[rmv]

names(dat3)


sort(unique(unlist(dat3[paste0("option_", letters[1:3])])))

# # add features
# files = list.files("data/raw", pattern = "POTATO", full.names = TRUE)
# 
# files
# 
# sp_feat = data.frame()
# 
# for (i in seq_along(files)) {
#   x = read_excel(files[i])
#   
#   names(x) = make_clean_names(names(x))
#   
#   x = as.data.frame(x)
#   
#   x = x[-c(1:3), ]
#   
#   rmv = !grepl("End|Coef|Start|S.|Average", as.vector(x[,1]))
#   
#   x = x[rmv, ]
#   
#   sp_feat = rbind(sp_feat, x)
# }
# 
# items = sort(unique(unlist(dat3[paste0("option_", letters[1:3])])))
# 
# items
# 
# for (i in seq_along(items)) {
#   rpl = grepl(items[i], sp_feat$test_id)
#   sp_feat$test_id[rpl] = items[i]
# }
# 
# items %in% sp_feat$test_id
# 
# 
# sp_feat = sp_feat[sp_feat$test_id %in% items, -c(2:3) ]

list.files("data/raw")
sp_feat = read_excel("data/raw/MUK_CIRAD Boiled SP KABALE ON FARM DSA DATA FINAL.xlsx", sheet = "mean per product")

names(sp_feat) = make_clean_names(names(sp_feat))

items = sort(unique(unlist(dat3[paste0("option_", letters[1:3])])))

items

names(sp_feat)[names(sp_feat) == "sample_code"] = "genotype"

for (i in seq_along(items)) {
  rpl = grepl(items[i], sp_feat$genotype)
  sp_feat$genotype[rpl] = items[i]
}

sp_feat$genotype[sp_feat$genotype=="MUGURUSI  KABALE ONFARM"] = "LOCAL CHECK"

sp_feat = sp_feat[, -1]

sp_feat = sp_feat[, -2]

write.csv(dat3, "data/sweetpotato-tricot-data.csv", row.names = FALSE)

write.csv(sp_feat, "data/sweetpotato-biophysical-features.csv", row.names = FALSE)

