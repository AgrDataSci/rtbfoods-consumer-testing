# .............................
# .............................
# Packages ####
library("PlackettLuce")
library("gosset")
library("ClimMobTools")
library("partykit")
library("gtools")
library("patchwork")
library("ggplot2")
library("tidytext")
library("tm")
library("wordcloud")

source("https://raw.githubusercontent.com/AgrDataSci/ClimMob-analysis/master/modules/01_functions.R")

# function to extract pladmm coeffs as data.frame 
pladmm_coeffs = function(object, ref = 1L, ...) {
  
  # Extract ids from terminal nodes
  node_id = partykit::nodeids(object, terminal = TRUE)
  
  # get models from each node
  nodes = list()
  for (i in seq_along(node_id)) {
    obj_i = object[[ node_id[i] ]]$node$info$object
    coefs = coef(obj_i)
    coefficients = matrix(NA, nrow = length(coefs), ncol = 4L, 
                           dimnames = list(names(coefs), c("Estimate", "Std. Error", 
                                                           "z value", "Pr(>|z|)")))
    coefficients[, 1L] = coefs
    se = sqrt(diag(vcov(obj_i)))
    coefficients[names(se), 2L] = se
    coefficients[, 3L] = coefficients[, 1L]/coefficients[, 2L]
    coefficients[, 4L] = 2L * pnorm(-abs(coefficients[, 3L]))
    
    coefficients = as.data.frame(coefficients)
    
    coefficients[, 5] = gtools::stars.pval(coefficients[, 4])
    
    coefficients[, 4] = formatC(coefficients[, 4], format = "e", digits = 2)
    
    coefficients[, 6] = node_id[i]
    
    coefficients[, 7] = rownames(coefficients)
    
    rownames(coefficients) = 1:nrow(coefficients)
    
    coefficients = coefficients[,c(6, 7, 1:5)]
    
    names(coefficients)[1] = "Node"
    
    names(coefficients)[c(2, 7)] = ""
    
    nodes[[i]] = coefficients
    
  }
  
  result = do.call("rbind", nodes)
  
  rownames(result) = 1:nrow(result)
  
  return(result)
  
}

PL_coeffs = function(object, ref = 1L, ...) {
  coefs <- coef(object, ref = ref)
  coefficients <- matrix(NA, nrow = length(coefs), ncol = 4L, 
                         dimnames = list(names(coefs), c("Estimate", "Std. Error", 
                                                         "z value", "Pr(>|z|)")))
  coefficients[, 1L] <- coefs
  se <- sqrt(diag(vcov(object, ref = ref, ...)))
  coefficients[names(se), 2L] <- se
  ref <- attr(coefs, "ref")
  if (length(ref) == 1L) 
    coefficients[attr(coefs, "ref"), 2L] <- NA
  coefficients[, 3L] <- coefficients[, 1L]/coefficients[, 2L]
  coefficients[, 4L] <- 2L * pnorm(-abs(coefficients[, 3L]))
  coefficients
}

# Read and organize data ####

# make PlackettLuce rankings
list.files("data", full.names = TRUE)

dat = read.csv("data/cassava-tricot-data.csv")

names(dat)

features = read.csv("data/blups-cassava.csv")

names(dat)

# replace names with readable names
items = unique(dat$option_a)

items

pack = c("option_a", "option_b", "option_c")

for (i in seq_along(pack)) {
  dat[pack[i]][dat[pack[i]] == "TMS13F1307P0016"] =	"TMS1"
  dat[pack[i]][dat[pack[i]] == "TMS13F1343P0044"] =	"TMS2"
  dat[pack[i]][dat[pack[i]] == "TMS14F1278P0003"] =	"TMS3"
  dat[pack[i]][dat[pack[i]] == "TMS13F1160P0004"] =	"Game Changer"
  dat[pack[i]][dat[pack[i]] == "TMS13F1343P0022"] =	"Obasanjo-2"
  dat[pack[i]][dat[pack[i]] == "TMS30572"] =	"TMS6"
  dat[pack[i]][dat[pack[i]] == "TMEB1_MS6"] =	"TMEB1"
  dat[pack[i]][dat[pack[i]] == "IITATMSIBA920326"] =	"TMSIBA"
}

items = unique(unlist(dat[pack]))

features$genotype[features$genotype == "TMS13F1307P0016"] =	"TMS1"
features$genotype[features$genotype == "TMS13F1343P0044"] =	"TMS2"
features$genotype[features$genotype == "TMS14F1278P0003"] =	"TMS3"
features$genotype[features$genotype == "TMS13F1160P0004"] =	"Game Changer"
features$genotype[features$genotype == "TMS13F1343P0022"] =	"Obasanjo-2"
features$genotype[features$genotype == "TMS30572"] =	"TMS6"
features$genotype[features$genotype == "TMEB1_MS6"] =	"TMEB1"
features$genotype[features$genotype == "IITATMSIBA920326"] =	"TMSIBA"

# check covariates
table(dat$gender)

dat$gender[dat$gender == "Female"] = "Woman"
dat$gender[dat$gender == "Male"] = "Man"

boxplot(dat$age ~ dat$gender)

plot(dat[,c("longitude", "latitude")], 
     pch = as.integer(factor(dat$gender)))

dat$consumption[dat$consumption == "Once A Month"] = "Once in a month"

dat$consumption = ClimMobTools:::.title_case(dat$consumption)

# check the rankings
table(unlist(dat[pack]), rep(dat$country, 3))

names(dat)[grepl("_pos", names(dat))]

# rmv = !grepl("odour|firmness", names(dat))
# dat = dat[rmv]

trait = c("colour", "smoothness", "mouldability", "stretchability", 
          "taste",  "overall", "odour", "firmness")

trait_list = getTraitList(dat, c("_pos", "_neg"),
                           trait.labels = ClimMobTools:::.title_case(trait))

trait_list

reference_trait_index = grep("overall", trait)

# put overall as the last trait in the list
o = rev(union(reference_trait_index, 1:length(trait_list)))

trait_list = trait_list[o]

trait = trait[o]

reference_trait_index = grep("overall", trait)

# make the PlackettLuce rankings
# since some of the traits don't have a fully network 
# this needs to be done using a loop instead of lapply()
R = list()
for(i in seq_along(trait)){
  R[[i]] = rank_tricot(dat[trait_list[[i]]$keep, ], 
                       pack, 
                       trait_list[[i]]$string)
}

# ....................................
# ....................................
# Kendall tau ####
kendall = lapply(trait_list[-reference_trait_index], function(x){
  # update the vector keep to match with dimensions from reference trait and 
  # the trait 'x' applied in this function 
  k = trait_list[[reference_trait_index]]$keep & x$keep
  
  r1 = rankTricot(dat[k, ],
                   items = pack,
                   input = trait_list[[reference_trait_index]]$string)
  
  r2 = rankTricot(dat[k, ],
                   items = pack,
                   input = x$string)
  
  kendall = kendallTau(r1, r2)
  
  kendall
  
})

kendall = do.call("rbind", kendall)

kendall$Trait = ClimMobTools:::.title_case(trait[-reference_trait_index])

kendall = kendall[rev(order(kendall$kendallTau)), ]

kendall

write.csv(kendall[-2], "output/cassava-kendall-cor.csv", row.names = FALSE)

# ..........................
# network representation #####
plot(network(R[[reference_trait_index]]))

pdf("output/cassava-experimental-network.pdf",
    width = 10,
    height = 10)
plot(network(R[[reference_trait_index]]))
dev.off()

# ..........................
# worth map ####
mod = lapply(R, PlackettLuce)

mod

tested = matrix(NA, nrow = length(items), ncol = length(mod), 
                dimnames = list(items, trait))

for(i in seq_along(mod)) {
  tested[names(coef(mod[[i]])), i] = 1
}

tested

worth = worth_map(mod, 
                   labels = ClimMobTools:::.title_case(trait),
                  na.replace = F) +
  labs(x = "Genotype", y = "Trait") +
  scale_fill_distiller(palette = "BrBG", 
                       direction = 1, 
                       na.value = "white", 
                       name = "")

worth

ggsave("output/cassava-worth-map.pdf",
       plot = worth,
       width = 15,
       height = 10,
       dpi = 400,
       units = "cm")

#..........................................................
# Reliability ####

# do it for all the data and also by country
# for some reason if I put the another vector for all the data 
# the reliability function breaks, since I don't have to 
# investigate I will duplicate the code
dat_split = list(all = rep(TRUE, nrow(dat)),
                 nig = dat$country == "Nigeria",
                 cam = dat$country == "Cameroon")


# run reliability over the different check varieties
items = sort(unique(unlist(dat[pack])))

checks = c("Akpu", "TMEB2", "TMEB3", "TMEB1", "Madame", "Sape")

all(checks %in% items)

trait = ClimMobTools:::.title_case(trait)

# run over the split list
for (s in seq_along(dat_split)) {
  
  # since some of the traits don't have a fully network 
  # this needs to be done using a loop instead of lapply()
  R_s = list()
  
  for(r in seq_along(trait)){
    
    k = trait_list[[r]]$keep & dat_split[[s]]
    
    if(sum(k) == 0) next
    
    R_s[[r]] = rank_tricot(dat[k, ], 
                           pack, 
                           trait_list[[r]]$string)
  }
  
  mod_s = lapply(R_s, function(x){
    try(PlackettLuce(x),
        silent = TRUE)
  })
  
  rel = data.frame()
  
  for (i in seq_along(checks)) {
    
    rel_i = lapply(mod_s, function(x){
      
      r = try(reliability(x, ref = checks[i]), silent = TRUE)
      
      if ("try-error" %in% class(r)) {
        return()
      }
      
      r$Check = checks[i]
      
      r
      
    })
    
    for(j in seq_along(trait)) {
      
      if (is.null(rel_i[[j]])) next
      
      rel_i[[j]]$Trait = trait[j]
      
    }
    
    rel_i = do.call("rbind", rel_i)
    
    rel = rbind(rel, rel_i)
    
  }
  
  # put traits in the right order
  rel$Trait = factor(rel$Trait, levels = trait)
  
  rel$improvement = round((rel$reliability / 0.5 - 1), 2)
  
  rel = rel[!rel$item %in% checks, ]
  
  relplot = 
    ggplot(data = rel,
           aes(x = improvement, 
               y = item, 
               group = Trait,
               fill = Trait))+
    geom_bar(stat = "identity",
             width = 0.7,
             position = "dodge", 
             show.legend = TRUE) +
    scale_fill_brewer(palette = 'BrBG') +
    geom_vline(xintercept = 0,
               colour = "#1f78b4",
               linewidth = 1) +
    facet_wrap(~ Check, strip.position = "top") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          strip.background =element_rect(fill="white"),
          text = element_text(color = "grey20"),
          strip.background.x = element_blank(),
          strip.placement = "outside",
          legend.position = "bottom",
          strip.text = element_text(size = 12, color = "grey20"),
          legend.text = element_text(size = 12, color = "grey20"),
          axis.text = element_text(size = 12, color = "grey20"),
          axis.title = element_text(size = 12, color = "grey20"),
          legend.title = element_blank()) +
    labs(x = "Probability of outperforming",
         y = "")
  
  ggsave(paste0("output/cassava-reliability-",
                names(dat_split)[s],
                ".pdf"),
                plot = relplot, 
                width = 22,
                height = 18,
                units = "cm",
                dpi = 400)
  
  rm(rel, relplot, R_s, mod_s)
  
}

#..........................................................
# PlackettLuce with aggregated rankings
# this put the rankings from all traits into a single
# grouped rankings to assess "overall technology performance"
# reference trait must be the first in this vector
names(trait_list) = trait

othertraits = union(trait[reference_trait_index],
                     trait[-reference_trait_index])

indicesbase = as.vector(which(trait_list[[reference_trait_index]]$keep))
resetindices = 1:length(indicesbase)

RG = list()

index = c()

for(i in seq_along(othertraits)) {

  trait_i = which(names(trait_list) %in% othertraits[i])

  # this should be combined with the baseline trait
  index_i = as.vector(which(trait_list[[trait_i]]$keep))

  keep_i = index_i %in% indicesbase

  index_i = index_i[keep_i]

  r_i = rankTricot(dat[index_i, ],
                    pack,
                    trait_list[[trait_i]]$string,
                    group = FALSE)

  # reset indices to match with grouped_rankings later
  index_i = resetindices[indicesbase %in% index_i]

  index = c(index, index_i)

  RG[[i]] = r_i

}

# make weights based on response
weight = as.vector(table(index))
weight = weight / max(weight)

# RG = do.call("rbind", RG)

RG = group(RG[[1]], index = indicesbase)

RG

# ...............................
# run the PLADMM model
names(features)

sel <- c("genotype", 'cohesiveness', 'adhesiveness_g_sec_1', 'hardness_g', 
         'gumminess', 'chewiness', "resilience_percent", "springiness_percent",
         "percent_swelling_power" , "percent_soluble", 
         'color_l', 'color_a', 'color_b')

features <- features[sel]

features

names(features)[-1] <- c("Cohesiveness", "Adhesiveness", "Hardness", 
                         "Gumminess", "Chewiness", "Resilience", "Springiness",
                        "SwellingPower", "Soluble", 
                         "L", "a", "b")

features

names(features)

# ..........................
# Fit the model using the PL tree 
names(dat)

covar = c("country", "consumption", "gender", "age")

covar = dat[, covar]

unique(covar$consumers_gender)
table(covar$consumers_gender)
hist(covar$age)

covar[c(1:3)] = lapply(covar[c(1:3)], as.factor)

str(covar)

names(covar) = c("Country", "Consumption", "Gender", "Age")

# get the ranking as grouped rankings
G = group(R[[6]], index = indicesbase)
pld = cbind(genotype = RG, covar, one = rep(1, length(indicesbase)))

head(pld)

# do the model combining different types of covariates
models <- list(inst_texture = c("Cohesiveness", "Adhesiveness", "Hardness", 
                                "Gumminess", "Chewiness", "Resilience", "Springiness"),
               pasting = c("SwellingPower", "Soluble"),
               color = c("L", "a", "b"))


models

# fit the pladmm using the group rankings of the three 
# traits assessed by participants
for (i in seq_along(models)) {
  
  form = formula(paste0("~ ", paste(models[[i]], collapse = " + ")))
  
  print(form)
  
  # model without participant covariates
  mod = pltree(genotype ~ one,
               worth = form,
               data = list(pld,
                           features))

  mod_coeff = pladmm_coeffs(mod)

  mod_coeff

  mod_coeff$Estimate[mod_coeff[,2] ==  "(Intercept)"] = 0
  
  # model with participant covariates
  mod3 = pltree(genotype ~ Country,
                alpha = 0.1,
                worth = form,
                data = list(pld,
                            features))
  
  mod3_coef = pladmm_coeffs(mod3)
  
  rules = node_rules(mod3)
  names(rules) = c("Node", "Group")
  
  mod3_coef = merge(rules, mod3_coef, by = "Node", all.y = TRUE)
  
  names(mod3_coef)[3] = ""
  names(mod3_coef)[8] = ""
  
  mod3_coef$Estimate[mod3_coef[, 3]  ==  "(Intercept)"] = 0
  
  write.csv(mod3_coef,
            file = paste0("output/cassava-pltree-consumer-choices-", 
                   names(models[i]), ".csv"),
            row.names = FALSE)
  
  write.csv(mod_coeff[, -1],
            paste0("output/cassava-pladmm-coeffs-overall-rankings-",
                   names(models[i]), ".csv"),
            row.names = FALSE)
  
}

# color has a computational issue in the three so I run it independently
form = formula(paste0("~ ", paste(models[[3]], collapse = " + ")))

print(form)

# for cameroon
itemsCam = dimnames(unclass(RCam))[[2]]

GCam = group(RCam, index = 1:nrow(RCam))

featuresCam = features[features$genotype %in% itemsCam, ]

pldCam = data.frame(genotype = GCam, one = 1)

modCam = pltree(genotype ~ one,
                worth = form,
                data = list(pldCam,
                            featuresCam))

mod_coeffCam = pladmm_coeffs(modCam)

mod_coeffCam

mod_coeffCam$Estimate[mod_coeffCam[,2] ==  "(Intercept)"] = 0

mod_coeffCam$Node = 2

# for nigeria
itemsNig = dimnames(unclass(RNig))[[2]]

GNig = group(RNig, index = 1:nrow(RNig))

featuresNig = features[features$genotype %in% itemsNig, ]

pldNig = data.frame(genotype = GNig, one = 1)

modNig = pltree(genotype ~ one,
                worth = form,
                data = list(pldNig,
                            featuresNig))

mod_coeffNig = pladmm_coeffs(modNig)

mod_coeffNig

mod_coeffNig$Estimate[mod_coeffNig[,2] ==  "(Intercept)"] = 0

mod_coeffNig$Node = 3

modCol = rbind(mod_coeffCam, mod_coeffNig)

write.csv(modCol,
          file = paste0("output/cassava-pltree-consumer-choices-", 
                        names(models[3]), ".csv"),
          row.names = FALSE)

# ...................................
# ...................................
# ...................................
# fit a simple plackett-luce tree
mod2 = pltree(genotype ~ .,
              alpha = 0.1,
              minsize = 50, 
              data = pld,
              gamma = TRUE,
              verbose = TRUE)

mod2

# extract the node models to build the tree using 
nodes = predict(mod2, type = "node")

table(nodes)

table(dat$country)

node_ids = sort(unique(nodes))

models = list()
nobs = integer()

for (i in seq_along(node_ids)) {
  x = pld[nodes  ==  node_ids[i], ]
  nobs = cbind(nobs, nrow(x))
  x = PlackettLuce(x$genotype)
  models[[i]] = x
}

branch = gosset:::build_tree_branches(mod2)
node = gosset:::build_tree_nodes(models,
                                  ref = reference,
                                  log = F,
                                  node.ids = node_ids, 
                                  multcomp = FALSE,
                                  n.obs = nobs,
                                  ci.level = 0.9)

tree = branch / node

tree

ggsave("output/cassava-plackett-luce-tree.png",
       width = 20,
       height = 15,
       units = "cm",
       dpi = 400)

ggsave("output/cassava-plackett-luce-tree.pdf",
       width = 20,
       height = 15,
       units = "cm",
       dpi = 400)




# # ..................................
# # text analysis
# like = tolower(dat$reason_like)
# 
# like = removeWords(like, stopwords("english"))
# 
# like = removeNumbers(like)
# 
# like = gsub("sour", "not-sour", like)
# 
# like = gsub(",", "", like)
# 
# corpus = stripWhitespace(like)
# corpus = iconv(corpus)
# corpus = Corpus(VectorSource(corpus))
# 
# tdm = TermDocumentMatrix(corpus)
# tdm = as.matrix(tdm)
# 
# w = sort(rowSums(tdm), decreasing = TRUE)
# set.seed(2112)
# wordcloud(words = names(w),
#           freq = w,
#           max.words = 200,
#           random.order = F,
#           min.freq = 5,
#           colors = brewer.pal(8, 'Dark2'),
#           scale = c(5, 0.3),
#           rot.per = 0.7)
# 
# w = sort(rowSums(tdm), decreasing = TRUE)
# set.seed(2112)
# png(filename = "output/cassava-reason-like.png",
#     width = 9, 
#     height = 9,
#     units = "cm",
#     res = 600)
# wordcloud(words = names(w),
#           freq = w,
#           max.words = 200,
#           random.order = F,
#           min.freq = 5,
#           colors = brewer.pal(8, 'Dark2'),
#           scale = c(5, 0.3),
#           rot.per = 0.7)
# dev.off()
# 
# # dislike
# dislike = tolower(dat$reason_dislike)
# 
# dislike = removeWords(dislike, stopwords("english"))
# 
# dislike = removeNumbers(dislike)
# 
# dislike = gsub(",", "", dislike)
# 
# corpus = stripWhitespace(dislike)
# corpus = iconv(corpus)
# corpus = Corpus(VectorSource(corpus))
# 
# tdm = TermDocumentMatrix(corpus)
# tdm = as.matrix(tdm)
# 
# w = sort(rowSums(tdm), decreasing = TRUE)
# set.seed(2112)
# wordcloud(words = names(w),
#           freq = w,
#           max.words = 200,
#           random.order = F,
#           min.freq = 5,
#           colors = brewer.pal(8, 'Dark2'),
#           scale = c(5, 0.3),
#           rot.per = 0.7)
# 
# w = sort(rowSums(tdm), decreasing = TRUE)
# set.seed(2112)
# png(filename = "output/cassava-reason-dislike.png",
#     width = 9, 
#     height = 9,
#     units = "cm",
#     res = 600)
# wordcloud(words = names(w),
#           freq = w,
#           max.words = 200,
#           random.order = F,
#           min.freq = 5,
#           colors = brewer.pal(8, 'Dark2'),
#           scale = c(5, 0.3),
#           rot.per = 0.7)
# dev.off()




