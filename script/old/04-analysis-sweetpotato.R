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
pladmm_coeffs = function(object, ...) {
  
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

# make PlackettLuce rankings
list.files("data", full.names = TRUE)

dat = read.csv("data/sweetpotato-tricot-data.csv")

names(dat)

features = read.csv("data/blups-sweetpotato.csv")

names(dat)

# check covariates
table(dat$gender)

dat$gender[dat$gender=="Female"] = "Woman"
dat$gender[dat$gender=="Male"] = "Man"

# ..................................
# text analysis
like = tolower(dat$reason_like)

like = removeWords(like, stopwords("english"))

like = removeNumbers(like)

like = gsub("sour", "not-sour", like)

like = gsub(",", "", like)

corpus = stripWhitespace(like)
corpus = iconv(corpus)
corpus = Corpus(VectorSource(corpus))

tdm = TermDocumentMatrix(corpus)
tdm = as.matrix(tdm)

w = sort(rowSums(tdm), decreasing = TRUE)
set.seed(2112)
wordcloud(words = names(w),
          freq = w,
          max.words = 200,
          random.order = F,
          min.freq = 5,
          colors = brewer.pal(8, 'Dark2'),
          scale = c(5, 0.3),
          rot.per = 0.7)

w = sort(rowSums(tdm), decreasing = TRUE)
set.seed(2112)
png(filename = "output/sweetpotato-reason-like.png",
    width = 9, 
    height = 9,
    units = "cm",
    res = 600)
wordcloud(words = names(w),
          freq = w,
          max.words = 200,
          random.order = F,
          min.freq = 5,
          colors = brewer.pal(8, 'Dark2'),
          scale = c(5, 0.3),
          rot.per = 0.7)
dev.off()

# dislike
dislike = tolower(dat$reason_dislike)

dislike = removeWords(dislike, stopwords("english"))

dislike = removeNumbers(dislike)

dislike = gsub(",", "", dislike)

corpus = stripWhitespace(dislike)
corpus = iconv(corpus)
corpus = Corpus(VectorSource(corpus))

tdm = TermDocumentMatrix(corpus)
tdm = as.matrix(tdm)

w = sort(rowSums(tdm), decreasing = TRUE)
set.seed(2112)
wordcloud(words = names(w),
          freq = w,
          max.words = 200,
          random.order = F,
          min.freq = 5,
          colors = brewer.pal(8, 'Dark2'),
          scale = c(5, 0.3),
          rot.per = 0.7)

w = sort(rowSums(tdm), decreasing = TRUE)
set.seed(2112)
png(filename = "output/sweetpotato-reason-dislike.png",
    width = 9, 
    height = 9,
    units = "cm",
    res = 600)
wordcloud(words = names(w),
          freq = w,
          max.words = 200,
          random.order = F,
          min.freq = 5,
          colors = brewer.pal(8, 'Dark2'),
          scale = c(5, 0.3),
          rot.per = 0.7)
dev.off()


# check the rankings
pack = c("option_a", "option_b", "option_c")

table(unlist(dat[pack]), rep(dat$country, 3))

names(dat)

trait = c("Sweetness", "Firmness", "PowderyTexture", "Overall")

trait_list = getTraitList(dat, c("_pos", "_neg"),
                           trait.labels = ClimMobTools:::.title_case(trait))


R = list()

for(i in seq_along(trait)){
  R[[i]] = rank_tricot(dat[trait_list[[i]]$keep, ], 
                        pack, 
                        trait_list[[i]]$string)
}

# compute kendall correlation
reference_trait_index = 4

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

kendall

write.csv(kendall[-2], "output/sweetpotato-kendall-cor.csv", row.names = FALSE)

# ..........................
# get the network representation
plot(network(R[[reference_trait_index]]))

png("output/sweetpotato-experimental-network.png",
    width = 20,
    height = 20,
    units = "cm",
    res = 400)
plot(network(R[[reference_trait_index]]))
dev.off()

# ..........................
# get worth map
reference = "NASPOT 8"

mod = lapply(R, PlackettLuce)

mod

worth = worth_map(mod, ref = reference,
                   labels = ClimMobTools:::.title_case(trait)) +
  labs(x = "Genotype", y = "Trait")

worth

ggsave("output/sweetpotato-worth-map.png",
       plot = worth,
       width = 15,
       height = 10,
       dpi = 400,
       units = "cm")

#..........................................................
# Reliability

# run reliability over the different check varieties
sort(unique(unlist(dat[pack])))

checks = c("NASPOT 8", "LOCAL CHECK", "NAROSPOT 1")

trait = ClimMobTools:::.title_case(trait)

rel = data.frame()

for (i in seq_along(checks)) {
  
  rel_i = lapply(mod, function(x){
    
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

# compute the improvement
rel$improvement = round(rel$reliability / 0.5 - 1, 2)

# remove the checks
rel = rel[!rel$item %in% checks, ]

# put traits in the right order
rel$Trait = factor(rel$Trait, levels = trait)

relplot = ggplot(data = rel,
                 aes(x = reliability, 
                     y = Trait,
                     color = Check,
                     shape = Check)) +
  geom_vline(xintercept = 0.5, 
             colour = "#de2d26", linewidth = 0.5) +
  scale_x_continuous(limits=c(0, 1)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap( ~ item) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        strip.background =element_rect(fill="white"),
        text = element_text(color = "grey20"),
        legend.position = "bottom") +
  labs(x = "Reliability",
       y = "Trait")

relplot

ggsave("output/sweetpotato-reliability.png",
       plot = relplot, 
       width = 20,
       height = 18,
       units = "cm",
       dpi = 400)

#..........................................................
# PlackettLuce with aggregated rankings
# this put the rankings from all traits into a single
# grouped rankings to assess "overall technology performance"
# reference trait must be the first in this vector
# names(trait_list) = trait
# 
# othertraits = union(trait[reference_trait_index],
#                      trait[-reference_trait_index])
# 
# indicesbase = as.vector(which(trait_list[[reference_trait_index]]$keep))
# resetindices = 1:length(indicesbase)
# 
# RG = list()
# 
# index = c()
# 
# for(i in seq_along(othertraits)) {
# 
#   trait_i = which(names(trait_list) %in% othertraits[i])
# 
#   # this should be combined with the baseline trait
#   index_i = as.vector(which(trait_list[[trait_i]]$keep))
# 
#   keep_i = index_i %in% indicesbase
# 
#   index_i = index_i[keep_i]
# 
#   r_i = rankTricot(dat[index_i, ],
#                     pack,
#                     trait_list[[trait_i]]$string,
#                     group = FALSE)
# 
#   # reset indices to match with grouped_rankings later
#   index_i = resetindices[indicesbase %in% index_i]
# 
#   index = c(index, index_i)
# 
#   RG[[i]] = r_i
# 
# }
# 
# # make weights based on response
# weight = as.vector(table(index))
# weight = weight / max(weight)
# 
# # RG = do.call("rbind", RG)
# RG = RG[[1]]

RG = group(R[[reference_trait_index]], index = 1:length(R[[reference_trait_index]]))

RG

# ...............................
# run the PLADMM model
names(features)




# ...............................
# run the PLADMM model
names(features)

sel <- c("genotype", "sweet_taste", 
         "crumbliness_mealiness_by_hand",
         "firmness_hardness", "fibrousness", "fibrous_appearance")

features <- features[sel]

features

names(features)[-1] <- c("SweetTaste", 
                         "Mealiness", "Firmness", "Fibrousness",
                         "FibrousAppearance")

features

names(features)

# ..........................
# Fit the model using the PL tree 
names(dat)

covar = c("gender", "district")

covar = dat[trait_list[[reference_trait_index]]$keep, covar]

head(covar)

table(covar$district)

covar$district[grepl("Adjumani", covar$district)] = "Adjumani"

covar[c(1:2)] = lapply(covar[c(1:2)], as.factor)

str(covar)

names(covar) = c("Gender", "District")

# get the ranking as grouped rankings
G = RG
pld = cbind(genotype = G, covar, one = rep(1, length(G)))

head(pld)

# do the model combining different types of covariates
models <- list(inst_texture = c("SweetTaste", 
                                "Mealiness", "Firmness", "Fibrousness",
                                "FibrousAppearance"))


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
  mod3 = pltree(genotype ~ .,
                alpha = 0.17,
                minsize = 50,
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
            file = paste0("output/sweetpotato-pltree-consumer-choices-", 
                          names(models[i]), ".csv"),
            row.names = FALSE)
  
  write.csv(mod_coeff[, -1],
            paste0("output/sweetpotato-pladmm-coeffs-overall-rankings-",
                   names(models[i]), ".csv"),
            row.names = FALSE)
  
}

# fit a simple plackett-luce tree
mod2 = pltree(genotype ~ District,
              alpha = 0.17,
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
                                 ci.level = 0.8)

tree = branch / node

tree

ggsave("output/sweetpotato-plackett-luce-tree.png",
       width = 20,
       height = 15,
       units = "cm",
       dpi = 400)

