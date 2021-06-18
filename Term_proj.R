#Term project
data<-read.delim("PLoS1MicrobiotaResponsetoDiet.txt")


#disregard the phylotype measurements and use the genotype level measurements
genotype_data<- data[1:156,1040:1169]
#trims the data down a lot, but seeing as the dataset was monstrously large
#to begin with, I think its reasonable

#proportional genotype data not working
#prop_genotype_data <- cbind(genotype_data)
#prop_genotype_data<- apply(prop_genotype_data,MARGIN=c(1,2), function(x) x^10)
#sums <- rowSums(prop_genotype_data)
#prop_genotype_data<- cbind(sums, prop_genotype_data)
#prop_genotype_data<- apply(prop_genotype_data, MARGIN=2, function(x) x/prop_genotype_data$sums)
#prop_genotype_data[,1:130] = apply(prop_genotype_data[,1:130],1,function(x){x/sum(x)})

#add the ID column
id <- data$id[1:156]
genotype_data<- cbind(id, genotype_data)
#add a study col 
Study <- substring(genotype_data$id, 1, 1)
genotype_data<-cbind(Study,genotype_data)

#add column and make pre/post datasets
preOrPost <- str_sub(genotype_data$id, start= -1)
genotype_data <- cbind(preOrPost, genotype_data)
split_geno_dat<- split(genotype_data, genotype_data$preOrPost)
pre_geno_dat<- as.data.frame(split_geno_dat[1])
post_geno_dat<- as.data.frame(split_geno_dat[2])

#anosim
ano1 = anosim(genotype_data[,4:132], genotype_data$preOrPost, 
              distance = "bray", permutations = 9999)
ano2 = anosim(pre_geno_dat[,4:132], pre_geno_dat$X1.Study, 
              distance = "bray", permutations = 9999)
ano3 = anosim(post_geno_dat[,4:132], post_geno_dat$X2.Study, 
              distance = "bray", permutations = 9999)
ano4 = anosim(pre_geno_dat[13:78,4:132], pre_geno_dat$X1.Study[13:78], 
              distance = "bray", permutations = 9999)
ano5 = anosim(post_geno_dat[13:78,4:132], post_geno_dat$X2.Study[13:78], 
              distance = "bray", permutations = 9999)

#NMDS
#community matrix
com1 = as.matrix(genotype_data[,4:133])
nmds1 = metaMDS(com1, distance = "bray")
plot(nmds1)
data.scores = as.data.frame(scores(nmds1))
data.scores$preOrpost = genotype_data$preOrPost
data.scores$study = genotype_data$Study

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = study, colour = preOrpost))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Pre or Post Intervention", y = "NMDS2", shape = "Study")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00")) 



#species difference
abund = genotype_data[ ,4:133]
PreOrPost = genotype_data$preOrPost 
inv = multipatt(abund, PreOrPost, func = "r.g", control = how(nperm=9999))

#C
abund2 = as.matrix(genotype_data[131:156,4:133])
PreOrPost2 = genotype_data$preOrPost[131:156]
inv2 = multipatt(abund2, PreOrPost2, func = "r.g", control = how(nperm=9999))

#AB
abund3 = as.matrix(genotype_data[27:130,4:133])
PreOrPost3 = genotype_data$preOrPost[27:130]
inv3 = multipatt(abund3, PreOrPost3, func = "r.g", control = how(nperm=9999))

#post
abund4 = as.matrix(post_geno_dat [,4:133])
study = post_geno_dat$X2.Study
inv4 = multipatt(abund4, study, func = "r.g", control = how(nperm=9999))



#pca's
pca_pre <- princomp(genotype_data[,4:132])

#nMDS w/ BC
vare.dis <- vegdist(genotype_data[4:132])
vare.mds0 <- isoMDS(vare.dis)
ordiplot(vare.mds0, type = "p", trace =TRUE)

#pre
NMDS_pre <- metaMDS(pre_geno_dat[3:132], k = 2,trace = F, 
                 autotransform = FALSE, distance="bray")
plot(NMDS_pre)
plot(NMDS_pre, display = "sites",type = "p")


pre_geno_dat[,3:132] %>%
  metaMDS(trace = F,distance = "bray") %>%
  ordiplot(type = "p")

#post
NMDS_post <- metaMDS(post_geno_dat[3:132], k = 2,trace = F, 
                    autotransform = FALSE, distance="bray")
plot(NMDS_post)
plot(NMDS_post, display = "sites",type = "p")

post_geno_dat[,3:132] %>%
  metaMDS(trace = F,distance = "bray") %>%
  ordiplot(type = "none") %>%
  text("sites") %>%
  plot()
  type = "p"

  