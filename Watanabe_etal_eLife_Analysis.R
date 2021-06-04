### << NOVEL NEUROANATOMICAL INTEGRATION AND SCALING DEFINE AVIAN BRAIN SHAPE EVOLUTION AND DEVELOPMENT >> ###
### << WATANABE ET AL. 2021 ELIFE >> ###
### << R CODE FOR PERFORMING ANALYSES >> ###
### << UPDATED 2 JUNE 2021 >> ###

rm(list=ls())

### << LIBRARIES >> ###

library(abind); library(ape); library(BTRTools); library(EMMLiv2); library(geiger); library(geomorph); library(geoscale); library(ggplot2); library(Morpho); library(paleomorph); library(phytools); library(plotrix); library(qgraph); library(RColorBrewer); library(RRPP)

### << FUNCTIONS >> ###

plotNetwork2 <- function (rhos, module_names = NULL, linecolour = "#7C7444",  # modified EMMLiv2::plotNetwork function ##
                          title = NULL, layout = NULL) 
{
  withins <- grep("Module*", colnames(rhos))
  nmodule <- length((grep("Module*", colnames(rhos))))
  rholist <- t(rhos)
  words <- strsplit(rownames(rholist), " ")
  plotcorr <- matrix(data = NA, nrow = nmodule, ncol = nmodule)
  modnums <- unlist(lapply(withins, function(x) strsplit(colnames(rhos)[x], 
                                                         " ")[[1]][2]))
  mods <- sapply(words[withins], function(x) paste(x, collapse = " "))
  mods <- gsub("Module", "M", mods)
  mods <- mods[order(sapply(mods, function(x) as.numeric(strsplit(x, 
                                                                  "M ")[[1]][[2]])))]
  colnames(plotcorr) <- rownames(plotcorr) <- mods
  for (i in 1:(length(words) - 1)) {
    if (length(words[[i]]) == 2) {
      module <- paste("M", words[[i]][2])
      plotcorr[module, module] <- rholist[i, "MaxL_p"]
    }
    if (length(words[[i]]) == 3) {
      from_module <- paste("M", (words[[i]][1]))
      to_module <- paste("M", words[[i]][3])
      plotcorr[from_module, to_module] <- rholist[i, "MaxL_p"]
      plotcorr[to_module, from_module] <- rholist[i, "MaxL_p"]
    }
  }
  within <- diag(plotcorr)
  between <- plotcorr
  mod.names <- module_names
  if (linecolour == "viridis") {
    cls <- viridis::viridis(100, direction = -1)
    vcols <- apply(between * 100, 1, function(x) cls[x])
    linecolour <- NULL
  }
  else {
    vcols <- NULL
  }
  qgraph::qgraph(between, shape = "circle", posCol = linecolour, 
                 edge.color = vcols, labels = mod.names, vsize = within * 
                   10, diag = FALSE, title = title, layout = layout, cut=mean(rhos[2,5:10]))
}

### << INITIALIZATION >> ###
{
setwd("...") # !! set working directory to where the coordinate data file is located !! #
left.lm <- c(1:54,109:137,168:170,173:175,178:180,183:185,188:190,193:195,198:200,203:205,208:210,213:215,218:220,222:224)
}

### << DATA >> ###
 ## < Coordinate Data > #
{
 slid.coord <- read.table("./Watanabe_etal_eLife_SIData.txt", row.names=1)
 slid.coord <- arrayspecs(slid.coord, ncol(slid.coord)/3, 3)
 alpha.order <- order(rownames(two.d.array(slid.coord)))
 slid.coord.alpha <- slid.coord[,,alpha.order]
}
 
 ## < Phylogenetic Tree > ##
 phy <- read.tree(text="((Alligator_Hatchling_3_endocast:0.0000001,Alligator_Hatchling_4_endocast:0.0000001,Alligator_Yearling_3_endocast:0.0000001,Alligator_Yearling_4_endocast:0.0000001,Alligator_Yearling_5_endocast:0.0000001,Alligator_Adult_endocast:0.0000001,Alligator_Hatchling_1_endocast:0.0000001,Alligator_Hatchling_2_endocast:0.0000001,Alligator_Yearling_1_endocast:0.0000001,Alligator_Yearling_2_endocast:0.0000001,Alligator_Juvenile_1_endocast:0.0000001,Alligator_Subadult_1_endocast:0.0000001,Alligator_Subadult_2_endocast:0.0000001,Alligator_Juvenile_2_endocast:0.0000001):251.2,(Alioramus_altai:129.2,((Incisivosaurus_gauthieri:28.2489,(Khaan_mckennai:42.07445,Citipati_osmolskae:42.07445):42.07445):28.2489,((Zanabazar_junior:47.8,Inga:47.8):47.8,(Archaeopteryx_lithographica:8.45,((Struthio_camelus:50.5891,(Rhea_americana:43.7645,(Apteryx_sp:41.6828,(((Casuarius_casuarius:8.355686,Casuarius_unappendiculatus:8.355686):0.983514,Dromaius_novaehollandiae:9.3392):29.8853,Crypturellus_tataupa:39.2245):2.4583):2.0817):6.8246):22.3181,(((Chicken_Adult_2_endocast:0.0000001,Chicken_1day_3_endocast:0.0000001,Chicken_1day_4_endocast:0.0000001,Chicken_8week_1_endocast:0.0000001,Chicken_8week_2_endocast:0.0000001,Chicken_1day_1_endocast:0.0000001,Chicken_1day_2_endocast:0.0000001,Chicken_1week_1_endocast:0.0000001,Chicken_1week_2_endocast:0.0000001,Chicken_3week_1_endocast:0.0000001,Chicken_3week_2_endocast:0.0000001,Chicken_6week_1_endocast:0.0000001,Chicken_6week_2_endocast:0.0000001,Chicken_Adult_1_endocast:0.0000001):54.6117,(Chauna_chavaria:46.2261,(Anas_platyrhynchos:5.810356,(Tachyeres_leucocephalus:0.661619,(Tachyeres_brachypterus:0.025335,Tachyeres_pteneres:0.025335):0.636284):5.148737):40.415744):8.3856):17.2449,(Chordeiles_minor:67.4195,((Ptilinopus_melanospilus:21.3604,(Raphus_cucullatus:15.1,Caloenas_nicobarica:15.1):6.2604):45.1467,(((Gallirallus_australis:6.987276,(Gallirallus_philippensis:1.444393,Gallirallus_rovianae:1.444393):5.542883):33.420824,(Podilymbus_podiceps:20.20405,Grus_canadensis:20.20405):20.20405):25.3318,(((Alca_torda:9.433,Pinguinus_impennis:9.433):54.5927,(Phaethon_rubricauda:63.1932,(Gavia_immer:62.2636,(Eudyptes_sp:61.3634,((Fregata_magnificens:26.2425705,Diomedea_sp:26.2425705):27.0041295,(Phalacrocorax_harrisi:30.604397,Phalacrocorax_penicillatus:30.604397):22.642303):8.1167):0.9002):0.9296):0.8325):1.0655,((Coragyps_atratus:56.8919,Haliaeetus_leucocephalus:56.8919):6.018,((Bucorvus_abyssinicus:56.9319,Melanerpes_aurifrons:56.9319):3.426,(Cariama_cristata:58.7846,(Nestor_meridionalis:32.7877,Brotogeris_chrysoptera:32.7877):25.9969):1.5733):2.552):2.1813):0.6487):0.7672):0.9124):4.4371):1.0506):83.340617):8.45):16.8):16.8):52.902183);")
 
 ## < Shape Data > ##
 {  
  labels <- dimnames(slid.coord.alpha)[[3]]
  genus.labels <- gsub("_[A-z,1-9]+","", labels)
  coord.data <- slid.coord.alpha[left.lm,,]
  gpa <- gpagen(coord.data, ProcD=F)
  shape.data <- gpa$coords
  cs <- gpa$Csize
  log.cs <- log(cs)
 }

### << TAXONOMIC ASSIGNMENT >> ###
{
dinos <- c(2,42,55:57,72)
aves <- c(1,17,18:26,41,43:54,58:71)
neornith <- c(1,17,18,20:26,41,43:54,58:71)
gator.dev <- 3:16
chick.dev <- 27:40
dev <- c(gator.dev, chick.dev)
no.dev <- c(3,dinos,aves,40)
phy <- drop.tip(phy, labels[-no.dev])

shape.dinos <- shape.data[,,dinos]
shape.aves <- shape.data[,,aves]
shape.neornith <- shape.data[,,neornith]
shape.gator <- shape.data[,,gator.dev]
shape.chick <- shape.data[,,chick.dev]
shape.nodev <- shape.data[,,no.dev]

cs.dinos <- cs[dinos]
cs.aves <- cs[aves]
cs.neornith <- cs[neornith]
cs.gator <- cs[gator.dev]
cs.chick <- cs[chick.dev]

tax.code <- rep(NA, dim(shape.data)[3])
tax.code[c(gator.dev,dinos)] <- "na.arch"  # non-avialan archosaur
tax.code[c(chick.dev, neornith)] <- "birds"
tax.code[19] <- "birds" # or remove from analysis

tax.dev.code <- rep(NA, dim(shape.data)[3])
tax.dev.code[dinos] <- "dinos"
tax.dev.code[neornith] <- "birds"
tax.dev.code[gator.dev] <- "gator.dev"
tax.dev.code[chick.dev] <- "chick.dev"
tax.dev.code[19] <- "archaeopteryx"

dev.code <- rep("no.dev", dim(shape.data)[3])
dev.code[c(chick.dev, gator.dev)] <- "dev"
}

### << REGIONAL ASSIGNMENT >> ###
{
cerebrum <- 1:54
opticlobe <- 55:83
cerebellum <- 84:101
medulla <- 102:119
part.cerebrum <- rep("1", length(cerebrum))
part.opticlobe <- rep("2", length(opticlobe))
part.cerebellum <- rep("3", length(cerebellum))
part.medulla <- rep("4", length(medulla))
partitions <- c(part.cerebrum, part.opticlobe, part.cerebellum, part.medulla)
}

### << REGIONAL ALIGNMENT >> ###
{
cerebrum.shape <- gpagen(coord.data[cerebrum,,])$coords
opticlobe.shape <- gpagen(coord.data[opticlobe,,])$coords
cerebellum.shape <- gpagen(coord.data[cerebellum,,])$coords
medulla.shape <- gpagen(coord.data[medulla,,])$coords

cerebrum.cs <- gpagen(coord.data[cerebrum,,])$Csize
opticlobe.cs <- gpagen(coord.data[opticlobe,,])$Csize
cerebellum.cs <- gpagen(coord.data[cerebellum,,])$Csize
medulla.cs <- gpagen(coord.data[medulla,,])$Csize
}

### << MORPHOSPACE >> ###
pca.shape <- prcomp(two.d.array(shape.data[,,]))
#pca.shape <- prcomp(two.d.array(shape.data[cerebrum,,])) # for regional PCA
pc1.shape <- pca.shape$x[,1]
pc2.shape <- pca.shape$x[,2]
pc3.shape <- pca.shape$x[,3]
{
palette <- brewer.pal(n=3, name="Dark2")
mspace.palette <- c("white", palette[3], NA, palette[2], NA)
mspace.ptshape <- c(5,21)
ptshape.palette <- c("black", "black", palette[3], "black", palette[1]) 
ptsize <- c(6,6,4,6,4)
p <- ggplot(data=NULL, aes(x=pc1.shape, y=pc2.shape))
p <- p + geom_point(size=ptsize[factor(tax.dev.code)], shape=mspace.ptshape[factor(dev.code)], colour=ptshape.palette[factor(tax.dev.code)], fill=mspace.palette[factor(tax.dev.code)], stroke=2)
p <- p + xlab("PC1 (53.0%)") + ylab("PC2 (7.9%)")
p <- p + xlim(-0.3, 0.4) + ylim(-0.2, 0.23)
p <- p + geom_text(aes(label=labels, hjust=-0.1), size=3, colour="black")
p <- p + theme_bw() + coord_fixed(ratio=1) 
p <- p + theme(axis.title=element_text(size=18), axis.text=element_text(size=12), panel.border = element_blank(), panel.grid.major = element_line(colour="grey"), panel.grid.minor = element_line(colour="grey"), axis.line = element_line(colour = "black"), legend.title=element_blank()) + guides(size=guide_legend(override.aes=list(shape=21)))
p
}
 ## Visualize PC Shape Changes ##
gm.pca <- gm.prcomp(shape.data)
pc1.max <- gm.pca$shapes$shapes.comp1$max
pc1.min<-gm.pca$shapes$shapes.comp1$min
pc2.max<-gm.pca$shapes$shapes.comp2$max
pc2.min<-gm.pca$shapes$shapes.comp2$min
 # PC 1 shape change #
plotRefToTarget(pc1.max, pc1.min, method="points")  # black points = pc1.min 
plotRefToTarget(pc2.max, pc2.min, method="points")


### << SHAPE DIFFERENCE BETWEEN DINOS AND BIRDS >> ###
y <- shape.nodev[,,-which(dimnames(shape.nodev)[[3]] == "Alligator_Adult_endocast")]  # remove alligator
x <- tax.code[no.dev]
x <- x[-1]  # remove alligator
summary(procD.lm(y~x))
summary(procD.pgls(y~x, phy=drop.tip(phy, tip="Alligator_Adult_endocast")))

### << ANCESTRAL RECONSTRUCTIONS >> ###
library(Rphylopars)
anc.reco <- anc.recon(trait_data=two.d.array(shape.data[,,no.dev]), tree=phy, CI=FALSE)  
anc.reco <- arrayspecs(anc.reco, ncol(anc.reco)/3, 3)
phy$node.label <- dimnames(anc.reco)[[3]]
plot.phylo(phy, show.node.label=TRUE, cex=1)
ancshape.neornith <- anc.reco[,,9]
ancshape.aves <- anc.reco[,,8]
ancshape.archo <- anc.reco[,,1]
ancshape.coeluro <- anc.reco[,,2]
ancshape.pennarapt <- anc.reco[,,3]
ancshape.paraves <- anc.reco[,,6]
ancshape.theropod <- anc.reco[,,c(2,3,6)]

### << CALCULATE ANGLE BETWEEN TRAJECTORIES >> ###
archo.traj <- ancshape.neornith - ancshape.coeluro
chick.traj <- shape.data[,,40] - shape.data[,,27]
gator.traj <- shape.data[,,3] - shape.data[,,4]
angleTest(chick.traj, gator.traj)
angleTest(archo.traj, chick.traj)
angleTest(archo.traj, gator.traj)

shape.data.wAnc <- abind(shape.data, anc.reco, along=3)
pca.shape.wAnc <- prcomp(two.d.array(shape.data.wAnc))
archo.pc.traj <- pca.shape.wAnc$x[81,] - pca.shape.wAnc$x[74,]
chick.pc.traj <- pca.shape.wAnc$x[40,] - pca.shape.wAnc$x[27,]
gator.pc.traj <- pca.shape.wAnc$x[3,] - pca.shape.wAnc$x[4,]
angleTest(chick.pc.traj, gator.pc.traj)  # p-value < 0.01 meaning that observed diff is > diff under assumption of normalized n-vectors represented by n-dimensional hypersphere.
angleTest(archo.pc.traj, chick.pc.traj)
angleTest(archo.pc.traj, gator.pc.traj) 

### << ALLOMETRY >> ###
 cerebrum.cs.nodev <- cerebrum.cs[no.dev]
 opticlobe.cs.nodev <- opticlobe.cs[no.dev]
 cerebellum.cs.nodev <- cerebellum.cs[no.dev]
 medulla.cs.nodev <- medulla.cs[no.dev]
  
 ## Archosaur ##
 {
 cs.no.dev <- cs[no.dev]
 allom <- procD.lm(shape.nodev~log(cs.no.dev))
 physig <- physignal(shape.nodev, phy)
 ev.allom <- procD.pgls(shape.nodev~log(cs.no.dev), phy)
 }
 ## Non-Avian Dinosaurs ##
 shape.dino <- shape.data[,,dinos]
 cs.dino <- cs[dinos]
 allom.dino <- procD.lm(shape.dino~log(cs.dino))
 physig.dino <- physignal(shape.dino, drop.tip(phy, labels[-dinos]))
 ev.allom.dino <- procD.pgls(shape.dino~log(cs.dino), drop.tip(phy, labels[-dinos]))
  # Region (Globally Aligned) #
  physig.dino <- physignal(shape.dino[medulla,,], drop.tip(phy, labels[-dinos]))
  allom.dino <- procD.lm(shape.dino[medulla,,]~log(cs.dino))
  ev.allom.dino <- procD.pgls(shape.dino[medulla,,]~log(cs.dino), drop.tip(phy, labels[-dinos]))
  # Region (Locally Aligned) #
  medulla.dino <- medulla.shape[,,dinos]
  medulla.cs.dino <- medulla.cs[dinos]
  physig.dino <- physignal(medulla.dino, drop.tip(phy, labels[-dinos]))
  allom.dino <- procD.lm(medulla.dino~log(medulla.cs.dino))
  ev.allom.dino <- procD.pgls(medulla.dino~log(medulla.cs.dino), drop.tip(phy, labels[-dinos]))
 
 ## Neornithes ##
 shape.neornith <- shape.data[,,neornith]
 cs.neornith <- cs[neornith]
 allom.neornith <- procD.lm(shape.neornith~log(cs.neornith))
 physig.neornith <- physignal(shape.neornith, drop.tip(phy, labels[-neornith]))
 ev.allom.neornith <- procD.pgls(shape.neornith~log(cs.neornith), drop.tip(phy, labels[-neornith]))
  # Region (Globally Aligned) #
  medulla.cs.neornith <- medulla.cs[neornith]
  physig.neornith <- physignal(shape.neornith[medulla,,], drop.tip(phy, labels[-neornith]))
  allom.neornith <- procD.lm(shape.neornith[medulla,,]~log(cs.neornith))
  ev.allom.neornith <- procD.pgls(shape.neornith[medulla,,]~log(cs.neornith), drop.tip(phy, labels[-neornith]))
  # Region (Locally Aligned) #
  medulla.neornith <- medulla.shape[,,neornith]
  physig.neornith <- physignal(medulla.neornith, drop.tip(phy, labels[-neornith]))
  allom.neornith <- procD.lm(medulla.neornith~log(medulla.cs.neornith))
  ev.allom.neornith <- procD.pgls(medulla.neornith~log(medulla.cs.neornith), drop.tip(phy, labels[-neornith]))
  
 ## Cerebrum ##
 allom.cerebrum <- procD.lm(shape.nodev[cerebrum,,]~log(cs.no.dev))
 physig.cerebrum <- physignal(shape.nodev[cerebrum,,], phy)
 ev.allom.cerebrum <- procD.pgls(shape.nodev[cerebrum,,]~log(cs.no.dev), phy)
 
 allom.cerebrum <- procD.lm(cerebrum.shape[,,no.dev]~log(cerebrum.cs.nodev))
 physig.cerebrum <- physignal(cerebrum.shape[,,no.dev], phy)
 ev.allom.cerebrum <- procD.pgls(cerebrum.shape[,,no.dev]~log(cerebrum.cs.nodev), phy)
 
 ## Optic Lobe ##
 allom.opticlobe <- procD.lm(shape.nodev[opticlobe,,]~log(cs.no.dev))
 physig.opticlobe <- physignal(shape.nodev[opticlobe,,], phy)
 ev.allom.opticlobe <- procD.pgls(shape.nodev[opticlobe,,]~log(cs.no.dev), phy)
 
 allom.opticlobe <- procD.lm(opticlobe.shape[,,no.dev]~log(opticlobe.cs.nodev))
 physig.opticlobe <- physignal(opticlobe.shape[,,no.dev], phy)
 ev.allom.opticlobe <- procD.pgls(opticlobe.shape[,,no.dev]~log(opticlobe.cs.nodev), phy)

 ## Cerebellum ##
 allom.cerebellum <- procD.lm(shape.nodev[cerebellum,,]~log(cs.no.dev))
 physig.cerebellum <- physignal(shape.nodev[cerebellum,,], phy)
 ev.allom.cerebellum <- procD.pgls(shape.nodev[cerebellum,,]~log(cs.no.dev), phy)
 
 allom.cerebellum <- procD.lm(cerebellum.shape[,,no.dev]~log(cerebellum.cs.nodev))
 physig.cerebellum <- physignal(cerebellum.shape[,,no.dev], phy)
 ev.allom.cerebellum <- procD.pgls(cerebellum.shape[,,no.dev]~log(cerebellum.cs.nodev), phy)

 ## Medulla ##
 allom.medulla <- procD.lm(shape.nodev[medulla,,]~log(cs.no.dev))
 physig.medulla <- physignal(shape.nodev[medulla,,], phy)
 ev.allom.medulla <- procD.pgls(shape.nodev[medulla,,]~log(cs.no.dev), phy)
 
 allom.medulla <- procD.lm(medulla.shape[,,no.dev]~log(medulla.cs.nodev))
 physig.medulla <- physignal(medulla.shape[,,no.dev], phy)
 ev.allom.medulla <- procD.pgls(medulla.shape[,,no.dev]~log(medulla.cs.nodev), phy)

 ## Alligator Development ##
 summary(procD.lm(shape.gator~log(cs.gator)))
 summary(procD.lm(shape.gator[cerebrum,,]~log(cs.gator)))
 summary(procD.lm(shape.gator[opticlobe,,]~log(cs.gator)))
 summary(procD.lm(shape.gator[cerebellum,,]~log(cs.gator)))
 summary(procD.lm(shape.gator[medulla,,]~log(cs.gator)))
 
 cerebrum.gator <- cerebrum.shape[,,gator.dev]
 summary(procD.lm(cerebrum.gator~log(cs.gator)))
 opticlobe.gator <- opticlobe.shape[,,gator.dev]
 summary(procD.lm(opticlobe.gator~log(cs.gator)))
 cerebellum.gator <- cerebellum.shape[,,gator.dev]
 summary(procD.lm(cerebellum.gator~log(cs.gator)))
 medulla.gator <- medulla.shape[,,gator.dev]
 summary(procD.lm(medulla.gator~log(cs.gator)))
 
 ## Chicken Development ##
 summary(procD.lm(shape.chick~log(cs.chick)))
 summary(procD.lm(shape.chick[cerebrum,,]~log(cs.chick)))
 summary(procD.lm(shape.chick[opticlobe,,]~log(cs.chick)))
 summary(procD.lm(shape.chick[cerebellum,,]~log(cs.chick)))
 summary(procD.lm(shape.chick[medulla,,]~log(cs.chick)))
 
 cerebrum.chick <- cerebrum.shape[,,chick.dev]
 summary(procD.lm(cerebrum.chick~log(cs.chick)))
 opticlobe.chick <- opticlobe.shape[,,chick.dev]
 summary(procD.lm(opticlobe.chick~log(cs.chick)))
 cerebellum.chick <- cerebellum.shape[,,chick.dev]
 summary(procD.lm(cerebellum.chick~log(cs.chick)))
 medulla.chick <- medulla.shape[,,chick.dev]
 summary(procD.lm(medulla.chick~log(cs.chick)))

 ### TEST DIFFERENCE BETWEEN DINOSAUR VS. BIRD ALLOMETRY ###
 shape.no.dev <- shape.data[,,no.dev]
 shape.no.gator <- shape.no.dev[,,-c(1,11)]  # remove alligator and archaeopteryx
 cs.no.dev <- cs[no.dev]
 cs.no.gator <- cs.no.dev[-c(1,11)]
 tax.code.no.dev <- tax.code[no.dev]
 tax.code.no.gator <- tax.code.no.dev[-c(1,11)]
 phy.no.gator <- drop.tip(phy, labels[-no.dev])
 phy.no.gator <- drop.tip(phy.no.gator, c("Alligator_Adult_endocast", "Archaeopteryx_lithographica"))
 fit1 <- procD.lm(shape.no.gator~log(cs.no.gator)+tax.code.no.gator, RRPP=TRUE)
 fit2 <- procD.lm(shape.no.gator~log(cs.no.gator)*tax.code.no.gator, RRPP=TRUE)

 ### TEST DIFFERENCE BETWEEN ALLIGATOR VS. CHICKEN ALLOMETRY ###
 shape.dev <- shape.data[,,dev]
 cs.dev <- cs[dev]
 tax.code.dev <- tax.code[dev]
 fit1 <- procD.lm(shape.dev~log(cs.dev)+tax.code.dev, RRPP=TRUE)
 fit2 <- procD.lm(shape.dev~log(cs.dev)*tax.code.dev, RRPP=TRUE)
 
### TEST DIFFERENCE BETWEEN CHICKEN VS. BIRDS ###
 shape.bird.chick <- shape.data[,,c(neornith,chick.dev)]
 cs.bird.chick <- cs[c(neornith, chick.dev)]
 code.bird.chick <- rep("bird", length(cs.bird.chick))
 code.bird.chick[38:51] <- "chicken"
 fit <- procD.lm(shape.bird.chick~log(cs.bird.chick)*code.bird.chick, RRPP=TRUE)
 
### TEST DIFFERENCE BETWEEN ALLIGATOR VS. NON-AVIAN DINOSAURS ###
 shape.dino.gator <- shape.data[,,c(dinos, gator.dev)]
 cs.dino.gator <- cs[c(dinos, gator.dev)]
 code.dino.gator <- rep("dino", length(cs.dino.gator))
 code.dino.gator[7:20] <- "gator"
 fit <- procD.lm(shape.dino.gator~log(cs.dino.gator)*code.dino.gator, RRPP=TRUE)
  
 ## Supplementary ##
{
  p <- ggplot(data=NULL, aes(x=cac, y=-rsc1))
  p <- p + geom_point(size=3, shape=21, colour="black", fill="white", stroke=1)
  p <- p + geom_text(aes(label=fig.label,hjust=-0.15), size=4.5, colour="black")
  p <- p + xlab("Common Allometric Component") + ylab("PC1 Residual Shape")
  p <- p + xlim(-0.15, 0.2)
  p <- p + theme_bw()
  p <- p + theme(axis.title=element_text(size=18), axis.text=element_text(size=12), panel.border = element_blank(), panel.grid.major = element_line(colour="grey"), panel.grid.minor = element_line(colour="grey"), axis.line = element_line(colour = "black"), legend.title=element_blank())
  p <- p + guides(size=guide_legend(override.aes=list(shape=21)))
  p
}

 # Test chicken vs. birds #
rsc1.bird.chick <- rsc1[c(neornith,chick.dev)]
cac.bird.chick <- cac[c(neornith, chick.dev)]
code.bird.chick <- rep("bird", length(cs.bird.chick))
code.bird.chick[38:51] <- "chicken"
fit <- glm(rsc1.bird.chick~cac.bird.chick*code.bird.chick)

 # Test alligator vs. dinos #
rsc1.dino.gator <- rsc1[c(dinos, gator.dev)]
cac.dino.gator <- cac[c(dinos, gator.dev)]
code.dino.gator <- rep("dino", length(cs.dino.gator))
code.dino.gator[7:20] <- "gator"
fit <- lm(rsc1.dino.gator~cac.dino.gator*code.dino.gator)

 # Cerebrum #
cac <- CAC(medulla.shape, medulla.cs, log=TRUE)$CACscores
rsc1 <- CAC(medulla.shape, medulla.cs, log=TRUE)$RSCscores[,1]
{
  p <- ggplot(data=NULL, aes(x=cac, y=-rsc1))
  p <- p + geom_point(size=ptsize[factor(tax.dev.code)], shape=mspace.ptshape[factor(dev.code)], colour=ptshape.palette[factor(tax.dev.code)], fill=mspace.palette[factor(tax.dev.code)], stroke=2)
  p <- p + xlab("Common Allometric Component") + ylab("PC1 Residual Shape")
  p <- p + xlim(-0.35, 0.35) + ylim(-0.15, 0.35)
  #p <- p + geom_text(aes(label=labels, hjust=-0.1), size=3, colour="black")
  p <- p + theme_bw()
  p <- p + theme(aspect.ratio=1, axis.title=element_text(size=18), axis.text=element_text(size=12), panel.border = element_blank(), panel.grid.major = element_line(colour="grey"), panel.grid.minor = element_line(colour="grey"), axis.line = element_line(colour = "black"), legend.title=element_blank()) + guides(size=guide_legend(override.aes=list(shape=21)))
  p
}

 ## Supplementary ##
cac <- CAC(medulla.shape, medulla.cs, log=TRUE)$CACscores
rsc1 <- CAC(medulla.shape, medulla.cs, log=TRUE)$RSCscores[,1]
{
  p <- ggplot(data=NULL, aes(x=cac, y=-rsc1))
  p <- p + geom_point(size=3, shape=21, colour="black", fill="white", stroke=1)
  p <- p + geom_text(aes(label=fig.label,hjust=-0.15), size=4, colour="black")
  p <- p + xlab("Common Allometric Component") + ylab("PC1 Residual Shape")
  p <- p + xlim(-0.35, 0.4)
  p <- p + theme_bw() #+ coord_fixed(ratio=1)
  p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=15), panel.border = element_blank(), panel.grid.major = element_line(colour="grey"), panel.grid.minor = element_line(colour="grey"), axis.line = element_line(colour = "black"), legend.title=element_blank())
  p <- p + guides(size=guide_legend(override.aes=list(shape=21)))
  p
}

### << MODULARITY & INTEGRATION (CR METHOD) >> ### 
edge.col <- "#7C7444"  # Nature Ecol Evol color
brain.lab <- c("cerebrum", "optic lobe", "cerebellum", "medulla")
local.shape.data <- abind(cerebrum.shape, opticlobe.shape, cerebellum.shape, medulla.shape, along=1)
 ## Archosauria ##
archo.pgls <- procD.pgls(shape.nodev~log(cs.no.dev), phy=phy, iter=1)$pgls.residuals
archo.pgls <- arrayspecs(archo.pgls, ncol(archo.pgls)/3, 3)
 # For locally aligned shape data #
archo.pgls <- procD.pgls(local.shape.data[,,no.dev]~log(cs.no.dev), phy=phy, iter=1)$pgls.residuals
archo.pgls <- arrayspecs(archo.pgls, ncol(archo.pgls)/3, 3)

cr.archo <- modularity.test(archo.pgls, partition.gp=partitions, iter=9999)
qgraph(cr.archo$CR.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(cr.archo$CR.mat), title="CR values")
pls.archo <- integration.test(archo.pgls, partition.gp=partitions, iter=9999)
phy.pls.archo <- phylo.integration(shape.nodev, phy=phy, partition.gp=partitions, iter=9999)
qgraph(pls.archo$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(pls.archo$r.pls.mat), title="r-PLS values")

 ## Non-Avian Dinosaurs ##
dino.pgls <- procD.pgls(shape.dinos~log(cs.dinos), phy=drop.tip(phy, labels[-dinos]), iter=1)$pgls.residuals
dino.pgls <- arrayspecs(dino.pgls, ncol(dino.pgls)/3, 3)
 # For locally aligned shape data #
dino.pgls <- procD.pgls(local.shape.data[,,dinos]~log(cs.dinos), phy=drop.tip(phy, labels[-dinos]), iter=1)$pgls.residuals
dino.pgls <- arrayspecs(dino.pgls, ncol(dino.pgls)/3, 3)
 #
cr.dino <- modularity.test(dino.pgls, partition.gp=partitions, iter=9999)
qgraph(cr.dino$CR.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=median(cr.dino$CR.mat), title="CR values")
phy.pls.dino <- phylo.integration(shape.dinos, phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions, iter=9999)
phy.pls.dino <- phylo.integration(local.shape.data[,,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions, iter=9999)
qgraph(pls.dino$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(pls.dino$r.pls.mat), title="r-PLS values")
qgraph(phy.pls.dino$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(phy.pls.dino$r.pls.mat), title="r-PLS values")

 ## Neornithes ##
neornith.pgls <- procD.pgls(shape.neornith~log(cs.neornith), phy=drop.tip(phy, labels[-neornith]), iter=1)$pgls.residuals
neornith.pgls <- arrayspecs(neornith.pgls, ncol(neornith.pgls)/3, 3)
 # For locally aligned shape data #
neornith.pgls <- procD.pgls(local.shape.data[,,neornith]~log(cs.neornith), phy=drop.tip(phy, labels[-neornith]), iter=1)$pgls.residuals
neornith.pgls <- arrayspecs(neornith.pgls, ncol(neornith.pgls)/3, 3)
 #
cr.neornith <- modularity.test(neornith.pgls, partition.gp=partitions, iter=9999)
qgraph(cr.neornith$CR.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=median(cr.neornith$CR.mat), title="CR values")
pls.neornith <- integration.test(neornith.pgls, partition.gp=partitions, iter=9999)
phy.pls.neornith <- phylo.integration(shape.neornith, phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions, iter=9999)
phy.pls.neornith <- phylo.integration(local.shape.data[,,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions, iter=9999)
qgraph(pls.neornith$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(pls.neornith$r.pls.mat), title="r-PLS values")
qgraph(phy.pls.neornith$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(phy.pls.neornith$r.pls.mat), title="r-PLS values")

 ## Alligator Development ##
cr.gator <- modularity.test(shape.gator, partition.gp=partitions, iter=9999)
 # for locally aligned data #
cr.gator <- modularity.test(local.shape.data[,,gator.dev], partition.gp=partitions, iter=9999)
 #
qgraph(cr.gator$CR.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=median(cr.gator$CR.mat), title="CR values")
pls.gator <- integration.test(shape.gator, partition.gp=partitions, iter=9999)
pls.gator <- integration.test(local.shape.data[,,gator.dev], partition.gp=partitions, iter=9999)
qgraph(pls.gator$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(pls.gator$r.pls.mat), title="r-PLS values")

 ## Chicken Development ##
cr.chick<- modularity.test(shape.data[,,chick.dev], partition.gp=partitions, iter=9999)
 # for locally aligned data #
cr.chick<- modularity.test(local.shape.data[,,chick.dev], partition.gp=partitions, iter=9999)
 #
qgraph(cr.chick$CR.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=median(cr.chick$CR.mat), title="CR values")
pls.chick <- integration.test(shape.data[,,chick.dev], partition.gp=partitions, iter=9999)
pls.chick <- integration.test(local.shape.data[,,chick.dev], partition.gp=partitions, iter=9999)
qgraph(pls.chick$r.pls.mat, shape="circle", posCol=edge.col, labels=brain.lab, diag=FALSE, cut=mean(pls.chick$r.pls.mat), title="r-PLS values")

### << MODULARITY & INTEGRATION (EMMLI) >> ### 
parts <- partitions
emmli.parts <- t(rbind(as.character(1:dim(shape.data)[1]), parts))
colnames(emmli.parts) <- c("LM","modules")
emmli.parts <- as.data.frame(emmli.parts)
for(i in 1:ncol(emmli.parts)) {
  emmli.parts[,i] <- factor(emmli.parts[,i])
}

 ## Archosaurs ##
shape.archo.emmli <- shape.data[,,no.dev]
rownames(shape.archo.emmli) <- NULL
colnames(shape.archo.emmli) <- NULL
emmli.archo <- phyloEmmli(shape.archo.emmli, phylo=phy, method="pgls", EMMLi=TRUE, N_sample=length(no.dev), mod=emmli.parts)
emmplotNetwork2(rhos=emmli.archo$EMMLi$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

 ## Non-Avian Dinosaurs ##
shape.dino.emmli <- shape.dino
rownames(shape.dino.emmli) <- NULL
colnames(shape.dino.emmli) <- NULL
EMMLi.dinos <- phyloEmmli(shape.dino.emmli, phylo=drop.tip(phy, labels[-dinos]), method="pgls", EMMLi=TRUE, N_sample=length(dinos), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.dinos$EMMLi$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)
 
 ## Neornithes ##
shape.neornith.emmli <- shape.neornith
rownames(shape.neornith.emmli) <- NULL
colnames(shape.neornith.emmli) <- NULL
EMMLi.neornith <- phyloEmmli(shape.neornith.emmli, phylo=drop.tip(phy, labels[-neornith]), method="pgls", EMMLi=TRUE, N_sample=length(neornith), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.neornith$EMMLi$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

 ## Alligator Development ##
vcm.gator <- dotcorr(shape.data[,,gator.dev])
EMMLi.gator <- EMMLi(corr=vcm.gator, N_sample=length(gator.dev), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.gator$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

 ## Chicken Development ##
vcm.chick <- dotcorr(shape.data[,,chick.dev])
EMMLi.chick <- EMMLi(corr=vcm.chick, N_sample=length(chick.dev), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.chick$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

 ## Local Shape Tests ##

## Archosaurs ##
shape.archo.emmli <- local.shape.data[,,no.dev]
rownames(shape.archo.emmli) <- NULL
colnames(shape.archo.emmli) <- NULL
emmli.archo <- phyloEmmli(shape.archo.emmli, phylo=phy, method="pgls", EMMLi=TRUE, N_sample=length(no.dev), mod=emmli.parts)
plotNetwork2(rhos=emmli.archo$EMMLi$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

## Non-Avian Dinosaurs ##
shape.dino.emmli <- local.shape.data[,,dinos]
rownames(shape.dino.emmli) <- NULL
colnames(shape.dino.emmli) <- NULL
EMMLi.dinos <- phyloEmmli(shape.dino.emmli, phylo=drop.tip(phy, labels[-dinos]), method="pgls", EMMLi=TRUE, N_sample=length(dinos), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.dinos$EMMLi$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

## Neornithes ##
shape.neornith.emmli <- local.shape.data[,,neornith]
rownames(shape.neornith.emmli) <- NULL
colnames(shape.neornith.emmli) <- NULL
EMMLi.neornith <- phyloEmmli(shape.neornith.emmli, phylo=drop.tip(phy, labels[-neornith]), method="pgls", EMMLi=TRUE, N_sample=length(neornith), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.neornith$EMMLi$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

## Alligator Development ##
vcm.gator <- dotcorr(local.shape.data[,,gator.dev])
EMMLi.gator <- EMMLi(corr=vcm.gator, N_sample=length(gator.dev), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.gator$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

## Chicken Development ##
vcm.chick <- dotcorr(local.shape.data[,,chick.dev])
EMMLi.chick <- EMMLi(corr=vcm.chick, N_sample=length(chick.dev), mod=emmli.parts)
plotNetwork2(rhos=EMMLi.chick$rho$`modules.sep.Mod + sep.between`, module_names=brain.lab)

### INTEGRATION TESTS ###

 ## Dinos vs. Birds ##
  # Overall #
  pls.dinos <- phylo.integration(shape.dino, phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions, iter=9999)
  pls.birds <- phylo.integration(shape.neornith, phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions, iter=9999)
  pls.dinos <- integration.test(dino.pgls, partition.gp=partitions, iter=9999)
  pls.birds <- integration.test(neornith.pgls, partition.gp=partitions, iter=9999)
  pls.dinos <- phylo.integration(local.shape.data[,,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions, iter=9999)
  pls.birds <- phylo.integration(local.shape.data[,,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions, iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  # Cerebrum - Optic Lobe #
  pls.dinos <- phylo.integration(shape.dino[c(cerebrum, opticlobe),,], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.birds <- phylo.integration(shape.neornith[c(cerebrum, opticlobe),,], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.dinos <- integration.test(dino.pgls[c(cerebrum, opticlobe),,], partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.birds <- integration.test(neornith.pgls[c(cerebrum, opticlobe),,], partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
    # local gpa #
  pls.dinos <- phylo.integration(local.shape.data[c(cerebrum, opticlobe),,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.birds <- phylo.integration(local.shape.data[c(cerebrum, opticlobe),,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  # Cerebrum - Cerebellum #
  pls.dinos <- phylo.integration(shape.dino[c(cerebrum, cerebellum),,], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.birds <- phylo.integration(shape.neornith[c(cerebrum, cerebellum),,], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.dinos <- integration.test(shape.dino[c(cerebrum, cerebellum),,], partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.birds <- integration.test(shape.neornith[c(cerebrum, cerebellum),,], partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
   # local gpa #
  pls.dinos <- phylo.integration(local.shape.data[c(cerebrum, cerebellum),,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.birds <- phylo.integration(local.shape.data[c(cerebrum, cerebellum),,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  # Cerebrum - Medulla #
  pls.dinos <- phylo.integration(shape.dino[c(cerebrum, medulla),,], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  pls.birds <- phylo.integration(shape.neornith[c(cerebrum, medulla),,], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  pls.dinos <- integration.test(shape.dino[c(cerebrum, medulla),,], partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  pls.birds <- integration.test(shape.neornith[c(cerebrum, medulla),,], partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
   # local gpa #
   pls.dinos <- phylo.integration(local.shape.data[c(cerebrum, medulla),,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
   pls.birds <- phylo.integration(local.shape.data[c(cerebrum, medulla),,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  # Optic Lobe - Cerebellum #
  pls.dinos <- phylo.integration(shape.dino[c(opticlobe, cerebellum),,], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.birds <- phylo.integration(shape.neornith[c(opticlobe, cerebellum),,], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.dinos <- integration.test(shape.dino[c(opticlobe, cerebellum),,], partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.birds <- integration.test(shape.neornith[c(opticlobe, cerebellum),,], partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
   # local gpa #
  pls.dinos <- phylo.integration(local.shape.data[c(opticlobe, cerebellum),,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.birds <- phylo.integration(local.shape.data[c(opticlobe, cerebellum),,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  # Optic Lobe - Medulla #
  pls.dinos <- phylo.integration(shape.dino[c(opticlobe, medulla),,], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  pls.birds <- phylo.integration(shape.neornith[c(opticlobe, medulla),,], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  pls.dinos <- integration.test(shape.dino[c(opticlobe, medulla),,], partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  pls.birds <- integration.test(shape.neornith[c(opticlobe, medulla),,], partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
   # local gpa #
   pls.dinos <- phylo.integration(local.shape.data[c(opticlobe, medulla),,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
   pls.birds <- phylo.integration(local.shape.data[c(opticlobe, medulla),,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  # Cerebellum - Medulla #
  pls.dinos <- phylo.integration(shape.dino[c(cerebellum, medulla),,], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.birds <- phylo.integration(shape.neornith[c(cerebellum, medulla),,], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.dinos <- integration.test(shape.dino[c(cerebellum, medulla),,], partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.birds <- integration.test(shape.neornith[c(cerebellum, medulla),,], partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  # local gpa #
  pls.dinos <- phylo.integration(local.shape.data[c(cerebellum, medulla),,dinos], phy=drop.tip(phy, labels[-dinos]), partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.birds <- phylo.integration(local.shape.data[c(cerebellum, medulla),,neornith], phy=drop.tip(phy, labels[-neornith]), partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  compare.pls(A=pls.dinos, A2=pls.birds, two.tailed=FALSE)
  
 ## Alligators vs. Chickens ##
  # Overall #
  pls.gator <- integration.test(shape.gator, partition.gp=partitions, iter=9999)
  pls.chick <- integration.test(shape.chick, partition.gp=partitions, iter=9999)
  local.shape.data <- abind(cerebrum.shape, opticlobe.shape, cerebellum.shape, medulla.shape, along=1)
  pls.gator <- integration.test(local.shape.data[,,gator.dev], partition.gp=partitions, iter=9999)
  pls.chick <- integration.test(local.shape.data[,,chick.dev], partition.gp=partitions, iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
  # Cerebrum - Optic Lobe #
  pls.gator <- integration.test(shape.gator[c(cerebrum, opticlobe),,], partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.chick <- integration.test(shape.chick[c(cerebrum, opticlobe),,], partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.gator <- integration.test(abind(cerebrum.shape[,,gator.dev], opticlobe.shape[,,gator.dev], along=1), partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  pls.chick <- integration.test(abind(cerebrum.shape[,,chick.dev], opticlobe.shape[,,chick.dev], along=1), partition.gp=partitions[c(cerebrum, opticlobe)], iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
  # Cerebrum - Cerebellum #
  pls.gator <- integration.test(shape.gator[c(cerebrum, cerebellum),,], partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.chick <- integration.test(shape.chick[c(cerebrum, cerebellum),,], partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.gator <- integration.test(abind(cerebrum.shape[,,gator.dev], cerebellum.shape[,,gator.dev], along=1), partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  pls.chick <- integration.test(abind(cerebrum.shape[,,chick.dev], cerebellum.shape[,,chick.dev], along=1), partition.gp=partitions[c(cerebrum, cerebellum)], iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
  
  # Cerebrum - Medulla #
  pls.gator <- integration.test(shape.gator[c(cerebrum, medulla),,], partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  pls.chick <- integration.test(shape.chick[c(cerebrum, medulla),,], partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  pls.gator <- integration.test(abind(cerebrum.shape[,,gator.dev], medulla.shape[,,gator.dev], along=1), partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  pls.chick <- integration.test(abind(cerebrum.shape[,,chick.dev], medulla.shape[,,chick.dev], along=1), partition.gp=partitions[c(cerebrum, medulla)], iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
  
  # Optic Lobe - Cerebellum #
  pls.gator <- integration.test(shape.gator[c(opticlobe, cerebellum),,], partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.chick <- integration.test(shape.chick[c(opticlobe, cerebellum),,], partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.gator <- integration.test(abind(opticlobe.shape[,,gator.dev], cerebellum.shape[,,gator.dev], along=1), partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  pls.chick <- integration.test(abind(opticlobe.shape[,,chick.dev], cerebellum.shape[,,chick.dev], along=1), partition.gp=partitions[c(opticlobe, cerebellum)], iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
  # Optic Lobe - Medulla #
  pls.gator <- integration.test(shape.gator[c(opticlobe, medulla),,], partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  pls.chick <- integration.test(shape.chick[c(opticlobe, medulla),,], partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  pls.gator <- integration.test(abind(opticlobe.shape[,,gator.dev], medulla.shape[,,gator.dev], along=1), partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  pls.chick <- integration.test(abind(opticlobe.shape[,,chick.dev], medulla.shape[,,chick.dev], along=1), partition.gp=partitions[c(opticlobe, medulla)], iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
  # Cerebellum - Medulla #
  pls.gator <- integration.test(shape.gator[c(cerebellum, medulla),,], partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.chick <- integration.test(shape.chick[c(cerebellum, medulla),,], partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.gator <- integration.test(abind(cerebellum.shape[,,gator.dev], medulla.shape[,,gator.dev], along=1), partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  pls.chick <- integration.test(abind(cerebellum.shape[,,chick.dev], medulla.shape[,,chick.dev], along=1), partition.gp=partitions[c(cerebellum, medulla)], iter=9999)
  compare.pls(A=pls.gator, A2=pls.chick, two.tailed=FALSE)
