library(vegan)
library(MKinfer)
library(ggplot2)
library(ggpubr)
library(rempsyc)

cohoscales <- read.csv("CohoGrowthAverages.csv")

custom.col <- c("#66bdb3", "#83919c")

cohoscales2 <- na.omit(cohoscales)
cohoscales2$MGRatio <- cohoscales2$MG2a / cohoscales2$MG1a

BigQ <- subset(cohoscales, cohoscales$Hatchery == "Big Qualicum")
Chilliwack <- subset(cohoscales, cohoscales$Hatchery == "Chilliwack")
Quinsam <- subset(cohoscales, cohoscales$Hatchery == "Quinsam")

cohoscalesnat <- subset(cohoscales, cohoscales$Origin == "Wild")

BigQnat <- subset(cohoscalesnat, cohoscalesnat$Hatchery == "Big Qualicum")
Chilliwacknat <- subset(cohoscalesnat, cohoscalesnat$Hatchery == "Chilliwack")
Quinsamnat <- subset(cohoscalesnat, cohoscalesnat$Hatchery == "Quinsam")

set.seed(1)
adonis2(cohoscales2$CirculiTotal ~ cohoscales2$Hatchery * cohoscales2$Region, by = "terms")




#Figure 5

tiff("Figure5.tiff", units = "in", width = 8, height = 8, res = 300, compression = "lzw")
ggplot(data = cohoscales2, aes(x = Hatchery, y = CirculiTotal, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("Number of Circuli")+
  theme_classic()
dev.off()

#Figure 6

a<- ggplot(data = cohoscales, aes(x = Hatchery, y = FWG1a, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("FWG1 μm/circulus")+
  theme_classic()

b<- ggplot(data = cohoscales, aes(x = Hatchery, y = FWG2a, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("FWG2 μm/circulus")+
  theme_classic()

c <- ggplot(data = cohoscales, aes(x = Hatchery, y = MG1a, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("EMG μm/circulus")+
  theme_classic()

d <- ggplot(data = cohoscales2, aes(x = Hatchery, y = MG2a, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("LMG μm/circulus")+
  theme_classic()

e <- ggplot(data = cohoscales2, aes(x = Hatchery, y = MGAvg, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("MG μm/circulus")+
  theme_classic()

f <- ggplot(data = cohoscales2, aes(x = Hatchery, y = MGRatio, fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("LMG / EMG")+
  theme_classic()

tiff("Figure6.tiff", units = "in", width = 8, height = 8, res = 300, compression = "lzw")
ggarrange(a, b, c, d, e, f, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
dev.off()


#period specific growth stats
set.seed(1)
adonis2(cohoscales$FWG1a ~ cohoscales$Hatchery + cohoscales$Region, by = "terms")

set.seed(1)
adonis2(cohoscales$FWG2a ~ cohoscales$Hatchery + cohoscales$Region, by = "terms")

set.seed(2)
adonis2(cohoscales$MG1a ~ cohoscales$Hatchery + cohoscales$Region, by = "terms")

set.seed(1)
adonis2(cohoscales2$MG2a ~ cohoscales2$Hatchery + cohoscales2$Region, by = "terms")

set.seed(1)
adonis2(cohoscales2$MGAvg ~ cohoscales2$Hatchery + cohoscales2$Region, by = "terms")

set.seed(2)
adonis2(cohoscales2$MGRatio ~ cohoscales2$Hatchery + cohoscales2$Region, by = "terms")

#Figure 7

tiff("Figure7.tiff", units = "in", width = 8, height = 8, res = 300, compression = "lzw")
ggplot(data = cohoscales, aes(x = Hatchery, y = FL.mm., fill = Region))+
  geom_boxplot()+
  scale_fill_manual(values = custom.col)+
  ylab("Fork Length (mm)")+
  theme_classic()
dev.off()

#FL stats
set.seed(1)
adonis2(cohoscales$FL.mm. ~ cohoscales$Region * cohoscales$Origin * cohoscales$Hatchery, by = "terms")

set.seed(1)
perm.t.test(Chilliwack$FL.mm. ~ Chilliwack$Region)
perm.t.test(BigQ$FL.mm. ~ BigQ$Region)
perm.t.test(Quinsam$FL.mm. ~ Quinsam$Region)

set.seed(1)
perm.t.test(Chilliwacknat$FL.mm. ~ Chilliwacknat$Region)
perm.t.test(BigQnat$FL.mm. ~ BigQnat$Region)
perm.t.test(Quinsamnat$FL.mm. ~ Quinsamnat$Region)
perm.t.test(cohoscalesnat$FL.mm. ~ cohoscalesnat$Region)


#growth region stats

cor.test(cohoscales$FWG1a, cohoscales$FWG2a, method = "kendall")
cor.test(cohoscales$FWG1a, cohoscales$MG1a, method = "kendall")
cor.test(cohoscales$FWG1a, cohoscales$MG2a, method = "kendall")
cor.test(cohoscales$FWG1a, cohoscales$FL.mm., method = "kendall")

cor.test(cohoscales$FWG2a, cohoscales$MG1a, method = "kendall")
cor.test(cohoscales$FWG2a, cohoscales$MG2a, method = "kendall")
cor.test(cohoscales$FWG2a, cohoscales$FL.mm., method = "kendall")

cor.test(cohoscales$MG1a, cohoscales$MG2a, method = "kendall")
cor.test(cohoscales$MG1a, cohoscales$FL.mm., method = "kendall")

cor.test(cohoscales$MG2a, cohoscales$FL.mm., method = "kendall")
