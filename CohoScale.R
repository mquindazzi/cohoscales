cohoscales <- read.csv("CohoGrowthAverages.csv")

custom.col <- c("#66bdb3", "#83919c")

adonis2(cohoscales2$TotalScaleLen ~ cohoscales2$FL.mm.)
adonis2(cohoscales2$TotalScaleLen ~ cohoscales2$Hatchery * cohoscales2$Region)
adonis2(cohoscales2$CirculiTotal ~ cohoscales2$Hatchery * cohoscales2$Region)
adonis2(cohoscales$FWG1a ~ cohoscales$Hatchery + cohoscales$Region)
adonis2(cohoscales$FWG1t ~ cohoscales$Hatchery * cohoscales$Region)
adonis2(cohoscales$FWG2a ~ cohoscales$Hatchery + cohoscales$Region)
adonis2(cohoscales$FWG2t ~ cohoscales$Hatchery * cohoscales$Region)
adonis2(cohoscales$MG1a ~ cohoscales$Hatchery + cohoscales$Region)
adonis2(cohoscales$MG1t ~ cohoscales$Hatchery + cohoscales$Region)
cohoscales2 <- na.omit(cohoscales)
adonis2(cohoscales2$MG2a ~ cohoscales2$Hatchery + cohoscales2$Region)
adonis2(cohoscales2$MG2t ~ cohoscales2$Hatchery + cohoscales2$Region)
adonis2(cohoscales2$MGAvg ~ cohoscales2$Hatchery + cohoscales2$Region)
adonis2(cohoscales2$MGRatio ~ cohoscales2$Hatchery + cohoscales2$Region)

cohoscales2$MGRatio <- cohoscales2$MG2a / cohoscales2$MG1a

lm1 <- lm(cohoscales2$TotalScaleLen ~ cohoscales2$FL.mm. * cohoscales2$Hatchery)

summary(lm1)

View(nice_assumptions(lm1))

ggplot(data = cohoscales, aes(x = Hatchery, y = , fill = Region))+
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

tiff("CohoGrowth.tiff", units = "in", width = 8, height = 8, res = 300, compression = "lzw")
ggarrange(a, b, c, d, e, f, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
dev.off()
