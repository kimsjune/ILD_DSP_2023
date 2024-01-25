# topTable() in this file is redundant with those in dge.R 

## fibroblast vs fibrosis

topT.fibroblast_v_fibrosis <-topTable(efit, coef=1, n=Inf)


#topTags(qlf)
#anov <- glmQLFTest(fit, contrast=con)
#topTags(anov)
#summary(decideTests(anov))
#plotMD(anov)


vol.fibroblast_v_fibrosis <- data.frame(Target.name = rownames(topT.fibroblast_v_fibrosis), cbind(logFC = topT.fibroblast_v_fibrosis$logFC, PValue = topT.fibroblast_v_fibrosis$adj.P.Val))


vol.fibroblast_v_fibrosis$de <- "NO"
vol.fibroblast_v_fibrosis$de[vol.fibroblast_v_fibrosis$logFC >=1 & vol.fibroblast_v_fibrosis$PValue < 0.05] <- "UP"
vol.fibroblast_v_fibrosis$de[vol.fibroblast_v_fibrosis$logFC <= -1 & vol.fibroblast_v_fibrosis$PValue < 0.05] <- "DN"
vol.fibroblast_v_fibrosis$deLab <- NA
vol.fibroblast_v_fibrosis$deLab[vol.fibroblast_v_fibrosis$de!="NO"] <- vol.fibroblast_v_fibrosis$Target.name[vol.fibroblast_v_fibrosis$de != "NO"]

vc.vol.fibroblast_v_fibrosis <- ggplot(data=vol.fibroblast_v_fibrosis,
                                       aes(x=logFC,
                                           y=-log10(PValue),
                                           col=de,
                                           label=deLab))+
  geom_point()+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    panel.border = element_rect(colour="black"),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    plot.margin=unit(c(1,1,1,1),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_vline(xintercept=c(-1,1), col="black",
             linetype=3)+
  geom_hline(yintercept=-log10(0.05), col="black",
             linetype=3)+
  geom_text_repel(
    segment.colour="black",
    segment.linetype=3,
    data=vol.fibroblast_v_fibrosis %>% 
      filter(logFC>=1|logFC<=-1),
    aes(label=deLab),
    position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
    min.segment.length=0,
    max.overlaps=4)+
  scale_color_manual(values=c("darkblue", "grey75", "red"),
                     breaks=c("DN","NA","UP"))+
  xlab(expression(log[2]~fold~change))+
  ylab(expression(-log[10]~italic(P)~value))+
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))



vc.vol.fibroblast_v_fibrosis.fixAxis <- vc.vol.fibroblast_v_fibrosis +
  xlim(-8,8)+
  ylim(0,13)

ggsave(plot=vc.vol.fibroblast_v_fibrosis.fixAxis,
       filename="vc.vol.fibroblast_v_fibrosis.fixAxis.png", bg="transparent",
       width=4, height=4, units=c("in"), dpi=2000)

#####################

## fibrosis vs neutral

topT.fibrosis_v_neutral <-topTable(efit, coef=2, n=Inf)


#topTags(qlf)
#anov <- glmQLFTest(fit, contrast=con)
#topTags(anov)
#summary(decideTests(anov))
#plotMD(anov)


vol.fibrosis_v_neutral <- data.frame(Target.name = rownames(topT.fibrosis_v_neutral), cbind(logFC = topT.fibrosis_v_neutral$logFC, PValue = topT.fibrosis_v_neutral$adj.P.Val))


vol.fibrosis_v_neutral$de <- "NO"
vol.fibrosis_v_neutral$de[vol.fibrosis_v_neutral$logFC >=1 & vol.fibrosis_v_neutral$PValue < 0.05] <- "UP"
vol.fibrosis_v_neutral$de[vol.fibrosis_v_neutral$logFC <= -1 & vol.fibrosis_v_neutral$PValue < 0.05] <- "DN"
vol.fibrosis_v_neutral$deLab <- NA
vol.fibrosis_v_neutral$deLab[vol.fibrosis_v_neutral$de!="NO"] <- vol.fibrosis_v_neutral$Target.name[vol.fibrosis_v_neutral$de != "NO"]

vc.vol.fibrosis_v_neutral <- ggplot(data=vol.fibrosis_v_neutral,
                                       aes(x=logFC,
                                           y=-log10(PValue),
                                           col=de,
                                           label=deLab))+
  geom_point()+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    panel.border = element_rect(colour="black"),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    plot.margin=unit(c(1,1,1,1),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_vline(xintercept=c(-1,1), col="black",
             linetype=3)+
  geom_hline(yintercept=-log10(0.05), col="black",
             linetype=3)+
  geom_text_repel(
    segment.colour="black",
    segment.linetype=3,
    data=vol.fibrosis_v_neutral %>% 
      filter(logFC>=2|logFC<=-1),
    aes(label=deLab),
    position=position_nudge_center(direction="radial",x=1.5,y=1.5, center_x=0, center_y=5),
    min.segment.length=0,
    max.overlaps=15)+
  scale_color_manual(values=c("darkblue", "grey75", "red"),
                     breaks=c("DN","NA","UP"))+
  xlab(expression(log[2]~fold~change))+
  ylab(expression(-log[10]~italic(P)~value))+
  scale_y_continuous(
    breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))),
    labels = label_number(accuracy = 0.1))+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))



vc.vol.fibrosis_v_neutral.fixAxis <- vc.vol.fibrosis_v_neutral +
  xlim(-8,8)+
  ylim(0,13)#+
  # scale_y_continuous(
  # #  breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))),
  #   labels = label_number(accuracy = 0.1))

ggsave(plot=vc.vol.fibrosis_v_neutral.fixAxis,
       filename="vc.vol.fibrosis_v_neutral.fixAxis.png", bg="transparent",
       width=4, height=4, units=c("in"), dpi=2000)

#####################

## fibroblast vs neutral

topT.fibroblast_v_neutral <-topTable(efit, coef=3, n=Inf)


#topTags(qlf)
#anov <- glmQLFTest(fit, contrast=con)
#topTags(anov)
#summary(decideTests(anov))
#plotMD(anov)


vol.fibroblast_v_neutral <- data.frame(Target.name = rownames(topT.fibroblast_v_neutral), cbind(logFC = topT.fibroblast_v_neutral$logFC, PValue = topT.fibroblast_v_neutral$adj.P.Val))


vol.fibroblast_v_neutral$de <- "NO"
vol.fibroblast_v_neutral$de[vol.fibroblast_v_neutral$logFC >=1 & vol.fibroblast_v_neutral$PValue < 0.05] <- "UP"
vol.fibroblast_v_neutral$de[vol.fibroblast_v_neutral$logFC <= -1 & vol.fibroblast_v_neutral$PValue < 0.05] <- "DN"
vol.fibroblast_v_neutral$deLab <- NA
vol.fibroblast_v_neutral$deLab[vol.fibroblast_v_neutral$de!="NO"] <- vol.fibroblast_v_neutral$Target.name[vol.fibroblast_v_neutral$de != "NO"]

vc.vol.fibroblast_v_neutral <- ggplot(data=vol.fibroblast_v_neutral,
                                    aes(x=logFC,
                                        y=-log10(PValue),
                                        col=de,
                                        label=deLab))+
  geom_point()+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    panel.border = element_rect(colour="black"),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    plot.margin=unit(c(1,1,1,1),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_vline(xintercept=c(-1,1), col="black",
             linetype=3)+
  geom_hline(yintercept=-log10(0.05), col="black",
             linetype=3)+
  geom_text_repel(
    segment.colour="black",
    segment.linetype=3,
    data=vol.fibroblast_v_neutral %>% 
      filter(logFC>=1.5|logFC<=-1.5),
    aes(label=deLab),
    position=position_nudge_center(direction="radial",x=1.5,y=1.5, center_x=0, center_y=2),
    min.segment.length=0,
    max.overlaps=4)+
  scale_color_manual(values=c("darkblue", "grey75", "red"),
                     breaks=c("DN","NA","UP"))+
  xlab(expression(log[2]~fold~change))+
  ylab(expression(-log[10]~italic(P)~value))+
  scale_y_continuous(breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))



vc.vol.fibroblast_v_neutral.fixAxis <- vc.vol.fibroblast_v_neutral+
  xlim(-8,8)+
  ylim(0,13)

ggsave(plot=vc.vol.fibroblast_v_neutral.fixAxis,
       filename="vc.vol.fibroblast_v_neutral.fixAxis.png", bg="transparent",
       width=4, height=4, units=c("in"), dpi=2000)

