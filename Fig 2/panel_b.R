## Pairwise comparisons
set.seed(0)
fit_contrast <- contrasts.fit(fit, contrasts = con)
efit <- eBayes(fit_contrast, robust = TRUE)


topGenes <- topTable(efit, coef=c(1:3), n=50, p.value=0.05, adjust.method="BH", lfc=1, sort.by="F")


lcpm_subset <- cbind(
 # id=rownames(assay(spe_ruv_subset,2)),
    assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno=="neutral"],
    assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno=="fibrosis"],
    assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno=="fibroblast"])

# scale
lcpm_subset_scale <- cbind(
  id=rownames(lcpm_subset),
  data.frame(t(scale(t(lcpm_subset))))
)

# subset topGenes only

lcpm_subset_scale_topGenes <- lcpm_subset_scale[rownames(topGenes),]

# make pairwise comparisons and annotate differentially expressed genes as either up, down, or NA
topT.fibroblast_v_fibrosis <-topTable(efit, coef=1, n=Inf)
vol.fibroblast_v_fibrosis <- data.frame(Target.name = rownames(topT.fibroblast_v_fibrosis), cbind(logFC = topT.fibroblast_v_fibrosis$logFC, PValue = topT.fibroblast_v_fibrosis$adj.P.Val))
vol.fibroblast_v_fibrosis$de <- "NO"
vol.fibroblast_v_fibrosis$de[vol.fibroblast_v_fibrosis$logFC >=1 & vol.fibroblast_v_fibrosis$PValue < 0.05] <- "UP"
vol.fibroblast_v_fibrosis$de[vol.fibroblast_v_fibrosis$logFC <= -1 & vol.fibroblast_v_fibrosis$PValue < 0.05] <- "DN"
vol.fibroblast_v_fibrosis$deLab <- NA
vol.fibroblast_v_fibrosis$deLab[vol.fibroblast_v_fibrosis$de!="NO"] <- vol.fibroblast_v_fibrosis$Target.name[vol.fibroblast_v_fibrosis$de != "NO"]
vol.fibroblast_v_fibrosis$seq <- seq_along(vol.fibroblast_v_fibrosis[,1])

topT.fibrosis_v_neutral <-topTable(efit, coef=2, n=Inf)
vol.fibrosis_v_neutral <- data.frame(Target.name = rownames(topT.fibrosis_v_neutral), cbind(logFC = topT.fibrosis_v_neutral$logFC, PValue = topT.fibrosis_v_neutral$adj.P.Val))
vol.fibrosis_v_neutral$de <- "NO"
vol.fibrosis_v_neutral$de[vol.fibrosis_v_neutral$logFC >=1 & vol.fibrosis_v_neutral$PValue < 0.05] <- "UP"
vol.fibrosis_v_neutral$de[vol.fibrosis_v_neutral$logFC <= -1 & vol.fibrosis_v_neutral$PValue < 0.05] <- "DN"
vol.fibrosis_v_neutral$deLab <- NA
vol.fibrosis_v_neutral$deLab[vol.fibrosis_v_neutral$de!="NO"] <- vol.fibrosis_v_neutral$Target.name[vol.fibrosis_v_neutral$de != "NO"]

topT.fibroblast_v_neutral <-topTable(efit, coef=3, n=Inf)
vol.fibroblast_v_neutral <- data.frame(Target.name = rownames(topT.fibroblast_v_neutral), cbind(logFC = topT.fibroblast_v_neutral$logFC, PValue = topT.fibroblast_v_neutral$adj.P.Val))
vol.fibroblast_v_neutral$de <- "NO"
vol.fibroblast_v_neutral$de[vol.fibroblast_v_neutral$logFC >=1 & vol.fibroblast_v_neutral$PValue < 0.05] <- "UP"
vol.fibroblast_v_neutral$de[vol.fibroblast_v_neutral$logFC <= -1 & vol.fibroblast_v_neutral$PValue < 0.05] <- "DN"
vol.fibroblast_v_neutral$deLab <- NA
vol.fibroblast_v_neutral$deLab[vol.fibroblast_v_neutral$de!="NO"] <- vol.fibroblast_v_neutral$Target.name[vol.fibroblast_v_neutral$de != "NO"]






## Volcano plots
vc.vol.fibroblast_v_fibrosis.new <- ggplot(data=vol.fibroblast_v_fibrosis,
       aes(y=logFC,
         x=-log10(PValue),
         col=de,
         label=deLab))+
  geom_point()+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    axis.line = element_line(linewidth=0.2, colour="black"),
    #panel.border = element_rect(colour="black"),
    panel.border = element_blank(),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    plot.margin=unit(c(1,0.1,1,0),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_vline(xintercept=-log10(0.05), col="black",
              linetype=3)+
  geom_hline(yintercept=c(-1,1), col="black",
             linetype=3)+
  geom_text_repel(
    segment.colour="black",
    segment.linetype=3,
    data=vol.fibroblast_v_fibrosis %>%
      filter(logFC>=1|logFC<=-1),
    aes(label=deLab),
    position=
      position_nudge_center(x=0, y=0.15,
                            # center_x=function(x){
                            #   quantile(x,
                            #            probs=1/4,
                            #            names=FALSE)
                            # },
                            center_y=0, center_x=5),
    # position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
    min.segment.length=0,
    max.overlaps=25)+
  scale_color_manual(values=c("darkblue", "grey75", "red"),
                     breaks=c("DN","NA","UP"))+
  ylab(expression(log[2]~fold~change))+
  #xlab(expression(-log[10]~italic(P)~value))+
  xlab("")+
  # scale_y_continuous(
  #   labels = scales::label_number(accuracy = 0.1),
  #   breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))




vc.vol.fibroblast_v_fibrosis.new.fixAxis <- vc.vol.fibroblast_v_fibrosis.new +
  ylim(-4,4)+
  xlim(0,10)

ggsave(plot=vc.vol.fibroblast_v_fibrosis.new.fixAxis,
       filename="vc.vol.fibroblast_v_fibrosis.new.fixAxis.png", bg="transparent",
       width=4, height=4, units=c("in"), dpi=2000)

####################################################################

vc.vol.fibrosis_v_neutral.new <- 
  ggplot(data=vol.fibrosis_v_neutral,
         aes(y=logFC,
             x=-log10(PValue),
             col=de,
             label=deLab))+
  geom_point()+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    #axis.ticks.y = element_blank(),
    axis.line = element_line(linewidth=0.2, colour="black"),
    #panel.border = element_rect(colour="black"),
    panel.border = element_blank(),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    axis.text.y=element_blank(),
    plot.margin=unit(c(1,0.1,1,0),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_vline(xintercept=-log10(0.05), col="black",
             linetype=3)+
  geom_hline(yintercept=c(-1,1), col="black",
             linetype=3)+
  geom_text_repel(
    segment.colour="black",
    segment.linetype=3,
    data=vol.fibrosis_v_neutral %>%
      filter(logFC>=1|logFC<=-1),
    aes(label=deLab),
    position=
      position_nudge_center(x=0, y=0.15,
                            # center_x=function(x){
                            #   quantile(x,
                            #            probs=1/4,
                            #            names=FALSE)
                            # },
                            center_y=0, center_x=5),
    # position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
    min.segment.length=0,
    max.overlaps=25)+
  scale_color_manual(values=c("darkblue", "grey75", "red"),
                     breaks=c("DN","NA","UP"))+
  #ylab(expression(log[2]~fold~change))+
  #xlab(expression(-log[10]~italic(P)~value))+
  ylab("")+
  xlab("")+
  # scale_y_continuous(
  #   labels = scales::label_number(accuracy = 0.1),
  #   breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))




vc.vol.fibrosis_v_neutral.new.fixAxis <- vc.vol.fibrosis_v_neutral.new +
  ylim(-4,4)+
  xlim(0,10)

ggsave(plot=vc.vol.fibrosis_v_neutral.new.fixAxis,
       filename="vc.vol.fibrosis_v_neutral.new.fixAxis.png", bg="transparent",
       width=4, height=4, units=c("in"), dpi=2000)

####################################################################

vc.vol.fibroblast_v_neutral.new <- 
  ggplot(data=vol.fibroblast_v_neutral,
                                           aes(y=logFC,
                                               x=-log10(PValue),
                                               col=de,
                                               label=deLab))+
  geom_point()+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    #axis.ticks.y = element_blank(),
    axis.line = element_line(linewidth=0.2, colour="black"),
    #panel.border = element_rect(colour="black"),
    panel.border = element_blank(),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    axis.text.y=element_blank(),
    plot.margin=unit(c(1,0.1,1,0),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_vline(xintercept=-log10(0.05), col="black",
             linetype=3)+
  geom_hline(yintercept=c(-1,1), col="black",
             linetype=3)+
  geom_text_repel(
    segment.colour="black",
    segment.linetype=3,
    data=vol.fibroblast_v_neutral %>%
      filter(logFC>=1|logFC<=-1),
    aes(label=deLab),
    position=
      position_nudge_center(x=0, y=0.15,
                            # center_x=function(x){
                            #   quantile(x,
                            #            probs=1/4,
                            #            names=FALSE)
                            # },
                            center_y=0, center_x=5),
    # position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
    min.segment.length=0,
    max.overlaps=25)+
  scale_color_manual(values=c("darkblue", "grey75", "red"),
                     breaks=c("DN","NA","UP"))+
  #ylab(expression(log[2]~fold~change))+
  #xlab(expression(-log[10]~italic(P)~value))+
  ylab("")+
  xlab("")+
  # scale_y_continuous(
  #   labels = scales::label_number(accuracy = 0.1),
  #   breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))




vc.vol.fibroblast_v_neutral.new.fixAxis <- vc.vol.fibroblast_v_neutral.new +
  ylim(-4,4)+
  xlim(0,10)

ggsave(plot=vc.vol.fibroblast_v_neutral.new.fixAxis,
       filename="vc.vol.fibroblast_v_neutral.new.fixAxis.png", bg="transparent",
       width=4, height=4, units=c("in"), dpi=2000)

#####

patch <- vc.vol.fibroblast_v_fibrosis.new.fixAxis + vc.vol.fibrosis_v_neutral.new.fixAxis +
  vc.vol.fibroblast_v_neutral.new.fixAxis

ggsave(plot=patch,
       filename="vc.patchwork.png", bg="transparent",
       width=11, height=4, units=c("in"), dpi=3000)

