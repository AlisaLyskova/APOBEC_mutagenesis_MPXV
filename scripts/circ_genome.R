library(circlize)
library(data.table)

gff <- read.csv("GCF_014621545.1_ASM1462154v1_genomic.gff", sep='\t', header=F, comment.char = "#")
gff <- data.table(gff)

genes <- gff[V3=="gene"]
genes$V3 <- as.factor(genes$V3)
genes[,strand := ifelse(V7=="+",1,0)]
genes[,strandColor := ifelse(strand==1,"red","green")]
genes[,tmp := tstrsplit(V9,"Name=")[2]]
genes[,geneName := tstrsplit(tmp,';')[1]]

genesUCSC <- read.csv("ucsc_early_late.txt",sep='\t')
genesUCSC <- data.table(genesUCSC)
genesUCSC[, chromStartFix := chromStart + 1]
genesUCSCmerge <- merge(genesUCSC,genes,by.x=c("name","chromStartFix"),by.y=c("geneName","V4"),all=T)
genesUCSCmerge[, chromEndFix := ifelse(is.na(V5),chromEnd,V5)]
genesUCSCmerge[,strandFixStr := ifelse(is.na(V7),strand.x,V7)]
genesUCSCmerge[,strandFix := ifelse(strandFixStr=="+",1,0)]
genesUCSCmerge[,stageFix := ifelse(is.na(stage),"unknown",stage)]
gfinal <- genesUCSCmerge[order(chromStartFix),.(name,chromStartFix,chromEndFix,strandFix,stageFix)]
gfinal[stageFix == "early",stageColor := "#8AB6F9"]
gfinal[stageFix == "intermediate",stageColor := "#00246B"]
gfinal[stageFix == "late",stageColor := "#993955"]
gfinal[stageFix == "unknown",stageColor := "grey"]

gbed <- gfinal[,.("genome",chromStartFix,chromEndFix,name)]

clair3rna <- read.csv("ERR10513574_clair3.txt",sep='\t')
clair3rna <- data.table(clair3rna)
clair3rna[,clair3rna := 1]

pos <- merge(clair3rna,clair3rna,by=c("V2","V4","V5","motif3","isAPOBEC","isAPOBECextended","vaf","altcnt"))
pos[isAPOBECextended==1, isAPOBECextended_color := "#993955"]
pos[isAPOBECextended==2, isAPOBECextended_color := "#00246B"]
pos[isAPOBECextended==0, isAPOBECextended_color := "#8AB6F9"]

pos[,antisense := 0]
for(i in 1:nrow(pos)){
  if(pos[i,isAPOBECextended] == 0)
    next
  gt <- genes[pos[i,V2] >= V4 & pos[i,V2] <= V5]
  if(nrow(gt) == 0 | length(unique(gt[,V7])) > 1){
   next
  }
  gstrand <- unique(gt[,V7])
  if((gstrand == "+" & pos[i,V4] == "G") | (gstrand == "-" & pos[i,V4] == "C")){
    pos[i,antisense := 1]
  }
}

pos[, antisense_pch := 16]
pos[antisense==1, antisense_pch := 18]

genome_length <- 197201

# Initialize circular plot
tiff("ERR10513574_clair3.tiff", width = 2000, height = 2000, res = 100)
circos.clear()
circos.par(start.degree = 90)
circos.initialize(sectors="genome",xlim = c(0, genome_length))

circos.track(ylim=c(0,1),bg.border=NA,track.height=0.1)
circos.genomicAxis(
  major.by = 10000,  # Tick interval (1 Mb in this example)
  labels.cex = 1.5     # Adjust the size of coordinate labels
)
circos.rect(xleft=gfinal$chromStartFix,xright=gfinal$chromEndFix,ybottom=gfinal$strandFix*0.5,ytop=0.5+gfinal$strandFix*0.5,sector.index = "genome", border=NA,col=gfinal$stageColor)
circos.genomicLabels(gbed,labels.column=4,cex=1)

circos.track(ylim=c(0,10),bg.border=NA)
circos.trackPoints(sectors=rep("genome",nrow(pos)),x=pos$V2,y=pos$vaf*10,col=pos$isAPOBECextended_color,pch=pos$antisense_pch,cex=2.5)
circos.trackLines(sectors=rep("genome",nrow(pos)),x=pos$V2,y=pos$vaf*10,type='h',col=pos$isAPOBECextended_color,lwd=2)
circos.yaxis()

circos.text(
  x = -0.5,  # X-coordinate outside the circular track
  y = 3,   # Centered vertically
  labels = "VAF", 
  facing = "reverse.clockwise", 
  adj = c(1, 0),  # Adjust alignment
  cex = 1.5  # Adjust text size
)
dev.off()
