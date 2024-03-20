#scRNA seq analysis of CNC-derived MSCs figure generation E12.5-E13.5 by ERIKA H. 
library(Seurat)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(viridis)
library(colourpicker)
library(dplyr)
library(VennDiagram)
library(ggVennDiagram)
library(EnhancedVolcano)
library(stringr) 
library(clusterProfiler)
library(org.Mm.eg.db)

#Figure1
seurat_all_UMAP<- readRDS(seurat_all_UMAP)

seurat_all_UMAP<- SetIdent(seurat_all_UMAP,value= "cell_type_1")

DefaultAssay(seurat_all_UMAP)<- "SCT"

#create subsets
seurat_mesenchymal<- SetIdent(seurat_mesenchymal,value= "Samplev2")
seurat_mesenchymal_ctrls<-subset(seurat_mesenchymal,idents = c("Control E12.5","Control E13.5"))
seurat_mesenchymal_ckos<-subset(seurat_mesenchymal,idents = c("Meis2 cKO E12.5","Meis2 cKO E13.5"))

clusters_colored_var2 <- c("mesenchymal" = "#d4838c", 
                           "neuroepithelial" = "#99c7e3",
                           "neuronal" = "#3d95c9",
                           "neuroepithelial +\nneuronal" = "#3d95c9",
                           "glial" = "#46ACC8",
                           "epithelial" = "#02818A",
                           "endothelial" = "#b5351b",
                           "erythroblasts" = "#7B1113")

###UMAPPLOT WITH MAJOR CELL TYPES BIG SEURAT
umapall<-UMAPPlot(seurat_all_UMAP, group.by = "cell_type_1", 
                  pt.size = 0.5,label=F)+ scale_color_manual(values = clusters_colored_var2)+
  theme_void() +  # Remove background and axis labels
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank())+NoLegend() +ggtitle(NULL)

#Donut chart from big seurat of cell counts for each cluster

# Step 1: Calculate the proportion of cells in each cluster
cluster_proportions <- prop.table(table(seurat_all_UMAP$cell_type_1))

# Step 2: Create a data frame
pie_data <- data.frame(Cluster = names(cluster_proportions), Proportion = as.numeric(cluster_proportions))
pie_data$Cluster<- factor(pie_data$Cluster, levels = c("mesenchymal", 
                                                       "neuroepithelial",
                                                       "neuronal",
                                                       "neuroepithelial +\nneuronal",
                                                       "glial",
                                                       "epithelial",
                                                       "endothelial",
                                                       "erythroblasts"))
# Step 3: Create the donut chart using ggplot
hsize <- 1.5
df <- df %>% mutate(x = hsize)

donut_chart<-ggplot(pie_data, aes(x = hsize, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +  xlim(c(0.1, hsize + 0.5))+ # This makes it a polar coordinate pie chart
  theme_void() + # Removes axis and labels
  scale_fill_manual(values=clusters_colored_var2)+NoLegend()

ggsave(filename= "donutchartall.tiff",path=ws, plot=donut_chart, device = 'tiff', dpi = 600, units = 'cm', width = 12, height = 12)

seurat_mesenchymal_ctrls<-readRDS(seurat_mesenchymal_ctrls)

DefaultAssay(seurat_mesenchymal_ctrls)<- "integrated"

seurat_mesenchymal_ctrls<- SetIdent(seurat_mesenchymal_ctrls,value= "seurat_clusters")

cluster_colors_var1 <- c("1" = "#3CB371",
                         "2" = "#FFD6DD", 
                         "3" = "#5E7AE1",
                         "4" = "#94B4E1",
                         "5" = "#F090A1",
                         "6" = "#2D54B6",
                         "7" = "#D76048",
                         "8" = "#CA77AA",
                         "9" = "#A9A9A9",
                         "10" = "#66CDAA",
                         "11" = "#87CEFF",
                         "12" = "#F8E687",
                         "13" = "#C6E2FF",
                         "14" = "#696969",
                         "15" = "#36648B",
                         "16" = "#F08080",
                         "17" = "#77B6E1",
                         "18" = "#FFA991",
                         "19" = "#A7D98F", 
                         "20" = "#B1439C",
                         "21" = "#EEAEEE")

#umap classic by cluster
seurat_mesenchymal_1<- readRDS(seurat_mesenchymal_ctrls)

mesenchymalumap<-UMAPPlot(seurat_mesenchymal_1, group.by = "seurat_clusters", 
                          pt.size = 0.1,label=F)+ scale_color_manual(values = cluster_colors_var1)+
  theme_void() +  #Remove background and axis labels
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank())+NoLegend() +ggtitle(NULL)

ggsave(filename= "UMAPplot_mesenchymal_colvar1.tiff",path=fig2, plot=mesenchymalumap, device = 'tiff', dpi = 600, units = 'cm', width = 16, height = 12)

#Manhttan plot of cell counts per cluster ctrl
cell_counts <- table(seurat_mesenchymal_ctrls$seurat_clusters)#get number of cells in each condition in each cluster

# Convert the table to a data frame
cell_counts_df <- as.data.frame(cell_counts)
manhattan_data <- data.frame(
  Cluster = as.character(names(cell_counts)),
  CellCount = as.numeric(cell_counts)
)
manhattan_data$Cluster <- factor(manhattan_data$Cluster, levels = c("1","2","3","4","5","6","7","8","9",
                                                                    "10","11","12","13","14","15","16","17",
                                                                    "18","19","20","21"))

manhattan <- ggplot(manhattan_data, aes(x = Cluster, y = CellCount, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Cell Count") +
  scale_fill_manual(values = cluster_colors_var1, labels = c("")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),  # Remove the border for the plot panel
    legend.text = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),axis.title = element_text(size = 11),
    legend.background = element_rect(fill = "white", color = NA),  # Remove the legend border
    legend.key = element_rect(fill = "white", color = NA)  # Remove the legend key border
  ) +NoLegend()
ggsave(filename= "cellcountsperclusterctrls.tiff",path=fig1, plot=manhattan, device = 'tiff', dpi = 600, units = 'cm', width = 14, height = 5)

#HEATMAP
DefaultAssay(seurat_mesenchymal_ctrls)<- "RNA"
d <- DoHeatmap(seurat_mesenchymal_ctrls, features = c("Itm2a", "Nr5a2", "Top2a", "Hist1h1a", "Hist1h2ae", "Hist1h1b", 
                                                      "Runx1", "Barx1", "Dlx6", "Tcf4", "Pdgfra", "Twist1", "Hoxa2", 
                                                      "Meis2", "Meis1", "Pax3", "Wnt5a", "Cxcl14", "Osr1", "Sfrp2", 
                                                      "Pcp4", "Meis1", "Foxp1", "Hoxa2", "Gsc", "Pax7", "Pax9", "Amot", 
                                                      "Scx", "Shox2", "Meox1", "Meox2", "Foxd1", "Tbx18", "Sox18", 
                                                      "Ptch1", "Foxc1", "Eya2", "Dkk2", "Crym", "Hey2", "Lef1", "Pitx2", 
                                                      "Gpc6", "Asb4", "Osr2", "Pax9", "Lhx6", "Lhx8", "Creb5", "Meg3", 
                                                      "Dcn", "Hic1", "Dlk1", "Tcf7l2", "Plagl1", "Foxf1", "Foxf2", 
                                                      "Hand2", "Meis2", "Hoxa2", "Dlx5", "Msx2", "Prrx2", "Alx3", "Msx1", 
                                                      "Alx1", "Foxp4", "Alcam", "Mgp", "Foxp2", "Foxp1", "Cxcl14", 
                                                      "Twist2", "Crabp1", "Irx1", "Pantr1", "Lix1", "Vim", "Alx1", "Zfhx3", 
                                                      "Stmn2", "Col9a1", "Sox9", "Col2a1", "Wwp2", "Col9a3", "Sp7", 
                                                      "Alpl", "Msx1", "Runx2"), group.bar = TRUE, group.colors = cluster_colors_var1,label = F, size = 6, draw.lines = T)+ theme(axis.text.y = element_text(size = 11, angle = 45, vjust = -0.5, hjust = 1,face ="bold")) +  # Adjust angle and vertical/horizontal justification
  scale_fill_gradient2( low = rev(c('#d1e5f0',"#a3d4f0",'#67a9cf',"#4c9fcf",'#2166ac')), high = rev(c("#ad1138",'#ba3847','#cf5f6c','#fab6be')), midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  theme(panel.background = element_blank())+NoLegend()

#Figure2
#Create module score from ctrls only for some clusters
seurat_mesenchymal_ctrls<-readRDS(seurat_mesenchymal_ctrls)

######### EXAMPLE for osteoblastic lineage##################
osteoblast_genes <-list(c("Runx2", "Msx1", "Alpl","Sp7","Dlx5"))

seurat_mesenchymal_ctrls<- AddModuleScore(object = seurat_mesenchymal_ctrls, 
                                          features = osteoblast_genes, ctrl = 100,
                                          name = "osteoblastic")

osteoblastic<-FeaturePlot(seurat_mesenchymal_ctrls,features = "osteoblastic1", pt.size = 0.1, min.cutoff = -0.2, max.cutoff = 1.2) +ggplot2::theme_void() +
  ggplot2::scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) +  # Remove background and axis labels
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank())+ggtitle(NULL)+theme(legend.text = element_text(size = 10),legend.position = c(0.9,0.8))

ggsave(filename= "osteoblastic.tiff",path=fig2, plot=osteoblastic, device = 'tiff', dpi = 600, units = 'cm', width = 16, height = 12)

#Figure3
#Feature plots creation
a<-FeaturePlot(seurat_mesenchymal_ctrls,features="gene of interest",pt.size = 0.1,min.cutoff = "q10", max.cutoff = "q90",cols = c("grey","#02239a"))+ggplot2::theme_void()+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
 axis.ticks = element_blank())+theme(legend.text = NULL,legend.position = c(0.9,0.8),
                                     plot.title = element_text(hjust = 0.5, size = 20)) 


#Figure4
seurat_mesenchymal<- readRDS(seurat_mesenchymal_final)

seurat_mesenchymal<- SetIdent(seurat_mesenchymal,value= "Samplev2")

DefaultAssay(seurat_mesenchymal)<- "integrated"
#umap all clusters all conditions
mesenchymalumapsplitbysample<-UMAPPlot(seurat_mesenchymal, group.by = "seurat_clusters",split.by="Samplev2",
                                       pt.size = 0.1,label=F,ncol=2)+ scale_color_manual(values = cluster_colors_var1)+
  theme_void() +  #Remove background and axis labels
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank(),strip.text = element_text(size = 13))+NoLegend()+ggtitle(NULL)

#Bargraph of cell counts per cluster per condition
cell_counts <- table(seurat_mesenchymal$Samplev2, seurat_mesenchymal$seurat_clusters) #get number of cells in each condition in each cluster

# Convert the table to a data frame
cell_counts_df <- as.data.frame(cell_counts)

proportions_df<- as.data.frame(cell_counts) #convert proportions to datafframe

colnames(proportions_df) <- c("Condition", "Cluster", "Proportion") #rename columns of dataframe

proportions_df$Condition <- factor(proportions_df$Condition, levels = c("Control E12.5", "Meis2 cKO E12.5", "Control E13.5", "Meis2 cKO E13.5"))
proportions_df$Condition_wrapped <- str_wrap(proportions_df$Condition, width = 10) 

cellproppercluster<-ggplot(proportions_df,aes(x = Condition, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position="fill",width = 0.8)+
  geom_text(aes(label = paste0(round(Proportion))), position = position_fill(vjust = 0.5), size = 3) + 
  labs(x = "Condition", y = "Proportion of cells", fill = "")+ scale_fill_manual(values=cluster_colors_var1,labels = c(""))+scale_x_discrete(labels=proportions_df$Condition_wrapped)+theme (
    panel.background = element_rect(fill = "white"),  # Set background color to white
    panel.grid.major = element_line(color = "gray90"),  # Customize major grid lines
    panel.grid.minor = element_line(color = NULL), legend.text = element_text(size = 11),axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size=11))+NoLegend()

# Meis2 violin plot
seurat_mesenchymal<- SetIdent(seurat_mesenchymal_ctrls,value= "seurat_clusters")
Meis2modulescore<-list("Meis2")
seurat_mesenchymal_ctrls<- AddModuleScore(object = seurat_mesenchymal_ctrls, 
                                          features =Meis2modulescore, ctrl = 100,
                                          name = "Meis2") ###create module score makes nicer floating violin

Meis2vln<- VlnPlot(seurat_mesenchymal_ctrls, features = "Meis21", cols = my_colors, pt.size = 0)+labs(title = "Meis2", x = "Cluster", y = "Expression Level")+
  NoLegend()+theme(
    axis.ticks = element_blank(),
    panel.background = element_blank())

#venn diagram DEGS 1 versus all

# Perform differential expression analysis for E12 Control vs all
degs_e12_ctrl<- FindMarkers(seurat_mesenchymal_final, ident.1 = "Control E12.5", ident.2 = NULL, assay = "SCT",only.pos = T)

# Define upregulated genes for E12
degs_e12_ctrl <- (degs_e12_ctrl[degs_e12_ctrl$avg_log2FC > 0 & degs_e12_ctrl$p_val_adj< 0.05, ])

write.xlsx(degs_e12_ctrl,"/Users/erika/Desktop/Dizertacna/scrna/integrated/tables/ControlE12vsall_upreg_genes_avglog_sct.xlsx", rowNames = TRUE)

# Perform differential expression analysis for E12 Mutant vs all 
degs_e12_mut <- FindMarkers(seurat_mesenchymal_final, ident.1 = "Meis2 cKO E12.5", ident.2 = NULL,only.pos = T)

degs_e12_mut <- (degs_e12_mut[degs_e12_mut$avg_log2FC > 0 & degs_e12_mut$p_val_adj< 0.05, ])

write.xlsx(degs_e12_mut,"/Users/erika/Desktop/Dizertacna/scrna/integrated/tables/MutantE12vsall_upreg_genes_avglog_sct.xlsx", rowNames = TRUE)

# Perform differential expression analysis for E13 Control vs all
degs_e13_ctrl <- FindMarkers(seurat_mesenchymal_final, ident.1 = "Control E13.5", ident.2 = NULL,only.pos = T)
degs_e13_ctrl <-(degs_e13_ctrl[degs_e13_ctrl$avg_log2FC> 0 & degs_e13_ctrl$p_val_adj < 0.05, ])

write.xlsx(degs_e13_ctrl,"/Users/erika/Desktop/Dizertacna/scrna/integrated/tables/ControlE13vsall_upreg_genes_avglog_sct.xlsx", rowNames = TRUE)

# Perform differential expression analysis for E13 Mutant vs all 
degs_e13_mut <- FindMarkers(seurat_mesenchymal_final, ident.1 = "Meis2 cKO E13.5", ident.2 = NULL,only.pos = T)

degs_e13_mut <- (degs_e13_mut[degs_e13_mut$avg_log2FC > 0 & degs_e13_mut$p_val_adj< 0.05, ])

write.xlsx(degs_e13_mut,"/Users/erika/Desktop/Dizertacna/scrna/integrated/tables/MutantE13vsall_upreg_genes_avglog_sct.xlsx", rowNames = TRUE)

# Create a list containing the vectors of upregulated genes
gene_lists <- list("Control E12.5"= degs_e12_ctrl,
                   "Meis2 cKO E12.5"= degs_e12_mut,
                   "Control E13.5" = degs_e13_ctrl,
                   "Meis2 cKO E13.5" = degs_e13_mut)

# Define colors for each set
venn_colors <- c("#5ca6d2", 
                 "#C65E6A",
                 "#99c7e3",
                 "#F8AFA8")

# Generate Venn diagram
vennplot <- venn.diagram(
  x=gene_lists,
  category.names = c("Control E12.5",
                     "Meis2 cKO E12.5",
                     "Control E13.5",
                     "Meis2 cKO E13.5"), filename = NULL,
  fill=venn_colors,
  output = F, 
  lwd=0,
  cex = 1.4,
  fontfamily = "Arial" )
# Display the plot
grid.draw(vennplot)

# example of volcano plot for E12
DEG_E12 <- FindMarkers(seurat_mesenchymal_final,  ident.1 = 'Control E12.5', ident.2 = 'Meis2 cKO E12.5',logfc.threshold = 0.00)
DEG_E12 <- DEG_E12[order(-DEG_E12$avg_log2FC),]

names <- rownames(DEG_E12)
rownames(DEG_E12) <- NULL
DEG_E12 <- cbind(names,DEG_E12)

filtered_DEG <- DEG_E12[!grepl("^HBB|^AA|^HIST|^MRPL|^MRPS|^RPL|^RPS|^MT-|^GM|Rik$",DEG_E12$names, ignore.case = TRUE), ]

keyvals.colour <- ifelse(
  filtered_DEG$avg_log2FC < -0.25 & filtered_DEG$p_val_adj < 0.05, '#F8AFA8',
  ifelse(filtered_DEG$avg_log2FC > 0.25 & filtered_DEG$p_val_adj <0.05, '#5ca6d2','#E7E7E7'))

names(keyvals.colour)[keyvals.colour == "#F8AFA8"] <- 'UP in Meis2 cKO'; names(keyvals.colour)[keyvals.colour == "#5ca6d2"] <- "UP in Control"; names(keyvals.colour)[keyvals.colour == '#E7E7E7'] <- 'insignificant'

volcano<- EnhancedVolcano(toptable = filtered_DEG,
                          lab = filtered_DEG$names,
                          x ='avg_log2FC',
                          y = 'p_val_adj',
                          axisLabSize = 15,
                          titleLabSize = 15,
                          title = "",
                          colAlpha = 0.7,
                          pCutoff = 0.05,
                          FCcutoff = 0.25,
                          colCustom = keyvals.colour,
                          subtitleLabSize = 15,
                          labSize = 4.8,
                          pointSize = 3.0,
                          drawConnectors = T,
                          max.overlaps = 30,
                          maxoverlapsConnectors = NULL,
                          arrowheads = F,
                          boxedLabels = F,
                          xlim = c(-1,1),
                          ylim = c(0,100),
                          selectLab = c("Sox9","Cald1", "Meis2","Runx2","Calr","Crabp1","Sp7","Xist","Usmg5","Tomm7","Prdx1","Sf3b2","Shox2","Six2", "Asb4","Alx1","Col27a1","Dlx1","Foxp1","Dnpep","Ran","C1qbp","Ost4","Rpa3"),
                          cutoffLineType = "blank") + theme_void()+ theme(legend.text = element_text(size = 15),axis.text.y = element_text(size = 15,margin = margin(r = -20, unit = "pt")), axis.text.x = element_text(size = 15),axis.title.x = element_text(size =15,face = "bold"),axis.title.y = element_text(size =15,face = "bold",angle = 90))+NoLegend()
#Figures 5 (Figure 6 was created in the same manner)
seurat_mesenchymal<-readRDS(seurat_mesenchymal_final)

condition_colors <- c("Control E12.5"= "#5ca6d2", 
                      "Meis2 cKO E12.5"= "#d49098",
                      "Control E13.5"= "#99c7e3",
                      "Meis2 cKO E13.5"= "#F8AFA8")

#Module score violin
oss_genes<-list(c("Runx2","Nrp2","Col12a1","Alx1","Col5a2","Igf2","Postn","Col5a1","Col1a1","Sp7")) ##expression of bone related collagens

seurat_mesenchymal<- AddModuleScore(object = seurat_mesenchymal, 
                                    features = oss_genes, ctrl = 100,
                                    name = "oss")

wrapped_labels <- str_wrap(seurat_mesenchymal$Samplev2, width = 10) 

violin <- VlnPlot(seurat_mesenchymal, features= "oss1", split.by = "Samplev2", cols = condition_colors, pt.size = 0) +labs(title = "Ossification markers", x = "Condition", y = "Expression Level")+
  NoLegend()+theme(
    axis.ticks = element_blank(),
    panel.background = element_blank())

new_labels <- c("Control\nE12.5", "Meis2 cKO\nE12.5", "Control\nE13.5", "Meis2 cKO\nE13.5")

# Modify the plot to use the new labels
violin <- violin + scale_x_discrete(labels = new_labels)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) 

#Module score for UMAP 
fp <- FeaturePlot(seurat_mesenchymal,features = "oss1", pt.size = 0.6, min.cutoff= "q10", max.cutoff= "q98",split.by = "Samplev2", cols = c("#e8c2be","#0b78b8")) #02239a   #09669c  #0282cf  #a32634 #0672c4 

fp <- fp & theme(
  axis.text = element_blank(), 
  axis.title = element_blank(), 
  axis.ticks = element_blank(), 
) & theme_void()&  
  theme(
    legend.position = c(0.9,0.8),  
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 9)  
  )
print(fp)

#Analysis E13 subset osteo
seurat_mesenchymal_E13 <- subset(seurat_mesenchymal, idents = c("Control E13.5", "Meis2 cKO E13.5"))
seurat_mesenchymal_E13_subsetosteo2<- subset(seurat_mesenchymal_E13, idents = c("1","5","7","16","19"))

#Analysis GO TERM
PrepSCTFindMarkers(seurat_mesenchymal_E13_subsetosteo)

DEG_E13_C1 <- FindMarkers(seurat_mesenchymal_E13_subsetosteo2,ident.1 ='1_MUT_E13', ident.2 = '1_WT_E13',logfc.threshold = 0.00)
DEG_E13_C1 <- DEG_E13_C1[order(-DEG_E13_C1$avg_log2FC),]

names <- rownames(DEG_E13_C1)
rownames(DEG_E13_C1) <- NULL
DEG_E13_C1 <- cbind(names,DEG_E13_C1)

filtered_DEG_E13_C1 <- DEG_E13_C1[!grepl("^HBB|^AA|^HIST|^MRPL|^MRPS|^RPL|^RPS|^MT-|^GM|Rik$",DEG_E13_C1$names, ignore.case = TRUE), ]

up1<- subset(filtered_DEG_E13_C1,filtered_DEG_E13_C1$avg_log2FC> 0) 

up_GO1 <- enrichGO(up1$names,OrgDb = "org.Mm.eg.db",  #run GO
                   ont = c("CC"), 
                   keyType = "SYMBOL")
#Gene enrichment bar charts
## Cluster 1 E13
go_terms<- c("mesenchymal cell differentiation","bone development","skeletal system morphogenesis","ossification","osteoblast differentiation",
             "regulation of osteoblast differentiation","bone morphogenesis","collagen fibril organization", "regulation of ossification",
             "bone mineralization", "osteoblast proliferation","phosphatase activity","collagen binding",
             "collagen-containing extracellular matrix", "collagen trimer")

p_values <- c(71.06E-39,3.39E-33,2.04522E-29,1.04416E-27,2.94603E-23, 
              2.17938E-17,3.26994E-17,9.29316E-12,8.79664E-11,
              2.65433E-09,9.30211E-08, 1.03784E-18,5.67161E-05,
              2.37E-13,8.50184E-06)# correct
custom_colors <- c(
  "#67a9cf",
  "#67a9cf",
  "#67a9cf",
  "#67a9cf", 
  "#67a9cf" ,
  "#67a9cf",
  "#67a9cf",
  "#67a9cf",
  "#67a9cf",
  "#67a9cf",
  "#67a9cf",
  "#fddbc7",
  "#fddbc7",
  "#faacb6",
  "#faacb6")

# Compute log10 of p-values
log_p_values <- -log10(p_values)

# Create a data frame
data <- data.frame(GO_Term = go_terms, Log_P_Value = log_p_values)

# Order the data in descending order of log10(p-values)
custom_order_reversed <- rev(go_terms)
data$GO_Term <- factor(data$GO_Term, levels = custom_order_reversed)

# Create a horizontal bar chart with descending order 
a <- ggplot(data, aes(x = Log_P_Value, y= (GO_Term))) + geom_bar(stat = "identity", fill=custom_colors) + 
labs(title = "Meis2 cKO E13.5 osteo", x = "-Log10(Padj)",y ="Go Analysis" )+ 
  theme_minimal() + theme(axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 18),
axis.title.y = element_text(size = 18),
axis.text.x = element_text(size = 18),panel.background=element_blank(),plot.title = element_text(size=18,hjust = -3,face = "bold"))

#Fluorescence staining analysis bar graphs
###NRP2 as an example

# Calculate means and SEMs for each condition from table
summary_df <- NRP2 %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Intensity),
    SEM = sd(Intensity)/sqrt(n()),
    SD = sd(Intensity)  # You can choose to use SEM or SD for the error bars
  )

# Plot
nrp2 <- ggplot(summary_df, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.5), width = 0.6) +
  geom_errorbar(
    aes(ymin = Mean - SEM, ymax = Mean + SEM), 
    width = 0.15, 
    position = position_dodge(0.5)
  ) +
  geom_jitter(
    data = NRP2, 
    aes(x = Condition, y = Intensity), 
    position = position_jitter(width = 0.25), size = 2, alpha = 0.6, color = "black"
  ) +
  labs(title = "",
       y = "Maximum Nrp2 ROI\nintensity", x = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 15),panel.background = element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.text= element_text(size = 15,colour = "black"),axis.title = element_text(size = 15))+NoLegend() +scale_fill_manual(values=condition_colors)+ylim(0, 1000)

# Split the data into two vectors based on the condition
intensity_ctrl <- NRP2$Intensity[NRP2$Condition == "Control E13.5"]
intensity_meis2 <- NRP2$Intensity[NRP2$Condition == "Meis2 cKO E13.5"]

# Perform the Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(intensity_ctrl, intensity_meis2, alternative = "two.sided")

# Print the test result
print(wilcox_test_result)

# data: intensity_ctrl and intensity_meis2
# p-value = 0.03147

