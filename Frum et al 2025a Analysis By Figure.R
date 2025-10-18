# Load required packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr)

library(ggrepel)
library(tibble)
library(Seurat)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db) # Replace with org.Mm.eg.db for mouse
library(KEGGREST)
library(ReactomePA)
library(msigdbr)
perform_enrichment_analysis <- function(seurat_obj, output_dir, organism = "hsa", pval_cutoff = 0.05, top_n_genes = 100) {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  clusters <- unique(Idents(seurat_obj))
  
  # Load Hallmark gene sets from MSigDB
  library(msigdbr)
  hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
  
  # Prepare trimmed features
  trimmed.features <- rownames(seurat_obj)
  var_regex <- '^MT|^RP'
  trimmed.features <- grep(var_regex, trimmed.features, invert = TRUE, value = TRUE)
  
  for (cluster in clusters) {
    message("Processing cluster: ", cluster)
    
    markers <- FindMarkers(
      seurat_obj,
      ident.1 = cluster,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      features = trimmed.features
    )
    
    write.csv(markers, file = file.path(output_dir, paste0("Cluster_", cluster, "_Markers.csv")), row.names = TRUE)
    
    # Filter genes with positive fold change only (keep all p-values)
    pos_genes <- markers[markers$avg_log2FC > 0, ]
    if (nrow(pos_genes) > 0) {
      gene_list <- pos_genes$avg_log2FC
      names(gene_list) <- rownames(pos_genes)
      gene_list <- sort(gene_list, decreasing = TRUE)
      top_gene_list <- head(gene_list, n = top_n_genes)
      
      # Run enrichment without pvalueCutoff
      gsea_go_bp <- enrichGO(gene = names(top_gene_list), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", readable = TRUE, pvalueCutoff = 1)
      gsea_go_cc <- enrichGO(gene = names(top_gene_list), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", readable = TRUE, pvalueCutoff = 1)
      gsea_go_mf <- enrichGO(gene = names(top_gene_list), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", readable = TRUE, pvalueCutoff = 1)
      
      # KEGG
      kegg_genes <- bitr(names(top_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      gsea_kegg <- if (!is.null(kegg_genes) && nrow(kegg_genes) > 0) {
        enrichKEGG(gene = kegg_genes$ENTREZID, organism = organism, keyType = "kegg", pvalueCutoff = 1)
      } else NULL
      
      # Reactome
      library(ReactomePA)
      reactome_results <- enrichPathway(gene = kegg_genes$ENTREZID, organism = "human", readable = TRUE, pvalueCutoff = 1)
      
      # Hallmark
      hallmark_results <- enricher(
        gene = names(top_gene_list),
        TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
        pvalueCutoff = 1
      )
      
      # Save enrichment
      save_enrichment_results <- function(result, prefix, category, top_n = 20, plot_width = 12, plot_height = 12) {
        df <- as.data.frame(result)
        if (nrow(df) > 0 && "p.adjust" %in% colnames(df)) {
          df <- df[order(df$p.adjust), ]
        }
        write.csv(df, file = file.path(output_dir, paste0("Cluster_", cluster, "_", prefix, ".csv")), row.names = FALSE)
        
        if (!is.null(result) && nrow(df) > 0) {
          p <- barplot(result, showCategory = min(top_n, nrow(df))) +
            ggtitle(paste("Cluster", cluster, category, "(Top ", min(top_n, nrow(df)), ")", sep = ""))
          ggsave(file.path(output_dir, paste0("Cluster_", cluster, "_", prefix, "_barplot.pdf")), plot = p, width = plot_width, height = plot_height)
        }
      }
      
      save_enrichment_results(gsea_go_bp, "GO_BP", "GO Enrichment (BP)")
      save_enrichment_results(gsea_go_cc, "GO_CC", "GO Enrichment (CC)")
      save_enrichment_results(gsea_go_mf, "GO_MF", "GO Enrichment (MF)")
      save_enrichment_results(gsea_kegg, "KEGG", "KEGG Pathway Enrichment")
      save_enrichment_results(reactome_results, "Reactome", "Reactome Pathways")
      save_enrichment_results(hallmark_results, "Hallmark", "Hallmark Gene Sets")
      
      # Extra GO_BP Top 5 barplot
      if (!is.null(gsea_go_bp) && nrow(as.data.frame(gsea_go_bp)) > 0) {
        p5 <- barplot(gsea_go_bp, showCategory = 5) +
          ggtitle(paste("Cluster", cluster, "Top 5 GO Enrichment (BP)")) +
          theme(plot.title = element_text(size = 8))
        ggsave(file.path(output_dir, paste0("Cluster_", cluster, "_GO_BP_Top5_barplot.pdf")), plot = p5, width = 6, height = 3)
      }
      
    } else {
      message("No positive fold-change genes found for cluster: ", cluster)
    }
  }
}

pediatric <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/RObjects/Pediatric Cell Atlas/FinalCellsProcessedandAnnotatedv2.RDS")
pediatric <- subset(pediatric, idents = "Basophils", invert = TRUE)
Idents(pediatric) <- "annotation_lvl3"
pediatric <- RenameIdents(pediatric, "CFTR+ AT2" = "AT2", "FMO5+ AT2" = "AT2", "Transitional 1" = "Transitional", "Transitional 2" = "Transitional", "AT1 1" = "AT1", "AT1 2" = "AT1", "AT1 3" = "AT1", "Terminal Bronchiole Cell (RASC or TRB)" = "Respiratory Bronchiole", "Respiratory Bronchiole Cell" = "Respiratory Bronchiole", "Goblet-Like Secretory" = "Respiratory Bronchiole")
pediatric <- AddMetaData(pediatric, col.name = "full.dotplot.annotation", Idents(pediatric))

Idents(pediatric) <- "annotation_lvl2"
pediatric <- RenameIdents(pediatric, c("AT2" = "Alveolar", "AT1" = "Alveolar", "Transitional" = "Alveolar", "Respiratory Bronchiole" = "Airway", "Neuroendocrine" = "Airway", "Multiciliated" = "Airway", "Basal" = "Airway"))
pediatric <- AddMetaData(pediatric, col.name = "lvl2.dotplot.annotation", Idents(pediatric))



color_palette <- c(
  "#E41A1C", "#F4A3A8", # Intense Red, Soft Red
  "#377EB8", "#A6CEE3", # Intense Blue, Soft Blue
  "#4DAF4A", "#B2DF8A", # Intense Green, Soft Green
  "#984EA3", "#CAB2D6", # Intense Purple, Soft Lavender
  "#A65628", "#E7C5B4", # Intense Brown, Soft Brown
  "#F781BF", "#FCCDE5", # Intense Pink, Soft Pink
  "#66C2A5", "#C7EAE5", # Intense Teal, Soft Teal
  "#8DA0CB", "#D1D1E0", # Intense Indigo, Soft Indigo
  "#E78AC3", "#F4CAE4", # Intense Magenta, Soft Magenta
  "#A6D854", "#D9F0A3", # Intense Lime, Soft Lime
  "#1F78B4", "#B3CDE3", # Intense Deep Blue, Soft Deep Blue
  "#33A02C", "#A1DAB4", # Intense Forest Green, Soft Forest Green
  "#6A3D9A", "#D4B9DA", # Intense Violet, Soft Violet
  "#B15928", "#E7BA8B", # Intense Rust, Soft Rust
  "#1B9E77", "#B2E2E2", # Intense Dark Teal, Soft Dark Teal
  "#FF69B4", "#FFC0CB", # Intense Hot Pink, Soft Blush Pink
  "#8C510A", "#C19A6B"  # Intense Earthy Brown, Soft Sand Brown
)

pdf(file.path("./", paste0("RNA Epithelium Louvain Figure 1B", ".pdf")), w=16, h=12.35)
DimPlot(pediatric,  label = FALSE, group.by = "annotation_lvl1", pt.size = 0.2) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("RNA Epithelium Louvain Figure 1C", ".pdf")), w=16, h=12.35)
DimPlot(pediatric,  label = FALSE, group.by = "lvl2.dotplot.annotation", pt.size = 0.2) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("RNA Epithelium Louvain Figure 1C Reference", ".pdf")), w=16, h=12.35)
DimPlot(pediatric,  label = TRUE, group.by = "lvl2.dotplot.annotation", pt.size = 0.2) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("RNA Epithelium Louvain Figure 1D", ".pdf")), w=16, h=12.35)
DimPlot(pediatric,  label = FALSE, pt.size = 0.2, cols = color_palette) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("RNA Epithelium Louvain Figure 1D Reference", ".pdf")), w=16, h=12.35)
DimPlot(pediatric,  label = TRUE, pt.size = 0.2) + NoLegend() + NoAxes()
dev.off()

Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 15, scale = TRUE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Average Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}

pediatric$full.dotplot.annotation <- factor(pediatric$full.dotplot.annotation, levels = c("Alveolar Fibroblast 1", "Alveolar Fibroblast 2","CCL2-Positive Fibroblast", "Fibroblast","Alveolar Myofibroblast",  "Adventitial-Like Mesenchyme",  "Pericyte", "NK Cell", "CD8+ T-Cell", "CD4+ T-Cell", "B-Cells", "Neutrophils",  "CD1c+ DC","pDC","Basophils", "Monocytes", "Macrophages", "Activated Macrophages", "Capillary Aerocyte","Capillary Endothelial",  "Arterial Endothelial", "Venous Endothelial", "Lymphatic Endothelial", "AT1", "Transitional", "AT2",  "Respiratory Bronchiole", "Basal", "Multiciliated", "Neuroendocrine"))
all.markers <- c(mesenchymal.dotplot.markers, immune.dotplot.markers, endothelial.dotplot.markers, epi.ann2.dotplot.markers)
all.markers <- make.unique((all.markers))
pdf(file.path("./", paste0("RNA All Annotation Level 3 DotPlot", ".pdf")), w=80, h=20)
Dotplot_Zhiwei_Version(pediatric, c(all.markers), "full.dotplot.annotation") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()
pdf(file.path("./", paste0("RNA All Annotation Level 3 Trimmed", ".pdf")), w=1, h=15)
Dotplot_Zhiwei_Version(pediatric, c(all.genes), "full.dotplot.annotation") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()

all.genes <- c( "RSPO2",  "CCL2", "KIAA1211", "MYH11", "PI15", "LAMC3",  "KLRF1", "CD8A", "ICOS", "BANK1", "FCN1", "NR4A3", "CLEC4C", "F13A1", "MARCO", "APOE", "EDNRB", "FCN3", "GJA5", "ACKR1", "RELN", "AGER","CLDN4","LAMP3", "NIBAN1", "TP63", "FOXJ1", "CHGA")

Idents(pediatric) <- "lvl2.dotplot.annotation"
trimmed.features = rownames(pediatric)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
pediatric.markers <- FindAllMarkers(pediatric, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
pediatric.markers %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(pediatric))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("All Top10  Markers per seurat cluster lvl2", ".pdf")), w=50, h=6)
Dotplot_Zhiwei_Version(pediatric, genes.for.dotplot, "lvl2.dotplot.annotation")
dev.off()
pdf(file.path("./", paste0("All Top10  Markers per seurat cluster full annotation", ".pdf")), w=50, h=15)
Dotplot_Zhiwei_Version(pediatric, genes.for.dotplot, "full.dotplot.annotation")
dev.off()
pediatric$lvl2.dotplot.annotation <- factor(pediatric$lvl2.dotplot.annotation, levels = c("Signaling Mesenchyme", "Structural Mesenchyme","Lymphoid", "Myeloid", "Capillary", "Large Vessel", "Alveolar", "Airway"))
lvl2.genes <- c("VIM", "PDGFRA", "PODN", "PTPRC", "ITK", "OLR1", "CLEC7A","PECAM1", "CDH5", "FENDRR" , "MMRN1",  "CDH1", "NKX2-1", "LAMA3",  "SOX4")
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 10, scale = TRUE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Average Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}



pdf(file.path("./", paste0("All Top10  Markers per lvl2 annotation Fig", ".pdf")), w=7.5, h=4.5)
Dotplot_Zhiwei_Version(pediatric, lvl2.genes, "lvl2.dotplot.annotation")
dev.off()



#Mesenchymal Subset




PDA.mes <- subset(PDA.fastmnn, idents = c(24, 7, 19, 24, 22, 26)) 
PDA.mes <- RunUMAP(PDA.mes, dims = 1:25,reduction = "mnn", return.model = TRUE)
PDA.mes <- FindNeighbors(PDA.mes, dims = 1:25,reduction = "mnn")
PDA.mes <- FindClusters(PDA.mes, resolution = 0.4, algorithm = 1) 
DimPlot(PDA.mes, label = TRUE, pt.size = 0.2)

DefaultAssay(PDA.mes) <- "RNA"
PDA.mes <- NormalizeData(PDA.mes)
PDA.mes$annotation_lvl3 <- factor(PDA.mes$annotation_lvl3, levels = c("Alveolar Fibroblast 1", "Alveolar Fibroblast 2","CCL2-Positive Fibroblast", "Fibroblast","Alveolar Myofibroblast",  "Adventitial-Like Mesenchyme",  "Pericyte"))

mesenchymal.dotplot.markers <- c("ITGA8", "PIEZO2","PDGFRA", "NPNT",  "RSPO2", "FGF7", "WNT2", "VEGFD",  "CD44", "CCL2","ZEB1", "SVIL", "FILIP1L", "FLRT2",  "KLHL13", "KIAA1211", "COL8A1", "ASPN", "MYH11", "ACTA2", "DACH2","FGF18","WNT5A","PLCB1", "ATP10A","PDGFD", "PI15",  "SCARA5", "SFRP2","EBF1","PDGFRB", "MYO1B", "LAMC3")
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
DotPlot(seurat_object, features = gene_list,  dot.scale = 14, scale = TRUE, group.by = grouping) + RotatedAxis() +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
scale_colour_distiller(palette = "YlGnBu", direction = 1) +
guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
labs(y=NULL, x= NULL) +
guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}

Idents(PDA.mes) <- "annotation_lvl2"
trimmed.features = rownames(PDA.mes)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
PDA.mes.markers <- FindAllMarkers(PDA.mes, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
PDA.mes.markers %>% group_by(cluster)  %>% top_n(n = 40, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(PDA.mes))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("MesenchymeTop10 Markers per seurat cluster lvl2", ".pdf")), w=50, h=6)
Dotplot_Zhiwei_Version(PDA.mes, genes.for.dotplot, "annotation_lvl3")
dev.off()


pdf(file.path("./", paste0("RNA Mesenchyme Annotation Level 3 DotPlot", ".pdf")), w=16, h=5.5)
Dotplot_Zhiwei_Version(PDA.mes, mesenchymal.dotplot.markers, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()
pdf(file.path("./", paste0("RNA Mesenchyme Louvain Labelled Lvl 3", ".pdf")), w=16, h=12.35)
DimPlot(PDA.mes, group.by = c("annotation_lvl3"), label = FALSE, pt.size = 0.5) + NoLegend() + NoAxes()
dev.off()

PDA.immune <- subset(PDA.fastmnn, idents = c(12, 28, 16, 4, 30, 10)) 
DefaultAssay(PDA.immune) <- "integrated"
PDA.immune <- ScaleData(PDA.immune)

PDA.immune <- RunPCA(PDA.immune, verbose = FALSE)
ElbowPlot(PDA.immune, ndims = 50)
PDA.immune <- RunUMAP(PDA.immune, reduction = "pca", dims = 1:16, return.model = TRUE)
PDA.immune <- FindNeighbors(PDA.immune, dims = 1:16)
PDA.immune <- FindClusters(PDA.immune, resolution = 0.3)

DimPlot(PDA.immune, label = TRUE, pt.size = 0.35)


PDA.immune <- RenameIdents(PDA.immune,  "12" = "pDC", "11" = "Macrophages","10" = "Basophils", "8" = "B-Cells", "7" = "CD1c+ DC", "6" = "Activated Macrophages", "5" = "Neutrophils", "4" = "Monocytes", "1" = "Macrophages", "0" = "Macrophages")
PDA.immune <- RenameIdents(PDA.immune, "13" = "Lymphoid", "12" = "Myeloid", "11" = "Myeloid", "10" = "Myeloid", "9" = "Lymphoid", "8" = "Lymphoid", "7" =  "Myeloid", "6" = "Myeloid", "5" = "Myeloid", "4" = "Myeloid", "3" = "Lymphoid", "2" = "Lymphoid", "1" = "Myeloid", "0" = "Myeloid")                                           
PDA.immune <- AddMetaData(PDA.immune, col.name = "annotation_lvl1", "Immune")
PDA.immune <- AddMetaData(PDA.immune, col.name = "annotation_lvl2", Idents(PDA.immune.named2))
PDA.immune <- AddMetaData(PDA.immune, col.name = "annotation_lvl3", c(Idents(PDA.fastmnn.rpca.tcell.3.named), Idents(PDA.immune.named)))
PDA.immune$annotation_lvl3 <- factor(PDA.immune$annotation_lvl3, levels = c("NK Cell", "CD8+ T-Cell", "CD4+ T-Cell", "B-Cells",  "Neutrophils",  "CD1c+ DC","pDC","Basophils", "Monocytes", "Macrophages", "Activated Macrophages")) <- RenameIdents(PDA.immune,  "12" = "pDC", "11" = "Macrophages","10" = "Basophils", "8" = "B-Cells", "7" = "CD1c+ DC", "6" = "Activated Macrophages", "5" = "Neutrophils", "4" = "Monocytes", "1" = "Macrophages", "0" = "Macrophages")
PDA.immune.named2 <- RenameIdents(PDA.immune, "13" = "Lymphoid", "12" = "Myeloid", "11" = "Myeloid", "10" = "Myeloid", "9" = "Lymphoid", "8" = "Lymphoid", "7" =  "Myeloid", "6" = "Myeloid", "5" = "Myeloid", "4" = "Myeloid", "3" = "Lymphoid", "2" = "Lymphoid", "1" = "Myeloid", "0" = "Myeloid")                                           
PDA.immune <- AddMetaData(PDA.immune, col.name = "annotation_lvl1", "Immune")
PDA.immune <- AddMetaData(PDA.immune, col.name = "annotation_lvl2", Idents(PDA.immune.named2))
PDA.immune <- AddMetaData(PDA.immune, col.name = "annotation_lvl3", c(Idents(PDA.fastmnn.rpca.tcell.3.named), Idents(PDA.immune.named)))
PDA.immune$annotation_lvl3 <- factor(PDA.immune$annotation_lvl3, levels = c("NK Cell", "CD8+ T-Cell", "CD4+ T-Cell", "B-Cells",  "Neutrophils",  "CD1c+ DC","pDC","Basophils", "Monocytes", "Macrophages", "Activated Macrophages"))
immune.dotplot.markers <- c("ETS1",  "CD96","CD247",  "KLRF1", "LINGO2", "KLRK1", "CCL5","ITGA1", "CD8A", "CD3G", "LEF1", "CD28", "ICOS",  "CD74", "PLXDC2","SLC11A1", "BANK1", "MS4A1", "PAX5", "CD19", "CSF3R",  "FCN1", "CD300E", "MXD1", "NR4A3", "GPAT3","CD1C",  "RUNX2", "CLEC4C","CUX2","CD4","FMN1", "RGL1", "F13A1", "MRC1", "CD163", "MARCO",  "APOE", "CST3", "NPC2", "TREM2")
Idents(PDA.immune.subset) <- "annotation_lvl3"
PDA.immune.subset <- subset(PDA.immune, idents = "Basophils", invert = TRUE)

DefaultAssay(PDA.immune) <- "RNA"
PDA.immune <- NormalizeData(PDA.immune)
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 10, scale = TRUE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}
pdf(file.path("./", paste0("RNA Immune Annotation Level 3 DotPlot", ".pdf")), w=16, h=5.5)
Dotplot_Zhiwei_Version(PDA.immune.subset, immune.dotplot.markers, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()
pdf(file.path("./", paste0("RNA Immune Louvain Labelled Lvl 3", ".pdf")), w=16, h=12.35)
DimPlot(PDA.immune.subset, group.by = c("annotation_lvl3"), label = FALSE, pt.size = 0.5) + NoLegend() + NoAxes()
dev.off()

activated.mes.markers <- c("C1QA","C1QC",  "CTSD", "CTSH","CTSA",  "LYZ", "S100A6","S100A4",  "S100A9")
activated.subset <- subset(PDA.immune.subset, idents = c("Activated Macrophages", "Macrophages"))
Dotplot_Zhiwei_Version_NoScale <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 10, scale = FALSE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Average Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}

Idents(PDA.immune) <- "annotation_lvl3"
activated.subset <- subset(PDA.immune, idents = c("Activated Macrophages", "Macrophages"))

pdf(file.path("./", paste0("RNA Activated v General MacrophagesDotPlot", ".pdf")), w=7, h=2)
Dotplot_Zhiwei_Version_NoScale(activated.subset, activated.mes.markers, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()


PDA.endothelial <- subset(PDA.fastmnn, idents = c(20, 25, 11, 14)) 
PDA.endothelial <- RunUMAP(PDA.endothelial, dims = 1:15,reduction = "mnn", return.model = TRUE)
PDA.endothelial <- FindNeighbors(PDA.endothelial, dims = 1:15, reduction = "mnn")
PDA.endothelial <- FindClusters(PDA.endothelial, resolution = 0.15, algorithm = 1) 
DimPlot(PDA.endothelial, label = TRUE, pt.size = 0.35)

PDA.endothelial.named <- RenameIdents(PDA.endothelial,   "3" = "Arterial and Venous Endothelial", "2" = "Lymphatic Endothelial", "1" = "Capillary Aerocyte", "0" = "Capillary Endothelial")
PDA.endothelial.named2 <- RenameIdents(PDA.endothelial,   "3" = "Large Vessel", "2" = "Large Vessel", "1" = "Capillary", "0" = "Capillary")
PDA.endothelial <- AddMetaData(PDA.endothelial, col.name = "annotation_lvl2", Idents(PDA.endothelial.named2))
PDA.endothelial <- AddMetaData(PDA.endothelial, col.name = "annotation_lvl1", "Endothelium")

PDA.endothelial <- AddMetaData(PDA.endothelial, col.name = "annotation_lvl3", c(Idents(PDA.fastmnn.artvein.named), Idents(PDA.endothelial.named)))
PDA.endothelial$annotation_lvl3 <- factor(PDA.endothelial$annotation_lvl3, levels = c("Capillary Aerocyte","Capillary Endothelial",  "Arterial Endothelial", "Venous Endothelial", "Lymphatic Endothelial"))

endothelial.dotplot.markers <- c("ADGRL2", "FENDRR", "CA4", "HPGD","EDNRB", "ROBO2", "APLN", "VWF", "NOSTRIN", "BMP6","IL7R", "FCN3", "EFNB2","DKK2","PCSK5", "GJA5", "IL1R1", "ACKR1","CLU","VCAM1","CCL21", "RELN", "PROX1",  "PDPN")

Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 18, scale = TRUE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}

pdf(file.path("./", paste0("RNA Endothelium Annotation Level 3 DotPlot", ".pdf")), w=16, h=5.5)
Dotplot_Zhiwei_Version(PDA.endothelial, endothelial.dotplot.markers, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()
pdf(file.path("./", paste0("RNA Endothelium Louvain Labelled Lvl 3", ".pdf")), w=16, h=12.35)
DimPlot(PDA.endothelial, group.by = c("annotation_lvl3"), label = FALSE, pt.size = 0.5) + NoLegend() + NoAxes()
dev.off()

Idents(PDA.endothelial) <- "annotation_lvl3"
trimmed.features = rownames(PDA.endothelial)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
PDA.endothelial.markers <- FindAllMarkers(PDA.endothelial, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
PDA.endothelial.markers %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(PDA.endothelial))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("Endothelial Top10  Markers per seurat cluster lvl3", ".pdf")), w=50, h=6)
Dotplot_Zhiwei_Version(PDA.endothelial, genes.for.dotplot, "annotation_lvl3")
dev.off()

PDA.epi <- subset(PDA.fastmnn, idents = c(0, 1, 6, 27, 21, 17, 9, 8, 18, 2, 3, 5, 29, 13, 23, 15))


PDA.epi <- RunUMAP(PDA.epi, dims = 1:20, reduction = "mnn", return.model = TRUE)
PDA.epi <- FindNeighbors(PDA.epi, dims = 1:20, reduction = "mnn")
PDA.epi <- FindClusters(PDA.epi, resolution = 0.8, algorithm = 1) 
DimPlot(PDA.epi, label = TRUE, pt.size = 0.2)

PDA.epi.named <- RenameIdents(PDA.epi,   "16" = "Neuroendocrine", "15" = "Basal", "14" = "Transitional 2", "13" = "Multiciliated", "12" = "Respiratory Bronchiole", "11"  = "AT2 2", "10" = "Transitional 1", "9" = "Multiciliated", "8" = "Multiciliated", "7" = "AT2 3","6" = "AT2 1",  "5" = "Respiratory Bronchiole", "4" = "AT1 2", "3" = "AT1 3", "2" = "AT2 2", "1" = "AT1 1", "0" = "AT2 1")
PDA.epi.named2 <- RenameIdents(PDA.epi,   "16" = "Neuroendocrine", "15" = "Basal", "14" = "Transitional", "13" = "Multiciliated", "12" = "Respiratory Bronchiole", "11"  = "AT2", "10" = "Transitional", "9" = "Multiciliated", "8" = "Multiciliated", "7" = "AT2","6" = "AT2",  "5" = "Respiratory Bronchiole", "4" = "AT1", "3" = "AT1", "2" = "AT2", "1" = "AT1", "0" = "AT2")
PDA.epi <- AddMetaData(PDA.epi, col.name = "annotation_lvl2", Idents(PDA.epi.named2))
PDA.epi <- AddMetaData(PDA.epi, col.name = "annotation_lvl1", "Epithelium")

PDA.epi <- AddMetaData(PDA.epi, col.name = "annotation_lvl3", c(Idents(PDA.fastmnn.rpca.trb.named), Idents(PDA.epi.named)))
PDA.epi$annotation_lvl3 <- factor(PDA.epi$annotation_lvl3, levels = c("Neuroendocrine", "Multiciliated", "Basal", "Goblet-Like Secretory", "Respiratory Bronchiole Cell", "Terminal Bronchiole Cell (RASC or TRB)", "AT2 1", "AT2 2", "AT2 3", "Transitional 2", "Transitional 1", "AT1 3", "AT1 2", "AT1 1"))
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 18, scale = TRUE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8)) 
}
pdf(file.path("./", paste0("RNA Epithelium Annotation Level 2 DotPlot", ".pdf")), w=16, h=5.5)
Dotplot_Zhiwei_Version(PDA.epi, epi.ann2.dotplot.markers, "annotation_lvl2") + theme(axis.text.x = element_text(face = "italic", size = 20)) #dotsize = 15
dev.off()
pdf(file.path("./", paste0("RNA Epithelium Louvain Labelled Lvl 2", ".pdf")), w=16, h=12.35)
DimPlot(PDA.epi, group.by = c("annotation_lvl2"), label = FALSE, pt.size = 0.2) + NoLegend() + NoAxes()
dev.off()
PDA.epi$annotation_lvl2 <- factor(PDA.epi$annotation_lvl2, levels = c("AT1", "AT2", "Transitional",  "Respiratory Bronchiole", "Basal", "Multiciliated", "Neuroendocrine"))
epi.ann2.dotplot.markers <- c("CAV1", "RTKN2", "AGER", "SFTPB", "SFTPC", "LAMP3","ETV5", "SCGB3A2","SCGB1A1", "CP", "NIBAN1",  "BMP6","GDF15", "RNASE1", "TP63", "KRT17", "KRT5",   "CFAP299", "FOXJ1","C6","MUC16",  "SCGN",  "CHGA", "GRP","GHRL")

Idents(PDA.epi) <- "annotation_lvl2"
trimmed.features = rownames(PDA.epi)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
PDA.epi.markers <- FindAllMarkers(PDA.epi, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
PDA.epi.markers %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(PDA.epi))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("Epithelial Top10  Markers per seurat cluster lvl3", ".pdf")), w=50, h=6)
Dotplot_Zhiwei_Version(PDA.epi, genes.for.dotplot, "annotation_lvl2")
dev.off()



PDA.epi <- FindClusters(PDA.epi, resolution = 1.5, algorithm = 1) 
DimPlot(PDA.epi, label = TRUE, pt.size = 0.2)

DefaultAssay(PDA.epi) <- "RNA"
PDA.epi <- JoinLayers(PDA.epi)
PDA.epi <- NormalizeData(PDA.epi)

Idents(PDA.epi) <- "seurat_clusters"
trimmed.features = rownames(PDA.epi)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
PDA.epi.markers <- FindAllMarkers(PDA.epi, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
PDA.epi.markers %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(PDA.epi))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("Epithelial Top10  Markers per seurat cluster res08", ".pdf")), w=50, h=6)
Dotplot_Zhiwei_Version(PDA.epi, genes.for.dotplot, "seurat_clusters")
dev.off()


pdf("scatter_plots_neuroendocrine.pdf", w = 11, h = 8.5)
# Extract the unique levels from the metadata category 'source.x.disease'

  # Subset the data for the current level
  subset_data <- subset(PDA.epi, annotation_lvl2 == "Neuroendocrine")
  
  # Create a scatter plot
  p <- FeatureScatter(PDA.epi.ne, pt.size = 3, group.by = "seurat_clusters", feature1 = "GRP", feature2 = "GHRL") +
    ggtitle(paste("Scatter Plot for Neuroendocrine")) + NoLegend()
  
  # Print the plot to the PDF
  print(p)


# Close the PDF file
dev.off()


marker.1 <- "GRP"
marker.2 <- "GHRL"

# Extract expression data for markers
expr_data <- FetchData(subset_data, vars = c(marker.1, marker.2))

# Define expression threshold (adjust if needed)
threshold <- 0  # Adjust if necessary

# Count cells in exclusive categories
num_grp_only <- sum(expr_data[[marker.1]] > threshold & expr_data[[marker.2]] <= threshold)
num_ghrl_only <- sum(expr_data[[marker.2]] > threshold & expr_data[[marker.1]] <= threshold)
num_both <- sum(expr_data[[marker.1]] > threshold & expr_data[[marker.2]] > threshold)
num_neither <- sum(expr_data[[marker.1]] <= threshold & expr_data[[marker.2]] <= threshold)

# Print results
cat("Cells expressing only GRP:", num_grp_only, "\n")
cat("Cells expressing only GHRL:", num_ghrl_only, "\n")
cat("Cells expressing both GRP and GHRL:", num_both, "\n")
cat("Cells expressing neither GRP nor GHRL:", num_neither, "\n")

#subcluster ne
Idents(PDA.epi) <- "annotation_lvl3"
PDA.epi.ne <- subset(PDA.epi, idents = "Neuroendocrine")
PDA.epi.ne <- RunUMAP(PDA.epi.ne, dims = 1:30, reduction = "mnn", return.model = TRUE)
PDA.epi.ne <- FindNeighbors(PDA.epi.ne, dims = 1:30, reduction = "mnn")
PDA.epi.ne <- FindClusters(PDA.epi.ne, resolution = 0.1, algorithm = 1) 


pdf(file.path("./", paste0("Pediatric Neuroendocrine", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi.ne, label = TRUE, pt.size = 3) + NoLegend() + NoAxes()
dev.off()


pdf(file.path("./", paste0("Pediatric Neuroendocrine GHRL Feature", ".pdf")), w=11, h=8.5)
FeaturePlot(PDA.epi.ne, features = "GHRL",  pt.size = 3, order = TRUE) & NoLegend() & NoAxes() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Pediatric Neuroendocrine GRP Feature", ".pdf")), w=11, h=8.5)
FeaturePlot(PDA.epi.ne, features = "GRP",  pt.size = 3, order = TRUE) & NoLegend() & NoAxes() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf("scatter_plots_neuroendocrine.pdf", w = 11, h = 8.5)
p <- FeatureScatter(PDA.epi.ne, pt.size = 3, group.by = "seurat_clusters", 
                    feature1 = "GRP", feature2 = "GHRL") +
  ggtitle("Scatter Plot for Neuroendocrine") + 
  NoLegend() + 
  theme(
    axis.text = element_text(size = 30))
print(p)
dev.off()

FeaturePlot(PDA.epi.ne, "GRP")
FeaturePlot(PDA.epi.ne, "GHRL")
PDA.epi.ghrl <- subset(PDA.epi.ne, idents = "1")
PDA.epi.grp <- subset(PDA.epi.ne, idents = "0")


##Figure S1L Endothelial
#signatures from Schupp et al., 2021 Circulation
aerocyte.sig <- list(c("SOSTDC1",
                       "EDNRB",
                       "HPGD",
                       "CYP3A5",
                       "PRKG1",
                       "TBX2",
                       "RTN1",
                       "RCSD1",
                       "CA4",
                       "EDA",
                       "EXPH5",
                       "B3GALNT1",
                       "ADGRL2",
                       "NCALD",
                       "CLEC4E",
                       "LINC02197",
                       "EMCN",
                       "CHRM2",
                       "HOPX",
                       "AFF3",
                       "S100A4",
                       "SPON2",
                       "ST5",
                       "TSPAN15",
                       "PDE1C",
                       "FOXP2",
                       "DTNA",
                       "KHDRBS2",
                       "PLA2G4C",
                       "EDIL3",
                       "AC138776.1",
                       "SPOCK2",
                       "CHL1",
                       "CPM",
                       "HLA-E",
                       "MALAT1",
                       "PRX",
                       "SLC14A1",
                       "CD83",
                       "ITGA1",
                       "CRPPA",
                       "HECW2",
                       "SGK1",
                       "CDH13",
                       "FMNL2",
                       "MYLK",
                       "PRICKLE2",
                       "DNAJC1",
                       "ZNF331",
                       "APP"))

gcapilary.sig <- list(c("FCN3",
                        "IL7R",
                        "CD36",
                        "RGCC",
                        "BTNL9",
                        "CA4",
                        "SLC6A4",
                        "AFF3",
                        "ADGRF5",
                        "NRXN3",
                        "TMEM100",
                        "NTRK2",
                        "EPAS1",
                        "KIAA1217",
                        "BTNL8",
                        "GPIHBP1",
                        "AL356124.1",
                        "F2RL3",
                        "TIMP3",
                        "ARHGAP6",
                        "ARHGAP18",
                        "MT2A",
                        "ADRB1",
                        "ADGRL2",
                        "AKAP12",
                        "VIPR1",
                        "GALNT18",
                        "IL18R1",
                        "MT1M",
                        "KIAA1549",
                        "NCKAP5",
                        "VWF",
                        "SORBS1",
                        "MT1E",
                        "HECW2",
                        "ADGRE5",
                        "PDK4",
                        "RNASE1",
                        "NOSTRIN",
                        "ITGA1",
                        "CAVIN2",
                        "TPM1",
                        "PTPRB",
                        "JAM3",
                        "CD14",
                        "IFNGR1",
                        "CD109",
                        "KHDRBS2",
                        "SH3D19",
                        "SPARC"
                      ))
fetal.endo <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/RObjects/Fetal Proximal to Distal Atlas/EndothelialSubsetAnnotated.RDS")
PDA.endothelial <- readRDS("~/Dropbox (University of Michigan)/0_Mac_scRNA_Seq_Processing/RObjects/Pediatric Cell Atlas/PDAVesselsAnnotated.RDS")
Idents(fetal.endo) <- "group"
fetal.endo.distal <- subset(fetal.endo, idents = "distal")
Idents(fetal.endo.distal) <- "annotation_lvl3"
fetal.endo.cap <- subset(fetal.endo.distal, idents = c("Capillary Endothelial"))
DefaultAssay(fetal.endo.cap) <- "integrated"
fetal.endo.cap <- ScaleData(fetal.endo.cap)
fetal.endo.cap <- RunPCA(fetal.endo.cap, verbose = FALSE)
ElbowPlot(fetal.endo.cap, ndims = 50)
fetal.endo.cap <- RunUMAP(fetal.endo.cap, reduction = "pca", dims = 1:12, return.model = TRUE)
fetal.endo.cap <- FindNeighbors(fetal.endo.cap, dims = 1:12)
fetal.endo.cap <- FindClusters(fetal.endo.cap, resolution = 0.1)

DimPlot(fetal.endo.cap, label = TRUE, pt.size = 0.35)


DefaultAssay(fetal.endo.cap) <- "RNA"
fetal.endo.cap <- NormalizeData(fetal.endo.cap)
Idents(fetal.endo.cap) <- "seurat_clusters"
trimmed.features = rownames(fetal.endo.cap)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
fetal.endo.cap.markers <- FindAllMarkers(fetal.endo.cap, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
fetal.endo.cap.markers %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(fetal.endo.cap))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("Fetal Endothelial Top10  Markers per seurat cluster res08", ".pdf")), w=50, h=6)
Dotplot_Zhiwei_Version(fetal.endo.cap, genes.for.dotplot, "seurat_clusters")
dev.off()

DimPlot(PDA.endothelial, group.by = "annotation_lvl3", label = TRUE) + NoLegend()
Idents(PDA.endothelial.cap) <- "annotation_lvl3"
PDA.endothelial.cap <- subset(PDA.endothelial, idents = c("Capillary Endothelial", "Capillary Aerocyte"))

DefaultAssay(PDA.endothelial.cap) <- "RNA"
PDA.endothelial.cap <- NormalizeData(PDA.endothelial.cap)
PDA.endothelial.cap <- AddModuleScore(PDA.endothelial.cap, features = aerocyte.sig, name = "aero", assay = "RNA")
PDA.endothelial.cap <- AddModuleScore(PDA.endothelial.cap, features = gcapilary.sig, name = "gcap", assay = "RNA")
pdf(file.path("./", paste0("Pediatric Capillary Aerocytes Score.pdf")), w=6, h=4)
VlnPlot(PDA.endothelial.cap, features = c("aero1", "gcap1"), group.by = NULL, pt.size = 0) + geom_boxplot(width = 0.25, fill = NA, color = "black", outlier.shape = NA) + NoLegend() + theme(
axis.title = element_text(size = 16),  # Increase axis title font size
axis.text = element_text(size = 20),   # Increase axis text font size
axis.text.x = element_blank(),       # Remove x-axis text labels
axis.title.x = element_blank(),
axis.line = element_line(size = 1)     # Increase axis line thickness
)
dev.off()

DefaultAssay(fetal.endo.cap) <- "RNA"
fetal.endo.cap <- NormalizeData(fetal.endo.cap)
fetal.endo.cap <- AddModuleScore(fetal.endo.cap, features = aerocyte.sig, name = "aero", assay = "RNA")
fetal.endo.cap <- AddModuleScore(fetal.endo.cap, features = gcapilary.sig, name = "gcap", assay = "RNA")
pdf(file.path("./", paste0("Fetal General Capillary Score.pdf")), w=3, h=4)
VlnPlot(fetal.endo.cap, features = c("aero1", "gcap1"), group.by = NULL, pt.size = 0) + geom_boxplot(width = 0.25, fill = NA, color = "black", outlier.shape = NA) + NoLegend() + theme(
  axis.title = element_text(size = 16),  # Increase axis title font size
  axis.text = element_text(size = 20),   # Increase axis text font size
  axis.text.x = element_blank(),       # Remove x-axis text labels
  axis.title.x = element_blank(),
  axis.line = element_line(size = 1)     # Increase axis line thickness
)
dev.off()

fetal.endo.cap <- AddMetaData(fetal.endo.cap, col.name = "title", "all fetal capillaries")
PDA.endothelial.cap <- AddMetaData(PDA.endothelial.cap, col.name = "title", "all pediatric capillaries")
endo.cap <- merge(fetal.endo.cap, PDA.endothelial.cap)
endo.cap <- AddModuleScore(endo.cap, features = aerocyte.sig, name = "aero", assay = "RNA")
endo.cap <- AddModuleScore(endo.cap, features = gcapilary.sig, name = "gcap", assay = "RNA")
pdf(file.path("./", paste0("Fetal and Pediatric Capillary Scores.pdf")), w=6, h=4)
VlnPlot(endo.cap, features = c("aero1", "gcap1"), group.by = "title", pt.size = 0) & geom_boxplot(width = 0.25, fill = NA, color = "black", outlier.shape = NA) & NoLegend() & theme(
  axis.title = element_text(size = 16),  # Increase axis title font size
  axis.text = element_text(size = 20),   # Increase axis text font size
  axis.text.x = element_blank(),       # Remove x-axis text labels
  axis.title.x = element_blank(),
  axis.line = element_line(size = 1)     # Increase axis line thickness
)
dev.off()

#Figure S1L (fetal v pediatric myeloid numbers/proportions)
fetal <- readRDS("~/University of Michigan Dropbox/Tristan Frum/Fetal Regional Aiway_snMultiomics/Tristan/RObjects/Figure Revisions/FigureS1AnnotationLVL1.RDS")
Idents(fetal) <- "group"
fetal.distal <- subset(fetal, idents = "distal")
Idents(fetal.distal) <- "annotation_lvl2"
fetal.myeloid <- subset(fetal.distal, idents = "Myeloid")
fetal.lymphoid <- subset(fetal.distal, idents = "Lymphoid")

Idents(pediatric) <- "annotation_lvl2"
PDA.myeloid <- subset(pediatric, idents = "Myeloid")
PDA.lymphoid <- subset(pediatric, idents = "Lymphoid")

fetal.myeloid
fetal.lymphoid
PDA.myeloid
PDA.lymphoid


#Figure S1O
PDA.myeloid <- AddMetaData(PDA.myeloid, col.name = "source", "pediatric")
fetal.myeloid <- AddMetaData(fetal.myeloid, col.name = "source", "fetal")
all.myeloid <- merge(fetal.myeloid, PDA.myeloid)

pdf(file.path("./", paste0("Fetal v Pediatric Myeloid Markers Dotplot", ".pdf")), w=5, h=3)
Dotplot_Zhiwei_Version_NoScale(all.myeloid, c("MSR1", "MRC1", "MARCO"), "source") + theme(axis.text.x = element_text(face = "italic", size = 20))
dev.off()
pdf(file.path("./", paste0("Fetal v Pediatric Myeloid Markers Dotplot Full", ".pdf")), w=5, h=5.5)
Dotplot_Zhiwei_Version_NoScale(all.myeloid, c("MSR1", "MRC1", "MARCO"), "source")
dev.off()

Idents(fetal) <- "group"
fetal.distal <- subset(fetal, idents = "distal")
Idents(fetal.distal) <- "annotation_lvl3"
fetal.mc <- subset(fetal.distal, idents = c("Ciliated 1", "Ciliated 2"))
fetal.mc <- AddMetaData(fetal.mc, col.name = "source", "Fetal")
Idents(PDA.epi) <- "annotation_lvl3"
PDA.mc <- subset(PDA.epi, idents = "Multiciliated")
PDA.mc <- AddMetaData(PDA.mc, col.name = "source", "Pediatric")
PDA.mc.merge <- merge(PDA.mc, fetal.mc)

PDA.mc.merge$source <- factor(PDA.mc.merge$source, levels = c("Fetal", "Pediatric"))


pdf(file.path("./", paste0("Fetal v Pediatric Multiciliated Dotplot Full", ".pdf")), w=4.5, h=3)
Dotplot_Zhiwei_Version_NoScale(PDA.mc.merge, c("C6", "MUC16"), "source") + theme(axis.text.x = element_text(face = "italic", size = 20))
dev.off()
pdf(file.path("./", paste0("Fetal v Pediatric Multiciliated Dotplot Full Size", ".pdf")), w=4.5, h=5)
Dotplot_Zhiwei_Version_NoScale(PDA.mc.merge, c("C6", "MUC16"), "source") + theme(axis.text.x = element_text(face = "italic", size = 20))
dev.off()
#Figure 2
pediatric <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/RObjects/Pediatric Cell Atlas/FinalCellsProcessedandAnnotatedv2.RDS")
Idents(pediatric) <- "annotation_lvl3"
AT1 <- subset(pediatric, idents = c("AT1 1", "AT1 2", "AT1 3"))
DefaultAssay(AT1) <- "RNA"
AT1 <- NormalizeData(AT1)
trimmed.features = rownames(AT1)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
output_dir <- "AT1 extraction"
perform_enrichment_analysis(AT1, output_dir)




AT12 <- subset(pediatric, idents = c("AT1 1", "AT1 2", "AT1 3", "Transitional 1", "Transitional 2", "FMO5+ AT2", "CFTR+ AT2"))
AT12 <- RenameIdents(AT12, "FMO5+ AT2" = "AT2 1", "CFTR+ AT2" = "AT2 2")
AT12 <- AddMetaData(AT12, col.name = "plot_annotation", Idents(AT12))

DefaultAssay(AT12) <- "RNA"
AT12 <- NormalizeData(AT12)
trimmed.features = rownames(AT12)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
goi.at12 <- c("SFTPC", "SFTPB", "SFTPA1", "ABCA3","LAMP3","NPC2","KRT8", "CLDN4", "KRT18","ACTG1","ACTB", "AGER","VEGFA", "CAV1","RTKN2")
AT12$plot_annotation <- factor(AT12$plot_annotation, levels = c("AT2 1", "AT2 2", "Transitional 2", "Transitional 1", "AT1 3", "AT1 2", "AT1 1"))

pdf(file.path("./", paste0("AT2 Transitional AT1 dotplot", ".pdf")), w=13, h=5.5)
Dotplot_Zhiwei_Version(AT12, goi.at12, "plot_annotation") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 26))
dev.off()
goi.at12.extended <- c("SFTPC", "SFTPB", "SFTPA1", "ABCA3","LAMP3","NPC2","KRT8", "CLDN4", "KRT18","VEGFA", "CAV1","AGER","RTKN2","MYRF", "ZSCAN31", "SPOCK2", "RNASE1", "SCGB3A2", "CCBE1", "RAP1GAP2", "ACTB", "ACTG1", "CCN1",  "LAMC2", "ITGA2", "CAVIN1", "EDN1", "HBEGF", "AQP5")
pdf(file.path("./", paste0("AT2 Transitional AT1 dotplot extended", ".pdf")), w=17, h=6)
Dotplot_Zhiwei_Version(AT12, goi.at12.extended, "plot_annotation") + theme(axis.text.x = element_text(face = "italic", size = 20))
dev.off()

DefaultAssay(AT12) <- "RNA"
AT12 <- NormalizeData(AT12)
trimmed.features = rownames(AT12)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
Idents(AT12) <- "plot_annotation"
markers <- FindAllMarkers(AT12, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
markers %>% group_by(cluster)  %>% top_n(n = 20, wt = avg_log2FC) -> top10.markers
top10.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top10.markers$cluster)
clusteridents <- levels(Idents(AT12))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top10.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
pdf(file.path("./", paste0("AT2 Transitional AT1 new genes dotplot", ".pdf")), w=50, h=10)
Dotplot_Zhiwei_Version(AT12, genes.for.dotplot, "plot_annotation")
dev.off()


output_dir <- "AT2transAT1 extraction"
perform_enrichment_analysis(AT12, output_dir)



library(tidyverse)
library(readr)
library(stringr)
library(ggplot2)


PDA.epi <- readRDS("~/Dropbox (University of Michigan)/0_Mac_scRNA_Seq_Processing/RObjects/Pediatric Cell Atlas/PDAEpitheliumAnnotated.RDS")

Idents(PDA.epi) <- "seurat_clusters"
PDA.epi.alv <- subset(PDA.epi, idents = c("6", "11", "0", "2", "7", "10", "14", "3", "4", "1"))
Idents(PDA.epi.alv) <- "annotation_lvl3"
PDA.epi.alv <- RenameIdents(PDA.epi.alv, "AT2 3" = "AT2 2")
PDA.epi.alv <-AddMetaData(PDA.epi.alv, col.name = "annotation_lvl3", Idents(PDA.epi.alv))
PDA.epi.alv <- RunUMAP(PDA.epi.alv, dims = 1:30, reduction = "mnn", return.model = TRUE)
PDA.epi.alv <- FindNeighbors(PDA.epi.alv, dims = 1:30, reduction = "mnn")
PDA.epi.alv <- FindClusters(PDA.epi.alv, resolution = 0.6, algorithm = 1)
DimPlot(PDA.epi.alv, label = TRUE, pt.size = 0.2)

pdf("Alveolar Subclustering No Label.pdf", width = 11, height = 8.5)
DimPlot(PDA.epi.alv, pt.size = 0.5, group.by = "annotation_lvl3", label = FALSE) + NoLegend() + NoAxes()
dev.off()
pdf("Alveolar Subclustering With Label.pdf", width = 11, height = 8.5)
DimPlot(PDA.epi.alv, pt.size = 0.5, group.by = "annotation_lvl3", label = TRUE) + NoLegend() + NoAxes()
dev.off()

DefaultAssay(PDA.epi.alv) <- "RNA"
PDA.epi.alv <- NormalizeData(PDA.epi.alv)



pdf(file.path("./", paste0("PDA ALV EPI SFTPC Feature", ".pdf")), w=11, h=8.5)
FeaturePlot(PDA.epi.alv , features = "SFTPC", pt.size = 0.5, order = TRUE)  & NoLegend() & NoAxes() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("PDA ALV EPI FMO5 Feature", ".pdf")), w=11, h=8.5)
FeaturePlot(PDA.epi.alv , features = "FMO5", pt.size = 0.5, order = TRUE)  & NoLegend() & NoAxes() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("PDA ALV EPI CFTR Feature", ".pdf")), w=11, h=8.5)
FeaturePlot(PDA.epi.alv , features = "CFTR", pt.size = 0.5, order = TRUE)  & NoLegend() & NoAxes() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

perform_enrichment_analysis(seurat_obj = PDA.epi.alv, output_dir = "all_terms_enrichment50", top_n_genes = 50)
perform_enrichment_analysis(seurat_obj = PDA.epi.alv, output_dir = "all_terms_enrichment200", top_n_genes = 200)


#Figure S2A
#########################################
# Define the list of IDs you want to plot
library(tidyverse)
library(readr)
library(stringr)
library(ggplot2)

# Define the list of IDs you want to plot
selected_ids <- c("GO:0052866", "GO:0000922", "GO:0030695","R-HSA-9012999","GO:1900180",  "R-HSA-373755","R-HSA-983170", "GO:0005925",
                  "R-HSA-5627117","R-HSA-446728",  "R-HSA-3000157", 
                  "GO:0003779", "GO:0042060","GO:0045785", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","R-HSA-2682334", "R-HSA-1474244", 
                  "R-HSA-3000171", "GO:0050919" )

# Define the desired order of clusters (must match the names extracted from file names)
cluster_order <- c("AT2 1", "AT2 2", "Transitional 2", 
                   "Transitional 1", "AT1 3", "AT1 2", "AT1 1")

# Function to read each file, filter by selected IDs, and compute logp and fold enrichment.
# This version forces p.adjust to be numeric and handles files where the description is in "Term".
extract_pvals <- function(file_path, selected_ids) {
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Force p.adjust to numeric (this fixes issues where p.adjust is read as character)
  df <- df %>% mutate(p.adjust = as.numeric(p.adjust))
  
  # If Description column is missing, check for Term and rename it.
  if (!"Description" %in% colnames(df)) {
    if ("Term" %in% colnames(df)) {
      df <- df %>% rename(Description = Term)
    }
  }
  
  # Extract the cluster name from the file name (e.g., "Cluster_AT1 3_GO_BP.csv" -> "AT1 3")
  cluster <- str_match(basename(file_path), "Cluster_(.*?)_")[,2]
  
  # Trim whitespace from the ID column if needed
  df <- df %>% mutate(ID = trimws(ID))
  
  # Filter rows based on the selected IDs
  df_filtered <- df %>% filter(ID %in% selected_ids)
  
  df_filtered <- df_filtered %>%
    mutate(
      cluster = cluster,
      logp = -log10(p.adjust)
    )
  
  # Determine fold enrichment from available columns
  if ("FoldEnrichment" %in% colnames(df_filtered)) {
    df_filtered <- df_filtered %>% mutate(fold_enrichment = as.numeric(FoldEnrichment))
  } else if (all(c("GeneRatio", "BgRatio") %in% colnames(df_filtered))) {
    parse_ratio <- function(ratio_str) {
      nums <- as.numeric(str_split(ratio_str, "/", simplify = TRUE))
      nums[1] / nums[2]
    }
    df_filtered <- df_filtered %>%
      mutate(
        fold_enrichment = mapply(function(g, b) parse_ratio(g) / parse_ratio(b), GeneRatio, BgRatio),
        fold_enrichment = as.numeric(fold_enrichment)
      )
  } else {
    stop(paste("Missing FoldEnrichment or ratio info in", file_path))
  }
  
  df_filtered %>% select(cluster, ID, Description, p.adjust, logp, fold_enrichment)
}
#########################################
# ========================
# 1. Read enrichment results and prepare plot data
# ========================
result_dir <- "~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 2/Alveolar Subclustering/all_terms_enrichment200"  # Update this path if needed
files <- list.files(result_dir, 
                    pattern = "Cluster_.*_(GO_BP|GO_MF|Reactome|KEGG|GO_CC|Hallmark)\\.csv", 
                    full.names = TRUE)
plot_data <- map_dfr(files, extract_pvals, selected_ids = selected_ids)

# Build a lookup table to create formatted term labels from the Description column.
labels_lookup <- plot_data %>%
  group_by(ID) %>%
  summarize(Description = dplyr::first(Description), .groups = "drop") %>%
  mutate(Description = coalesce(Description, ID),
         term_label = case_when(
           str_detect(ID, "^GO")        ~ paste0("GO:\n", toupper(Description)),
           str_detect(ID, "^R-HSA-")     ~ paste0("REACTOME:\n", toupper(Description)),
           str_detect(ID, "^hsa")        ~ paste0("KEGG:\n", toupper(Description)),
           str_detect(ID, "^HALLMARK")   ~ paste0("HALLMARK:\n", toupper(Description)),
           TRUE                         ~ toupper(Description)
         ))

# Ensure the term labels appear in the order of selected_ids.
term_labels_ordered <- labels_lookup %>%
  filter(ID %in% selected_ids) %>%
  arrange(match(ID, selected_ids)) %>%
  pull(term_label)

# Merge the formatted term labels into plot_data and factor them appropriately.
plot_data <- left_join(plot_data, labels_lookup %>% select(ID, term_label), by = "ID") %>%
  mutate(
    term_label = factor(term_label, levels = term_labels_ordered),
    cluster = factor(cluster, levels = cluster_order)
  )

# Create a full grid of clusters and selected IDs (with their term labels) to account for missing combinations.
all_combos <- expand_grid(
  cluster = factor(cluster_order, levels = cluster_order),
  ID = selected_ids
) %>%
  left_join(labels_lookup, by = "ID") %>%
  mutate(term_label = factor(term_label, levels = term_labels_ordered))

# Merge the complete grid with the data and fill missing values with zeros.
plot_data_full <- all_combos %>%
  left_join(plot_data, by = c("cluster", "ID", "term_label")) %>%
  mutate(
    logp = replace_na(logp, 0),
    fold_enrichment = replace_na(fold_enrichment, 0),
    cluster = factor(cluster, levels = cluster_order)
  )

# ========================
# 2. Generate the dot plot and save as PDF
# ========================
pdf("Enrichment_DotPlot200.pdf", width = 8, height = 4.5)
ggplot(plot_data_full, aes(x = term_label, y = cluster, size = fold_enrichment, fill = logp)) +
  geom_point(shape = 21, color = "black", stroke = 0.4) +
  scale_size_continuous(range = c(2, 9), name = "Fold Enrichment") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, na.value = "grey90",
                       name = "-log10(adj p-value)",
                       guide = guide_colourbar(ticks = TRUE, frame.colour = "black")) +
  # Wrap x-axis labels every 40 characters (adjust width as needed)
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.line.y.left = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Enrichment of Selected IDs Across Clusters"
  )
dev.off()
# ========================


pediatric.at2 <- subset(pediatric, idents = "AT2")

Idents(pediatric.at2) <- "annotation_lvl3"
DefaultAssay(pediatric.at2) <- "RNA"
pediatric.at2 <- JoinLayers(pediatric.at2)
pediatric.at2 <- NormalizeData(pediatric.at2)
Idents(pediatric.at2) <- "annotation_lvl3"
trimmed.features = rownames(pediatric.at2)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)

a.markers <- FindMarkers(pediatric.at2, ident.1 = "FMO5+ AT2", ident.2 = "CFTR+ AT2", 
                         min.pct = 0.05, logfc.threshold = 0.05, features = trimmed.features)
a.markers$Significant <- ifelse(a.markers$p_val_adj < 0.05, "Significant", "Not Significant")
genes_to_label <- c("FMO5", "FASN", "PLA2G4","HMGCS1", "SCD", "PLCXD1", "CD36", "FGFR2", "SYBU", "SLC04C1","ICAM1", "CFTR", "CD44", "RND1", "SERPINA1", "LCN2", "DDX21", "MAP2K3")
write.csv(a.markers, "FMO5vCFTRThresholdedGenes.csv")

a.volcano_data <- a.markers %>%
  rownames_to_column(var = "gene") %>%
  mutate(logP = -log10(p_val_adj))

a.volcano_data$label <- ifelse(a.volcano_data$gene %in% genes_to_label, a.volcano_data$gene, NA)
pdf(file.path("./", paste0("FMO5+ v CFTR+ AT2 Volcano", ".pdf")), w=8, h=8)
ggplot(a.volcano_data, aes(x = avg_log2FC, y = logP, color = Significant)) +
  
  # Plot the points
  geom_point() +
  
  # Scale for significant and non-significant points
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  
  # Labels for specific genes with italicized text
  geom_label_repel(aes(label = label),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.size = 0.5,
                   direction = "y",
                   nudge_y = 2,
                   force = 5,
                   max.overlaps = Inf,
                   fontface = "italic") +  # Make gene names italic
  
  # Remove plot title and add axis labels
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  
  # Add a black border around the plot
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.0),  # Add black border
    axis.text = element_text(size = 12),  # Axis text: 8 pt Arial
    axis.title = element_text(size = 12),  # Axis labels: 8 pt Arial
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust plot margins
  )
dev.off()

perform_enrichment_analysis(pediatric.at2, output_dir = "/")
#######FMO5GeneSet#######
fmo5geneset <- list(c("FMO5",
                      "AC242426.2",
                      "CD36",
                      "IRX3",
                      "TOB2",
                      "SYBU",
                      "SCD",
                      "SLCO4C1",
                      "HHIP",
                      "TMEM268",
                      "CEP170",
                      "LIFR",
                      "GUCY1A1",
                      "HMGCS1",
                      "AC026992.2",
                      "CREBRF",
                      "PLCXD1",
                      "ID4",
                      "DEPTOR",
                      "HLF",
                      "KLHL24",
                      "ARID4A",
                      "SACM1L",
                      "SELENOP",
                      "IDI1",
                      "ASAH1",
                      "ACSS2",
                      "CA2",
                      "BNIP3L",
                      "PXMP4",
                      "ACAT2",
                      "VSIG10",
                      "FGFR2",
                      "INSIG1",
                      "GPR160",
                      "OSBPL11",
                      "SLC12A6",
                      "TBC1D2B",
                      "CNNM2",
                      "ACSS3",
                      "PDCD4",
                      "AC098588.3",
                      "C16orf89",
                      "RBL2",
                      "TSC22D3",
                      "AXIN2",
                      "SECISBP2L",
                      "FBXL20",
                      "FAM117B",
                      "RMST"))
############################
#########CFTRGeneSet########
cftrgeneset <- list(c("MAP2",
                      "ITGA2",
                      "CSF3",
                      "AC011511.2",
                      "ICAM1",
                      "NABP1",
                      "SERPINA1",
                      "TXNRD1",
                      "RND1",
                      "SERPINB1",
                      "TNFAIP3",
                      "TNFRSF12A",
                      "SGPP2",
                      "TM4SF1",
                      "TNC",
                      "SLC6A14",
                      "ANXA2",
                      "AL133268.4",
                      "LAMB3",
                      "TNFRSF21",
                      "BIRC3",
                      "IFNAR2",
                      "CXCL2",
                      "AC096564.1",
                      "CD83",
                      "LYN",
                      "JDP2",
                      "ZC3H12A",
                      "S100A10",
                      "PELI1",
                      "CTSB",
                      "TNIP1",
                      "KDM6B",
                      "ITPKC",
                      "CLIC5",
                      "HDAC7",
                      "AC003991.1",
                      "RELB",
                      "DDX21",
                      "NFKB2",
                      "SGMS2",
                      "EMP1",
                      "TGIF1",
                      "MAP2K3",
                      "MIR222HG",
                      "SLC7A1",
                      "SOD2",
                      "MED24",
                      "ANXA5",
                      "HIVEP2"))
############################





PDA.epi.alv <- AddModuleScore(PDA.epi.alv, features = cftrgeneset, name = "cftr", assay = "RNA")
PDA.epi.alv <- AddModuleScore(PDA.epi.alv, features = fmo5geneset, name = "fmo5", assay = "RNA")
pdf(file.path("./", paste0("PDA EPI ALV FMO5-50 Score Feature.pdf")), w=11, h=8.5)
FeaturePlot(PDA.epi.alv, features = "fmo51", order = TRUE, pt.size = 0.5) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("PDA EPI CFTR-50 Score Score Feature.pdf")), w=11, h=8.5) 
FeaturePlot(PDA.epi.alv, features = "cftr1", order = TRUE, pt.size = 0.5) & NoAxes() & NoLegend() &scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()



Idents(PDA.epi.alv) <- "annotation_lvl3"
PDA.at2 <- subset(PDA.epi.alv, idents = c("AT2 1", "AT2 2"))

DefaultAssay(PDA.at2) <- "RNA"
PDA.at2 <- NormalizeData(PDA.at2)
#AT2 markers
at2_markers <- c("SFTPC", "SFTPB", "ABCA3", "LAMP3")

# AT2 subtype 1 - Non-lipogenic
at2_1_nonlipogenic <- c("SECISBP2L", "ZBTB16", "FGFR2", "LIFR", "FBXL20", "TSC22D3", "FMO5", "HHIP", "KLHL24", "BNIP3L")

# AT2 subtype 1 - Lipogenic
at2_1_lipogenic <- c("CD36", "ACSS3", "ACSS2", "HMGCS1", "PLA2G4F", "FASN",  "ASAH1", "INSIG1", "SCD", "ACAT2")

# AT2 subtype 2 - Non-NF-kB/TNF
at2_2_nonnfkb <- c("MYOF", "EML4", "TEAD1", "VMP1", "SVIL", "FAM20A", "HIVEP2", "SGPP2", "SGMS2", "ELF3")

# AT2 subtype 2 - NF-kB/TNF regulated
at2_2_nfkb <- c("NAMPT","ITPKC","GPRC5A", "SERPINB1", "ICAM1", "CXCL2", "SLC6A14",   "CFTR", "RND1", "CSF3")

StandardizedDotPlot_NoScale <- function(seurat_obj, genes, grouping) {
  DotPlot(seurat_obj, features = genes, dot.scale = 20, scale = FALSE, group.by = grouping) +
    RotatedAxis() +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.4) +
    scale_color_distiller(palette = "YlGnBu", direction = 1, limits = c(0, 4)) +
    scale_size(range = c(3, 20), limits = c(0, 100)) +
    guides(
      size = guide_legend(title = "Percent Expressed", override.aes = list(shape = 21, colour = "black", fill = "black")),
      color = guide_colourbar(title = "Average Expression", ticks = TRUE, frame.colour = "black")
    ) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    )
}

PDA.at2$annotation_lvl3 <- factor(PDA.at2$annotation_lvl3, levels = c("AT2 1", "AT2 2"))
pdf(file.path("./", paste0("AT2 1 v AT2 2 DotPlot", ".pdf")), w=25, h=3.5)
StandardizedDotPlot_NoScale(PDA.at2, c(at2_1_nonlipogenic, at2_1_lipogenic, at2_2_nonnfkb, at2_2_nfkb), "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 16)) + NoLegend()
dev.off()



pdf(file.path("./", paste0("AT2 1 v AT2 2 DotPlot With Scale", ".pdf")), w=25, h=6.5)
StandardizedDotPlot_NoScale(PDA.at2, c(at2_1_nonlipogenic, at2_1_lipogenic, at2_2_nonnfkb, at2_2_nfkb), "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 16)) 
dev.off()

pdf(file.path("./", paste0("AT2 1 v AT2 2 DotPlot for DOIM", ".pdf")), w=28, h=3.5)
Dotplot_Zhiwei_Version_NoScale(PDA.at2, c(at2_markers, at2_1_nonlipogenic, at2_1_lipogenic, at2_2_nonnfkb, at2_2_nfkb), "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 16)) 
dev.off()
pdf(file.path("./", paste0("AT2 1 v AT2 2 DotPlot for DOIM Large", ".pdf")), w=28, h=6.5)
Dotplot_Zhiwei_Version_NoScale(PDA.at2, c(at2_markers, at2_1_nonlipogenic, at2_1_lipogenic, at2_2_nonnfkb, at2_2_nfkb), "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 16)) 
dev.off()

at2_2_core_nfkb <- c("LYN", "NFKB1", "RELB", "BCL3", "TNIP1", "TNFAIP3","NFKB2", "DDX21", "MAP2K3")
at2_2_core_tnfa <- c("MAP3K14","TNFRSF21", "BIRC3", "TNFRSF12A")



pdf(file.path("./", paste0("AT2 1 v AT2 2 NFKbTNFA", ".pdf")), w=9, h=3.5)
StandardizedDotPlot_NoScale(PDA.at2, c(at2_2_core_nfkb, at2_2_core_tnfa), "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 16)) + NoLegend()
dev.off()


pdf(file.path("./", paste0("AT2 1 v AT2 2 NFKbTNFA with Scale", ".pdf")), w=9, h=6.5)
StandardizedDotPlot_NoScale(PDA.at2, c(at2_2_core_nfkb, at2_2_core_tnfa), "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 26)) + theme(axis.text.y = element_text(size = 16)) 
dev.off()


#Figure 2H
Idents(PDA.mes) <- "at2class"
PDA.mes <- subset(PDA.mes, idents = "unknown", invert = TRUE)

pdf(file.path("./", paste0("PDA FastMNN Mesenchyme FastMNN Louvain by AT2 class", ".pdf")), w=33, h=8.5)
DimPlot(PDA.mes, label = TRUE, pt.size = 0.5, split.by = "at2class", group.by = "annotation_lvl3")
dev.off()

pdf(file.path("./", paste0("PDA FastMNN Mesenchyme LouvainAT2 class One", ".pdf")), w=11, h=8.5)
DimPlot(PDA.mes, label = FALSE, pt.size = 0.75, cells.highlight = WhichCells(subset(PDA.mes, idents = "one")), cols.highlight = "#F8766D", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Mesenchyme FastMNN LouvainAT2 class Two", ".pdf")), w=11, h=8.5)
DimPlot(PDA.mes, label = FALSE, pt.size = 0.75, cells.highlight = WhichCells(subset(PDA.mes, idents = "two")), cols.highlight = "#00A9FF", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Mesenchyme FastMNN LouvainAT2 class Three", ".pdf")), w=11, h=8.5)
DimPlot(PDA.mes, label = FALSE, pt.size = 0.75, cells.highlight = WhichCells(subset(PDA.mes, idents = "three")), cols.highlight = "#0CB702", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
PDA.mes.plot <- PDA.mes
PDA.mes.plot$annotation_lvl3 <- factor(PDA.mes.plot$annotation_lvl3, levels = c("Alveolar Fibroblast 1", "Alveolar Fibroblast 2","CCL2-Positive Fibroblast", "Fibroblast","Alveolar Myofibroblast",  "Adventitial-Like Mesenchyme",  "Pericyte"))
DefaultAssay(PDA.mes.plot) <- "RNA"
PDA.mes.plot <- NormalizeData(PDA.mes.plot)
PDA.mes.genes <- c("TNC","RUNX1", "CCL2", "ICAM1","FGF7", "FGF2", "FGF10","RSPO2","NFATC2","LDLR","PTGIR",  "IL6R", "TEAD4", "CCN1","RELB",  "CTHRC1", "COL1A1", "SFRP2", "SFRP4", "CXCL14")
#dotsize = 15
pdf(file.path("./", paste0("Mesenchyme dotplot", ".pdf")), w=12, h=4.0)
Dotplot_Zhiwei_Version(PDA.mes.plot, PDA.mes.genes, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 18))
dev.off()
pdf(file.path("./", paste0("Mesenchyme dotplot Large", ".pdf")), w=12, h=7.0)
Dotplot_Zhiwei_Version(PDA.mes.plot, PDA.mes.genes, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 18))
dev.off()

Idents(PDA.epi.alv) <- "at2class"
pdf(file.path("./", paste0("PDA FastMNN Alv Epi LouvainAT2 class One", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi.alv, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.epi.alv, idents = "one")), cols.highlight = "#F8766D", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Alv Epi  FastMNN LouvainAT2 class Two", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi.alv, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.epi.alv, idents = "two")), cols.highlight = "#00A9FF", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Alv Epi  FastMNN LouvainAT2 class Three", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi.alv, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.epi.alv, idents = "three")), cols.highlight = "#0CB702", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()

Idents(PDA.endothelial) <- "at2class"
pdf(file.path("./", paste0("PDA FastMNN Endothelial LouvainAT2 class One", ".pdf")), w=11, h=8.5)
DimPlot(PDA.endothelial, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.endothelial, idents = "one")), cols.highlight = "#F8766D", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Endothelial  FastMNN LouvainAT2 class Two", ".pdf")), w=11, h=8.5)
DimPlot(PDA.endothelial, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.endothelial, idents = "two")), cols.highlight = "#00A9FF", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Endothelial  FastMNN LouvainAT2 class Three", ".pdf")), w=11, h=8.5)
DimPlot(PDA.endothelial, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.endothelial, idents = "three")), cols.highlight = "#0CB702", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()

Idents(PDA.immune.subset) <- "annotation_lvl3"
PDA.immune.subset <- subset(PDA.immune, idents = "Basophils", invert = TRUE)
Idents(PDA.immune) <- "at2class"
pdf(file.path("./", paste0("PDA FastMNN Immune LouvainAT2 class One", ".pdf")), w=11, h=8.5)
DimPlot(PDA.immune, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.immune, idents = "one")), cols.highlight = "#F8766D", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Immune  FastMNN LouvainAT2 class Two", ".pdf")), w=11, h=8.5)
DimPlot(PDA.immune, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.immune, idents = "two")), cols.highlight = "#00A9FF", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Immune  FastMNN LouvainAT2 class Three", ".pdf")), w=11, h=8.5)
DimPlot(PDA.immune, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.immune, idents = "three")), cols.highlight = "#0CB702", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()


Idents(PDA.epi) <- "annotation_lvl3"
PDA.epi <- RenameIdents(PDA.epi, "AT2 3" = "AT2 2", "Terminal Bronchiole Cell (RASC or TRB)" = "Respiratory Bronchiole Cell","Goblet-Like Secretory" = "Respiratory Bronchiole Cell")
PDA.epi <- AddMetaData(PDA.epi, col.name = "annotation_lvl3", Idents(PDA.epi))
pdf("PDA Epithelium Mixed Resolution Annotation No Label.pdf", width = 11, height = 8.5)
DimPlot(PDA.epi, pt.size = 0.2, group.by = "annotation_lvl3", label = FALSE) + NoLegend() + NoAxes()
dev.off()


Idents(PDA.epi) <- "at2class"
pdf(file.path("./", paste0("PDA FastMNN Epi LouvainAT2 class One", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.epi, idents = "one")), cols.highlight = "#F8766D", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Epi  FastMNN LouvainAT2 class Two", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.epi, idents = "two")), cols.highlight = "#00A9FF", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("PDA FastMNN Epi  FastMNN LouvainAT2 class Three", ".pdf")), w=11, h=8.5)
DimPlot(PDA.epi, label = FALSE, pt.size = 0.5, cells.highlight = WhichCells(subset(PDA.epi, idents = "three")), cols.highlight = "#0CB702", sizes.highlight = 0.5) + NoLegend() + NoAxes()
dev.off()


#Figure 3

#Fig. S3A
library(geneBasisR)
celltype_mapping = get_celltype_mapping(Pediatric.sce ,celltype.id = "annotation_lvl3", genes.selection = xenium.pediatric.trimmed, batch = "ID", return.stat = F)
desired_order <- c("FMO5+ AT2", "CFTR+ AT2", "Transitional 1", "Transitional 2", "AT1 1", "AT1 2", "AT1 3", "Respiratory Bronchiole", "Basal",
  "Multiciliated", "Neuroendocrine", "Macrophages", "TREM2+ Macrophages",
  "Monocytes", "CD1c+ DC", "pDC", "NK Cell", "CD4+ T-Cell", "CD8+ T-Cell",
  "B-Cells", "Neutrophils", "Lymphatic Endothelial", "Arterial Endothelial",
  "Venous Endothelial", "Capillary Endothelial", "Capillary Aerocyte",
  "Alveolar Fibroblast 1", "Alveolar Fibroblast 2" ,"Alveolar Fibroblast 3",
  "Alveolar Myofibroblast", "Adventitial Mesenchyme",
  "Pericyte")

pdf(file.path("./", paste0("Confusion Matrix", ".pdf")), w=12, h=12)
p <- p +
scale_x_discrete(limits = desired_order) +
scale_y_discrete(limits = rev(desired_order)) +  # reverse for heatmap readability
theme(
axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
axis.text.y = element_text(size = 14, face = "bold"),
axis.title.x = element_text(size = 16, face = "bold"),
axis.title.y = element_text(size = 16, face = "bold"),
plot.title  = element_text(size = 18, face = "bold", hjust = 0.5)
)
p
dev.off()


#####load Xenium samples#####
path <- "~/Dropbox (University of Michigan)/0_Xenium/11155-TF/output-XETG00077__0040306__11155-TF-1_ROI_A2__20240621__173354"
data <- ReadXenium(
data.dir = path,
type = c("centroids", "segmentations"),
)
assay <- "Xenium"
segmentations.data <- list(
"centroids" = CreateCentroids(data$centroids),
"segmentation" = CreateSegmentation(data$segmentations)
)
coords <- CreateFOV(
coords = segmentations.data,
type = c("segmentation", "centroids"),
molecules = data$microns,
assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
#xenium.obj <- AddMetaData(xenium.obj, col.name = "oUMI", Idents(xenium.obj))
#xenium.obj <- RenameCells(xenium.obj, add.cell.id = "A2")
xenium.obj.2 <- xenium.obj
features_to_keep <- intersect(desired_features, rownames(xenium.obj))
xenium.obj <- subset(xenium.obj, features = features_to_keep)
xenium.obj[["Xenium_All"]] <- xenium.obj.2[["Xenium"]]
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
fov <- "test"
xenium.obj[[fov]] <- coords
TF11155.1A2 <- xenium.obj
path <- "~/Dropbox (University of Michigan)/0_Xenium/11155-TF/output-XETG00077__0040306__11155-TF-1_ROI_A1__20240621__173354"
data <- ReadXenium(
data.dir = path,
type = c("centroids", "segmentations"),
)
assay <- "Xenium"
segmentations.data <- list(
"centroids" = CreateCentroids(data$centroids),
"segmentation" = CreateSegmentation(data$segmentations)
)
coords <- CreateFOV(
coords = segmentations.data,
type = c("segmentation", "centroids"),
molecules = data$microns,
assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
#xenium.obj <- AddMetaData(xenium.obj, col.name = "oUMI", Idents(xenium.obj))
#xenium.obj <- RenameCells(xenium.obj, add.cell.id = "A1")
xenium.obj.2 <- xenium.obj
features_to_keep <- intersect(desired_features, rownames(xenium.obj))
xenium.obj <- subset(xenium.obj, features = features_to_keep)
xenium.obj[["Xenium_All"]] <- xenium.obj.2[["Xenium"]]
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
fov <- "test"
xenium.obj[[fov]] <- coords
TF11155.1A1 <- xenium.obj
path <- "~/Dropbox (University of Michigan)/0_Xenium/11155-TF/output-XETG00077__0040306__11155-TF-1_ROI_B1__20240621__173354"
data <- ReadXenium(
data.dir = path,
type = c("centroids", "segmentations"),
)
assay <- "Xenium"
segmentations.data <- list(
"centroids" = CreateCentroids(data$centroids),
"segmentation" = CreateSegmentation(data$segmentations)
)
coords <- CreateFOV(
coords = segmentations.data,
type = c("segmentation", "centroids"),
molecules = data$microns,
assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
#xenium.obj <- AddMetaData(xenium.obj, col.name = "oUMI", Idents(xenium.obj))
#xenium.obj <- RenameCells(xenium.obj, add.cell.id = "B1")
xenium.obj.2 <- xenium.obj
features_to_keep <- intersect(desired_features, rownames(xenium.obj))
xenium.obj <- subset(xenium.obj, features = features_to_keep)
xenium.obj[["Xenium_All"]] <- xenium.obj.2[["Xenium"]]
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
fov <- "test"
xenium.obj[[fov]] <- coords
TF11155.1B1 <- xenium.obj
path <- "~/Dropbox (University of Michigan)/0_Xenium/11155-TF/output-XETG00077__0040306__11155-TF-1_ROI_B2__20240621__173354"
data <- ReadXenium(
data.dir = path,
type = c("centroids", "segmentations"),
)
assay <- "Xenium"
segmentations.data <- list(
"centroids" = CreateCentroids(data$centroids),
"segmentation" = CreateSegmentation(data$segmentations)
)
coords <- CreateFOV(
coords = segmentations.data,
type = c("segmentation", "centroids"),
molecules = data$microns,
assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
#xenium.obj <- AddMetaData(xenium.obj, col.name = "oUMI", Idents(xenium.obj))
#xenium.obj <- RenameCells(xenium.obj, add.cell.id = "B2")
xenium.obj.2 <- xenium.obj
features_to_keep <- intersect(desired_features, rownames(xenium.obj))
xenium.obj <- subset(xenium.obj, features = features_to_keep)
xenium.obj[["Xenium_All"]] <- xenium.obj.2[["Xenium"]]
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
fov <- "test"
xenium.obj[[fov]] <- coords
TF11155.1B2 <- xenium.obj
path <- "~/Dropbox (University of Michigan)/0_Xenium/11155-TF/output-XETG00077__0040306__11155-TF-1_ROI_C__20240621__173354"
data <- ReadXenium(
data.dir = path,
type = c("centroids", "segmentations"),
)
assay <- "Xenium"
segmentations.data <- list(
"centroids" = CreateCentroids(data$centroids),
"segmentation" = CreateSegmentation(data$segmentations)
)
coords <- CreateFOV(
coords = segmentations.data,
type = c("segmentation", "centroids"),
molecules = data$microns,
assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
#xenium.obj <- AddMetaData(xenium.obj, col.name = "oUMI", Idents(xenium.obj))
#xenium.obj <- RenameCells(xenium.obj, add.cell.id = "C1")
xenium.obj.2 <- xenium.obj
features_to_keep <- intersect(desired_features, rownames(xenium.obj))
xenium.obj <- subset(xenium.obj, features = features_to_keep)
xenium.obj[["Xenium_All"]] <- xenium.obj.2[["Xenium"]]
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
fov <- "test"
xenium.obj[[fov]] <- coords


TF11155.1C1 <- xenium.obj
cell_names <- colnames(TF11155.1A1)
TF11155.1A1 <- AddMetaData(TF11155.1A1, metadata = cell_names, col.name = "cell_names")
TF11155.1A1 <- AddMetaData(TF11155.1A1, metadata = "1A1", col.name = "ID")
TF11155.1A1 <- AddMetaData(TF11155.1A1, metadata = "S016", col.name = "Sample")
TF11155.1A1 <- RenameCells(TF11155.1A1, add.cell.id = "Obj1")
cell_names <- colnames(TF11155.1A2)
TF11155.1A2 <- AddMetaData(TF11155.1A2, metadata = cell_names, col.name = "cell_names")
TF11155.1A2 <- AddMetaData(TF11155.1A2, metadata = "1A2", col.name = "ID")
TF11155.1A2 <- AddMetaData(TF11155.1A2, metadata = "S022", col.name = "Sample")
TF11155.1A2 <- RenameCells(TF11155.1A2, add.cell.id = "Obj2")

cell_names <- colnames(TF11155.1B1)
TF11155.1B1 <- AddMetaData(TF11155.1B1, metadata = cell_names, col.name = "cell_names")
TF11155.1B1 <- AddMetaData(TF11155.1B1, metadata = "1B1", col.name = "ID")
TF11155.1B1 <- AddMetaData(TF11155.1B1, metadata = "S010", col.name = "Sample")
TF11155.1B1 <- RenameCells(TF11155.1B1, add.cell.id = "Obj3")

cell_names <- colnames(TF11155.1B2)
TF11155.1B2 <- AddMetaData(TF11155.1B2, metadata = cell_names, col.name = "cell_names")
TF11155.1B2 <- AddMetaData(TF11155.1B2, metadata = "1B2", col.name = "ID")
TF11155.1B2 <- AddMetaData(TF11155.1B2, metadata = "S020", col.name = "Sample")
TF11155.1B2 <- RenameCells(TF11155.1B2, add.cell.id = "Obj4")

cell_names <- colnames(TF11155.1C1)
TF11155.1C1 <- AddMetaData(TF11155.1C1, metadata = cell_names, col.name = "cell_names")
TF11155.1C1 <- AddMetaData(TF11155.1C1, metadata = "1C1", col.name = "ID")
TF11155.1C1 <- AddMetaData(TF11155.1C1, metadata = "S013", col.name = "Sample")
TF11155.1C1 <- RenameCells(TF11155.1C1, add.cell.id = "Obj5")

options(future.globals.maxSize = 8000 * 1024^3)
xenium.obj <- merge(TF11155.1A1, c(TF11155.1A2, TF11155.1B1, TF11155.1B2, TF11155.1C1))
DefaultAssay(xenium.obj) <- "Xenium_All"
xenium.obj <- JoinLayers(xenium.obj)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 25) #used 50 on the first iteration - eliminated too many AT1 cells.
#xenium.obj=xenium.obj[,unname(which(colSums(GetAssayData(xenium.obj))!=0))] #remove cells in which counts log_umi = 0
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium_All")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)
pdf(file.path("./", paste0("Louvain 0.3", ".pdf")), w=11, h=8.5)
DimPlot(xenium.obj, label = TRUE) + NoLegend()
dev.off()

pdf(file.path("./", paste0("Louvain Split 0.3", ".pdf")), w=55, h=8.5)
DimPlot(xenium.obj, label = TRUE, split.by = "ID") + NoLegend()
dev.off()
xenium.obj <- FindClusters(xenium.obj, resolution = 0.8)

pdf(file.path("./", paste0("Louvain 0.8", ".pdf")), w=11, h=8.5)
DimPlot(xenium.obj, label = TRUE) + NoLegend()
dev.off()

pdf(file.path("./", paste0("Louvain Split 0.8", ".pdf")), w=55, h=8.5)
DimPlot(xenium.obj, label = TRUE, split.by = "ID") + NoLegend()
dev.off()




Idents(xenium.obj) <- "seurat_clusters"
seurat.obj.markers <- FindAllMarkers(xenium.obj, min.pct = 0.25, logfc.threshold = 0.25)
xenium.obj <- ScaleData(xenium.obj)

top20.markers <- seurat.obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top20.markers$cluster <- gsub("/", "-", top20.markers$cluster)
genes.for.dotplot <- make.unique(top20.markers$gene)

pdf(file.path("./", paste0("FullData_Top10_Markers_Per_Cluster Res 0.8.pdf")), w=100, h=10)
print(Dotplot_Zhiwei_Version(xenium.obj, genes.for.dotplot, "seurat_clusters"))
dev.off()


#####Figure S3A######
xenium.at2 <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/RObjects/XeniumAT2HandAnnotated.RDS")
xenium.at2 <- RenameIdents(xenium.at2, c("Unspecified AT2" = "Unassigned AT2", "SCGB3A2 High AT2" = "Unassigned AT2", "DMBT1 High AT2" = "Unassigned AT2"))
xenium.at2 <- AddMetaData(xenium.at2, col.name = "Hand_Annotation_S", Idents(xenium.at2))
xenium.at2$Hand_Annotation_S <- factor(xenium.at2$Hand_Annotation_S, levels = c("CFTR+ AT2", "FMO5+ AT2", "Unassigned AT2"))
pdf(file.path("./", paste0("AT2 Xenium Dotplot.pdf")), w=9, h=3.5)
Dotplot_Zhiwei_Version(xenium.at2, at2.xenium.markers, "Hand_Annotation_S") + theme(axis.text.x = element_text(face = "italic", size = 18))
dev.off()
pdf(file.path("./", paste0("AT2 Xenium Dotplot Large.pdf")), w=9, h=4.5)
Dotplot_Zhiwei_Version(xenium.at2, at2.xenium.markers, "Hand_Annotation_S") + theme(axis.text.x = element_text(face = "italic", size = 18))
dev.off()
write.csv(table(xenium.at2$Hand_Annotation_S), "HandAnnotatedAT2Simplified.csv")
xenium.at2$Hand_Annotation_S  <- as.character(xenium.at2$Hand_Annotation_S )
xenium.mes <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/RObjects/XeniumMesHandAnnotated.RDS")
xenium.mes <- RenameIdents(xenium.mes, c("SCARA5+ Pericyte" = "Pericyte", "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme"))
xenium.mes <- AddMetaData(xenium.mes, col.name = "Hand_Annotation_S", Idents(xenium.mes))
xenium.mes$Hand_Annotation_S <- factor(xenium.mes$Hand_Annotation_S, levels = c("Alveolar Fibroblast", "Alveolar Myofibroblast", "Pericyte", "Adventitial", "CCL2+ Mesenchyme"))

mes.xenium.markers <- c("ITGA8", "BMP5", "MYH11", "ACTA2", "PDGFRB", "COL4A1","C3", "OGN", "CCL2", "FOSB", "JUNB", "CCN1")

pdf(file.path("./", paste0("Mes Xenium Dotplot.pdf")), w=9, h=4.0)
Dotplot_Zhiwei_Version(xenium.mes, mes.xenium.markers, "Hand_Annotation_S") + theme(axis.text.x = element_text(face = "italic", size = 18))
dev.off()
write.csv(table(xenium.mes$Hand_Annotation_S), "HandAnnotatedMesSimplified.csv")
xenium.mes$Hand_Annotation_S  <- as.character(xenium.mes$Hand_Annotation_S )


xenium.obj <- AddMetaData(xenium.obj, col.name = "Hand_S", c(xenium.at2$Hand_Annotation_S, xenium.mes$Hand_Annotation_S))

xenium.obj.s <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/RObjects/PediatricXenium1to5SCTTransformedHandAnnotatedFixed.RDS")
Idents(xenium.obj.s) <- "All_Hand_Annotation"
xenium.obj.s <- RenameIdents(xenium.obj.s, c("SCARA5+ Pericyte" = "Pericyte", "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme", "Unspecified AT2" = "Unassigned AT2", "SCGB3A2 High AT2" = "Unassigned AT2", "DMBT1 High AT2" = "Unassigned AT2"))
xenium.obj.s.s <- subset(xenium.obj.s, idents = c("Alveolar Fibroblast", "Alveolar Myofibroblast", "Pericyte", "Adventitial", "CCL2+ Mesenchyme", "CFTR+ AT2", "FMO5+ AT2", "Unassigned AT2"))
xenium.obj.s <- AddMetaData(xenium.obj.s, col.name = "Hand_Annotation_S", Idents(xenium.obj.s.s))
distinct_palette <- c(
  "#5ed6e6", "#c9262e", "#328d36", "#eec072",
  "#7163a7", "#56bd85", "#b8b5f2", "#6b2620", "lightgrey"
)
pdf(file.path("./", paste0("All Xenium DimPlot with Hand Annotation AT2 Mes No Label.pdf")), w=11.5, h=8.5)
DimPlot(xenium.obj.s, group.by = "Hand_Annotation_S", label = FALSE, pt.size = 0.2, raster = FALSE, cols = distinct_palette) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("All Xenium DimPlot with Hand Annotation AT2 Mes Label.pdf")), w=11.5, h=8.5)
DimPlot(xenium.obj.s, group.by = "Hand_Annotation_S", label = TRUE, pt.size = 0.2, raster = FALSE, cols = distinct_palette) + NoLegend() + NoAxes()
dev.off()

#####Figure 3B######
#Load Xenium File, do this for each
path <- "path/to/file"
data <- ReadXenium(
  data.dir = path,
  type = c("centroids", "segmentations"),
)
assay <- "Xenium"
segmentations.data <- list(
  "centroids" = CreateCentroids(data$centroids),
  "segmentation" = CreateSegmentation(data$segmentations)
)

coords <- CreateFOV(
  coords = segmentations.data,
  type = c("segmentation", "centroids"),
  molecules = data$microns,
  assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)


xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])

data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])

fov <- "test"
xenium.obj[[fov]] <- coords
TF11155.1A1 <- xenium.obj
Pediatric <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/RObjects/Pediatric Cell Atlas/FullPDAObjectFixedforPaperRound1.RDS")

library(dplyr)

Pediatric.ref <- SCTransform(Pediatric, ncells = 3000, verbose = FALSE) %>%
       RunPCA(verbose = FALSE) %>%
       RunUMAP(dims = 1:30)
Pediatric.ref <- RenameIdents(Pediatric.ref, "Alveolar Fibroblast 1" = "Alveolar Fibroblast", "Alveolar Fibroblast 2" = "Alveolar Fibroblast", "Terminal Bronchiole Cell (RASC or TRB)" = "Respiratory Bronchiole", "Respiratory Bronchiole Cell" = "Respiratory Bronchiole", "Goblet-Like Secretory" = "Respiratory Bronchiole")
Pediatric.ref <- AddMetaData(Pediatric.ref, col.name = "annotation_lvl3", Idents(Pediatric.ref))
object_names <- c("TF11155.1A1", "TF11155.1A2", "TF11155.1B1", "TF11155.1B2", "TF11155.1C1")

# Run the function on each object
for (obj_name in object_names) {
  seurat_obj <- get(obj_name)
  process_seurat.object.integrated.ref(seurat_obj, prefix = obj_name)
}

process_seurat.object.integrated.ref <- function(seurat.obj, prefix) {
  message("🔁 Starting processing for: ", prefix)
  
  # Backup the original assay
  
  # Basic filtering on total counts
  seurat.obj <- subset(seurat.obj, subset = nCount_Xenium > 25)
  seurat.obj <- seurat.obj[, unname(which(colSums(GetAssayData(seurat.obj)) != 0))]
  
  # Drop transcripts with count < 3 per cell
  
  
  DefaultAssay(seurat.obj) <- "Xenium_All"
  
  # Violin plot for QC
  pdf(file.path("./", paste0(prefix, "_ViolinPlots.pdf")), w = 11, h = 8.5)
  print(VlnPlot(seurat.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0))
  dev.off()
  
  # SCTransform with top 500 variable genes
  seurat.obj <- SCTransform(seurat.obj, assay = "Xenium_All", variable.features.n = 250)
  
  # Remove common high-expression junk genes
  var.features <- VariableFeatures(seurat.obj)
  bad.pattern <- "^MALAT1$|^MT-|^RPS|^RPL|^GAPDH$|^ACTB$"
  var.features <- var.features[!grepl(bad.pattern, var.features)]
  
  # Also remove genes expressed in >95% of cells
  expr.frac <- rowMeans(GetAssayData(seurat.obj, slot = "counts") > 0)
  informative.genes <- names(expr.frac[expr.frac < 0.95])
  var.features <- intersect(var.features, informative.genes)
  
  # Finalize filtered feature set
  VariableFeatures(seurat.obj) <- var.features
  
  # PCA + UMAP + Clustering
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(seurat.obj))
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:30)
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
  
  pdf(file.path("./", paste0(prefix, "_Louvain.pdf")), w = 11, h = 8.5)
  print(DimPlot(seurat.obj, label = TRUE) + NoLegend())
  dev.off()
  
  # Normalize for marker plotting
  seurat.obj <- NormalizeData(seurat.obj)
  
  # Label Transfer
  anchors <- FindTransferAnchors(
    reference = Pediatric.ref,
    query = seurat.obj,
    normalization.method = "SCT",
    dims = 1:20,
    features = VariableFeatures(seurat.obj),
    k.filter = 200
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = Pediatric.ref$annotation_lvl3,
    prediction.assay = FALSE,
    weight.reduction = seurat.obj[["pca"]],
    dims = 1:20
  )
  
  seurat.obj$predicted.id <- predictions$predicted.id
  
  # Visualize prediction results
  pdf(file.path("./", paste0(prefix, "_PredictedID_Grouped.pdf")), w = 11, h = 8.5)
  print(DimPlot(seurat.obj, group.by = "predicted.id", label = TRUE) + NoLegend())
  dev.off()
  
  pdf(file.path("./", paste0(prefix, "_PredictedID_Split.pdf")), w = 100, h = 3)
  print(DimPlot(seurat.obj, split.by = "predicted.id", group.by = "seurat_clusters", label = TRUE) + NoLegend())
  dev.off()
  
  pdf(file.path("./", paste0(prefix, "_ImageDimPlot_PredictedID_Split.pdf")), w = 110, h = 85)
  print(ImageDimPlot(seurat.obj, cols = "polychrome", size = 2, split.by = "predicted.id", group.by = "predicted.id"))
  dev.off()
  
  # Marker analysis
  Idents(seurat.obj) <- "predicted.id"
  seurat.obj.markers <- FindAllMarkers(seurat.obj, min.pct = 0.25, logfc.threshold = 0.25)
  seurat.obj <- ScaleData(seurat.obj)
  
  top20.markers <- seurat.obj.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  top20.markers$cluster <- gsub("/", "-", top20.markers$cluster)
  genes.for.dotplot <- make.unique(top20.markers$gene)
  
  pdf(file.path("./", paste0(prefix, "_Top20_Markers_Per_Cluster.pdf")), w = 100, h = 10)
  print(Dotplot_Zhiwei_Version(seurat.obj, genes.for.dotplot, "predicted.id"))
  dev.off()
  
  # Export for Xenium Explorer
  new_colnames <- seurat.obj@meta.data$cell_names
  colnames(seurat.obj) <- new_colnames
  
  cell_data <- data.frame(
    cell_id = names(seurat.obj$predicted.id),
    group = as.character(seurat.obj$predicted.id)
  )
  
  write.csv(cell_data, paste0(prefix, "CellIdentitiesForXeniumExplorer.csv"), row.names = FALSE)
  saveRDS(seurat.obj, paste0(prefix, "LabelTransfer_FilteredFeatures_CleanedCounts.RDS"))
}

xenium.obj.1 <- RenameCells(xenium.obj.1, add.cell.id = "Obj1")
xenium.obj.2 <- RenameCells(xenium.obj.2, add.cell.id = "Obj2")
xenium.obj.3 <- RenameCells(xenium.obj.3, add.cell.id = "Obj3")
xenium.obj.4 <- RenameCells(xenium.obj.4, add.cell.id = "Obj4")
xenium.obj.5 <- RenameCells(xenium.obj.5, add.cell.id = "Obj5")

xenium.obj <- AddMetaData(xenium.obj, col.name = "predicted.id", c(xenium.obj.1$predicted.id, xenium.obj.2$predicted.id, xenium.obj.3$predicted.id, xenium.obj.4$predicted.id, xenium.obj.5$predicted.id))

pdf(file.path("./", paste0("Louvain Label Transfered Identities No Label", ".pdf")), w=11, h=8.5)
DimPlot(xenium.obj, group.by = "predicted.id", label = FALSE) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("Louvain Label Transfered Identities Label", ".pdf")), w=11, h=8.5)
DimPlot(xenium.obj, group.by = "predicted.id", label = TRUE) + NoLegend() + NoAxes()
dev.off()
####Niche Analysis######
xenium.obj.5 <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-01/Variable Features Label Transfer Simplified Annotation CellCut 25/TF11155.1C1LabelTransfer_FilteredFeatures_CleanedCounts.RDS")
xenium.obj.4 <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-01/Variable Features Label Transfer Simplified Annotation CellCut 25/TF11155.1B2LabelTransfer_FilteredFeatures_CleanedCounts.RDS")
xenium.obj.3 <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-01/Variable Features Label Transfer Simplified Annotation CellCut 25/TF11155.1B1LabelTransfer_FilteredFeatures_CleanedCounts.RDS")
xenium.obj.2 <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-01/Variable Features Label Transfer Simplified Annotation CellCut 25/TF11155.1A2LabelTransfer_FilteredFeatures_CleanedCounts.RDS")
xenium.obj.1 <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-01/Variable Features Label Transfer Simplified Annotation CellCut 25/TF11155.1A1LabelTransfer_FilteredFeatures_CleanedCounts.RDS")

xenium.obj.1 <- TF11155.1A1
xenium.obj.2 <- TF11155.1A2
xenium.obj.3 <- TF11155.1B1
xenium.obj.4 <- TF11155.1B2
xenium.obj.5 <- TF11155.1C1


xenium.obj.1 <- BuildNicheAssay(object = xenium.obj.1, , group.by = "predicted.id", niches.k = 4, neighbors.k = 50, fov = "test")
xenium.obj.2 <- BuildNicheAssay(object = xenium.obj.2, , group.by = "predicted.id", niches.k = 4, neighbors.k = 50, fov = "test")
xenium.obj.3 <- BuildNicheAssay(object = xenium.obj.3, , group.by = "predicted.id", niches.k = 4, neighbors.k = 50, fov = "test")
xenium.obj.4 <- BuildNicheAssay(object = xenium.obj.4, , group.by = "predicted.id", niches.k = 4, neighbors.k = 50, fov = "test")
xenium.obj.5 <- BuildNicheAssay(object = xenium.obj.5, , group.by = "predicted.id", niches.k = 4, neighbors.k = 50, fov = "test")
write.csv(table(xenium.obj.1$predicted.id, xenium.obj.1$niches), "obj_1NicheCompositionK=4.csv")
write.csv(table(xenium.obj.2$predicted.id, xenium.obj.2$niches), "obj_2NicheCompositionK=4.csv")
write.csv(table(xenium.obj.3$predicted.id, xenium.obj.3$niches), "obj_3NicheCompositionK=4.csv")
write.csv(table(xenium.obj.4$predicted.id, xenium.obj.4$niches), "obj_4NicheCompositionK=4.csv")
write.csv(table(xenium.obj.5$predicted.id, xenium.obj.5$niches), "obj_5NicheCompositionK=4.csv")
pdf(file.path("./", paste0("obj_1NichePlotK=4", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.1, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_2NichePlotK=4", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.2, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_3NichePlotK=4", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.3, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_4NichePlotK=4", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.4, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_5NichePlotK=4", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.5, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
seurat_objects <- list(xenium.obj.1, xenium.obj.2, xenium.obj.3, xenium.obj.4, xenium.obj.5)

# Function to extract, tabulate, and prefix column names
extract_and_tabulate <- function(seurat_obj, prefix) {
  # Extract metadata
  metadata <- seurat_obj@meta.data[, c("predicted.id", "niches")]
  # Create contingency table
  tab <- table(metadata$predicted.id, metadata$niche)
  # Convert to data frame
  df <- as.data.frame.matrix(tab)
  # Add 'predicted.id' as a column
  df$predicted.id <- rownames(df)
  # Reorder columns to have 'predicted.id' first
  df <- df[, c("predicted.id", setdiff(names(df), "predicted.id"))]
  # Rename columns
  colnames(df)[-1] <- paste0(prefix, "_", colnames(df)[-1])
  return(df)
}

# Apply the function to each Seurat object with appropriate prefix
df_list <- map2(seurat_objects, 1:5, extract_and_tabulate)

# Merge all data frames by 'predicted.id'
merged_df <- reduce(df_list, full_join, by = "predicted.id")

# Replace NAs with 0
merged_df[is.na(merged_df)] <- 0

# View the merged table
print(merged_df)

write.csv(merged_df, "merged_Seurat_Niche_CompositionK=4.csv", row.names = FALSE)
rownames(merged_df) <- merged_df[,1]
merged_df <- merged_df[,-1]
Feature <- as.matrix(merged_df)
cor_matrix <- round(cor(Feature), digits = 3)
pheatmap(cor_matrix, display_numbers = TRUE)
pdf(file.path("./", paste0("NicheCorrelationHeatmapK=4", ".pdf")), w=11, h=8.5)
p1 <- pheatmap(cor_matrix, display_numbers = TRUE)
p1
dev.off()
#Niches aggregated in Excel based on Pearson's correlation plot (uploaded .csv below)
######Alveolar Niche 1 and 2 comparison######
# Load required libraries
library(tidyverse)
library(ggtext)

# --- [Read and Normalize Data] ---
df <- read.csv("2025-05-21 Variable Features SCT CellCut 25 K = 4 Neighbors = 50 For Upload to R.csv", row.names = 1)

# Convert counts to percent composition
df_percent <- sweep(df, 2, colSums(df), FUN = "/") * 100

# Define F and C column groups
f_cols <- grep("^F", colnames(df_percent), value = TRUE)
c_cols <- grep("^C", colnames(df_percent), value = TRUE)

# --- [Calculate P-values Separately] ---
pval_df <- data.frame(
  CellType = rownames(df_percent),
  p_value = apply(df_percent, 1, function(x) {
    f_vals <- as.numeric(x[f_cols])
    c_vals <- as.numeric(x[c_cols])
    t.test(f_vals, c_vals, alternative = "two.sided", var.equal = TRUE)$p.value
  }),
  row.names = NULL
)

# --- [Summary Statistics by Cell Type and Group] ---
summary_df <- df_percent %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-CellType, names_to = "Group", values_to = "Percent") %>%
  mutate(GroupType = ifelse(str_starts(Group, "F"), "F", "C")) %>%
  group_by(CellType, GroupType) %>%
  summarise(mean = mean(Percent), sd = sd(Percent), .groups = "drop") %>%
  pivot_wider(names_from = GroupType, values_from = c(mean, sd))

# --- [Custom Cell Type Order] ---
desired_order <- c("FMO5+ AT2", "Alveolar Fibroblast", "CFTR+ AT2", "CCL2-Positive Fibroblast", 
                   "Transitional 1", "Transitional 2", "AT1 1", "AT1 2", "AT1 3", 
                   "Respiratory Bronchiole", "Basal", "Multiciliated", "Neuroendocrine", 
                   "Macrophages", "Activated Macrophages", "Monocytes", "CD1c+ DC", "pDC", 
                   "NK Cell", "CD4+ T-Cell", "CD8+ T-Cell", "B-Cells", "Neutrophils", 
                   "Lymphatic Endothelial", "Arterial Endothelial", "Venous Endothelial", 
                   "Capillary Endothelial", "Capillary Aerocyte", "Fibroblast", 
                   "Alveolar Myofibroblast", "Adventitial Mesenchyme", "Pericyte")

# --- [Apply Order, Join P-values, Compute Significance Label] ---
summary_df <- summary_df %>%
  filter(CellType %in% desired_order) %>%
  left_join(pval_df, by = "CellType") %>%
  mutate(
    CellType = factor(CellType, levels = desired_order),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label_y = pmax(mean_F + sd_F, mean_C + sd_C) + 2
  ) %>%
  arrange(CellType)

# --- [Plotting] ---
pdf("Niche Analysis.pdf", width = 12, height = 6)
ggplot(summary_df, aes(x = CellType)) +
  geom_bar(aes(y = mean_F), stat = "identity", fill = "#F8766D", color = "black", size = 0.2, width = 0.4,
           position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(ymin = mean_F - sd_F, ymax = mean_F + sd_F), width = 0.2, size = 0.2,
                position = position_nudge(x = -0.2)) +
  geom_bar(aes(y = mean_C), stat = "identity", fill = "#ABA300", color = "black", size = 0.2, width = 0.4,
           position = position_nudge(x = 0.2)) +
  geom_errorbar(aes(ymin = mean_C - sd_C, ymax = mean_C + sd_C), width = 0.2, size = 0.2,
                position = position_nudge(x = 0.2)) +
  geom_text(aes(y = label_y, label = significance), size = 6, fontface = "bold") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.line = element_line(size = 0.5, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_markdown(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    y = "Percent Composition",
    x = "",
    title = "<span style='color:#F8766D'>Alveolar Niche 1</span> v <span style='color:#ABA300'>Alveolar Niche 2</span>"
  ) +
  expand_limits(y = 0)
dev.off()

#######Figure S3E - G#########
xenium.obj.1.hand <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-07/1A1HandAnnotatedAll.RDS")
xenium.obj.2.hand <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-07/1A2HandAnnotatedAll.RDS")
xenium.obj.3.hand <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-07/1B1HandAnnotatedAll.RDS")
xenium.obj.4.hand <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-07/1B2HandAnnotatedAll.RDS")
xenium.obj.5.hand <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/2025-04-07/1C1HandAnnotatedAll.RDS")

Idents(xenium.obj.1.hand) <- "All_Hand_Annotation"
xenium.obj.1.hand.subset <- subset(xenium.obj.1.hand, idents = c("FMO5+ AT2", "CFTR+ AT2", "CCL2+ Mesenchyme", "CCL2+ Mesenchyme IL6", "Alveolar Fibroblast"))
Idents(xenium.obj.1.hand.subset) <- "All_Hand_Annotation"
xenium.obj.1.hand.subset <- RenameIdents(xenium.obj.1.hand.subset, "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme")
xenium.obj.1.hand.subset <- AddMetaData(xenium.obj.1.hand.subset, col.name = "All_Hand_Annotation", Idents(xenium.obj.1.hand.subset))

Idents(xenium.obj.2.hand) <- "All_Hand_Annotation"
xenium.obj.2.hand.subset <- subset(xenium.obj.2.hand, idents = c("FMO5+ AT2", "CFTR+ AT2", "CCL2+ Mesenchyme", "CCL2+ Mesenchyme IL6", "Alveolar Fibroblast"))
Idents(xenium.obj.2.hand.subset) <- "All_Hand_Annotation"
xenium.obj.2.hand.subset <- RenameIdents(xenium.obj.2.hand.subset, "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme")
xenium.obj.2.hand.subset <- AddMetaData(xenium.obj.2.hand.subset, col.name = "All_Hand_Annotation", Idents(xenium.obj.2.hand.subset))

Idents(xenium.obj.3.hand) <- "All_Hand_Annotation"
xenium.obj.3.hand.subset <- subset(xenium.obj.3.hand, idents = c("FMO5+ AT2", "CFTR+ AT2", "CCL2+ Mesenchyme", "CCL2+ Mesenchyme IL6", "Alveolar Fibroblast"))
Idents(xenium.obj.3.hand.subset) <- "All_Hand_Annotation"
xenium.obj.3.hand.subset <- RenameIdents(xenium.obj.3.hand.subset, "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme")
xenium.obj.3.hand.subset <- AddMetaData(xenium.obj.3.hand.subset, col.name = "All_Hand_Annotation", Idents(xenium.obj.3.hand.subset))

Idents(xenium.obj.4.hand) <- "All_Hand_Annotation"
xenium.obj.4.hand.subset <- subset(xenium.obj.4.hand, idents = c("FMO5+ AT2", "CFTR+ AT2", "CCL2+ Mesenchyme", "CCL2+ Mesenchyme IL6", "Alveolar Fibroblast"))
Idents(xenium.obj.4.hand.subset) <- "All_Hand_Annotation"
xenium.obj.4.hand.subset <- RenameIdents(xenium.obj.4.hand.subset, "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme")
xenium.obj.4.hand.subset <- AddMetaData(xenium.obj.4.hand.subset, col.name = "All_Hand_Annotation", Idents(xenium.obj.4.hand.subset))

Idents(xenium.obj.5.hand) <- "All_Hand_Annotation"
xenium.obj.5.hand.subset <- subset(xenium.obj.5.hand, idents = c("FMO5+ AT2", "CFTR+ AT2", "CCL2+ Mesenchyme", "CCL2+ Mesenchyme IL6", "Alveolar Fibroblast"))
Idents(xenium.obj.5.hand.subset) <- "All_Hand_Annotation"
xenium.obj.5.hand.subset <- RenameIdents(xenium.obj.5.hand.subset, "CCL2+ Mesenchyme IL6" = "CCL2+ Mesenchyme")
xenium.obj.5.hand.subset <- AddMetaData(xenium.obj.5.hand.subset, col.name = "All_Hand_Annotation", Idents(xenium.obj.5.hand.subset))
xenium.obj.1.hand.subset <- BuildNicheAssay(object = xenium.obj.1.hand.subset, , group.by = "All_Hand_Annotation", niches.k = 2, neighbors.k = 10, fov = "test")
xenium.obj.2.hand.subset <- BuildNicheAssay(object = xenium.obj.2.hand.subset, , group.by = "All_Hand_Annotation", niches.k = 2, neighbors.k = 10, fov = "test.2")
xenium.obj.3.hand.subset <- BuildNicheAssay(object = xenium.obj.3.hand.subset, , group.by = "All_Hand_Annotation", niches.k = 2, neighbors.k = 10, fov = "test.3")
xenium.obj.4.hand.subset <- BuildNicheAssay(object = xenium.obj.4.hand.subset, , group.by = "All_Hand_Annotation", niches.k = 2, neighbors.k = 10, fov = "test.4")
xenium.obj.5.hand.subset <- BuildNicheAssay(object = xenium.obj.5.hand.subset, , group.by = "All_Hand_Annotation", niches.k = 2, neighbors.k = 10, fov = "test.5")
write.csv(table(xenium.obj.1.hand.subset$All_Hand_Annotation, xenium.obj.1.hand.subset$niches), "obj_1NicheCompositionK=2.csv")
write.csv(table(xenium.obj.2.hand.subset$All_Hand_Annotation, xenium.obj.2.hand.subset$niches), "obj_2NicheCompositionK=2.csv")
write.csv(table(xenium.obj.3.hand.subset$All_Hand_Annotation, xenium.obj.3.hand.subset$niches), "obj_3NicheCompositionK=2.csv")
write.csv(table(xenium.obj.4.hand.subset$All_Hand_Annotation, xenium.obj.4.hand.subset$niches), "obj_4NicheCompositionK=2.csv")
write.csv(table(xenium.obj.5.hand.subset$All_Hand_Annotation, xenium.obj.5.hand.subset$niches), "obj_5NicheCompositionK=2.csv")

pdf(file.path("./", paste0("obj_1NichePlotK=2", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.1.hand.subset, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_2NichePlotK=2", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.2.hand.subset, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_3NichePlotK=2", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.3.hand.subset, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_4NichePlotK=2", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.4.hand.subset, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
pdf(file.path("./", paste0("obj_5NichePlotK=2", ".pdf")), w=11, h=8.5)
print(ImageDimPlot(xenium.obj.5.hand.subset, group.by = "niches", size = 0.75, dark.background = F))
dev.off()
seurat_objects <- list(xenium.obj.1.hand.subset, xenium.obj.2.hand.subset, xenium.obj.3.hand.subset, xenium.obj.4.hand.subset, xenium.obj.5.hand.subset)

# Function to extract, tabulate, and prefix column names
extract_and_tabulate <- function(seurat_obj, prefix) {
  # Extract metadata
  metadata <- seurat_obj@meta.data[, c("All_Hand_Annotation", "niches")]
  # Create contingency table
  tab <- table(metadata$All_Hand_Annotation, metadata$niche)
  # Convert to data frame
  df <- as.data.frame.matrix(tab)
  # Add 'All_Hand_Annotation' as a column
  df$All_Hand_Annotation <- rownames(df)
  # Reorder columns to have 'All_Hand_Annotation' first
  df <- df[, c("All_Hand_Annotation", setdiff(names(df), "All_Hand_Annotation"))]
  # Rename columns
  colnames(df)[-1] <- paste0(prefix, "_", colnames(df)[-1])
  return(df)
}

# Apply the function to each Seurat object with appropriate prefix
df_list <- map2(seurat_objects, 1:5, extract_and_tabulate)

# Merge all data frames by 'All_Hand_Annotation'
merged_df <- reduce(df_list, full_join, by = "All_Hand_Annotation")

# Replace NAs with 0
merged_df[is.na(merged_df)] <- 0

# View the merged table
print(merged_df)

write.csv(merged_df, "merged_Seurat_Niche_CompositionK=2.csv", row.names = FALSE)
rownames(merged_df) <- merged_df[,1]
merged_df <- merged_df[,-1]
Feature <- as.matrix(merged_df)
cor_matrix <- round(cor(Feature), digits = 3)
pdf(file.path("./", paste0("NicheCorrelationHeatmapK=2", ".pdf")), w=11, h=8.5)
p1 <- pheatmap(cor_matrix, display_numbers = TRUE)
p1
dev.off()

######Figure S3G######
write.csv(merged_df, "merged_Seurat_Niche_CompositionK=2.csv", row.names = FALSE)
#####Make Pie Charts in Prism##### ~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/Simple Test Case Merge IL6 Mes/2025-05-29/Simple Niche Accounting.xlsx
#~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Output/2025-01-28 Frum 2025/Figure 3/Xenium Label Transfer/Simple Test Case Merge IL6 Mes/2025-05-29/Simple niche Accounting (S3GH).prism

##Figure 5

identities <- c("AT2", "Macrophages", "Transitional", "Multiciliated",
                "Respiratory Bronchiole", "TREM2+ Macrophages", "Pericyte", "Lymphatic Endothelial",
                "AT1", "Monocytes", "CD1c+ DC", "Basal",
                "Capillary Endothelial", "Venous Endothelial", "Alveolar Fibroblast", "Neuroendocrine",
                "Arterial Endothelial", "B-Cells", "Alveolar Myofibroblast", "Basophils",
                "CCL2-Positive Fibroblast", "Neutrophils", "Capillary Aerocyte", "CD8+ T-Cell",
                "Fibroblast", "Adventitial Mesenchyme", "NK Cell", "CD4+ T-Cell",
                "pDC")
colors <- c(
  "#E41A1C", "#F4A3A8", # Intense Red, Soft Red
  "#377EB8", "#7FC97F", # Intense Blue, Soft Green-Blue
  "#4DAF4A", "#B2DF8A", # Intense Green, Soft Green
  "#984EA3", "#CAB2D6", # Intense Purple, Soft Lavender
  "#A65628", "#E7C5B4", # Intense Brown, Soft Brown
  "#F781BF", "#FCCDE5", # Intense Pink, Soft Pink
  "#66C2A5", "#C7EAE5", # Intense Teal, Soft Teal
  "#FF7F00", "#FDBF6F", # Intense Orange, Soft Orange
  "#6A3D9A", "#D4B9DA", # Intense Violet, Soft Violet
  "#A6D854", "#D9F0A3", # Intense Lime, Soft Lime
  "#1F78B4", "#A6CEE3", # Intense Medium Blue, Soft Medium Blue
  "#33A02C", "#A1DAB4", # Intense Forest Green, Soft Forest Green
  "#8E0152", "#C51B7D", # Intense Burgundy, Soft Rose
  "#B15928", "#E7BA8B", # Intense Rust, Soft Rust
  "#1B9E77", "#B2E2E2", # Intense Dark Teal, Soft Dark Teal
  "#FF69B4", "#FFC0CB", # Intense Hot Pink, Soft Blush Pink
  "#8C510A", "#C19A6B"  # Intense Earthy Brown, Soft Sand Brown
)
color_map <- setNames(colors, identities)


#Fig 5A

###RPCA Process Pediatric for Label Transfering
pediatric <- readRDS("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/RObjects/Pediatric Cell Atlas/FullPDAObjectFixedforPaperRound1.RDS")
PDA.RPCA <- pediatric
PDA.RPCA <- SplitObject(pediatric, split.by = "ID")


i = 1
for (i in seq_along(PDA.RPCA)) {
  PDA.RPCA[[i]] <- NormalizeData(PDA.RPCA[[i]]) %>% FindVariableFeatures() %>% CellCycleScoring(g2m.features = g2m.genes, s.features = s.genes)
}
features <- SelectIntegrationFeatures(PDA.RPCA, nfeatures = 3000)
for (i in seq_along(along.with = PDA.RPCA)) {
  PDA.RPCA[[i]] <- ScaleData(PDA.RPCA[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(PDA.RPCA, anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 3)
PDA.RPCA.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(PDA.RPCA.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

PDA.RPCA.integrated <- ScaleData(PDA.RPCA.integrated)
PDA.RPCA.integrated <- RunPCA(PDA.RPCA.integrated)
ElbowPlot(PDA.RPCA.integrated, ndims = 50)
PDA.RPCA.integrated <- RunUMAP(PDA.RPCA.integrated, dims = 1:30, reduction.name = "umap", return.model = TRUE)
PDA.RPCA.integrated <- FindNeighbors(PDA.RPCA.integrated, dims = 1:30)
PDA.RPCA.integrated <- FindClusters(PDA.RPCA.integrated, resolution = 1.0)
pdf(file.path("./", paste0("Full Data Pediatric RPCA Louvain", ".pdf")), w=11, h=8.5)
DimPlot(PDA.RPCA.integrated, label = TRUE, pt.size = 0.5)
dev.off() 

pediatric <- PDA.RPCA.integrated
# Ensure both datasets have been preprocessed (normalized, variable features, PCA)


Idents(pediatric) <- "annotation_lvl3"
pediatric <- RenameIdents(pediatric, "CFTR+ AT2" = "AT2", "FMO5+ AT2" = "AT2", "Transitional 1" = "Transitional", "Transitional 2" = "Transitional", "AT1 1" = "AT1", "AT1 2" = "AT1", "AT1 3" = "AT1", "Terminal Bronchiole Cell (RASC or TRB)" = "Respiratory Bronchiole", "Respiratory Bronchiole Cell" = "Respiratory Bronchiole", "Goblet-Like Secretory" = "Respiratory Bronchiole", "Alveolar Fibroblast 1" = "Alveolar Fibroblast", "Alveolar Fibroblast 2" = "Alveolar Fibroblast", "Activated Macrophages" = "TREM2+ Macrophages", "Adventitial-Like Mesenchyme" = "Adventitial Mesenchyme") ##merge subclusters to single annotation
pediatric <- AddMetaData(pediatric, col.name = "full.dotplot.annotation", Idents(pediatric))



##load and QC filter PIG samples

##RPCA Integrate PIG
DefaultAssay(PIG.fastmnn.cleaned) <- "RNA"
PIG.fastmnn.cleaned <- SplitObject(PIG.fastmnn.cleaned, split.by = "ID")
i = 1
for (i in seq_along(PIG.fastmnn.cleaned)) {
  PIG.fastmnn.cleaned[[i]] <- NormalizeData(PIG.fastmnn.cleaned[[i]]) %>% FindVariableFeatures() %>% CellCycleScoring(g2m.features = g2m.genes, s.features = s.genes)
}
features <- SelectIntegrationFeatures(PIG.fastmnn.cleaned, nfeatures = 3000)
for (i in seq_along(along.with = PIG.fastmnn.cleaned)) {
  PIG.fastmnn.cleaned[[i]] <- ScaleData(PIG.fastmnn.cleaned[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(PIG.fastmnn.cleaned, anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 3)
PIG.fastmnn.cleaned.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(PIG.fastmnn.cleaned.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

PIG.fastmnn.cleaned.integrated <- ScaleData(PIG.fastmnn.cleaned.integrated)
PIG.fastmnn.cleaned.integrated <- RunPCA(PIG.fastmnn.cleaned.integrated)
ElbowPlot(PIG.fastmnn.cleaned.integrated, ndims = 50)
PIG.fastmnn.cleaned.integrated <- RunUMAP(PIG.fastmnn.cleaned.integrated, dims = 1:30, reduction.name = "umap", return.model = TRUE)
PIG.fastmnn.cleaned.integrated <- FindNeighbors(PIG.fastmnn.cleaned.integrated, dims = 1:30)
PIG.fastmnn.cleaned.integrated <- FindClusters(PIG.fastmnn.cleaned.integrated, resolution = 0.8)

#Transfer labels from pediatric
anchors <- FindTransferAnchors(
  reference = pediatric,
  query = PIG.fastmnn.cleaned.integrated,
  dims = 1:30,
  reference.reduction = "pca",
  k.filter = 200, 
  reference.assay = "integrated",
  query.assay = "integrated"
)

# Transfer Labels
PIG.transferred <- TransferData(
  anchorset = anchors,
  refdata = pediatric$full.dotplot.annotation,
  dims = 1:30
)
PIG.fastmnn.cleaned.integrated$predicted_lvl3 <- PIG.transferred$predicted.id
DimPlot(PIG.fastmnn.cleaned.integrated, label = TRUE, group.by = "predicted_lvl3", pt.size = 0.2) + NoLegend()
#Fig 5A
pdf(file.path("./", paste0("PIG Lvl 3 Pediatric Label Transfer No Label", ".pdf")), w=11, h=8.5)
DimPlot(PIG.fastmnn.cleaned.integrated, label = FALSE, group.by = "predicted_lvl3", pt.size = 0.5, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()

##RPCA Integrate BPD
BPD.fastmnn.cleaned <- SplitObject(BPD.fastmnn.cleaned, split.by = "ID")


i = 1
for (i in seq_along(BPD.fastmnn.cleaned)) {
  BPD.fastmnn.cleaned[[i]] <- NormalizeData(BPD.fastmnn.cleaned[[i]]) %>% FindVariableFeatures() %>% CellCycleScoring(g2m.features = g2m.genes, s.features = s.genes)
}
features <- SelectIntegrationFeatures(BPD.fastmnn.cleaned, nfeatures = 3000)
for (i in seq_along(along.with = BPD.fastmnn.cleaned)) {
  BPD.fastmnn.cleaned[[i]] <- ScaleData(BPD.fastmnn.cleaned[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(BPD.fastmnn.cleaned, anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 3)
BPD.fastmnn.cleaned.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(BPD.fastmnn.cleaned.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

BPD.fastmnn.cleaned.integrated <- ScaleData(BPD.fastmnn.cleaned.integrated)
BPD.fastmnn.cleaned.integrated <- RunPCA(BPD.fastmnn.cleaned.integrated)
ElbowPlot(BPD.fastmnn.cleaned.integrated, ndims = 50)
BPD.fastmnn.cleaned.integrated <- RunUMAP(BPD.fastmnn.cleaned.integrated, dims = 1:30, reduction.name = "umap", return.model = TRUE)
BPD.fastmnn.cleaned.integrated <- FindNeighbors(BPD.fastmnn.cleaned.integrated, dims = 1:30)
BPD.fastmnn.cleaned.integrated <- FindClusters(BPD.fastmnn.cleaned.integrated, resolution = 0.8)

#Transfer labels from Pediatric onto BPD
anchors <- FindTransferAnchors(
  reference = pediatric,
  query = BPD.fastmnn.cleaned.integrated,
  dims = 1:30,
  reference.reduction = "pca",
  k.filter = 200, 
  reference.assay = "integrated",
  query.assay = "integrated"
)

# Transfer Labels
BPD.transferred <- TransferData(
  anchorset = anchors,
  refdata = pediatric$full.dotplot.annotation,
  dims = 1:30
)
BPD.fastmnn.cleaned.integrated$predicted_lvl3 <- BPD.transferred$predicted.id
pdf(file.path("./", paste0("BPD RPCA Lvl 3 Pediatric Label Transfer", ".pdf")), w=11, h=8.5)
DimPlot(BPD.fastmnn.cleaned.integrated, label = TRUE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE) + NoLegend() + NoAxes()
dev.off()

#Fig 5B
pdf(file.path("./", paste0("BPD RPCA Lvl 3 Pediatric Label Transfer by Condition No Label", ".pdf")), w=22, h=8.5)
DimPlot(BPD.fastmnn.cleaned.integrated, label = FALSE, group.by = "predicted_lvl3", split.by = "disease", pt.size = 0.5, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("BPD Lvl 3 Pediatric Label Transfer No Label", ".pdf")), w=11, h=8.5)
DimPlot(BPD.fastmnn.cleaned.integrated, label = FALSE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()



#Fig 5SA
library(tidyr)
group_colors <- c(
  "Pediatric"        = "#66C2A5",  # will show as “Non-BPD/PIG” in your legend via labels
  "BPD - Organizing" = "#FC8D62",
  "BPD - Severe"     = "#E78AC3",
  "PIG"              = "#8DA0CB"
)
pediatric <- AddMetaData(pediatric, col.name = "predicted_lvl3", pediatric$full.dotplot.annotation)
desired_order <- c("AT2", "Transitional", "AT1", "Respiratory Bronchiole", "Basal", 
                   "Multiciliated", "Neuroendocrine", "Macrophages", "TREM2+ Macrophages", 
                   "Monocytes", "CD1c+ DC", "pDC", "NK Cell", "CD4+ T-Cell", "CD8+ T-Cell", 
                    "B-Cells", "Neutrophils", "Lymphatic Endothelial", "Arterial Endothelial",
                   "Venous Endothelial", "Capillary Endothelial", "Capillary Aerocyte", 
                   "Alveolar Fibroblast", "CCL2-Positive Fibroblast", "Fibroblast", 
                    "Alveolar Myofibroblast", "Adventitial Mesenchyme", 
                   "Pericyte")

# Initialize collection list
individual_level_list <- list()

# Define proper group levels
original_group_levels <- c("Pediatric", "BPD - Organizing", "BPD - Severe", "PIG")
updated_group_names <- c("Non-BPD/PIG", "BPD - Organizing", "BPD - Severe", "PIG")
names(updated_group_names) <- original_group_levels

# Collect sample counts for group labeling later
group_sample_counts <- list()

seurat_list <- list(
  pediatric = pediatric, 
  BPD_fastmnn_cleaned_integrated = BPD.fastmnn.cleaned.integrated, 
  PIG_fastmnn_cleaned_integrated = PIG.fastmnn.cleaned.integrated
)
# Loop through Seurat objects
for (dataset_name in names(seurat_list)) {
  seurat_obj <- seurat_list[[dataset_name]]
  metadata <- seurat_obj@meta.data
  
  if (!inherits(seurat_obj, "Seurat")) next
  if (!all(c("predicted_lvl3", "ID") %in% colnames(metadata))) {
    warning(paste("Skipping", dataset_name, "- missing required columns"))
    next
  }
  
  cat("Processing:", dataset_name, "\n")
  
  # Initial % per identity per ID
  df <- metadata %>%
    dplyr::select(predicted_lvl3, ID) %>%
    group_by(ID, predicted_lvl3) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ID) %>%
    mutate(percentage = n / sum(n) * 100,
           dataset = dataset_name) %>%
    ungroup()
  
  # ---- GROUP assignment ----
  if (grepl("pediatric", dataset_name, ignore.case = TRUE)) {
    df <- df %>% mutate(group = "Pediatric")
  } else if (grepl("PIG", dataset_name)) {
    df <- df %>% mutate(group = "PIG")
  } else if (grepl("BPD", dataset_name)) {
    if (!"disease" %in% colnames(metadata)) {
      stop("BPD dataset is missing 'disease' column.")
    }
    disease_info <- metadata %>% dplyr::select(ID, disease) %>% distinct()
    cat("Disease labels in BPD:", paste(unique(disease_info$disease), collapse = ", "), "\n")
    df <- df %>%
      left_join(disease_info, by = "ID") %>%
      mutate(group = case_when(
        disease == "organizing" ~ "BPD - Organizing",
        disease == "severe" ~ "BPD - Severe",
        TRUE ~ NA_character_
      ))
  }
  
  if (!"group" %in% colnames(df)) {
    warning("No group assigned in dataset:", dataset_name)
    next
  }
  
  df <- df %>% filter(!is.na(group))
  group_sample_counts[[dataset_name]] <- table(unique(df[, c("ID", "group")])$group)
  individual_level_list[[dataset_name]] <- df
}

# Combine all individuals
individual_level_df <- bind_rows(individual_level_list)
individual_level_df$group <- factor(individual_level_df$group, levels = original_group_levels)
individual_level_df$predicted_lvl3 <- factor(individual_level_df$predicted_lvl3, levels = desired_order)

# Save per-individual data
write.csv(individual_level_df, "individual_identity_composition.csv", row.names = FALSE)

# Generate updated group labels with n counts
group_ns <- individual_level_df %>% 
  dplyr::select(ID, group) %>% 
  distinct() %>% 
  group_by(group) %>% 
  summarise(n = n(), .groups = "drop")
group_levels_labeled <- paste0(updated_group_names[as.character(group_ns$group)], " (n=", group_ns$n, ")")
names(group_levels_labeled) <- as.character(group_ns$group)

# ---- Summary stats and stats testing ----
grouped_summary <- individual_level_df %>%
  group_by(predicted_lvl3, group) %>%
  summarise(mean_pct = mean(percentage),
            sd_pct = sd(percentage), .groups = "drop")

# Pairwise stats: Pediatric vs each group
pairwise_results <- individual_level_df %>%
  group_by(predicted_lvl3) %>%
  summarise(
    mean_Pediatric = mean(percentage[group == "Pediatric"]),
    sd_Pediatric = sd(percentage[group == "Pediatric"]),
    mean_BPD_Org = mean(percentage[group == "BPD - Organizing"]),
    sd_BPD_Org = sd(percentage[group == "BPD - Organizing"]),
    p_BPD_Org = tryCatch(wilcox.test(percentage[group == "Pediatric"], percentage[group == "BPD - Organizing"])$p.value, error = function(e) NA),
    mean_BPD_Sev = mean(percentage[group == "BPD - Severe"]),
    sd_BPD_Sev = sd(percentage[group == "BPD - Severe"]),
    p_BPD_Sev = tryCatch(wilcox.test(percentage[group == "Pediatric"], percentage[group == "BPD - Severe"])$p.value, error = function(e) NA),
    mean_PIG = mean(percentage[group == "PIG"]),
    sd_PIG = sd(percentage[group == "PIG"]),
    p_PIG = tryCatch(wilcox.test(percentage[group == "Pediatric"], percentage[group == "PIG"])$p.value, error = function(e) NA),
    .groups = "drop"
  )
write.csv(pairwise_results, "group_composition_stats.csv", row.names = FALSE)

# ---- Plot with significance ----
positioner <- position_dodge(width = 0.9)

p <- ggplot(grouped_summary, aes(x = predicted_lvl3, y = mean_pct, fill = group)) +
  geom_bar(stat = "identity", position = positioner, color = "black") +
  geom_errorbar(
    data = grouped_summary,
    aes(x = predicted_lvl3, ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct, group = group),
    position = positioner, width = 0.4, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_colors, labels = group_levels_labeled) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Cell Type Composition Across Groups",
       x = "Cell Type", y = "Mean % Composition") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "top"
  )

# Use raw p-values and convert to * labels
annotation_df <- pairwise_results %>%
  pivot_longer(cols = c(p_BPD_Org, p_BPD_Sev, p_PIG),
               names_to = "comparison", values_to = "p_val") %>%
  mutate(
    group = case_when(
      comparison == "p_BPD_Org" ~ "BPD - Organizing",
      comparison == "p_BPD_Sev" ~ "BPD - Severe",
      comparison == "p_PIG" ~ "PIG"
    ),
    label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(label != "") %>%
  left_join(grouped_summary, by = c("predicted_lvl3", "group")) %>%
  mutate(
    y_position = mean_pct + sd_pct + 5,
    color = group_colors[group],
    x_nudge = as.numeric(predicted_lvl3) +
      (as.numeric(factor(group, levels = original_group_levels)) - 2.5) * 0.32 +
      ifelse(group == "BPD - Organizing", 0.22, 0)  # Shift orange *
  )

p <- p +
  geom_text(data = annotation_df,
            aes(x = x_nudge, y = y_position, label = label, color = group),
            size = 9, fontface = "bold", angle = 90, vjust = 0.5, inherit.aes = FALSE) +
  scale_color_manual(values = group_colors, guide = "none")

# Save plot
ggsave("grouped_interleaved_composition_plot_with_significance.pdf", p, width = 14, height = 7)

######Figure S5B######

identities <- c("AT2", "Macrophages", "Transitional", "Multiciliated",
                "Respiratory Bronchiole", "TREM2+ Macrophages", "Pericyte", "Lymphatic Endothelial",
                "AT1", "Monocytes", "CD1c+ DC", "Basal",
                "Capillary Endothelial", "Venous Endothelial", "Alveolar Fibroblast", "Neuroendocrine",
                "Arterial Endothelial", "B-Cells", "Alveolar Myofibroblast", "Basophils",
                "CCL2-Positive Fibroblast", "Neutrophils", "Capillary Aerocyte", "CD8+ T-Cell",
                "Fibroblast", "Adventitial Mesenchyme", "NK Cell", "CD4+ T-Cell",
                "pDC")
colors <- c(
  "#E41A1C", "#F4A3A8", # Intense Red, Soft Red
  "#377EB8", "#7FC97F", # Intense Blue, Soft Green-Blue
  "#4DAF4A", "#B2DF8A", # Intense Green, Soft Green
  "#984EA3", "#CAB2D6", # Intense Purple, Soft Lavender
  "#A65628", "#E7C5B4", # Intense Brown, Soft Brown
  "#F781BF", "#FCCDE5", # Intense Pink, Soft Pink
  "#66C2A5", "#C7EAE5", # Intense Teal, Soft Teal
  "#FF7F00", "#FDBF6F", # Intense Orange, Soft Orange
  "#6A3D9A", "#D4B9DA", # Intense Violet, Soft Violet
  "#A6D854", "#D9F0A3", # Intense Lime, Soft Lime
  "#1F78B4", "#A6CEE3", # Intense Medium Blue, Soft Medium Blue
  "#33A02C", "#A1DAB4", # Intense Forest Green, Soft Forest Green
  "#8E0152", "#C51B7D", # Intense Burgundy, Soft Rose
  "#B15928", "#E7BA8B", # Intense Rust, Soft Rust
  "#1B9E77", "#B2E2E2", # Intense Dark Teal, Soft Dark Teal
  "#FF69B4", "#FFC0CB", # Intense Hot Pink, Soft Blush Pink
  "#8C510A", "#C19A6B"  # Intense Earthy Brown, Soft Sand Brown
)
color_map <- setNames(colors, identities)
pdf(file.path("./", paste0("Fetal 16 week Distal and Small RPCA Lvl 3 Pediatric Label Transfer No Label", ".pdf")), w=11, h=8.5)
DimPlot(fetal.distal, label = FALSE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("16, 21, 23 Pediatric Label Transfer No Label", ".pdf")), w=11, h=8.5)
DimPlot(intermediate.timepoints.integrated, label = FALSE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("Pediatric FastMNN Lvl 3 Pediatric Label Transfer No Label", ".pdf")), w=11, h=8.5)
DimPlot(pediatric.og, label = FALSE, group.by = "full.dotplot.annotation", pt.size = 0.2, repel = TRUE,  cols = color_map) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("Fetal 16 week Distal and Small RPCA Lvl 3 Pediatric Label Transfer  Label", ".pdf")), w=11, h=8.5)
DimPlot(fetal.distal, label = TRUE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("16, 21, 23 Pediatric Label Transfer Label", ".pdf")), w=11, h=8.5)
DimPlot(intermediate.timepoints.integrated, label = TRUE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE, cols = color_map) + NoLegend() + NoAxes()
dev.off()
pdf(file.path("./", paste0("Pediatric FastMNN Lvl 3 Pediatric Label Transfer  Label", ".pdf")), w=11, h=8.5)
DimPlot(pediatric.og, label = TRUE, group.by = "full.dotplot.annotation", pt.size = 0.2, repel = TRUE,  cols = color_map) + NoLegend() + NoAxes()
dev.off()
seurat_list <- list(
  pediatric = pediatric, 
  fetal_distal = fetal.distal, 
  intermediate_timepoints_integrated = intermediate.timepoints.integrated
  
)

# Ensure all elements in the list are valid Seurat objects
if (!all(sapply(seurat_list, inherits, "Seurat"))) {
  stop("Error: One or more elements in seurat_list are not valid Seurat objects.")
}


# Function to calculate cell percentages and plot bar chart
color_mapping <- setNames(colors, identities)

# Extract the ordered levels of "predicted_lvl3" from pediatric
ordered_levels <- unique(as.character(pediatric@meta.data$predicted_lvl3))  # Ensure it's a vector

# Function to calculate cell percentages and plot bar chart


summary_list <- list()
plot_identity_distribution <- function(seurat_obj, dataset_name) {
  if (!inherits(seurat_obj, "Seurat")) {
    warning(paste("Skipping", dataset_name, ": Not a valid Seurat object"))
    return(NULL)
  }
  
  cat("Processing dataset:", dataset_name, "\n")  # Debugging print
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Check if the metadata column exists
  if (!"predicted_lvl3" %in% colnames(metadata)) {
    warning(paste("Skipping", dataset_name, ": 'predicted_lvl3' not found in metadata"))
    return(NULL)
  }
  
  # Convert predicted_lvl3 to a character vector
  metadata$predicted_lvl3 <- as.character(metadata$predicted_lvl3)
  
  # Count occurrences of each identity using table()
  identity_counts <- as.data.frame(table(metadata$predicted_lvl3))
  colnames(identity_counts) <- c("predicted_lvl3", "n")
  
  # Convert counts to percentages
  identity_counts <- identity_counts %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # Ensure all identities are represented
  identity_counts <- identity_counts %>%
    complete(predicted_lvl3 = ordered_levels, fill = list(n = 0, percentage = 0)) 
  
  # Ensure label exists before using it in geom_text()
  identity_counts <- identity_counts %>%
    mutate(label = ifelse(percentage == 0, "Not Detected", ""))
  
  # Define the desired column order
  desired_order <- c("AT2", "Transitional", "AT1", "Respiratory Bronchiole", "Basal", 
                     "Multiciliated", "Neuroendocrine", "Macrophages", "TREM2+ Macrophages", 
                     "Monocytes", "CD1c+ DC", "pDC", "NK Cell", "CD4+ T-Cell", "CD8+ T-Cell", 
                     "B-Cells", "Neutrophils", "Lymphatic Endothelial", "Arterial Endothelial",
                     "Venous Endothelial", "Capillary Endothelial", "Capillary Aerocyte", 
                     "Alveolar Fibroblast", "CCL2-Positive Fibroblast", "Fibroblast", 
                     "Alveolar Myofibroblast", "Adventitial Mesenchyme", 
                     "Pericyte")
  
  # Ensure predicted_lvl3 follows the desired order
  identity_counts$predicted_lvl3 <- factor(identity_counts$predicted_lvl3, levels = desired_order)
  
  # Export to CSV
  write.csv(identity_counts, file = paste0(dataset_name, "_identity_distribution.csv"), row.names = FALSE)
  
  # Create bar plot
  p <- ggplot(identity_counts, aes(x = predicted_lvl3, y = percentage, fill = predicted_lvl3)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = label, 
                  y = ifelse(label == "Not Detected", 10, percentage + 2)),  # Move "Not Detected" to y = 10
              angle = 90,
              size = 6,
              fontface = "bold") +
    scale_fill_manual(values = color_mapping, na.value = "gray50") +
    labs(title = paste("Cell Type Distribution in", dataset_name), 
         x = "Cell Identity", 
         y = "Percentage of Cells") +
    ylim(0, 50) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(20, 20, 20, 50, "pt"),
      axis.line.y = element_line(color = "black", size = 1.2),
      axis.line.x = element_line(color = "black", size = 1.2)
    )
  # Save plot as PDF
  ggsave(filename = paste0(dataset_name, "_identity_distribution.pdf"), plot = p, width = 11, height = 8.5)
  return(identity_counts %>% mutate(dataset = dataset_name))
  
}

for (dataset_name in names(seurat_list)) {
  identity_df <- plot_identity_distribution(seurat_list[[dataset_name]], dataset_name)
  if (!is.null(identity_df)) {
    summary_list[[dataset_name]] <- identity_df
  }
}

# Combine all datasets into a single summary table
summary_table <- bind_rows(summary_list)

# Reorder columns
summary_table <- summary_table %>%
  select(dataset, predicted_lvl3, n, percentage)

# Save the master summary CSV
write.csv(summary_table, "combined_identity_distribution_summary.csv", row.names = FALSE)



#Figure 5F
###Merge all ITGA8 positive fibroblast populations together#####
Idents(pediatric) <- "predicted_lvl3"
pediatric <- RenameIdents(pediatric, "CCL2-Positive Fibroblast" = "Alveolar Fibroblast", "Fibroblast" = "Alveolar Fibroblast")
pediatric <- AddMetaData(pediatric, col.name = "predicted_lvl3", Idents(pediatric))
Idents(PIG.fastmnn.cleaned.integrated) <- "predicted_lvl3"
PIG.fastmnn.cleaned.integrated <- RenameIdents(PIG.fastmnn.cleaned.integrated, "CCL2-Positive Fibroblast" = "Alveolar Fibroblast", "Fibroblast" = "Alveolar Fibroblast")
PIG.fastmnn.cleaned.integrated <- AddMetaData(PIG.fastmnn.cleaned.integrated, col.name = "predicted_lvl3", Idents(PIG.fastmnn.cleaned.integrated))

Idents(BPD.fastmnn.cleaned.integrated) <- "predicted_lvl3"

BPD.fastmnn.cleaned.integrated <- RenameIdents(BPD.fastmnn.cleaned.integrated, "CCL2-Positive Fibroblast" = "Alveolar Fibroblast",  "Fibroblast" = "Alveolar Fibroblast")
BPD.fastmnn.cleaned.integrated <- AddMetaData(BPD.fastmnn.cleaned.integrated, col.name = "predicted_lvl3", Idents(BPD.fastmnn.cleaned.integrated))

individual_level_list <- list()
####RERUN Stats####
# Define proper group levels
original_group_levels <- c("Pediatric", "BPD - Organizing", "BPD - Severe", "PIG")
updated_group_names <- c("Non-BPD/PIG", "BPD - Organizing", "BPD - Severe", "PIG")
names(updated_group_names) <- original_group_levels

# Collect sample counts for group labeling later
group_sample_counts <- list()

seurat_list <- list(
  pediatric = pediatric, 
  BPD_fastmnn_cleaned_integrated = BPD.fastmnn.cleaned.integrated, 
  PIG_fastmnn_cleaned_integrated = PIG.fastmnn.cleaned.integrated
)
# Loop through Seurat objects
for (dataset_name in names(seurat_list)) {
  seurat_obj <- seurat_list[[dataset_name]]
  metadata <- seurat_obj@meta.data
  
  if (!inherits(seurat_obj, "Seurat")) next
  if (!all(c("predicted_lvl3", "ID") %in% colnames(metadata))) {
    warning(paste("Skipping", dataset_name, "- missing required columns"))
    next
  }
  
  cat("Processing:", dataset_name, "\n")
  
  # Initial % per identity per ID
  df <- metadata %>%
    dplyr::select(predicted_lvl3, ID) %>%
    group_by(ID, predicted_lvl3) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ID) %>%
    mutate(percentage = n / sum(n) * 100,
           dataset = dataset_name) %>%
    ungroup()
  
  # ---- GROUP assignment ----
  if (grepl("pediatric", dataset_name, ignore.case = TRUE)) {
    df <- df %>% mutate(group = "Pediatric")
  } else if (grepl("PIG", dataset_name)) {
    df <- df %>% mutate(group = "PIG")
  } else if (grepl("BPD", dataset_name)) {
    if (!"disease" %in% colnames(metadata)) {
      stop("BPD dataset is missing 'disease' column.")
    }
    disease_info <- metadata %>% dplyr::select(ID, disease) %>% distinct()
    cat("Disease labels in BPD:", paste(unique(disease_info$disease), collapse = ", "), "\n")
    df <- df %>%
      left_join(disease_info, by = "ID") %>%
      mutate(group = case_when(
        disease == "organizing" ~ "BPD - Organizing",
        disease == "severe" ~ "BPD - Severe",
        TRUE ~ NA_character_
      ))
  }
  
  if (!"group" %in% colnames(df)) {
    warning("No group assigned in dataset:", dataset_name)
    next
  }
  
  df <- df %>% filter(!is.na(group))
  group_sample_counts[[dataset_name]] <- table(unique(df[, c("ID", "group")])$group)
  individual_level_list[[dataset_name]] <- df
}

# Combine all individuals
individual_level_df <- bind_rows(individual_level_list)
individual_level_df$group <- factor(individual_level_df$group, levels = original_group_levels)
individual_level_df$predicted_lvl3 <- factor(individual_level_df$predicted_lvl3, levels = desired_order)

# Save per-individual data
write.csv(individual_level_df, "individual_identity_composition.csv", row.names = FALSE)

# Generate updated group labels with n counts
group_ns <- individual_level_df %>% 
  dplyr::select(ID, group) %>% 
  distinct() %>% 
  group_by(group) %>% 
  summarise(n = n(), .groups = "drop")
group_levels_labeled <- paste0(updated_group_names[as.character(group_ns$group)], " (n=", group_ns$n, ")")
names(group_levels_labeled) <- as.character(group_ns$group)

grouped_summary <- individual_level_df %>%
  group_by(predicted_lvl3, group) %>%
  summarise(mean_pct = mean(percentage),
            sd_pct = sd(percentage), .groups = "drop")

# Pairwise stats: Pediatric vs each group
pairwise_results <- individual_level_df %>%
  group_by(predicted_lvl3) %>%
  summarise(
    mean_Pediatric = mean(percentage[group == "Pediatric"]),
    sd_Pediatric = sd(percentage[group == "Pediatric"]),
    mean_BPD_Org = mean(percentage[group == "BPD - Organizing"]),
    sd_BPD_Org = sd(percentage[group == "BPD - Organizing"]),
    p_BPD_Org = tryCatch(wilcox.test(percentage[group == "Pediatric"], percentage[group == "BPD - Organizing"])$p.value, error = function(e) NA),
    mean_BPD_Sev = mean(percentage[group == "BPD - Severe"]),
    sd_BPD_Sev = sd(percentage[group == "BPD - Severe"]),
    p_BPD_Sev = tryCatch(wilcox.test(percentage[group == "Pediatric"], percentage[group == "BPD - Severe"])$p.value, error = function(e) NA),
    mean_PIG = mean(percentage[group == "PIG"]),
    sd_PIG = sd(percentage[group == "PIG"]),
    p_PIG = tryCatch(wilcox.test(percentage[group == "Pediatric"], percentage[group == "PIG"])$p.value, error = function(e) NA),
    .groups = "drop"
  )
write.csv(pairwise_results, "group_composition_stats.csv", row.names = FALSE)





#####Maturation Scoring########
##Align Fetal Distal to Pediatric
Idents(fetal) <- "group"
fetal.distal <- subset(fetal, idents = c("distal", "small"))


DefaultAssay(fetal.distal) <- "integrated"
fetal.distal <- RunPCA(fetal.distal, verbose = FALSE)
ElbowPlot(fetal.distal, ndims = 50)
fetal.distal <- RunUMAP(fetal.distal, reduction = "pca", dims = 1:15, return.model = TRUE)
fetal.distal <- FindNeighbors(fetal.distal, dims = 1:15)
fetal.distal <- FindClusters(fetal.distal, resolution = 0.4)
pdf(file.path("./", paste0("Fetal Distal Louvain", ".pdf")), w=11, h=8.5)
DimPlot(fetal.distal, label = TRUE, pt.size = 0.5) + NoLegend() + NoAxes()
dev.off()

# Ensure both datasets have been preprocessed (normalized, variable features, PCA)

DefaultAssay(fetal.distal) <- "integrated"
DefaultAssay(pediatric) <- "integrated"
Idents(pediatric) <- "annotation_lvl3"
pediatric <- RenameIdents(pediatric, "CFTR+ AT2" = "AT2", "FMO5+ AT2" = "AT2", "Transitional 1" = "Transitional", "Transitional 2" = "Transitional", "AT1 1" = "AT1", "AT1 2" = "AT1", "AT1 3" = "AT1", "Terminal Bronchiole Cell (RASC or TRB)" = "Respiratory Bronchiole", "Respiratory Bronchiole Cell" = "Respiratory Bronchiole", "Goblet-Like Secretory" = "Respiratory Bronchiole", "Alveolar Fibroblast 1" = "Alveolar Fibroblast", "Alveolar Fibroblast 2" = "Alveolar Fibroblast")
pediatric <- AddMetaData(pediatric, col.name = "full.dotplot.annotation", Idents(pediatric))

anchors <- FindTransferAnchors(
  reference = pediatric,
  query = fetal.distal,
  dims = 1:30,
  reference.reduction = "pca",
  k.filter = 200, 
  reference.assay = "integrated",
  query.assay = "integrated"
)

# Transfer Labels
fetal.distal.transferred <- TransferData(
  anchorset = anchors,
  refdata = pediatric$full.dotplot.annotation,
  dims = 1:30
)

# Add predicted annotations to metadata
fetal.distal$predicted_lvl3 <- fetal.distal.transferred$predicted.id
pdf(file.path("./", paste0("Fetal 16 week Distal and Small RPCA Lvl 3 Pediatric Label Transfer", ".pdf")), w=11, h=8.5)
DimPlot(fetal.distal, label = TRUE, group.by = "predicted_lvl3", pt.size = 0.2, repel = TRUE) + NoLegend() + NoAxes()
dev.off()


######Determine Maturation Gene Sets######
# Iterate over shared cell types
cell_types <- intersect(unique(pediatric$full.dotplot.annotation ), unique(fetal.distal$predicted_lvl3))

maturation_genes_on <- list()
maturation_genes_off <- list()

for (cell_type in cell_types) {
  fetal.distal_subset <- subset(fetal.distal, predicted_lvl3 == cell_type)
  pediatric_subset <- subset(pediatric, full.dotplot.annotation == cell_type)
  fetal.distal_subset <- AddMetaData(fetal.distal_subset, col.name = "stage", "fetal.distal")
  pediatric_subset <- AddMetaData(pediatric_subset, col.name = "stage", "pediatric")
  
  # Merge fetal.distal and pediatric subsets for direct comparison
  combined <- merge(fetal.distal_subset, pediatric_subset)
  DefaultAssay(combined) <- "RNA"
  combined <- NormalizeData(combined, verbose = FALSE)
  # Create a metadata column to differentiate groups
  Idents(combined) <- "stage"
  # Find markers that are higher in pediatric vs. fetal.distal
  markers <- FindMarkers(
    combined, ident.1 = "pediatric", ident.2 = "fetal.distal",
    logfc.threshold = log(1.5), min.pct = 0.25, test.use = "wilcox"
  )
  top_on <- markers %>%
    filter(avg_log2FC > log2(1.5)) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 100) %>%
    rownames()
  
  # Filter top 100 genes that turn OFF (logFC < 0)
  top_off <- markers %>%
    filter(avg_log2FC < -log2(1.5)) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = 100) %>%
    rownames()
  
  # Store in the maturation gene lists
  maturation_genes_on[[cell_type]] <- top_on
  maturation_genes_off[[cell_type]] <- top_off
  
}

# Save results
write.csv(maturation_genes_on, "maturation_genes_on.csv") # doesn't work but gene list appears OK
write.csv(maturation_genes_off, "maturation_genes_off.csv") # doesn't work but gene list appears OK


#####Load and Process Independent 16 week, 21 and 23 week cells#####
counts.pd16 <- ReadCB_h5("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Multiome/Pediatric Cell Atlas/Sample_6730-TF-1/ Cellbender v03/output_filtered.h5")
pd16 <- CreateSeuratObject(
  counts = counts.pd16$`Gene Expression`, project = "PD16",
  assay = "RNA", min.cells = 3, min.features = 200
)
pd16 <- AddMetaData(pd16, "distal", col.name = "group")
pd16 <- AddMetaData(pd16, "PD16", col.name = "batch")
pd16 <- RenameCells(pd16, add.cell.id = "PD16")
pd16[["percent.mt"]] <- PercentageFeatureSet(pd16, pattern = "^MT-")
pd16[["percent.ribo"]] <- PercentageFeatureSet(pd16, pattern = "^RP[SL][[:digit:]]")
pd16 <- AddMetaData(pd16, "112", "age")
pd16 <- AddMetaData(pd16, col.name = "ID", "16 Week")

counts.pd23 <- ReadCB_h5("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Multiome/Pediatric Cell Atlas/Sample_6730-TF-2/Cellbender v03/output_filtered.h5")
pd23 <- CreateSeuratObject(
  counts = counts.pd23$`Gene Expression`, project = "pd23",
  assay = "RNA", min.cells = 3, min.features = 200
)
pd23 <- AddMetaData(pd23, "distal", col.name = "group")
pd23 <- AddMetaData(pd23, "pd23", col.name = "batch")
pd23 <- RenameCells(pd23, add.cell.id = "pd23")
pd23[["percent.mt"]] <- PercentageFeatureSet(pd23, pattern = "^MT-")
pd23[["percent.ribo"]] <- PercentageFeatureSet(pd23, pattern = "^RP[SL][[:digit:]]")
pd23 <- AddMetaData(pd23, "163", "age")
pd23 <- AddMetaData(pd23, col.name = "ID", "23 Week")



setwd("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Nuclei/HT650 150d distal lung/Cellbender v03")
HT650.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)
HT650 <- CreateSeuratObject(counts = HT650.data, project = "HT650", min.cells = 3, min.features = 200)
HT650 <- RenameCells(HT650, add.cell.id = "HT650")
HT650[["percent.mt"]] <- PercentageFeatureSet(HT650, pattern = "^MT-")
HT650[["percent.ribo"]] <- PercentageFeatureSet(HT650, pattern = "^RP[SL][[:digit:]]")


HT650 <- AddMetaData(HT650, "150", "age")
HT650 <- AddMetaData(HT650, "HT650", "batch")
HT650 <- AddMetaData(HT650, "HT650", "name")
HT650 <- AddMetaData(HT650, "distal", "group")
HT650 <- AddMetaData(HT650, col.name = "ID", "21 week")

intermediate.timepoints <- merge(HT650, c(pd16, pd23))
intermediate.timepoints <- JoinLayers(intermediate.timepoints)
intermediate.timepoints <- subset(intermediate.timepoints, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 5 & percent.ribo < 7.5 & nCount_RNA > 1000)

intermediate.timepoints <- SplitObject(intermediate.timepoints, split.by = "ID")


i = 1
for (i in seq_along(intermediate.timepoints)) {
  intermediate.timepoints[[i]] <- NormalizeData(intermediate.timepoints[[i]]) %>% FindVariableFeatures() %>% CellCycleScoring(g2m.features = g2m.genes, s.features = s.genes)
}
features <- SelectIntegrationFeatures(intermediate.timepoints, nfeatures = 3000)
for (i in seq_along(along.with = intermediate.timepoints)) {
  intermediate.timepoints[[i]] <- ScaleData(intermediate.timepoints[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(intermediate.timepoints, anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 3)
intermediate.timepoints.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(intermediate.timepoints.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

intermediate.timepoints.integrated <- ScaleData(intermediate.timepoints.integrated)
intermediate.timepoints.integrated <- RunPCA(intermediate.timepoints.integrated)
ElbowPlot(intermediate.timepoints.integrated, ndims = 50)
intermediate.timepoints.integrated <- RunUMAP(intermediate.timepoints.integrated, dims = 1:30, reduction.name = "umap", return.model = TRUE)
intermediate.timepoints.integrated <- FindNeighbors(intermediate.timepoints.integrated, dims = 1:30)
intermediate.timepoints.integrated <- FindClusters(intermediate.timepoints.integrated, resolution = 1.0)


anchors <- FindTransferAnchors(
  reference = pediatric,
  query = intermediate.timepoints.integrated,
  dims = 1:30,
  reference.reduction = "pca",
  k.filter = 200, 
  reference.assay = "integrated",
  query.assay = "integrated"
)

# Transfer Labels
intermediate.timepoints.transferred <- TransferData(
  anchorset = anchors,
  refdata = pediatric$full.dotplot.annotation,
  dims = 1:30
)
intermediate.timepoints.integrated$predicted_lvl3 <- intermediate.timepoints.transferred$predicted.id


####Score Datasets Using Maturation Gene Sets: Intermediate Timepoints ######
DefaultAssay(intermediate.timepoints.integrated) <- "RNA"
intermediate.timepoints.integrated <- JoinLayers(intermediate.timepoints.integrated)
intermediate.timepoints.integrated <- NormalizeData(intermediate.timepoints.integrated)

intermediate.timepoints.integrated_cell_types <- unique(intermediate.timepoints.integrated$predicted_lvl3)
expr_matrix <- GetAssayData(intermediate.timepoints.integrated, slot = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix)

for (cell_type in intermediate.timepoints.integrated_cell_types) {
  intermediate.timepoints.integrated[[paste0(cell_type, "_Maturation_AUC_Score")]] <- NA
}

for (cell_type in intermediate.timepoints.integrated_cell_types) {
  if (!cell_type %in% names(maturation_genes_on)) next
  
  genes_on <- maturation_genes_on[[cell_type]]
  genes_off <- maturation_genes_off[[cell_type]]
  
  if (length(genes_on) == 0 | length(genes_off) == 0) next
  
  intermediate.timepoints.integrated_subset <- subset(intermediate.timepoints.integrated, predicted_lvl3 == cell_type)
  
  genes_on_set <- GeneSet(genes_on, setName = paste0(cell_type, "_Maturation_On"))
  genes_off_set <- GeneSet(genes_off, setName = paste0(cell_type, "_Maturation_Off"))
  
  gene_sets <- GeneSetCollection(list(genes_on_set, genes_off_set))
  
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings[, colnames(intermediate.timepoints.integrated_subset)])
  
  intermediate.timepoints.integrated_subset[[paste0(cell_type, "_Maturation_AUC_Score")]] <- 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_On"), ] - 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_Off"), ]
  
  matching_cells <- Cells(intermediate.timepoints.integrated_subset)
  intermediate.timepoints.integrated@meta.data[matching_cells, paste0(cell_type, "_Maturation_AUC_Score")] <- 
    intermediate.timepoints.integrated_subset[[paste0(cell_type, "_Maturation_AUC_Score")]]
}

####Score Datasets Using Maturation Gene Sets: PIG ######
DefaultAssay(PIG.fastmnn.cleaned.integrated) <- "RNA"
PIG.fastmnn.cleaned.integrated <- JoinLayers(PIG.fastmnn.cleaned.integrated)

PIG.fastmnn.cleaned.integrated <- NormalizeData(PIG.fastmnn.cleaned.integrated)


# Get unique PIG.fastmnn.cleaned.integrated cell types
PIG.fastmnn.cleaned.integrated_cell_types <- unique(PIG.fastmnn.cleaned.integrated$predicted_lvl3)

# Get expression matrix
expr_matrix <- GetAssayData(PIG.fastmnn.cleaned.integrated, slot = "counts")

# Compute AUCell rankings once (since it's dataset-wide)
cells_rankings <- AUCell_buildRankings(expr_matrix)

# Ensure metadata columns exist before assignment
for (cell_type in PIG.fastmnn.cleaned.integrated_cell_types) {
  PIG.fastmnn.cleaned.integrated[[paste0(cell_type, "_Maturation_AUC_Score")]] <- NA  # Initialize empty column
}

# Iterate through each PIG.fastmnn.cleaned.integrated cell type
for (cell_type in PIG.fastmnn.cleaned.integrated_cell_types) {
  if (!cell_type %in% names(maturation_genes_on)) next  # Skip if no matching maturation genes
  
  genes_on <- maturation_genes_on[[cell_type]]
  genes_off <- maturation_genes_off[[cell_type]]
  
  if (length(genes_on) == 0 | length(genes_off) == 0) next  # Skip if gene sets are empty
  
  PIG.fastmnn.cleaned.integrated_subset <- subset(PIG.fastmnn.cleaned.integrated, predicted_lvl3 == cell_type)
  
  # AUCell Scoring
  genes_on_set <- GeneSet(genes_on, setName = paste0(cell_type, "_Maturation_On"))
  genes_off_set <- GeneSet(genes_off, setName = paste0(cell_type, "_Maturation_Off"))
  
  gene_sets <- GeneSetCollection(list(genes_on_set, genes_off_set))
  
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings[, colnames(PIG.fastmnn.cleaned.integrated_subset)])
  
  PIG.fastmnn.cleaned.integrated_subset[[paste0(cell_type, "_Maturation_AUC_Score")]] <- 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_On"), ] - 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_Off"), ]
  
  matching_cells <- Cells(PIG.fastmnn.cleaned.integrated_subset)  
  PIG.fastmnn.cleaned.integrated@meta.data[matching_cells, paste0(cell_type, "_Maturation_AUC_Score")] <- 
    PIG.fastmnn.cleaned.integrated_subset[[paste0(cell_type, "_Maturation_AUC_Score")]]
}
####Score Datasets Using Maturation Gene Sets: BPD ######

DefaultAssay(BPD.fastmnn.cleaned.integrated) <- "RNA"
BPD.fastmnn.cleaned.integrated <- JoinLayers(BPD.fastmnn.cleaned.integrated)
BPD.fastmnn.cleaned.integrated <- NormalizeData(BPD.fastmnn.cleaned.integrated)


# Get unique BPD.fastmnn.cleaned.integrated cell types
BPD.fastmnn.cleaned.integrated_cell_types <- unique(BPD.fastmnn.cleaned.integrated$predicted_lvl3)

# Get expression matrix
expr_matrix <- GetAssayData(BPD.fastmnn.cleaned.integrated, slot = "counts")

# Compute AUCell rankings once (since it's dataset-wide)
cells_rankings <- AUCell_buildRankings(expr_matrix)

# Ensure metadata columns exist before assignment
for (cell_type in BPD.fastmnn.cleaned.integrated_cell_types) {
  BPD.fastmnn.cleaned.integrated[[paste0(cell_type, "_Maturation_AUC_Score")]] <- NA  # Initialize empty column
}

# Iterate through each BPD.fastmnn.cleaned.integrated cell type
for (cell_type in BPD.fastmnn.cleaned.integrated_cell_types) {
  if (!cell_type %in% names(maturation_genes_on)) next  # Skip if no matching maturation genes
  
  genes_on <- maturation_genes_on[[cell_type]]
  genes_off <- maturation_genes_off[[cell_type]]
  
  if (length(genes_on) == 0 | length(genes_off) == 0) next  # Skip if gene sets are empty
  
  BPD.fastmnn.cleaned.integrated_subset <- subset(BPD.fastmnn.cleaned.integrated, predicted_lvl3 == cell_type)
  
  # AUCell Scoring
  genes_on_set <- GeneSet(genes_on, setName = paste0(cell_type, "_Maturation_On"))
  genes_off_set <- GeneSet(genes_off, setName = paste0(cell_type, "_Maturation_Off"))
  
  gene_sets <- GeneSetCollection(list(genes_on_set, genes_off_set))
  
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings[, colnames(BPD.fastmnn.cleaned.integrated_subset)])
  
  BPD.fastmnn.cleaned.integrated_subset[[paste0(cell_type, "_Maturation_AUC_Score")]] <- 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_On"), ] - 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_Off"), ]
  
  matching_cells <- Cells(BPD.fastmnn.cleaned.integrated_subset)  
  BPD.fastmnn.cleaned.integrated@meta.data[matching_cells, paste0(cell_type, "_Maturation_AUC_Score")] <- 
    BPD.fastmnn.cleaned.integrated_subset[[paste0(cell_type, "_Maturation_AUC_Score")]]
}

####Score Datasets Using Maturation Gene Sets: Pediatric ######

#pediatric <- JoinLayers(pediatric)
DefaultAssay(pediatric) <- "RNA"
pediatric <- NormalizeData(pediatric)

pediatric_cell_types <- unique(pediatric$full.dotplot.annotation)
expr_matrix <- GetAssayData(pediatric, slot = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix)

for (cell_type in pediatric_cell_types) {
  pediatric[[paste0(cell_type, "_Maturation_AUC_Score")]] <- NA
}

for (cell_type in pediatric_cell_types) {
  if (!cell_type %in% names(maturation_genes_on)) next
  
  genes_on <- maturation_genes_on[[cell_type]]
  genes_off <- maturation_genes_off[[cell_type]]
  
  if (length(genes_on) == 0 | length(genes_off) == 0) next
  
  pediatric_subset <- subset(pediatric, full.dotplot.annotation == cell_type)
  
  genes_on_set <- GeneSet(genes_on, setName = paste0(cell_type, "_Maturation_On"))
  genes_off_set <- GeneSet(genes_off, setName = paste0(cell_type, "_Maturation_Off"))
  
  gene_sets <- GeneSetCollection(list(genes_on_set, genes_off_set))
  
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings[, colnames(pediatric_subset)])
  
  pediatric_subset[[paste0(cell_type, "_Maturation_AUC_Score")]] <- 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_On"), ] - 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_Off"), ]
  
  matching_cells <- Cells(pediatric_subset)
  pediatric@meta.data[matching_cells, paste0(cell_type, "_Maturation_AUC_Score")] <- 
    pediatric_subset[[paste0(cell_type, "_Maturation_AUC_Score")]]
}

####Score Datasets Using Maturation Gene Sets: Fetal Distal ######
fetal.distal <- JoinLayers(fetal.distal)
DefaultAssay(fetal.distal) <- "RNA"
fetal.distal <- NormalizeData(fetal.distal)

fetal.distal_cell_types <- unique(fetal.distal$predicted_lvl3)
expr_matrix <- GetAssayData(fetal.distal, slot = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix)

for (cell_type in fetal.distal_cell_types) {
  fetal.distal[[paste0(cell_type, "_Maturation_AUC_Score")]] <- NA
}

for (cell_type in fetal.distal_cell_types) {
  if (!cell_type %in% names(maturation_genes_on)) next
  
  genes_on <- maturation_genes_on[[cell_type]]
  genes_off <- maturation_genes_off[[cell_type]]
  
  if (length(genes_on) == 0 | length(genes_off) == 0) next
  
  fetal.distal_subset <- subset(fetal.distal, predicted_lvl3 == cell_type)
  
  genes_on_set <- GeneSet(genes_on, setName = paste0(cell_type, "_Maturation_On"))
  genes_off_set <- GeneSet(genes_off, setName = paste0(cell_type, "_Maturation_Off"))
  
  gene_sets <- GeneSetCollection(list(genes_on_set, genes_off_set))
  
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings[, colnames(fetal.distal_subset)])
  
  fetal.distal_subset[[paste0(cell_type, "_Maturation_AUC_Score")]] <- 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_On"), ] - 
    getAUC(cells_AUC)[paste0(cell_type, "_Maturation_Off"), ]
  
  matching_cells <- Cells(fetal.distal_subset)
  fetal.distal@meta.data[matching_cells, paste0(cell_type, "_Maturation_AUC_Score")] <- 
    fetal.distal_subset[[paste0(cell_type, "_Maturation_AUC_Score")]]
}

saveRDS(fetal.distal, "fetal_distal_with_maturation_scores.rds")
saveRDS(pediatric, "pediatric_with_maturation_scores.rds")
saveRDS(PIG.fastmnn.cleaned.integrated, "PIG_with_maturation_scores.rds")
saveRDS(BPD.fastmnn.cleaned.integrated, "BPD_with_maturation_scores.rds")
saveRDS(intermediate.timepoints, "intermediate_timepoints_with_maturation_scores.rds")

####Prepare Data for Violin Plots####

BPD.fastmnn.cleaned.integrated <- AddMetaData(BPD.fastmnn.cleaned.integrated, col.name = "harmonized_annotation", "BPD")
pediatric <- AddMetaData(pediatric, col.name = "harmonized_annotation", "Pediatric Older Than 6 months")
fetal.distal <- AddMetaData(fetal.distal, col.name = "harmonized_annotation", "Fetal")
PIG.fastmnn.cleaned.integrated <- AddMetaData(PIG.fastmnn.cleaned.integrated, col.name = "harmonized_annotation", "P.I.G.")
intermediate.timepoints.integrated <- AddMetaData(intermediate.timepoints.integrated, col.name = "harmonized_annotation", intermediate.timepoints.integrated$ID)
intermediate.timepoints.integrated <- AddMetaData(intermediate.timepoints.integrated, col.name = "dataset", intermediate.timepoints.integrated$ID)

pediatric$dataset <- "Pediatric Older Than 6 months"
fetal.distal$dataset <- "Fetal Distal"
PIG.fastmnn.cleaned.integrated$dataset <- "PIG"
BPD.fastmnn.cleaned.integrated$dataset <- "BPD"

#subset scored pediatric for 6 month and less comparitor
Idents(pediatric) <- "age2"
pediatric.age.match <- subset(pediatric, idents = c(6, 0, 3, 5))
pediatric.age.match <- AddMetaData(pediatric.age.match, col.name = "dataset", "Pediatric 6 Months and Younger")
pediatric.age.match <- AddMetaData(pediatric.age.match, col.name = "harmonized_annotation", "Pediatric 6 Months and Younger")
pediatric <- AddMetaData(pediatric, col.name = "dataset", pediatric.age.match$dataset)
pediatric <- AddMetaData(pediatric, col.name = "harmonized_annotation", pediatric.age.match$harmonized_annotation)

Idents(pediatric) <- "age2"
pediatric.birth <- subset(pediatric, idents = c(0))
pediatric.birth <- AddMetaData(pediatric.birth, col.name = "dataset", "Pediatric Birth")
pediatric.birth <- AddMetaData(pediatric.birth, col.name = "harmonized_annotation", "Pediatric Birth")
pediatric <- AddMetaData(pediatric, col.name = "dataset", pediatric.birth$dataset)
pediatric <- AddMetaData(pediatric, col.name = "harmonized_annotation", pediatric.birth$harmonized_annotation)


intermediate.timepoints.integrated <- AddMetaData(intermediate.timepoints.integrated, col.name = "dataset", intermediate.timepoints.integrated$ID)



merged_disease <- merge(pediatric, y = list(fetal.distal, PIG.fastmnn.cleaned.integrated, BPD.fastmnn.cleaned.integrated, intermediate.timepoints.integrated), add.cell.ids = c("Pediatric", "FetalDistal", "PIG", "BPD", "Int"))
pediatric <- RenameCells(pediatric, add.cell.id = "Pediatric")
merged_disease <- AddMetaData(merged_disease, col.name = "harmonized_annotation", pediatric$harmonized_annotation )
merged_disease <- AddMetaData(merged_disease, col.name = "dataset", merged_disease$harmonized_annotation )

merged_disease$harmonized_annotation <- factor(merged_disease$harmonized_annotation, levels = c("Fetal", "16 Week", "21 week", "23 Week", "P.I.G.", "BPD", "Pediatric Birth", "Pediatric 6 Months and Younger", "Pediatric Older Than 6 months"))
merged_disease$dataset <- factor(merged_disease$dataset, levels = c("Fetal Distal", "16 Week", "21 week", "23 Week", "PIG", "BPD", "Pediatric Birth", "Pediatric 6 Months and Younger", "Pediatric Older Than 6 months"))

#####Violin Plots of Maturation Scores#####


vln_plots <- list()

# Iterate through each cell type
for (cell_type in pediatric_cell_types) {
  cell_type_2 <- gsub("[ +\\-]", ".", cell_type)
  
  if (!cell_type %in% names(maturation_genes_on)) {
    print(paste("Skipping", cell_type, "- No maturation genes found"))
    next
  }
  
  auc_score_col <- paste0(cell_type_2, "_Maturation_AUC_Score")
  
  # Ensure the feature exists in the merged dataset
  if (!(auc_score_col %in% colnames(merged_disease@meta.data))) {
    print(paste("Skipping", cell_type, "- No AUC score in metadata"))
    next
  }
  
  # Remove NA values
  merged_filtered <- subset(merged_disease, cells = which(
    !is.na(merged_disease[[auc_score_col]])
  ))
  
  if (ncol(merged_filtered) == 0) {
    print(paste("Skipping", cell_type, "- No non-NA scores found"))
    next
  }
  
  # Create violin plot with the title now only being the cell type
  p <- VlnPlot(
    merged_filtered, 
    features = c(auc_score_col), 
    group.by = "harmonized_annotation", 
    split.by = "dataset", 
    pt.size = 0
  ) + ggtitle(cell_type)
  
  # Overlay a box plot, increase font sizes, and remove the legend
  p <- p + 
    geom_boxplot(width = 0.3, outlier.shape = NA, color = "black") + 
    theme(
      text = element_text(size = 36),       # increases all text elements
      axis.text = element_text(size = 36),    # adjust axis text if needed
      axis.title = element_text(size = 0),    # adjust axis titles if needed
      plot.title = element_text(size = 40),   # adjust plot title size
      legend.position = "none"                # removes the legend
    )
  
  vln_plots[[cell_type]] <- p
}

# Debugging output
print(paste("Total plots generated:", length(vln_plots)))
if (length(vln_plots) < length(pediatric_cell_types)) {
  print("Warning: Not all cell types have plots! Check logs above.")
}

# Reorder the plots based on the desired order
identities <- c("AT2", "Transitional", "AT1", "Respiratory Bronchiole", "Basal",
                "Multiciliated", "Neuroendocrine", "Macrophages", "TREM2+ Macrophages", 
                "Monocytes", "CD1c+ DC", "pDC", "NK Cell", "CD4+ T-Cell", "CD8+ T-Cell", 
                "Basophils", "B-Cells", "Neutrophils", "Lymphatic Endothelial", "Arterial Endothelial",
                "Venous Endothelial", "Capillary Endothelial", "Capillary Aerocyte", 
                "Alveolar Fibroblast", "CCL2-Positive Fibroblast", "Fibroblast", 
                "Alveolar Myofibroblast", "Adventitial Mesenchyme", 
                "Pericyte")

ordered_plots <- list()
for (cell in identities) {
  if (cell %in% names(vln_plots)) {
    ordered_plots[[cell]] <- vln_plots[[cell]]
  }
}

# Save all plots to a single PDF with dimensions 100 x 30 inches.
pdf("Maturation_AUC_Scores_Violin_Plots_with_intermediate_timepoints.pdf", width = 66, height = 51)
if (length(ordered_plots) > 0) {
  print(wrap_plots(ordered_plots, nrow = 4))  # Arrange plots in 2 rows
} else {
  print("No valid plots were generated. Check if AUC scores exist in metadata.")
}
dev.off()



#####FORCE ALL ALVEOLAR FIBROBLASTS TO ALIGN####
Idents(fetal.distal) <- "predicted_lvl3"
fetal.distal <- RenameIdents(fetal.distal, "CCL2-Positive Fibroblast" = "Alveolar Fibroblast", "Fibroblast" = "Alveolar Fibroblast")
fetal.distal.alveolar.mes <- subset(fetal.distal, idents = "Alveolar Fibroblast")
Idents(pediatric) <- "full.dotplot.annotation"
pediatric.alveolar.mes <- subset(pediatric, idents = "Alveolar Fibroblast")
fetal.distal.alveolar.mes <- AddMetaData(fetal.distal.alveolar.mes, col.name = "Timing", "Fetal")
pediatric.alveolar.mes <- AddMetaData(pediatric.alveolar.mes, col.name = "Timing", "Pediatric")

alveolar.mes.merge <- merge(fetal.distal.alveolar.mes, pediatric.alveolar.mes)
DefaultAssay(alveolar.mes.merge) <- "RNA"
#alveolar.mes.merge <- JoinLayers(alveolar.mes.merge)
alveolar.mes.merge <- NormalizeData(alveolar.mes.merge)
trimmed.features = rownames(alveolar.mes.merge)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
Idents(alveolar.mes.merge) <- "Timing"
markers <- FindMarkers(alveolar.mes.merge, ident.1 = "Pediatric", ident.2 = "Fetal", features = trimmed.features, min.pct = 0.25, logfc.threshold = 0.25)


top_on <- markers %>%
  filter(avg_log2FC > 1.5) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100) %>%
  rownames()

# Filter top 100 genes that turn OFF (logFC < 0)
top_off <- markers %>%
  filter(avg_log2FC < -1.5) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 100) %>%
  rownames()

# Store in the maturation gene lists
alv.fib.maturation_genes_on <- top_on
alv.fib.maturation_genes_off<- top_off



Idents(intermediate.timepoints.integrated) <- "predicted_lvl3"
intermediate.timepoints.alveolar.fibroblast <- subset(intermediate.timepoints.integrated, idents = c("Alveolar Fibroblast", "Fibroblast"))
intermediate.timepoints.alveolar.fibroblast <- AddMetaData(intermediate.timepoints.alveolar.fibroblast, col.name = "Timing", intermediate.timepoints.integrated$ID)
Idents(BPD.fastmnn.cleaned.integrated) <- "predicted_lvl3"
BPD.alveolar.fibroblast <- subset(BPD.fastmnn.cleaned.integrated, idents = c("Alveolar Fibroblast", "CCL2-Positive Fibroblast", "Fibroblast"))
BPD.alveolar.fibroblast <- AddMetaData(BPD.alveolar.fibroblast, col.name = "Timing", "BPD")

Idents(PIG.fastmnn.cleaned.integrated) <- "predicted_lvl3"
PIG.alveolar.fibroblast <- subset(PIG.fastmnn.cleaned.integrated, idents = c("Alveolar Fibroblast", "CCL2-Positive Fibroblast", "Fibroblast"))
PIG.alveolar.fibroblast <- AddMetaData(PIG.alveolar.fibroblast, col.name = "Timing", "PIG")


alveolar.all.mes.merge <- merge(PIG.alveolar.fibroblast, c(BPD.alveolar.fibroblast, intermediate.timepoints.alveolar.fibroblast, pediatric.alveolar.mes, fetal.distal.alveolar.mes))
DefaultAssay(alveolar.all.mes.merge) <- "RNA"
alveolar.all.mes.merge <- JoinLayers(alveolar.all.mes.merge)
alveolar.all.mes.merge <- NormalizeData(alveolar.all.mes.merge)

alveolar.all.mes.merge[["Maturation_AUC_Score"]] <- NA

expr_matrix <- GetAssayData(alveolar.all.mes.merge, slot = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix)
if (exists("alv.fib.maturation_genes_on") && exists("alv.fib.maturation_genes_off")) {
  genes_on <- alv.fib.maturation_genes_on
  genes_off <- alv.fib.maturation_genes_off
  
  if (length(genes_on) > 0 && length(genes_off) > 0) {
    genes_on_set <- GeneSet(genes_on, setName = "Maturation_On")
    genes_off_set <- GeneSet(genes_off, setName = "Maturation_Off")
    
    gene_sets <- GeneSetCollection(list(genes_on_set, genes_off_set))
    
    # Compute AUC scores
    cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings[, colnames(alveolar.all.mes.merge)])
    
    # Assign AUC scores to metadata
    alveolar.all.mes.merge[["Maturation_AUC_Score"]] <- 
      getAUC(cells_AUC)["Maturation_On", ] - 
      getAUC(cells_AUC)["Maturation_Off", ]
  }
}
pdf("Alveolar Fibroblast Forced Score VlnPlot.pdf", width = 5, height = 4)
VlnPlot(alveolar.all.mes.merge, "Maturation_AUC_Score", split.by = "Timing", group.by = "Timing", pt.size = 0) + NoLegend()
dev.off()

filtered_markers <- markers[markers$p_val_adj < 0.05, ]

# Get positive fold change (upregulated in Pediatric)
positive_genes <- rownames(filtered_markers[filtered_markers$avg_log2FC > 0, ])

# Get negative fold change (downregulated in Pediatric)
negative_genes <- rownames(filtered_markers[filtered_markers$avg_log2FC < 0, ])

NormalvPIG.alv.fib <- merge(pediatric.alveolar.mes, PIG.alveolar.fibroblast)
DefaultAssay(NormalvPIG.alv.fib) <- "RNA"
NormalvPIG.alv.fib <- JoinLayers(NormalvPIG.alv.fib)
NormalvPIG.alv.fib <- NormalizeData(NormalvPIG.alv.fib)
Idents(NormalvPIG.alv.fib) <- "Timing"
perform_enrichment_analysis_two_clusters(NormalvPIG.alv.fib, "PIG", "Pediatric", output_dir = "PIGvPediatricEnrichment")

NormalvBPD.alv.fib <- merge(pediatric.alveolar.mes, BPD.alveolar.fibroblast)
DefaultAssay(NormalvBPD.alv.fib) <- "RNA"
NormalvBPD.alv.fib <- JoinLayers(NormalvBPD.alv.fib)
NormalvBPD.alv.fib <- NormalizeData(NormalvBPD.alv.fib)
Idents(NormalvBPD.alv.fib) <- "Timing"
perform_enrichment_analysis_two_clusters(NormalvBPD.alv.fib, "BPD", "Pediatric", output_dir = "BPDvPediatricEnrichment")