#My R functions

dist2pcoaplotAnno = function(distobj, metadata, groupColumn) {
  # 20170920 - Patrick Munk  
  # Takes (1) a distance / dissim object, (2) sample metadata dataframe, (3) metadata column used for color
  # Returns a PCoA plot with color-annotated sample groups
  pcoa = cmdscale(distobj, eig=T)
  #Separate the different objects
  points = pcoa$points
  eig = pcoa$eig
  #Calculate the variation explained by PCoA1 and 2
  eig_1_2 = eig[1:2]/sum(eig)*100
  eig_1 = paste("PCoA1", round(eig_1_2[1], digits=2),"% variance")
  eig_2 = paste("PCoA2", round(eig_1_2[2], digits=2),"% variance")
  #Construct a dataset for plotting and add color variable
  data = data.frame(sample=rownames(points),points)
  metadata = data.frame(sample=metadata[,1], Group=metadata[,groupColumn])
  data = merge(data, metadata, by="sample", all.x = T)
  #Calculate a reasonable offset for text label
  yoffset = (abs(max(data[,3])) + abs(min(data[,3]))) / 50
  #Produce plot
  p1 = ggplot(data, aes(X1,X2)) + geom_point(size=4, aes(col=Group)) + geom_text(nudge_y=yoffset, aes(label=sample)) + labs(x=eig_1, y=eig_2) + scale_color_brewer(palette = "Set1")
  return(p1)
}


geneHeatmap = function(data, metadata, annotateFeature) {
  #20170920 - Patrick Munk
  #Produces a gene heatmap with gene rows (BrayCurtis clustered) and sample columns (Corr. clustered)
  #Input: (1) abundance dataframe, (2) metadata (sample x attributes), (3) annotation attribute (colname in metadata)
  #Return: A sample-annotated heatmap with BC-clustered columns
  #
  #Make sure the required packages are loaded
  require(pheatmap)
  require(vegan)
  require(RColorBrewer)
  #Calculate sample-distance matrix with vegan (Bray-Curtis)
  #Note, that this is done on the full set, not just the shown genes
  distmatrix <- vegdist(decostand(t(data), method="hellinger"), method="bray")
  #Subset the data, so a maximum of 50 features will be shown (most abundant) or most variable
  if (nrow(data) > 50) {
    #data <- data[order(rowSums(data), decreasing = T),]
    data <- data[order(apply(data,1,median), decreasing = T),]
    #data <- data[order(apply(data,1,sd) / apply(data,1,mean), decreasing = T),]
    data <- data[1:50,]
  }
  #Log-transform the data to be plotted
  data <- log(data+1)
  #Make annotation table for the heatmap
  colannodf <- data.frame(metadata[,annotateFeature], row.names = metadata[,1])
  colnames(colannodf) <- annotateFeature
  #Create the colors used for the heatmap sample annotation
  ann_colors <- brewer.pal(n = length(unique(colannodf[,annotateFeature])),name = "Set1")
  names(ann_colors) <- sort(unique(colannodf[,annotateFeature]))
  #Create the color list required for the heatmap
  ann_colors = list(annotateFeature = ann_colors)
  names(ann_colors) = annotateFeature
  #Draw the heatmap
  pheatmap(data, margins=c(8,8), treeheight_row = 100, treeheight_col = 100, scale="none", clustering_distance_cols = distmatrix, clustering_distance_rows = "correlation", annotation_col = colannodf, annotation_colors = ann_colors)
}
