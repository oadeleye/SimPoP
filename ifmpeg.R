if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
#generate the full graph

g <- graph.edgelist(as.matrix(data[,c(1,2)]),directed=F)
write.table(g,"g.txt",sep = ",")
E(g)$time <- data[,3]
#generate a cool palette for the graph (darker colors = older nodes)
YlOrBr.pal <- colorRampPalette(brewer.pal(8,"YlOrRd"))
#colors for the nodes are chosen from the very beginning
V(g)$color <- rev(YlOrBr.pal(vcount(g)))[as.numeric(V(g)$name)]
ti<-1
gt <- delete_edges(g,which(E(g)$time > ti))
layout.old <- norm_coords(layout.graphopt(gt), xmin = -1, xmax = 1, ymin = -1, ymax = 1)
total_time <- max(E(g)$time)
dt <- 1
png(file="example%03d.png", width=800,height=450)
#Time loop starts
for(time in seq(1,total_time,dt)){
  
  gt <- delete_edges(g,which(E(g)$time > time))
  layout.new <- layout_with_fr(gt,coords=layout.old,niter=10,start.temp=0.05,grid="nogrid")
  plot(gt,layout=layout.new,vertex.label="",vertex.size=1+2*log(degree(gt)),vertex.frame.color=V(g)$color,edge.width=1.5,asp=9/16,margin=-0.15)
  #use the new layout in the next round
  fp<-degree(gt)/mean(degree(gt))
  c=fp*f
  print(c)
  layout.old <- layout.new
}
dev.off()
  
  