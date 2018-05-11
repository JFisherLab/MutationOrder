library(dplyr) # Version 0.7.0 or above
library(reshape2)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(tidyr)

#------Visualisation Initialisation------

edge_width = 3.5 #thickness of edges
arrow_size = 0.5 #size of arrow heads
arrow_width = 3 #width of arrow head
label_dist = 0 #distance of edge label along arrow
legend_cex = 2.7 #size of text in legends
legend_point_size = 4 #size of symbol in legend for phenotype vertices
legend_line_size = 6 #size of symbol in legend for edges
edge_label_size = 2.1 #size of text in edge labels

#------Functions------

blue.pal <- brewer.pal(n = 7, "Blues")
red.pal <- brewer.pal(n = 7, "YlOrRd")

edge.pal <- brewer.pal(8, "Dark2")[3:8] #remove first colour as green is too confusing

Edge.col <- function(edge.list, pal = edge.pal) {
  # Add colours to edges based on node which is mutated on traversing that edge
  #
  # Args: 
  #      edge.list: List of edges with columns from, to, and change with the last noting the node which mutated and the former being the vertices
  #      pal: palette from presets for RColorBrewer to use
  #
  # Returns: list of colours to apply to edges, to be used for edge.color in igraph::plot
  
  edge.index <- as.factor(edge.list$change)
  edge.pal <- pal
  edge.col <- edge.pal[edge.index]
  return(edge.col)
}

to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  # From https://github.com/dbarneche/2014-07-14-Dalhousie/
  # Takes a plotting function and outputs the plot. Can do any output type: png, svg, pdf
  # Args: 
  #      expr: A function that plots something, which alone evaluates within R
  #      dev: Type of output e.g. pdf, no need for quotes
  #      filename: name under which to save file
  #      ... : takes arguments for dev e.g. width and height
  # Returns: 
  #         Status of output
  
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

make.verts.table <- function(summary.line) {
  # Make a unique vertex id in the correct format for `igraph::graph_from_data_frame`
  # 
  # Args:
  #      List of dataframes with the node info for each cell line e.g. `expanded.summary.by.mutation.attractor`
  # Returns:
  #         Same list of dataframes with FIRST column being the VertID
  
  verts.table <- mutate(summary.line, VertID = paste(Mutation, Attractor, sep = "_AT_")) %>% 
    select(VertID, everything()) #put id column first, req for igraph::graph_from_data_frame
  return(verts.table)
}

make.edge.list <- function(verts){
  # Dependencies: library(dplyr) of version 0.7.0 or later 
  # Added mutation level as well as node to `Change`
  # 
  # Take a list of nodes and make all the edges for a powerset
  # 
  # Args: node.table: Data.frame with:
  ## unique node id: VertID, a column for every netw Node, 
  ## Column for all BMA nodes mutated Node1, Node2, etc
  # 
  # Returns: 
  #         List of edges for the visualisation, as a edge list not an adjacency matrix
  
  node.cols <- verts %>% 
    select(contains("Node")) %>% 
    colnames # makes a vector of all the Node columns
  
  mut.cols <- verts %>% 
    select(contains("Mut")) %>% 
    select(-Mutation) %>%
    select(-Number.Muts) %>%
    colnames # makes a vector of all the Node columns
  
  #Bind together Node and Mut for reporting as final overall change
  for( i in 1:length(node.cols)){
    verts <- unite(verts, !!paste0("NodeMut", i), node.cols[i], mut.cols[i], remove = FALSE, sep = "@") #!! is to enquote paste for unite in dplyr fashion
  }
  
  nodemut.cols <- verts %>% 
    select(contains("NodeMut")) %>% 
    colnames
  
  verts <- verts %>% unite(NodeMuts, nodemut.cols, remove = FALSE, sep = "_") # Unite the Nodes into one string (but preserve dropping of the Attractor component from VertID)
  
  set <- verts[["NodeMuts"]] %>% as.list 
  names(set) <- verts[["VertID"]] # Convert to names list for use of set functions of `dplyr`
  set <- lapply(set, strsplit, split = "_") # Split nodes into a vector of nodes
  set <- lapply(set, function(x) {x %>% unlist %>% setdiff("NA@NA")}) # Remove `NA`
  set <- set[order(sapply(set, length), decreasing=FALSE)] # Put in order of number of mutations from least to most
  
  edge.list <- data.frame(from = character(), to = character(), change = character()) # Initialise
  
  
  # Look through pairs of mutatations, from top to bottom as arranged in ascending order of length
  # If the number of mutations of the latter is 1 greater than the former, consider for edge. Otherwise we would be losing a mutation e.g. going from A to B necesitates losing A, which we assume does not occur (losing a gene counts as a new mutation)
  # Then check for intersection, only allow edge is the intersection contains all old genes, so cannot go from A to BC, as this again entails losing mutation A. 
  for(i in 1:length(set)){
    for (j in i:length(set)){
      
      # Check if only one node more between vertices
      if(length(set[[i]]) == length(set[[j]])-1){ 
        
        sect <- intersect(set[[j]], set[[i]]) # Nodes shared between vertices
        diff <- setdiff(set[[j]], set[[i]]) # Nodes which differ between vertices. Note, order matters for this function
        
        # Check if only one node has CHANGED between vertices
        if(length(sect) == length(set[[j]])-1){
          
          edge = data.frame(from = names(set[i]), to = names(set[j]), change = diff) # Make a new row for output data frame
          edge.list <- rbind(edge.list, edge) # Write to output frame
        }
      }
    }
  }
  return(edge.list)
}

plot.rel <- function(line, pheno, inet, vertices, label = TRUE, v.size = 15, e.label.size = edge_label_size){
  # Function to plot relative graphs 
  # 
  # Args: 
  #      line = cell line to plot
  #      pheno = phenotype to plot, either apop or prolif
  #      inet = list of igraph net objects
  #      vertices = list of vertices for above objects
  #      label = Plot edge labels or not
  #      v.size = vertex size
  #      e.label.size = edge label size
  # 
  # Returns a plot
  
  data <- inet[[line]]
  
  if(pheno == "apop"){
    plot(data
         , vertex.color = V(data)$color.apop.rel
         #, layout = layout_as_tree(net$`SKBR3`, root = "EGF__1__Wnt1__1_AT_0")
         , layout = layout.sugiyama(data, hgap = 2, vgap = 10)$layout
         , vertex.label = NA
         , vertex.size = v.size
         , edge.label = if (label == TRUE) {E(data)$change} else {NA}
         , edge.label.family = "sans"
         , edge.width = edge_width
         , edge.arrow.width = arrow_width
         , edge.arrow.size = arrow_size
         , edge.label.color = "black"
         , edge.label.cex = e.label.size
         , edge.label.dist = label_dist
         , edge.color = E(data)$colour
         #, main = paste0("Apoptosis in ", line)
    )
    lgnd = legend("topleft"
           , title = "Mean\n Apoptosis"
           , bty = "n"
           , pch=21
           , cex = legend_cex
           , pt.cex = legend_point_size
           , pt.bg = unique(as.vector(arrange(vertices[[line]], mean_Apoptosis))$color.apop.rel) #Colours to fill legend with
           , legend=unique(as.vector(arrange(vertices[[line]], mean_Apoptosis))$mean_Apoptosis)
    ) 
    params = lgnd$rect
    rect(xleft = params[['left']],
         ybottom = params[['top']] - params[['h']],
         xright = params[['left']] + params[['w']],
         ytop = params[['top']] + 0.15) 
    legend(
      "topright"
      , title = "Mutation"
      ##, bty = "n"
      , pch = "-"
      , pt.cex = legend_line_size
      , cex = legend_cex
      , pt.bg = unique(E(data)$colour) #add colour to non-line symbol in legend
      , col = unique(E(data)$colour) #add colour to line symbol in legend
      , legend = gsub("[\n]", "", unique(E(data)$change))
    )
  }
  
  if(pheno == "prolif"){
    plot(data
         , vertex.color = V(data)$color.prolif.rel
         #, layout = layout_as_tree(net$`SKBR3`, root = "EGF__1__Wnt1__1_AT_0")
         , layout = layout.sugiyama(data, hgap = 2, vgap = 10)$layout
         , vertex.label = NA
         , vertex.size = v.size
         , edge.label = if (label == TRUE) {E(data)$change} else {NA}
         , edge.label.family = "sans"
         , edge.width = edge_width
         , edge.arrow.width = arrow_width
         , edge.arrow.size = arrow_size
         , edge.label.color = "black"
         , edge.label.cex = e.label.size
         , edge.label.dist = label_dist
         , edge.color = E(data)$colour
         #, main = paste0("Proliferation in ", line)
    )
    lgnd = legend("topleft"
           , title = "Mean\n Proliferation"
           , bty = "n"
           , pch=21            
           , cex = legend_cex
           , pt.cex = legend_point_size
           , pt.bg = unique(as.vector(arrange(vertices[[line]], mean_Proliferation))$color.prolif.rel) #Colorus to fill legend with
           , legend=unique(as.vector(arrange(vertices[[line]], mean_Proliferation))$mean_Proliferation)
    ) 
    params = lgnd$rect
    rect(xleft = params[['left']],
         ybottom = params[['top']] - params[['h']],
         xright = params[['left']] + params[['w']],
         ytop = params[['top']] + 0.15) 
    legend(
      "topright"
      , title = "Mutation"
      #, bty = "n"
      , pch = "-"
      , col = unique(E(data)$colour) #add colour to line symbol in legend
      , pt.cex = legend_line_size
      , cex = legend_cex
      , pt.bg = unique(E(data)$colour)
      , legend = gsub("[\n]", "", unique(E(data)$change))
    )
  }
}

plot.abs <- function(line, pheno, inet, vertices, label = TRUE, v.size = 15, e.label.size = edge_label_size){
  # Function to plot absolute graphs 
  # 
  # Args: 
  #      line = cell line to plot
  #      pheno = phenotype to plot, either apop or prolif
  #      inet = list of igraph net objects
  #      vertices = list of vertices for above objects
  #      label = Plot edge labels or not
  #      v.size = vertex size
  #      e.label.size = edge label size
  # 
  # Returns a plot
  
  data <- inet[[line]]
  
  if(pheno == "apop"){
    plot(data
         , vertex.color = V(data)$color.apop.absolute
         , layout = layout.sugiyama(data, hgap = 2, vgap = 10)$layout
         , vertex.label = NA
         , vertex.size = v.size
         , edge.label = if (label == TRUE) {E(data)$change} else {NA}
         , edge.label.family = "sans"
         , edge.width = edge_width
         , edge.arrow.width = arrow_width
         , edge.arrow.size = arrow_size
         , edge.label.color = "black"
         , edge.label.cex = e.label.size
         , edge.label.dist = label_dist
         , edge.color = E(data)$colour
    )
    lgnd = legend("topleft"
           , title = "Mean\n Apoptosis"
           , bty = "n"
           , pch=21            
           , cex = legend_cex
           , pt.cex = legend_point_size
           , pt.bg = unique(as.vector(arrange(vertices[[line]], mean_Apoptosis))$color.apop.absolute) #Colors to fill legend with
           , legend=unique(as.vector(arrange(vertices[[line]], mean_Apoptosis))$mean_Apoptosis)
    ) 
    params = lgnd$rect
    rect(xleft = params[['left']],
         ybottom = params[['top']] - params[['h']],
         xright = params[['left']] + params[['w']],
         ytop = params[['top']] + 0.15) 
    legend(
      "topright"
      , title = "Mutation"
      , pch = "-"
      , col = unique(E(data)$colour) #add colour to line symbol in legend
      , pt.cex = legend_line_size
      , cex = legend_cex
      , pt.bg = unique(E(data)$colour)
      , legend = gsub("[\n]", "", unique(E(data)$change))
    )
  }
  
  if(pheno == "prolif"){
    plot(data
         , vertex.color = V(data)$color.prolif.absolute
         , layout = layout.sugiyama(data, hgap = 2, vgap = 10)$layout
         , vertex.label = NA
         , vertex.size = v.size
         , edge.label = if (label == TRUE) {E(data)$change} else {NA}
         , edge.label.family = "sans"
         , edge.width = edge_width
         , edge.arrow.width = arrow_width
         , edge.arrow.size = arrow_size
         , edge.label.color = "black"
         , edge.label.cex = e.label.size 
         , edge.label.dist = label_dist
         , edge.color = E(data)$colour
    )
    lgnd = legend("topleft"
           , title = "Mean\n Proliferation"
           , pch=21  
           , bty = "n"
           , cex = legend_cex
           , pt.cex = legend_point_size
           , pt.bg = unique(as.vector(arrange(vertices[[line]], mean_Proliferation))$color.prolif.absolute) #Colors to fill legend with
           , legend=unique(as.vector(arrange(vertices[[line]], mean_Proliferation))$mean_Proliferation)
    ) 
    params = lgnd$rect
    rect(xleft = params[['left']],
         ybottom = params[['top']] - params[['h']],
         xright = params[['left']] + params[['w']],
         ytop = params[['top']] + 0.15) 
    legend(
      "topright"
      , title = "Mutation"
      , pch = "-"
      , col = unique(E(data)$colour) #add colour to line symbol in legend
      , pt.cex = legend_line_size
      , cex = legend_cex
      , pt.bg = unique(E(data)$colour)
      , legend = gsub("[\n]", "", unique(E(data)$change))
    )
  }
}

scale.for.color <- function(df, col, round.to = 2){
  # Scale values to integers so can be used to colour vertices or edges
  # 
  # Args: 
  #      df:  dataframe to be acted on
  #      col: column to scale
  #      round.to: need to round decimals off, as cannot multiply up a recurring number. This is the number of decimal places to round to
  # Return:
  #       data frame with one column rounded to `round.to` decimal places and then multiplied up to be an integer. Also add 1 as this is needed to use as index for a vector of colours 
  
  new.col <- paste0("scaled_", col)
  df[[new.col]] <- (10**round.to)*round(df[[col]], digits = round.to)+1
  return(df)
}

visualise <- function(expanded.summary.by.mutation.attractor, cell.line.names) {
  verts <- lapply(expanded.summary.by.mutation.attractor, make.verts.table) #make vert table for all cell lines
  edges <- lapply(verts, make.edge.list) # Make edge list for all cell lines
  edges <- lapply(edges, function(x) {x$colour <- Edge.col(x, pal = edge.pal); return(x)})
  # Replace "@" symbol
  edges <- lapply(edges, function(x) {x$change <-  gsub("@2", "\n Act.", x$change); return(x)})
  edges <- lapply(edges, function(x) {x$change <-  gsub("@0", "\n Inact.", x$change); return(x)})
  verts <- lapply(verts, scale.for.color, col = "mean_Proliferation", round.to = 2) #Scale list of verts for colour for proliferation
  verts <- lapply(verts, scale.for.color, col = "mean_Apoptosis", round.to = 2) #Scale list of verts for colour for apoptosis
  verts.max.prolif <- lapply(verts, function(x) {max(x[["scaled_mean_Proliferation"]])}) #Find maximum value of scaled colours for each df of vertices
  verts.max.apop <- lapply(verts, function(x) {max(x[["scaled_mean_Apoptosis"]])}) #Find maximum value of scaled colours for each df of vertices
  pal.prolif.rel <- lapply(verts.max.prolif, function(x) {pal <- colorRampPalette(c("white", blue.pal))(n = x); return(pal)}) #Makes a colour palette where red is the highest value of apoptosis for THAT cell line
  pal.prolif.absolute <- colorRampPalette(c("white", blue.pal))(n = max(unlist(verts.max.prolif))) #Makes a colour palette where red is the highest value of apoptosis for ANY cell line
  # Make Apoptosis colour palettes
  pal.apop.rel <- lapply(verts.max.apop, function(x) {pal <- colorRampPalette(c("white", red.pal))(n = x); return(pal)}) #Makes a colour palette where red is the highest value of apoptosis for THAT cell line
  pal.apop.absolute <- colorRampPalette(c("white", red.pal))(n = max(unlist(verts.max.apop))) #Makes a colour palette where red is the highest value of apoptosis for ANY cell line
  verts <- mapply(function(x,y) {x$color.prolif.rel <- y[x$scaled_mean_Proliferation]; return(x)}, x = verts, y = pal.prolif.rel, SIMPLIFY = FALSE) #Must apply BEFORE conversion to igraph object, as cannot loop over it easily due to needing to access nodes with function V(net)
  verts <- lapply(verts, function(x) {x$color.prolif.absolute <- pal.prolif.absolute[x$scaled_mean_Proliferation]; return(x)}) #Must apply BEFORE conversion to igraph object, as cannot loop over it easily due to needing to access nodes with function V(net)
  verts <- mapply(function(x,y) {x$color.apop.rel <- y[x$scaled_mean_Apoptosis]; return(x)}, x = verts, y = pal.apop.rel, SIMPLIFY = FALSE) #Must apply BEFORE conversion to igraph object, as cannot loop over it easily due to needing to access nodes with function V(net)
  verts <- lapply(verts, function(x) {x$color.apop.absolute <- pal.apop.absolute[x$scaled_mean_Apoptosis]; return(x)}) #Must apply BEFORE conversion to igraph object, as cannot loop over it easily due to needing to access nodes with function V(net)
  
  # Make a list of igraph network objects for each cell line
  net <- mapply(graph_from_data_frame, d = edges, vertices = verts, list(directed = TRUE), SIMPLIFY = FALSE)
  
  cell.line.names.proper <- cell.line.names [! cell.line.names %in% "Background"] # Remove background as not a network to plot, only one mutation
  
  heightf = 7 # Factor to multiply resolution by for height
  widthf = 11 # Factor to multiply resolution by for height
  
  lapply(cell.line.names.proper, function(x) {to.dev(plot.rel(x, "apop", net, verts), png, file.path(paste0(x, "_apop_rel.png")), width = 300*widthf, height = 300*heightf, res = 300, pointsize = 8)})
  lapply(cell.line.names.proper, function(x) {to.dev(plot.rel(x, "prolif", net, verts), png, file.path(paste0(x, "_prolif_rel.png")), width = 300*widthf, height = 300*heightf, res = 300, pointsize = 8)})
  lapply(cell.line.names.proper, function(x) {to.dev(plot.abs(x, "apop", net, verts), png, file.path(paste0(x, "_apop_abs.png")), width = 300*widthf, height = 300*heightf, res = 300, pointsize = 8)})
  lapply(cell.line.names.proper, function(x) {to.dev(plot.abs(x, "prolif", net, verts), png, file.path(paste0(x, "_prolif_abs.png")), width = 300*widthf, height = 300*heightf, res = 300, pointsize = 8)})
  
  # Write out pdf plots
  lapply(cell.line.names.proper, function(x) {to.dev(plot.rel(x, "apop", net, verts), pdf, file.path(paste0(x, "_apop_rel.pdf")), width = widthf, height = heightf, pointsize = 8)})
  lapply(cell.line.names.proper, function(x) {to.dev(plot.rel(x, "prolif", net, verts), pdf, file.path(paste0(x, "_prolif_rel.pdf")), width = widthf, height = heightf, pointsize = 8)})
  lapply(cell.line.names.proper, function(x) {to.dev(plot.abs(x, "apop", net, verts), pdf, file.path(paste0(x, "_apop_abs.pdf")), width = widthf, height = heightf, pointsize = 8)})
  lapply(cell.line.names.proper, function(x) {to.dev(plot.abs(x, "prolif", net, verts), pdf, file.path(paste0(x, "_prolif_abs.pdf")), width = widthf, height = heightf, pointsize = 8)})
}
