library(dplyr) # Version 0.7.0 or above
library(reshape2)
library(purrr)

pull.pheno <- function(df, keep){
    # Pull out only certain behaviours from data frame
    # 
    # Args
    #     df: data frame where you want to subset out certain columns
    #     keep: columns to keep, of the format "A|B|C"
    # Returns: 
    #         df_short: data frame with only the columns that have been specified in keep
  
    if(length(grep(keep, colnames(df))) == 0){
        return(df)
    } else {
        df.short <- df[, grep(keep, colnames(df))]
        return(df.short)
    }
}

read.direc <- function(cell) {
    # Read all the csv files in a directory into a list of data.frames
    # Args:
    #      path = path to the root directory containing directories for each cell line. Likekly `results.direc`
    #      cell = cell line for which to read results e.g. `"SKBR3"`
    #      head = whether or not to use headers
    # Return: 
    #      list of data frames, one for each csv file for the cell line                         
    
    file.names.direc <- list.files(path=file.path(cell), pattern="*(Attractor|Fixpoints)")
    if (length(file.names.direc) == 0) {
        stop(paste0("Error: Zero results files in directory: \n", file.path(path, cell),
                    "\n Have you specified the directory correctly?"))
    }
    data <- lapply(file.path(cell, file.names.direc), read.csv)
    names(data) <- gsub("*.csv$", "", file.names.direc)
    data
}

read.cell.line.results <- function(cell.lines) {
    # Read the results for all the cell lines, and store as a list of list of data frames
    # Args:
    #      cell.lines = list of the cell lines for which there are directories of results
    #      path = path to directory of results
    # Returns:
    #      list of list of data frames of results, root list is of cell lines, leaves are of the stable states for each of those cell lines
    
    data <- lapply(cell.lines, read.direc)
    names(data) <- cell.lines
    data
}

name.results <- function(l, colname) {
    # Takes the name of a variable and turns it into the value of a column, with specified columnme, to allow for later merging without losing that info
    # Args:
    #      l = a list of data frames with each data frame named such that names(l) is meaningful
    #      colname = name for column containing name data
    # Returns: 
    #         same list as input but each data frame has a column `colname` which contains the name of each data frame
    
    mapply(function(x, y) {x[colname] <- y; x}, x=l, y=names(l), SIMPLIFY=F)
}

expand.range <- function(x) {
    ## if x is a character, interpreted as a range. Otherwise, as a scalar
    ## Expand ranges, e.g. "[1; 2; 3]", to a vector of numbers
    if (is.character(x) && grepl("\\[", x)) {
        nums <- strtoi(trimws(unlist(strsplit(substr(x, 2, nchar(x)-1), ";"))))
        return(nums)
    }
    return(as.integer(x))
}

expand.ranges <- function(df) {
    lists <- apply(df, 1, function(r) lapply(r, expand.range))
    dfs <- lapply(lists, cross_df)
    bind_rows(dfs)
}

expand.fixpoints <- function(l) {
    fixpoints <- grepl("Fixpoints", names(l))
    results <- l[!fixpoints]
    num.loops <- length(results)
  
    for (f in names(l)[fixpoints]) {
        expansion <- expand.ranges(l[[f]])
        expansion <- split(expansion, seq(nrow(expansion)))
        names(expansion) <- sapply(num.loops - 1 + seq(length(expansion)),
                                   function(i) paste(f, i, sep=""))
        results <- append(results, expansion)
    }
    names(results) <- gsub("Fixpoints", "Attractor", names(results))
    results
}

expand.mutations <- function(df){
  # Takes a data frame with Mutation as a single identifier and splits. 
  # Dependacies: dplyr
  # 
  # Args:
  #      df: data frame, one in which the Mutation variable is in a column, and has already had Attractors split off, e.g. summary.by.mutation.attractor
  # Returns:
  #        data frame with a set of columns for each node and mutation value, as well as the number of mutations applied, for each attractor. Should be able to be used to make both the node and edge table for visualisation
  
    out <- strsplit(as.character(df$Mutation), split = "__")
    number.muts <- unlist(lapply(out, function(x) length(x)/2)) #Number of mutations applied, will be needed to build edge list
    out <- lapply(out, FUN=function(x) c(unlist(x), rep(NA, max(lengths(out))-length(x))))
    out <- as.data.frame(cbind(do.call(rbind, out)))
  
    n.newcol <- ncol(out)
    Node.names <- paste0("Node", seq(1:(n.newcol/2)))
    Mut.names <- paste0("Mut", seq(1:(n.newcol/2)))
    new.colnames <- c(rbind(Node.names, Mut.names))
    colnames(out) <- new.colnames
  
    df <- cbind(out, unclass(df)) #Needed to remove effects of `summarise` which prevents `cbind` and `dply::arrange` from working
    df$Number.Muts <- number.muts
    df <- select(df, Cell_Line, Mutation, Attractor, Number.Muts, everything()) 
    df <- arrange(df, Number.Muts, Attractor)
    df
}

analyse = function(results) {
    short.results <- lapply(results, function(x) lapply(x, expand.ranges)) %>%
                     lapply(function (x) lapply(x, pull.pheno, "Proliferation|Apoptosis")) #Pull out the nodes that affect how well the mutations are going

    short.results.muts <- lapply(short.results, name.results, colname = "Full_Mutation_Name") # Add mutations from name of data frames to a column
    short.results.muts.unlist <- lapply(short.results.muts, function(x) Reduce(function(df1, df2) merge(df1, df2, by = c("Proliferation", "Apoptosis", "Full_Mutation_Name"), all = TRUE), x)) #
    short.results.named <- name.results(short.results.muts.unlist, colname = "Cell_Line") %>% # Add cell_line from name of data frames to a column
                           lapply(mutate, Full_Mutation_Name = gsub("^Mut", Full_Mutation_Name, replacement = "")) #Remove the `"Mut"` prefix only from the very begining of a string on column `Full_Mutation_Name`
    # Split Mutations into mutation and attractor. so before here we need to filter
    short.results.named.split <- lapply(short.results.named, function(x) {y <- cbind(x, colsplit(x[["Full_Mutation_Name"]], pattern = "Attractor", c("Mutation", "Attractor"))); return(y)})

    ## # Group by mutations AND THEN by attractor
    by.mutation.attractor <- lapply(short.results.named.split, group_by, Mutation, Attractor)
    summary.by.mutation.attractor <- lapply(by.mutation.attractor, summarise, 
                                            min_Proliferation = min(Proliferation), 
                                            max_Proliferation = max(Proliferation), 
                                            mean_Proliferation = mean(Proliferation), 
                                            range_Proliferation = max(Proliferation)-min(Proliferation), 
                                            min_Apoptosis = min(Apoptosis), 
                                            max_Apoptosis = max(Apoptosis), 
                                            mean_Apoptosis = mean(Apoptosis), 
                                            range_Apoptosis = max(Apoptosis)-min(Apoptosis),
                                            loop_length = n()) %>%
                                     name.results(colname = "Cell_Line") #need to add cell_line name as lost on the `summarise` process

    expanded.summary.by.mutation.attractor <- lapply(summary.by.mutation.attractor, expand.mutations)
    expanded.summary.by.mutation.attractor
}
