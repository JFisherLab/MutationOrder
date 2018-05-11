library(magrittr)
library(reshape2)
library(jsonlite)

remove.na.lines <- function(df) df[rowSums(is.na(df)) != ncol(df), ]

load.cell.lines <- function(mutations.filename) {
    read.csv(file=mutations.filename, check.names=F, stringsAsFactors=F, na.strings="") %>%
    remove.na.lines
}

space.warn <- function(name) {
    if (grepl(" ", name)) stop("Error: node names cannot contain spaces")
    if (grepl("__", name)) stop("Error: node names cannot contain double underscore")
}

pullID <- function(filename) {
    variables <- fromJSON(readLines(filename))$Model$Variables
    idnames <- variables[c("Id", "Name")]
    sapply(idnames$Name, space.warn)
    idnames
}

cell.line.w.range <- function(input.row, id) {
    # Converts data frame of all the cell lines with their mutations into list of data frames of the nodes to be mutated for each, with id, ko and range
    # Args: data frame of rows of cell lines and columns of nodes
    # Returns: list of data frames

    pert <- melt(input.row) 
    pert$Name <- rownames(pert)
    colnames(pert) <- c("Val", "Name") 
    na.omit(pert) %>% 
    merge(id, by="Name") #Add corresponding Id names matched with Node names by merging idnames
}

comb.upto.n <- function(collapse) {
    function(x) {
        comb <- function(i) combn(x, i) %>% apply(2, function(x) paste(x, collapse=collapse))
        sapply(1:length(x), comb)
    }
}

copy.background <- function(cell.line) {
    background.files <- list.files("Background", "(Attractor|Fixpoints)", full.names = TRUE)
    file.copy(background.files, cell.line)
}

load.all.attractors <- function(cell.line) {
    function(mutation.name) {
        file.names <- list.files(path=cell.line,
                                 pattern=paste0("Mut", mutation.name, "(Attractor|Fixpoints)"))
        dfs <- lapply(file.path(cell.line, file.names), read.csv)
        do.call(rbind, dfs) %>% unique
    }
}

select.odd.elements <- function(v) v[c(T,F)]
extract.mutation <- function(x) paste("-ko", trimws(x$Id), trimws(x$Val))
extract.mutation.name <- function(x) paste(trimws(x$Name), trimws(x$Val), sep="__")
extract.genes <- function(x) strsplit(x, "__") %>% unlist %>% select.odd.elements

generate <- function(network.filename, mutations.filename, temp.bma.path) {
    idnames <- pullID(network.filename)
    cell.lines <- load.cell.lines(mutations.filename)
    cell.lines.range <- apply(cell.lines, 1, cell.line.w.range, id=idnames)
    names(cell.lines.range) <- cell.lines[["Cell_Line"]]
    mutations <- lapply(cell.lines.range, extract.mutation)
    mutation.names <- lapply(cell.lines.range, extract.mutation.name)

    background.mutations <- paste0(mutations$Background, collapse=" ")
    mutations$Background <- NULL
    background.mutation.names <- paste0(mutation.names$Background, collapse="__")
    mutation.names$Background <- NULL
    mutations <- lapply(mutations, comb.upto.n(" ")) %>%
        lapply(function(l) lapply(l, function(x) paste(background.mutations, x)))

    mutation.names <- lapply(mutation.names, comb.upto.n("__")) %>%
        lapply(function(l) lapply(l, function(x) paste(background.mutation.names, x, sep="__")))
    
    dir.create("Background/", showWarnings=F)
    system(paste0(temp.bma.path, " -model ", network.filename,
                  " -engine ATTRACTORS ", is.async, " -out Background/Mut", background.mutation.names,
                  " ", background.mutations), intern=T) %>% print
    temp.csv <- "temp.csv"
    
    for (cell.line in names(mutations)) {
        dir.create(cell.line, showWarnings=F)
        copy.background(cell.line)
        write.csv(load.all.attractors("Background")(background.mutation.names), temp.csv, row.names=F)

        muts <- mutations[[cell.line]][[1]]
        names <- mutation.names[[cell.line]][[1]]
        
        for (j in 1:length(muts)) {
            system(paste0(temp.bma.path, " -model ", network.filename,
                          " -engine ATTRACTORS ", is.async, " -initial ", temp.csv,
                          " -out ", cell.line, "/Mut", names[j],
                          " ", muts[j]), intern=T) %>% print
        }

        for (i in 2:length(mutations[[cell.line]])) {
            muts <- mutations[[cell.line]][[i]]
            names <- mutation.names[[cell.line]][[i]]
            for (j in 1:length(muts)) {
                genes = extract.genes(names[j])
                idx <- mutation.names[[cell.line]][[i-1]] %>%
                       sapply(function(x) all(extract.genes(x) %in% genes))
                parents <- mutation.names[[cell.line]][[i-1]][idx]
                do.call(rbind, lapply(parents, load.all.attractors(cell.line))) %>%
                    write.csv(temp.csv, row.names=F)
                system(paste0(temp.bma.path, " -model ", network.filename,
                              " -engine ATTRACTORS ", is.async, " -initial ", temp.csv,
                              " -out ", cell.line, "/Mut", names[j],
                              " ", muts[j]), intern=T) %>% print
            }
        }
    }
}
