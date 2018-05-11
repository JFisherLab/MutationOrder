bma.path <- "C:\\Users\\t-stwood\\Desktop\\BMApush\\BioModelAnalyzer\\src\\BioCheckConsole\\bin\\x64\\Release\\BioCheckConsole.exe"

source("generation.R")
source("analysis.R")
source("visualisation.R")
library(magrittr)

args <- commandArgs(trailingOnly=T)
if (length(args) < 2) {
    stop("Error: too few arguments.")
}
if (length(args) > 3) {
    stop("Error: too many arguments.")
}

network.filename <- args[1]
mutations.filename <- args[2]
is.async <- ""
if (length(args) == 3 && args[3] == "-async") {
    is.async <- "-async"
}

generate(network.filename, mutations.filename, bma.path)

cell.lines <- load.cell.lines(mutations.filename)[["Cell_Line"]]

results <- read.cell.line.results(cell.lines) %>%
           lapply(expand.fixpoints)

summary <- analyse(results)

visualise(summary, cell.line.names = cell.lines)
