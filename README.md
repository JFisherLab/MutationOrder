# MutationOrder

MutationOrder is a tool to explore how the order in which mutations are acquired in an evolving cancer is constrained by the change in cell phenotypes that these mutations cause. 
It takes a Qualitative Network (QN) model of the gene regulatory network of a cell, as built in [Bio Model Analyzer](http://biomodelanalyzer.org/) (BMA) and a list of mutations observed in cancer, and returns an analysis of the phenotypes expected at each stage of the acquisition of those mutations, in every possible sequence. 

Details of the methodology can be found in [Clarke et al., _Automated Reasoning for Systems Biology and Medicine_, 2019](https://doi.org/10.1007/978-3-030-17297-8_5).

## R scripts

MutationOrder is implemented as a set of R scripts which take as input a QN model file and a mutations csv file. The scripts call the command line interface of Bio Model Analyzer to compute model attractors, and thent outputs csv files for each of these attractors, along with images charting the predicted phenotypes at each step of cancer evolution.

### Dependencies

To run MutationOrder, you must install [R](https://www.r-project.org/) and build a local installation of [BMA](https://github.com/Microsoft/BioModelAnalyzer). BMA currently only builds on Windows.

You will also require the following R libraries:

 - tidyr
 - dplyr (Version 0.7.0 or above)
 - reshape2
 - RColorBrewer
 - igraph
 - ggplot2
 - magrittr
 - jsonlite
 - purrr
 
### Network file

Networks are represented in the `.json` format of BMA. These can be built and exported using the [BMA web interface](http://biomodelanalyzer.org/). An example is given in `example_network.json`.
 
### Mutation file format

The general format is shown below, and a concrete example is given in `example_mutations.csv`.
Columns correspond to mutated genes and rows to cancer genotypes.
The file must contain a row labelled "Background", which gives the configuration of the "healthy" initial network. Other rows should not mutate the same genes as the background row.

**Cell_Line**|**Gene1**|**Gene2**|**Gene3**|**Gene4**|**Gene5**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
Background|1|1| | |
Cancer1| | |2|3|
Cancer2| | |1 | |3
⋮

### Running MutationOrder

After installing R and the required libraries, adding Rscript.exe to your path, and building BMA (see [here](https://github.com/hallba/BioModelAnalyzer#build-and-test))*, you are ready to run the R scripts.

Firstly, edit line 1 of ``order.R`` to point to your BMA executable:

```
bma.path <- "C:\\MyPath\\BioModelAnalyzer\\src\\BioCheckConsole\\bin\\x64\\Release\\BioCheckConsole.exe"
```

Now, assuming ``Rscript.exe`` is in your PATH, you can run MutationOrder:

```
Rscript.exe order.R example_network.json example_mutations.csv
```

Adding the ``-async`` command line flag will run the model with asynchronous rather than synchronous semantics:

```
Rscript.exe order.R example_network.json example_mutations.csv -async
```

After the script has finished executing, you will find a folder for each row in ``example_mutations.csv``, containing ``.csv`` files for the attractors of the model under every combination of mutations. You will also find ``.png`` and ``.pdf`` images for each row, showing these attractors (and the paths from background to fully mutated cell phenotypes) visually.

\* Note that for most applications one can install BMA with one click ([here](http://biomodelanalyzer.org/installer/bma.install.msi)), but to use MutationOrder requires the underlying [ATTRACTOR engine](https://github.com/swoodhouse/attractors), which can only be used after a manual compilation of BMA from source ([here](https://github.com/hallba/BioModelAnalyzer#build-and-test)). 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
