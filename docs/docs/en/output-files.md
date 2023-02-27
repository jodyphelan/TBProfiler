# Output files

## Output formats


## Generating summary files

The results from numerous runs can be collated into one table using the following command:

```
tb-profiler collate
```

This will automatically create a number of colled result files from all the individual result files in the result directory. If you would like to generate this file for a subset of the runs you can provide a list with the run sames using the `--samples` flag. The prefix for the output files is tbprofiler by default but this can be changed with the `--prefix` flag.

### Writing your own summary scripts

The collate function extracts the drug-resistance mutations and lineage, however you may want to extract more features that are present in the individual json result files. I have created a little tutorial on how to do this [here](https://jodyphelan.gitbook.io/tb-profiler/writing-a-custom-collate-script).

### iTOL files

Several files are produced by the `tb-profile collate` function. Among these are several config files that can be used with [iTOL](http://itol.embl.de/). to annotate phylogenetic trees. A small tree and config files have been placed in the example_data directory. To use navigate to the iTOL website and upload the tbprofiler.tree file using the upload button on the navigation bar. Once this has been uploaded you will be taken to a visualisation of the tree. To add the annotation, drag and drop the file onto the tree in the browser. You should now see a figure similar to the one below. The following annotations are included:

* Lineage
* Drug resistance classes (Sensitive, drug-resistant, MDR, XDR)
* Drug resistance calls for individual drugs, were filled circles represent resistance.

<img href="https://github.com/jodyphelan/TBProfiler/raw/docs/docs/docs/assets/images/itol_example.png">