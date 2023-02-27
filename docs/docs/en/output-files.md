# Output files

## Output formats

### Default output 

By default a `.json` formatted output file is produced in a directory called `results`. This file contains information including mutations found, lineage as well as some QC metrics. This format is perfect to load into script for downstream processing but it isn't very human-readable. 

### Text outputs

For a more human-readlable format you can use the `--csv` and `--txt` flags to generate outputs in csv and text format too. These files will contain several tables listing similar information to the json format. For advanced users, it is possible to customise this format by providing a template file. This template should be written in the [jinja templating language](https://jinja.palletsprojects.com/en/3.1.x/). A variable named `d` will be available within the template. This variable contains the exact same structure as the json file. For example, a simple template could be:

```
Custom TB-Profiler report

Lineage: {{d['sublin']}}
Drug-resistance: {{d['drtype']}}
```

### Docx output

It is also possible produce a nice-looking docx formatted report that can be viewed in Word or converted into pdf. The advantage of this is that it can contain images, text formatting, etc. To create this report, a template file must be provided. As with the text custom format, the docx template should have jinja variables defined that will be filled in with sample data when a report is generated. Details on the available variabled for the templating engine can be found [here](https://github.com/jodyphelan/TBProfiler/blob/master/tbprofiler/docx.py).

![example_report](https://github.com/jodyphelan/TBProfiler/raw/docs/docs/assets/img/report_example.png)

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