# Mutation library

TB-Profiler ships with a default database. The development of the mutation library is hosted on the [tbdb repository](https://github.com/jodyphelan/tbdb). Please visit this repo if you would like to get involved in the database or would like to modify and create your own.

## Why is there a separate github repository?

With analysis pipelines pretty much standardised, it is evident that accuracy of prediction is affected mostly by the underlying library of mutations. As new evidence for the inclusion or exclusion of mutations is generated, there is a constant need to update and re-evaluate the mutation library. Moreover, it is important for the control of the library to be put in the hands of the end-users. By hosting the library on a separate repository (rather than buried deep in the profiling tool code) it makes it easier to find out exactly which mutations are present. Additionally, github has a number of useful features which can be utilised:

* All changes to the library are tracked and can be easily be reverted to previous versions.
* Users can raise concerns or discuss the library using the [Issues](https://github.com/jodyphelan/tbdb/issues) section of the github repository.