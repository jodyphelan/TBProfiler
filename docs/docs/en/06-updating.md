# Updating

TB-Profiler is under constant rapid development. If you plan to use the program in your work please make sure you are using the most up to date version! Similarly, the database is not static and is continuously being improved so make sure you are using the most latest version. If you use TBProfiler in your work please state the version of both the tool and the database as they are deveoped independantly from each other.

### Updating the database

New mutations/genes are periodically added to the database. Run the following to make sure you are up to date.

```
tb-profiler update_tbdb
```

### Quick re-profiling

If you have a new mutation database but no new genes have been added you can quickly re-profile your samples by running the following.

```
tb-profiler reprofile /path/to/result.json
```

This can be useful when you have added in a few mutations yourself or you are sure that no new genes have been added in the update. If you are not sure then it is safest to run the full profiling step again.
