# Change log

## Version 7.0.0

### Added
- Added rules to make enable use with the WHO v3 catalogue
    - Intrinsic resistance based on lineage
    - Update katG coding variants confidence (from Uncertain to Assoc w R - Interim) when ahpC compensatory mutations are present
- Added new interpretation section to the json file
    - Contains R, R (low), R (high) or W (Warning) for each drug based on the interpretation rules
- Adds IS6110 intertion calling (can confer resistance)

### Changed
- Updated text report generation to include drug resistance interpretation.


### DB
Slightly customised version of WHO v3 catalogue for TB drug resistance interpretation.