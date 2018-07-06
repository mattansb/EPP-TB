# EPP-TB NEWS

## Version 3.1.2 +

| Function | Update | Notes |
|---------:|:-------|:------|
|`epp_plotchannels`| NEW | Plot ERPs per-channel in either a grid, or according to channel locations (Essentially, this is a stripped down version of `eeglab`'s `plottopo` function).|
|`epp_filter_by` | NEW | Filter data in study by some variable in `study.IDs` table.|
|`epp_matchsubjects` | NEW | Use to remove subjects that are missing data in some conditions.|
|`epp_combineconds` | IMPROVEMENT | Supports combining conds based on REGEXP.|
|`epp_plotbutterfly` | IMPROVEMENT | The `trace` parameter replaces `all` parameter.|
|`epp_get*` | IMPROVEMENT | All `epp_get*` functions now export the identification column with the correct name (*ID*).|
|`epp_plot*` | IMPROVEMENT | All plot functions are faster and more efficiant when exporting plot data to `.csv`.|
|`f_WaveletConv`|FIX| When creating wavelets via `f_WaveletConv` (formerly `suppWaveletConv`), baseline correcting no-longer stops on error.|
`epp_getTF` | FIX| Fixed bug in `epp_getTF` that caused error when conditions were not of equal size N.|
