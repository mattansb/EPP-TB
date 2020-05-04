# EPP-TB NEWS

This doc details user-facing changes only.

## Version 3.1.x

| Function | Update | Notes |
|---------:|:-------|:------|
| `epp_exportResults` | NEW | New function to export results produced by any `epp_get*` |
| `f_jackknife` | IMPROVEMENT | Jackknifing procedure re-writen to allow for weighting and better centering. See `help epp_getAmp` or `help epp_getLat` for more details |
|`f_makeColormap`| IMPROVEMENT | Produces better colors  |
|`epp_getLat` | IMPROVEMENT | in criterion based latency measures, it is now possible to measure not only the first crossing of the the criterion (the onset) but also the last crossing (of the offset) |
| `m_peak` | IMPROVEMENT | Replace tsmovavg with movmean for computing moving average | 

## Version 3.1.3

| Function | Update | Notes |
|---------:|:-------|:------|
|`epp_plotchannels`| NEW | Plot ERPs per-channel in either a grid, or according to channel locations (Essentially, this is a stripped down version of `eeglab`'s `plottopo` function).|
|`epp_filter_by` | NEW | Filter data in study by some variable in `study.IDs` table.|
|`epp_matchsubjects` | NEW | Use to remove subjects that are missing data in some conditions.|
|`epp_extractIDs` | NEW | Extract (and save) data from `study.IDs` tables. |  
|`epp_combineconds` | IMPROVEMENT | Supports combining conds based on REGEXP.|
|`epp_plotbutterfly` | IMPROVEMENT | The `trace` parameter replaces `all` parameter.|
|`epp_get*` | IMPROVEMENT | All `epp_get*` functions now export the identification column with the correct name (*ID*).|
|`epp_plot*` | IMPROVEMENT | All plot functions are faster and more efficient when exporting plot data to `.csv`.|
|`f_WaveletConv`|FIX| When creating wavelets via `f_WaveletConv` (formerly `suppWaveletConv`), baseline correcting no-longer stops on error.|
|`epp_getTF` | FIX | Fixed bug in `epp_getTF` that caused error when conditions were not of equal size N.|
