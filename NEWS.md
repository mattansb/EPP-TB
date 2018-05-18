# EPP-TB 3.1.2 +

* `epp_plotbutterfly` - the `trace` parameter replaces `all` parameter.
* `epp_plotchannels` - New function to plot ERPs per-channel in either a grid, or according to channel locations (Essentially, this is a stripped down version of `eeglab`'s `plottopo` function).
* `epp_get*` functions export the identification column with the correct name (*ID*).
* `epp_plot*` functions are faster and more efficiant when exporting plot data to `.csv`.