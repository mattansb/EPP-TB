%{
make new epp_loadeeglab

input is the address of multi .set files - if empty, pop gui input (select
set? or just file selector).

don't load, instead use:(?)
S = whos('-file',fullfile(File.Path, File.Name{s}));
find onw with mose channels, then interpolate if needed?
(or just assume all are the same - reorder if needed?
or import with chanlocs?)

id from EEG.subject
condition from EEG.condition

if nargout==2, also give split with rels?
option for ERP, ERSP+ITC, or both - to be used with wavelet_conv2?

[ERPs, ERPsR] = epp_loadeeglab(setlist\ALLEEG(\STUDY?),'ERP',true,'WAVE',true)
%}