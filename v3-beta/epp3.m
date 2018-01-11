%{
import
epp_loadeeglab(filelist,'ERP',vararg,'WAVE',vararg)
epp_loaderplab?
epp_loadegimat(filelist,..?)

data reduction functions:
epp_*(set,conditions,varargin)

plotting: (all should be executable to R)
epp_plot*(set,type,conditions,channels,varargin)
- grand
    - ERP - errors
    - WAVE (then will plot  square, no errors)
- butterfly
    - ERP
    - WAVE
- topo (+ chanlocs)
    - ERP
    - WAVE - also freq band(s?)

measuring: (all should be executable to R). all should also have a pop up
when nargin < 3?
epp_measure*(set,type,...?)
%}