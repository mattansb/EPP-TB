% PURPOSE:  exports epp-tb data to ERPLAB's ALLERP.
%
% FORMAT
% ------
% epp_erplab_export(ERPs)
% 
%
% INPUTS
% ------
% ERPs      - structure built by epp_load OR epp_erplab_import.
%
%
% See also epp_load & epp_erplab_import
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
19-01-2017  Added chanlocs support for GSNHydroCel129 nets
25-11-2016  New function (written in MATLAB R2015a)
%}


function epp_erplab_export(ERPs)

%% Find for each subject where they have conditions:

% make one subject list
IDlist = ERPs(1).IDs;
for c = 2:length(ERPs)
    IDlist = outerjoin(IDlist,ERPs(c).IDs,'MergeKeys',true);
end

% for each find ind in each condition
for c = 1:length(ERPs)
    for s = 1:height(IDlist)
        try
            [~,~,temp.(ERPs(c).Condition)(s)] = innerjoin(IDlist(s,1),ERPs(c).IDs);
        catch
            temp.(ERPs(c).Condition)(s) = nan;
        end
    end
    IDlist = [IDlist, array2table(temp.(ERPs(c).Condition)', 'VariableNames',{ERPs(c).Condition})];
end


IDlist.(ERPs(c).Condition)(s)
%% Build the file

ALLERP = evalin('base', 'ALLERP');

for s = 1:height(IDlist)
    erpset.erpname      = IDlist.ID{s};
    erpset.filename     = '';
    erpset.filepath     = '';
    erpset.workfiles    = {}; % <<<<<<<<<<<<<< okay to leave empty?
    erpset.subject      = [];
    erpset.nchan        = size(ERPs(1).Data,1); 
    erpset.nbin         = sum(~isnan(table2array(IDlist(s,2:end))));
    erpset.pnts         = size(ERPs(1).Data,2);
    erpset.srate        = ERPs(1).samplingRate;
    erpset.xmin         = []; % <<<<<<<<<<<<<<
    erpset.xmax         = []; % <<<<<<<<<<<<<<
    erpset.times        = ERPs(1).timeLine;
    
    condsB = ~isnan(table2array(IDlist(s,2:end)));
    for c = 1:width(IDlist)-1
        if condsB(c) && c==1
            erpset.bindata(:,:,1) = ERPs(c).Data(:,:,IDlist.(ERPs(c).Condition)(s)); % arranged: (electrode,time,CONDITION)
        elseif condsB(c)
            erpset.bindata(:,:,[end+1]) = ERPs(c).Data(:,:,IDlist.(ERPs(c).Condition)(s)); % arranged: (electrode,time,CONDITION)
        end
    end
    erpset.xmin = min(erpset.bindata(1:end));
    erpset.xmax = max(erpset.bindata(1:end));
    
    
    
    erpset.binerror = [];
    erpset.datatype = 'ERP';
    
    % ntrials
    erpset.ntrials.accepted = ones(1,erpset.nbin);
    erpset.ntrials.rejected = zeros(1,erpset.nbin);
    erpset.ntrials.invalid  = zeros(1,erpset.nbin);
    erpset.ntrials.arflags  = zeros(erpset.nbin,8);
    
    erpset.pexcluded    = 0;
    erpset.isfilt       = 0;
    
    % chanlocs
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    XYZ = load([Rpath '\GSNHydroCel129.mat'], '-mat');
    XYZ = XYZ.GSNHydroCel129;
    for e = 1:erpset.nchan
        erpset.chanlocs(e,1).labels       = ['e' num2str(e)];
        erpset.chanlocs(e,1).ref          = [];
        erpset.chanlocs(e,1).theta        = [];
        erpset.chanlocs(e,1).radius       = [];
        erpset.chanlocs(e,1).X            = XYZ(e,1);
        erpset.chanlocs(e,1).Y            = XYZ(e,2);
        erpset.chanlocs(e,1).Z            = XYZ(e,3);
        erpset.chanlocs(e,1).sph_theta    = [];
        erpset.chanlocs(e,1).sph_phi      = [];
        erpset.chanlocs(e,1).sph_radius   = [];
        erpset.chanlocs(e,1).type         = [];
        erpset.chanlocs(e,1).urchan       = [];
    end
    
    if erpset.nchan<129 % make last channel location CZ
        erpset.chanlocs(end,1).X = XYZ(end,1);
        erpset.chanlocs(end,1).Y = XYZ(end,2);
        erpset.chanlocs(end,1).Z = XYZ(end,3);
    end
    
    erpset.ref          = 'common';
    erpset.bindescr     = {ERPs(~isnan(table2array(IDlist(s,2:end)))).Condition}; % cell list, one cell per condition
    erpset.saved        = 'no';
    erpset.history      = '';
    erpset.version      = ''; % get ERPLAB version from somewhere
    erpset.splinefile   = '';
    
    % EVENTLIST
    erpset.EVENTLIST.setname        = '';
    erpset.EVENTLIST.report         = '';
    erpset.EVENTLIST.bdfname        = '';
    erpset.EVENTLIST.nbin           = erpset.nbin; 
    erpset.EVENTLIST.version        = ''; % get ERPLAB version from somewhere
    erpset.EVENTLIST.account        = '';
    erpset.EVENTLIST.username       = '';
    erpset.EVENTLIST.trialsperbin   = ones(1,erpset.nbin); % similar to ntrials: integer for each condition
    erpset.EVENTLIST.elname         = ''; % <<<<<<<<<<<<<< whats this?
    for c = 1:erpset.nbin % bdf
        erpset.EVENTLIST.bdf(c).description = erpset.bindescr{c};
        erpset.EVENTLIST.bdf(c).namebin     = sprintf('BIN%i',c);
    end
    erpset.EVENTLIST.eldate = datestr(datetime);
    % eventinfo % for each event - neccacery?
    erpset.EVENTLIST.eventinfo.item         = [];
    erpset.EVENTLIST.eventinfo.code         = [];
    erpset.EVENTLIST.eventinfo.binlabel     = '""';
    erpset.EVENTLIST.eventinfo.codelabel    = '""';
    erpset.EVENTLIST.eventinfo.time         = [];
    erpset.EVENTLIST.eventinfo.spoint       = [];
    erpset.EVENTLIST.eventinfo.dura         = [];
    erpset.EVENTLIST.eventinfo.flag         = [];
    erpset.EVENTLIST.eventinfo.enable       = [];
    erpset.EVENTLIST.eventinfo.bini         = [];
    erpset.EVENTLIST.eventinfo.bepoch       = [];
    
    if isempty(ALLERP)
        ALLERP = erpset;
    else
        ALLERP(end+1) = erpset;
    end
    clear erpset
end

assignin('base','ALLERP',ALLERP);
updatemenuerp(ALLERP)
end