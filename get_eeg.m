% Created on 4 Mar 2015
% Primary author: Kiefer Forseth
%                 kjforseth@gmail.com
%                 602.531.6430
%
% Summary
%   Get EEG trace from any nkdata organization.
% 
% Input
%   D: nkdata struct
%   task: task corresponding to nkdata struct (needed for audio buffer)
%   chs: logical vector (or list) of which channels to grab
%   trs: logical vector of trials to use (from get_trials)
%   fb: filter buffer, generally chosen to be 500 (samples)
%   varargin: 1st pair for traces, 2nd pair for baseline
%       center: alignment point between trials
%       win: time window relative to center
%
% Output
%   eeg: EEG traces with size (1) x (# samples)
%   varargout: 1st pair for traces, 2nd pair for baseline
%       traces & base of size (# trials) x (# samples)
%       t & bt of size (1) x (# samples)
% 
% Example
%   Complete EEG trace for channel 32:
%       eeg = get_eeg(D,'auditory',32,trs,500);
%   Aligned EEG traces for channel 8, center pulse_on, time win [-250 500]:
%       [traces,t] = get_eeg(D,'auditory',8,trs,500,'pulse_on',[-250 500]);
%   Aligned EEG traces & baseline:
%       [traces,t,base,bt] = get_eeg(D,'auditory',8,trs,500,'pulse_on',[-250 500],'pulse_on',[-750 -250]);
% 
% Change log
%   02 Apr 2015: fixed bug for trace without baseline input
%   21 Sep 2016: code doesn't break if audio_buffer field isn't present

function [varargout] = get_eeg(D,task,chs,trs,fb,varargin)

% handle chs
if islogical(chs), chs = find(chs); end

% handle varargin
if nargin == 5
    opts.tr   = false;
    opts.base = false;
elseif nargin == 7
    center.tr = varargin{1};
    time.tr   = varargin{2};
    opts.tr   = true;
    opts.base = false;
elseif nargin == 9
    center.tr   = varargin{1};
    time.tr     = varargin{2};
    center.base = varargin{3};
    time.base   = varargin{4};
    opts.tr     = true;
    opts.base   = true;
else
    error('Unexpected number of input arguments')
end

% read, align, & store eeg
Ne = length(chs);
if opts.tr,   traces = cell(Ne,1); end
if opts.base, base = cell(Ne,1); end
if ~opts.tr && ~opts.base, unaligned = cell(Ne,1); end
for e = 1:Ne
    elec = chs(e);
    
    % get eeg from all file types
    if ~isinteger(D.eeg)
        paths = get_paths;
        fname = sprintf('%s%s',paths.vol.S1,D.eeg(18:end));
        fid = fopen(fname);
        fseek(fid,D.HeaderOffset,'bof');
        fseek(fid,2*(elec-1),'cof');
        eeg = fread(fid,[1 D.num_dpoints],...
            '1*short=>short',2*(D.nchannels - 1));
        eeg = double(eeg(:)') - double(D.ref_tseries);
        fclose all;
    elseif isfield(D,'ref_tseries')
        eeg = double(D.eeg(elec,:)) - ...
                double(D.ref_tseries);
    else
        eeg = double(D.eeg(elec,:));
    end
    eeg = eeg/D.multiplier;

    % generate time and sample vectors
    T = D.ms_per_sample;

    if opts.tr
        t = (time.tr(1)     :T:time.tr(2)     );
        s = (time.tr(1)-fb*T:T:time.tr(2)+fb*T)/T;
        c = D.(center.tr);
    else
        unaligned{e} = eeg;
        t = T*(1:length(eeg));
        s = true(1,length(eeg));
    end
    Ns  = length(s);
    
    trsI = find(trs);

    if opts.base
        if strcmpi(center.base,'fixed')
            bt = (time.base(1)     :T:time.base(2)     );
            bs = (time.base(1)-fb*T:T:time.base(2)+fb*T)/T;
            bc = zeros(size(trs));
            Nbs = length(bs);
        else
            bt = (time.base(1)     :T:time.base(2)     );
            bs = (time.base(1)-fb*T:T:time.base(2)+fb*T)/T;
            bc = D.(center.base);
            Nbs = length(bs);
        end
    end

    % adjust centers for auditory buffer
    if strcmp(task,'ccep')
        ab = 0;
    elseif isfield(D,'audio_buffer')
        ab = D.audio_buffer;
    else
        ab = 0;
    end
    
    if any(strcmp(task,{'auditory' 'audscramble'})) % for auditory task, realign everything
        if opts.tr, c = c - ab; end
        if opts.base, bc = bc - ab; end
    elseif any(strcmp(task,{'perdata3_tone','perdata35'})) % for hickok white noise tasks, realign pulse_on & off 
        if opts.tr   && any(strcmp(center.tr,{'pulse_on','pulse_off','tone'})), c = c - ab; end
        if opts.base && any(strcmp(center.base,{'pulse_on','pulse_off','tone'})), bc = bc - ab; end
    elseif any(strcmp(task,{'A1localizerB','A1localizerHP','A1localizerHT' ,'A1localizerT'})) % for A1localizers, realign pulse_on
        if opts.tr   && any(strcmp(center.tr,{'pulse_on'})), c = c - ab; end
        if opts.base && any(strcmp(center.base,{'pulse_on'})), bc = bc - ab; end        
    elseif strcmp(task,'ccep')
        % do nothing
    else % by default, only realign articulation
        if opts.tr   && strcmp(center.tr,'articulation'), c = c - ab; end
        if opts.base && strcmp(center.base,'articulation'), bc = bc - ab; end
    end

    % align traces
    if opts.tr
        traces{e} = zeros(length(trsI),Ns);
        for m = 1:length(trsI)
            win = s + round(c(trsI(m)))/T;
            traces{e}(m,:) = eeg(win);
        end
        
        % detrend
        traces{e} = detrend(traces{e}','constant')';
    end

    if opts.base
        base{e} = zeros(length(trsI),Nbs);
        for m = 1:length(trsI)
            win = bs + round(bc(trsI(m)))/T;
            base{e}(m,:) = eeg(win);
        end
        
        % detrend
        base{e} = detrend(base{e}','constant')';
    end
end

% handle varargout
if opts.tr
    if (Ne == 1), traces = cell2mat(traces); end
    varargout{1} = traces;
    varargout{2} = t;
end
if opts.base
    if (Ne == 1), base = cell2mat(base); end
    varargout{3} = base;
    varargout{4} = bt;
end
if ~opts.tr && ~opts.base
    if (Ne == 1), unaligned = cell2mat(unaligned); end
    varargout{1} = unaligned;
end