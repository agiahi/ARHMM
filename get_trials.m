% Created on 4 Feb 2015
% Primary author: Kiefer Forseth
%                 kjforseth@gmail.com
%                 602.531.6430
%
% Summary
%   Get trials that fit a certain set of characteristics.
% 
% Input
%   D: nkdata struct
%   varargin: +/- determines whether to keep/remove trs by a given marker
%       all
%       correct
%       incorrect
%       unscrambled
%       scrambled
%       noise
%       rxn
%       tech
%
%   BACKWARD COMPATIBILITY NOTE:
%   varargin: list everything you do NOT want to DELETE
%       correct
%       incorrect
%       unscrambled
%       scrambled
%       noise
%       rxn
%       tech
%
% Output
%   trs: logical vector with true for good trials
% 
% Example
%   Correct only (no tech, scramble, noise, long rxns);
%       trs = get_trials(D,'+correct','-scrambled','-noise','-rxn','-tech')
%   All trials
%       trs = get_trials(D,'+all')
%   All except noise & tech
%       trs = get_trials(D,'+all','-noise','-tech')
%
%   BACKWARD COMPATIBILITY NOTE:
%   Correct only (no tech, scramble, noise, long rxns):
%       trs = get_trials(D,'correct','unscrambled');
%   Correct & incorrect:
%       trs = get_trials(D,'correct','incorrect','unscrambled');
% 
% Change log
%   % 21 Sep 2016: fundamental change to input structure
%       Previously, varargin was a list of everything to NOT DELETE
%       Afterwards, each item on the list is preceded by +/-
%       Backwards compatibility: lack of +/- assumes +all and -varargin

function trs = get_trials(D,varargin)

% handle inputs
if ~iscell(D), D = {D}; end

if iscell(varargin) && iscell(varargin{1}), varargin = varargin{1}; end

% backwards compatibility mode
if ~strcmp('-',varargin{1}(1)) && ~strcmp('+',varargin{1}(1))
%     warning('Backwards compatibility mode for get_trials')

    if any(strcmp('correct',varargin)), opts.correct = true; else opts.correct = false; end
    if any(strcmp('incorrect',varargin)), opts.incorrect = true; else opts.incorrect = false; end
    if any(strcmp('unscrambled',varargin)), opts.unscrambled = true; else opts.unscrambled = false; end
    if any(strcmp('scrambled',varargin)), opts.scrambled = true; else opts.scrambled = false; end
    if any(strcmp('noise',varargin)), opts.noise = true; else opts.noise = false; end
    if any(strcmp('rxn',varargin)), opts.rxn = true; else opts.rxn = false; end
    if any(strcmp('tech',varargin)), opts.tech = true; else opts.tech = false; end

    % calculate logical vector
    trs = cell(size(D));
    for p = 1:size(D,1)
        for t = 1:size(D,2)
            correct     = (~opts.correct     &  D{p,t}.accuracy(:));
            incorrect   = (~opts.incorrect   & ~D{p,t}.accuracy(:));
            unscrambled = (~opts.unscrambled & D{p,t}.scramble(:));
            scrambled   = (~opts.scrambled   & ~D{p,t}.scramble(:));
            noise       = (~opts.noise       & ~D{p,t}.noise(:));
            rxn         = (~opts.rxn         &  (D{p,t}.rxn_time(:) > 2000));
            tech        = (~opts.tech        & ~D{p,t}.tech(:));

            trs{p,t} = ~(correct | incorrect | unscrambled | scrambled | noise | rxn | tech);
        end
    end
else
    trs = cell(size(D));
    for p = 1:size(D,1)
        for t = 1:size(D,2)
            add = zeros(length(D{p,t}.pulse_on),1);
            rem = zeros(length(D{p,t}.pulse_on),1);
            for n = 1:length(varargin)
                switch varargin{n}(2:end)
                    case 'all'
                        val = true(length(D{p,t}.pulse_on),1);
                    case 'correct'
                        val = D{p,t}.accuracy(:);
                    case 'incorrect'
                        val = ~D{p,t}.accuracy(:);
                    case 'unscrambled'
                        val = D{p,t}.scramble(:);
                    case 'scrambled'
                        val = ~D{p,t}.scramble(:);
                    case 'noise'
                        val = ~D{p,t}.noise(:);
                    case 'rxn'
%                         if isfield(D{p,t},'task')
%                             if strcmp(D{p,t}.task(1:6),'common')
%                                 val = (D{p,t}.articulation(:) - D{p,t}.pulse_on(:)) > 2000;
%                             elseif strcmp(D{p,t}.task(1:8),'auditory') || ...
%                                    strcmp(D{p,t}.task(1:11),'audscramble')
%                                 val = (D{p,t}.articulation(:) - D{p,t}.pulse_off(:)) > 2000;
%                             else
%                                 warning('kjf: get_trials does not recognize task, using rxn_times')
%                                 val = D{p,t}.rxn_time(:) > 2000;
%                             end
%                         else
%                             warning('kjf: get_trials does not recognize task, using rxn_times')
%                             val = D{p,t}.rxn_time(:) > 2000;
%                         end
                        val = D{p,t}.rxn_time(:) > 2000;
                    case 'tech'
                        val = ~D{p,t}.tech(:);
                    case 'stim'
                        val = D{p,t}.behav.stimulus(:);
                    otherwise
                        error('Unrecognized trial characteristic')
                end

                if strcmp('+',varargin{n}(1))
                    add = add | val;
                elseif strcmp('-',varargin{n}(1))
                    rem = rem | val;
                else
                    error('Unrecognized trial add/remove symbol: use + or -')
                end
            end
            trs{p,t} = add & ~rem;
        end
    end
end

% handle output
if ~any(size(trs) ~= 1)
    trs = cell2mat(trs);
end