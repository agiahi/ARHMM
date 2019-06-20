%
% Probabilistic Inference for EEG/ECoG data
% created by Aram Giahi
%
% to run principal component analysis, set dopca=true

dw=10; % frames per quasistationary time window for AMVAR
cluster=1; % k-means clustering of AMVAR estimates
nq=6; % # EM recursions (set to 1 for crossvalidation) 
tau=10; % dyn. time scale (ms)
subject = 'sim'; % subject ('sim' for simulated data)
dopca=false; % performing principal component analysis
shownet=false; % show network of interactions
scrambled=0; % condition: scrambled / unscrambled pictures
z2=0; % trial normalization
u=.5; % state residency prior ('stickiness') with time constant (u-1)*tau
crossval=0; alpha=.8; % crossvalidation (alpha: crossval ratio)
adaptcut=0; % cut-off relative to articulation onset

clear seq;
if ~strcmp(subject,'sim')
    ecogdata=1;
    if ~exist('psub','var')
        paths = get_paths;
        fprintf('%s %s %s\n','loading', subject, '...');
        D = get_nkdata(paths,subject,'common');
    else
        if ~strcmp(psub,subject)
            fprintf('%s %s %s\n','loading', subject, '...');
            D = get_nkdata(paths,subject,'common');
        else
            fprintf('subject: %s\n',subject)
        end
    end
    psub=subject;
    
    if 0
        trs = get_trials(D,'+correct','+unscrambled','+scrambled','-noise','-tech');
    else
        if scrambled
            trs = get_trials(D,'+correct','-unscrambled','-noise','-tech');
        else
            trs = get_trials(D,'+correct','-scrambled','-noise','-tech');
        end
    end
    trialind=find(trs>0); % indices of selected trials
    rxntimes=D.articulation(trialind)-D.pulse_on(trialind);
    time = [-200 3500]; offset=500+time(1); cutoff=3500-time(2);
    baseline = [-750 -250];
    freq = [61 119];
    
    % electrode selection for each subject
    switch subject
        case 'ta421' % right sided
            ns=4; % #states
            p=2; % model order
            % SO6/7/8 PST2/3 LF13 LF21 LF23 LTO2
            channels={'*PST4-AV', '*LF13-AV', '*LF15-AV', '*LF18-AV', '*LTO1-AV'};
            clear rxntimesoff;
        case 'ta436'
            u=4;
            ns=4; % #states
            p=3; % model order
            channels={'*SubO3-AV', '*PST4-AV', '*LF10-AV', '*LF13-AV', '*LF5-AV', '*LT30-AV'};
            clear rxntimesoff;
            %channels={'*LF13-AV', '*LF7-AV','*LT30-AV'};
        case 'ta356'
            %EVC: LSO5 LSO6
            %mFusi: LPST5
            %pTri: LF9 LF10
            %pOp: LF12 LF4
            %SCG: LF7 LF6 LT21
            %STG: LT17 LT16 LT15
            ns=4; % #states
            p=3; % model order
            channels={'*LSO5-AV', '*LPST5-AV', '*LF9-AV', '*LF12-AV', '*LF6-AV', '*LT17-AV'};
            %channels={'*LF12-AV', '*LF6-AV'};
            clear rxntimesoff;
        case 'ta510'
            u=.5;
            ns=4; % #states
            p=3; % model order
            %pOrb: LT25/26 pTri: LF9/10/11 pOp: LF4/12/28 M1: LF5/13/29
            channels={'OP6','PST7','LF19','LF4','LF5','LT30'};
            %channels={'PST7','LF19','LF4','LF5','LT30'};
            %channels={'LF4','LF5','LT30'};
            if scrambled||0
                h=load(sprintf('%s/Volumes/Sector1/ta510_lfp_files/ta510_common_supplement.mat',paths.root));
                rxntimesoff=h.aoff(trialind); % pulse_on already substracted!
            else
                load rxntimes510;
            end
            if crossval
                trials_train=1:floor(alpha*size(rxntimesoff,1));
                trials_test=ceil(alpha*size(rxntimesoff,1)):size(rxntimesoff,1);
            end
        case 'ta505'
            ns=4; % #states
            p=3; % model order
            channels={'LO8','PST8','LF14','LF5','LF6','LT34'};
            %channels={'OP2','PST8','LF12','LF4','LF6','LT34'};
            %channels={'LO8','PST8','LF12','LF4','LF6','LT34'};
            %channels={'LF4','LF5','LF6'};
        case 'ts016'
            ns=4; % #states
            p=2; % model order
            %pOrb: LF2/3
            %pTri: LF4/5
            %pOp: LF6
            %M1: LF7/15/28
            %channels={'*PST6-AV','*LF12-AV','*LF14-AV','*LF7-AV','*LT21-AV'};%,'*LO8-AV'}'*LO4-AV''*PST2-AV';
            channels={'*LO4-AV','*MST4-AV','*LF12-AV','*LF14-AV','*LF7-AV','*LT21-AV'};%,'*LO8-AV'}'*LO4-AV''*PST2-AV';
            %channels={'*LF2-AV','*LF4-AV','*LF6-AV','*LF7-AV'};
        case 'ts081' % right sided
            u=.5;
            ns=4; % #states
            p=3; % model order
            %EVC: PSO6 (or: PSO5, LT30)
            channels={'LT30','ASO8','LF21','LF23','LP50','LT5'};%,'*LO8-AV'}'*LO4-AV''*PST2-AV';'LT30','ASO8',
            h=load(sprintf('%s/Volumes/Sector1/ts081_lfp_files/ts081_common_supplement.mat',paths.root));
            rxntimesoff=h.artic_off(trialind)-D.pulse_on(trialind);
        case 'ts085'
            u=.5;
            ns=4; % #states
            p=3; % model order
            channels={'PST5','MST5','LF38','LF25','LF12','LT58'};
            %channels={'PST6','PSTG5','LF38','LF25','LF12','LT58'};%,'*LO8-AV'}'*LO4-AV''*PST2-AV';'PST6','MST4',
            %EVC: PST6 (or: PST5, LO1)
            h=load(sprintf('%s/Volumes/Sector1/ts085_lfp_files/ts085_common_supplement.mat',paths.root));
            rxntimesoff=h.artic_off(trialind)-D.pulse_on(trialind);
    end
    if ~exist('rxntimesoff','var'),rxntimesoff=rxntimes;end
    rxntimesoff=min(rxntimesoff,3500);
    % extract data for selected electrodes
    % ****************************************************************************************************
    elect=sfind(cellstr(D.ch_names),channels); indxs=sum(elect>0);
    if indxs>0
        elect=elect(find(elect>0));
        yt=zeros(indxs,size(trialind,1),diff(time)*D.sampHz/1000+1);
    end
    

    for j=1:size(elect,1)
        if elect(j)>0
            
            [traces,t,base,bt] = get_eeg(D,'common',elect(j),trs,500,'pulse_on',time,'pulse_on',baseline);
            [ftraces, fbase] = filt_sigmFFT(freq,false,D.sampHz,500,traces,base);            
            %[f2, fbase] = filt_sigmFFT([2 15],false,D.sampHz,500,traces,base);

            % get power and normalize to baseline
            if 1
                Ptraces = abs(ftraces).^2;
                Pbtraces = abs(fbase).^2;
                bPm = mean(Pbtraces(:));
                yt(j,:,:) = bsxfun(@rdivide,bsxfun(@minus,Ptraces,bPm),bPm);
            else
                % imaginary part of signal to recover phase dependencies
                % (negative for 180 degree phase coupling)
                %[ftraces, fbase] = filt_sigmFFT([2 15],false,D.sampHz,500,traces,base);
                yt(j,:,:) = imag(ftraces);
                %yt(j+nd,:,:) = real(ftraces);
                if 0
                    [ftraces, fbase] = filt_sigmFFT(freq,false,D.sampHz,500,traces,base);
                    Ptraces = abs(ftraces).^2;
                    Pbtraces = abs(fbase).^2;
                    bPm = mean(Pbtraces(:));
                    yt(j+nd,:,:) = bsxfun(@rdivide,bsxfun(@minus,Ptraces,bPm),bPm);
                end
            end
        end
    end
    % ****************************************************************************************************
else % simulated data
    ecogdata=0;
    ns=3;p=1;
    paths = get_paths;
end
% run ARHMM inference
addpath(paths.code); 
fprintf('aHMM8 (ns=%i, p=%i)',ns,p); 
if dopca, fprintf('\n principal component analysis');Elpca;end
aHMM8v; if shownet, allgr; end;