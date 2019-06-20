% calculating principal components for groups of electrodes

% get channels
switch lower(subject)
    case 'ta510'
%         chs = [86 103 19 4 5 62];
%         chs = [85 103 19 4 5 62]; % used for SfN
        chs = {
            85:88
            [90:93 103]
            [9:11 18:19 59]
            [4 12 20]% 60]
            [5:6 13:14]
            [56 62]%[56 61:62]
            };
    case 'ts085'
        chs = {
            [64:65 233:234]
            [56:58 219:222 225:226]%[56:58 219:221]%
            [72:73 85:86 99:100 209]
            [75 88 102:103 208:210]
            [77:78 90:92 104:106 213:214]
            [28 34 193:194 201:204]
            };
    otherwise
        error('kjf: unrecognized patient')
end

% constants
fs = D.sampHz;
Ntr = sum(trs);
Nch = length(chs);
Nt = diff(time)*D.sampHz/1000+1;

% eeg features
Xo = zeros(Nch,Ntr,Nt);
for n = 1:Nch
    Xraw = zeros(length(chs{n}),Ntr,Nt);
    for ch = 1:length(chs{n})
        [traces,t,base,bt] = get_eeg(D,'common',chs{n}(ch),trs,500,'pulse_on',time,'pulse_on',baseline);
        [ftraces,fbase] = filt_sigmFFT(freq,false,D.sampHz,500,traces,base);  
        Ptraces = abs(ftraces).^2;
        Pbtraces = abs(fbase).^2;
        bPm = mean(Pbtraces(:));
        
        Xraw(ch,:,:) = bsxfun(@rdivide,bsxfun(@minus,Ptraces,bPm),bPm);
    end
    
    for m = 1:Ntr
        [coeff,score,latent,tsq,expp] = pca(squeeze(Xraw(:,m,:))','algorithm','svd');
        Xo(n,m,:) = score(:,1);
    end
end
yt=Xo;