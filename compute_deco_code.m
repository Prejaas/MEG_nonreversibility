function [FowRev, Hierarchy, xcorr_f, xcorr_r, FowRev_reg] = compute_deco_code(Y,no_lag)
% Tewarie 2022 v2, compute non-reversibility for MEG data (Non-reversibility outperforms functional connectivity in characterisation 
% of brain states in MEG data, 2023 Neuroimage)
%
% code partly based on code from Gustavo Deco, 
% which can be found at: https://github.com/decolab/insideout (paper:
% The INSIDEOUT framework provides precise signatures of the balance of intrinsic 
% and extrinsic dynamics in brain states)

% input: 
%           Y           amplitude envelope data/raw signals
%           no_lag      number of lags (in samples)
%
% output    FowRev      whole brain non-reversibility  
%           FowRev_reg  non-reversibility for every region
%           Hierarchy   hierarchy measure (see https://github.com/decolab/insideout)
%           xcorr_f     cross correlations (forward)
%           xcorr_r     cross correlations (reversed)

tic
Y = transpose(Y);                       % input can be raw signals, in our case amplitude envelopes
Tm = size(Y,1);                         % number of samples
p1 = size(Y,2);                         % number of signals 
xcorr_f    = zeros(p1,p1,no_lag);
xcorr_r    = zeros(p1,p1,no_lag);
FowRev    = zeros(no_lag,1);
% AsymFow   = zeros(no_lag,1);
% AsymRev   = zeros(no_lag,1);
Hierarchy = zeros(no_lag,1);
% Itaufout  = zeros(no_lag,p1);
% Itaufin   = zeros(no_lag,p1);
FowRev_reg = zeros(no_lag,p1);
 
% compute lagged correlations
tic
for Tau = 1: no_lag

    % cor forward
    Y1              = Y(1:Tm-Tau,:);
    Y2              = Y(Tau+1:Tm,:);
    Y1              = bsxfun(@minus,Y1,sum(Y1,1)/(Tm-Tau));  % Remove mean
    Y2              = bsxfun(@minus,Y2,sum(Y2,1)/(Tm-Tau));  % Remove mean
    coef            = Y1'*Y2/(Tm-Tau-1);
    d1              = std(Y1,[],1);
    d2              = std(Y2,[],1);
    xcorr_f(:,:,Tau) = bsxfun(@rdivide,coef,(d1'*d2)).*~eye(p1);
	FCtf=squeeze(xcorr_f(:,:,Tau));
    Itauf=-0.5*log(1-FCtf.*FCtf);

    % cor reverse
    Y1              = Y(Tm:-1:Tau+1,:);
    Y2              = Y(Tm-Tau:-1:1,:);
    Y1              = bsxfun(@minus,Y1,sum(Y1,1)/(Tm-Tau));  % Remove mean
    Y2              = bsxfun(@minus,Y2,sum(Y2,1)/(Tm-Tau));  % Remove mean
    coef            = Y1'*Y2/(Tm-Tau-1);
    d1              = std(Y1,[],1);
    d2              = std(Y2,[],1);
    xcorr_r(:,:,Tau) = bsxfun(@rdivide,coef,(d1'*d2)).*~eye(p1);
    FCtr=squeeze(xcorr_r(:,:,Tau));
    Itaur=-0.5*log(1-FCtr.*FCtr);
    
    Reference=((Itauf(:)-Itaur(:)).^2)';
    index=find(Reference>quantile(Reference,0.95));
    FowRev(Tau)= mean(Reference(index),'omitnan');
%     AsymFow(Tau)=mean(mean(abs(Itauf-Itauf')));
%     AsymRev(Tau)=mean(mean(abs(Itaur-Itaur')));
    Hierarchy(Tau) = std(Reference(index),'omitnan'); 
%     Itaufout(Tau,:) = sum(Itauf,1);
%     Itaufin(Tau,:)  = sum(Itauf,2)';

    for reg = 1: size(Itauf,1)
        Reference_reg=((Itauf(reg,:)-Itaur(reg,:)).^2)';
        index_r=find(Reference_reg>quantile(Reference_reg,0.95));
        FowRev_reg(Tau,reg)= mean(Reference_reg(index_r),'omitnan');
    end


end
time = toc;
fprintf('computed insideout metrics in %f seconds \n',time)

end
