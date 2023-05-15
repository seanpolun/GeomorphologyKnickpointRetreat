function [K, K_unc, K_out] = SolveK_celerity(Samples,DrainStruct)

numSamps = length(Samples); 

m = DrainStruct.m; 
n = DrainStruct.n; 

x = DrainStruct.MDist; 
Z = DrainStruct.Elev;
xKP = DrainStruct.xKP;
Area = DrainStruct.FlowArea;

sampWind = 1;
[~,finder] = min(abs(x - xKP));

xds = x(finder:length(x));

zds = Z(finder:length(Z));
Ads = Area(finder:length(Area));
cf = fit(xds',zds','power1');
S = abs(differentiate(cf,xds'));
power = (Ads.^m).*(S.^n);
nsim = 4000;

figure; 
yyaxis left
plot(xds,zds,'-k')
yyaxis right
% hold on
plot(xds, Ads,'--b')
% plot(xds, S,'--r')

% plot(xds, power,'-g')
yyaxis left
hold on
K_out = zeros(nsim,numSamps);
%K_unc_out = zeros(nsim,numSamps);
for i = 1:numSamps
    U = Samples(i).U ;
    U = U ./ 1000;
    U_unc = Samples(i).U_unc; 
    U_unc = U_unc ./ 1000;
    Dist = Samples(i).Dist;
    xDist = max(x) - Dist;
    
    
    
    [~, sampleFind] = min(abs(xds - xDist));
    plot(xds(sampleFind),zds(sampleFind),'or')
    text(xds(sampleFind),zds(sampleFind)+20,num2str(i))
    sampL = sampleFind - sampWind;
    sampH = sampleFind + sampWind;
    A_samp = mean(Ads(1,sampL:sampH));
    A_std = std(Ads(1,sampL:sampH));
    S_samp = mean(S(sampL:sampH,1));
    S_std = std(S(sampL:sampH,1));
    %K_vec = zeros(nsim,1);
    
    for j = 1:nsim
        U1 = normrnd(U,U_unc);
        A = normrnd(A_samp,A_std);
        S1 = normrnd(S_samp,S_std);
        P2 = (A.^m).*(S1'.^n);
        %P2 = nanmean(PowVec);
        K_out(j,i) = U1./P2;
    end
    
end
K = nanmean(K_out(:));
K_unc = nanstd(K_out(:)); 
end

