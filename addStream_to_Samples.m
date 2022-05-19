function Samples = addStream_to_Samples(Samples, DrainStruct)
numSamps = length(Samples); 

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

for i = 1:numSamps
    Dist = Samples(i).Dist;
    xDist = max(x) - Dist;
    [~, sampleFind] = min(abs(xds - xDist));
    sampL = sampleFind - sampWind;
    sampH = sampleFind + sampWind;
    A_samp = mean(Ads(1,sampL:sampH));
    A_std = std(Ads(1,sampL:sampH));
    S_samp = mean(S(sampL:sampH,1));
    S_std = std(S(sampL:sampH,1));
    %K_vec = zeros(nsim,1);
    Samples(i).A_mean = A_samp;
    Samples(i).A_std = A_std; 
    Samples(i).S_mean = S_samp;
    Samples(i).S_std = S_std;   
    Samples(i).A_norm = A_samp / max(Ads);
    
end

end
