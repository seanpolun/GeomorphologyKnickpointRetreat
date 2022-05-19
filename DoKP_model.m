function DrainStruct = DoKP_model(DrainStruct)

K = 1.383e-6; % 1/yr
K_unc = 8.622e-7; 
%K = 1.500e-6;
%K_unc = 7.056e-7; % This is 
% K_unc = 5.552e-7;
%K = 6.07e-6; 
%K_unc = 7.41e-7; 

%U = 0.653e-3 % [m/yr]
%U_unc = 0.056e-3 % [m/yr]s?
%DrainStruct = thetaFit(DrainStruct); 


DrainStruct = FindCelerity(DrainStruct,K,K_unc); 

%mChi = DrainStruct.mChi2; 

if DrainStruct.m ~= DrainStruct.Theta
    DrainStruct.m = DrainStruct.Theta;
end
m = DrainStruct.m; 
n = DrainStruct.n;  
A= DrainStruct.A_ds; 
S = DrainStruct.S; 
%unc_perc = K_unc / K; 
%K = K_unc.*randn(size(A))+K; 
%UpliftRate = (mChi^n)*(K*(A_0^m)); 
UpliftRates = K.*(A.^m).*(S'.^n); 
%UpliftUnc = UpliftRate*unc_perc; 
lowerWind = UpliftRates((length(UpliftRates) - 4):length(UpliftRates)); 

UpliftRate = min(UpliftRates);
UpliftUnc = UpliftRate .* (K_unc./K); 
if UpliftRate < (nanmean(lowerWind) - nanstd(lowerWind)) 
    UpliftRate = nanmean(lowerWind); 
    UpliftUnc = sqrt(((nanstd(lowerWind)./UpliftRate).^2) + ((K_unc./K).^2)) .* UpliftRate; 
end

%UpliftUnc = std(UpliftRates); 
DrainStruct.UpliftRate = UpliftRate*1000; 
DrainStruct.UpliftUnc = UpliftUnc*1000; 