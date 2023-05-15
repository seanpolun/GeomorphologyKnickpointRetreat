function DrainStruct = DoKP_model(DrainStruct)

% Refined K using quasi-bayesian analysis...
K = 1.4335E-6; 
K_unc = 2.1587E-7; 


% Using only values from Hann... 
% 2/10/23
% K = 1.9170E-6; 
% K_unc = 5.2249E-7; 

% % As of 2/9/23
% K = 1.66e-6; % 1/yr
% K_unc = 6.42E-7; % 1/yr

% K = 1.67e-6; % 1/yr
% K_unc = 7.18e-7; 
%K = 1.500e-6;
%K_unc = 7.056e-7; % This is 
% K_unc = 5.552e-7;
%K = 6.07e-6; 
%K_unc = 7.41e-7; 

%U = 0.653e-3 % [m/yr]
%U_unc = 0.056e-3 % [m/yr]s?
%DrainStruct = thetaFit(DrainStruct); 


DrainStruct = FindCelerity(DrainStruct,K,K_unc); 
DrainStruct.K_used = K; 
DrainStruct.K_unc_used = K_unc; 
%mChi = DrainStruct.mChi2; 

if DrainStruct.m ~= DrainStruct.Theta
    DrainStruct.m = DrainStruct.Theta;
end
x = DrainStruct.MDist; 
m = DrainStruct.m; 
n = DrainStruct.n;  
A= DrainStruct.A_ds; 
S = DrainStruct.S; 
A_pix_size = (x(2) - x(1))^2; 

% Initialize default (uncorrelated) uncertainties
def_m_unc = 0.01;  % 10% 
def_n_unc = 0.005; % 0.5% 
def_S_unc = 0.05; % 5%
def_A_unc = 3 * A_pix_size; 
%unc_perc = K_unc / K; 
%K = K_unc.*randn(size(A))+K; 
%UpliftRate = (mChi^n)*(K*(A_0^m)); 
n_sim = 4000; 
UpliftOut = zeros(n_sim, length(A)); 
for i=1:n_sim
    K1 = normrnd(K, K_unc); 
    m1 = normrnd(m, def_m_unc); 
    n1 = normrnd(n,def_n_unc);
    S1 = abs(S + randn(length(S),1).*std(S).*def_S_unc)';
    A_err = (randn(length(A),1) .* def_A_unc)'; 
    A_ds1 = abs(A + A_err);
    UpliftRates = K1.*(A_ds1.^m1).*(S1.^n1); 
    UpliftOut(i, :) = UpliftRates; 
end
UVec = nanmean(UpliftOut); 
U_unc_vec = nanstd(UpliftOut); 
%UpliftUnc = UpliftRate*unc_perc; 
lowerWind = UVec((length(UVec) - 4):length(UVec)); 
lowerWind_unc = U_unc_vec((length(U_unc_vec) - 4):length(U_unc_vec)); 

[UpliftRate, I] = min(UVec);
% UpliftUnc = UpliftRate .* (K_unc./K); 
UpliftUnc = U_unc_vec(I); 
if UpliftRate < (nanmean(lowerWind) - nanstd(lowerWind)) 
    UpliftRate = nanmean(lowerWind); 
    % UpliftUnc = sqrt(((nanstd(lowerWind)./UpliftRate).^2) + ((K_unc./K).^2)) .* UpliftRate; 
    UpliftUnc = nanmean(lowerWind_unc); 
end

%UpliftUnc = std(UpliftRates); 
DrainStruct.UpliftRate = UpliftRate*1000; 
DrainStruct.UpliftUnc = UpliftUnc*1000; 