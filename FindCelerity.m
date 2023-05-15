function DrainStruct = FindCelerity(DrainStruct,K,K_unc)

% DrainStruct = FindCelerity(DrainStruct,K,K_unc)
% Determines n and celerity rate from stream parameters
% Inputs: 
% DrainStruct: A MATLAB structure containing xKP (downstream location of
% knickpoint (m)), Theta (concavity index), mChi (slope in Chi space), h.
% These are determined using ChiFit.m and StreamPowerLaws.m. 
% m: assumed value of the m exponent for the stream power law. This is
% assumed to be constant across drainage networks for a region. 
% K: Regional erodability constant (1/yr). 
% Outputs: 
% Ce (mm/yr): Rate at which a knickpoint retreats 

% if DrainStruct.Theta > 0.7 || DrainStruct.Theta < 0.3
%     theta = 0.42; 
% else 
%     theta = DrainStruct.Theta; 
% end

% if DrainStruct.Theta > 0.6 || DrainStruct.Theta < 0.4
%     theta = 0.5; 
% else 
%     theta = DrainStruct.Theta; 
% end

%n = m / DrainStruct.Theta; 
theta = DrainStruct.Theta; 
n = 1; 
m = theta ./ n; 
DrainStruct.m = m; 
m_unc = DrainStruct.sigmaTheta ./ n; 
DrainStruct.n = n; 
x = DrainStruct.MDist; 
Z = DrainStruct.Elev; 
xKP = DrainStruct.xKP; 
Area = DrainStruct.FlowArea; 
A_pix_size = (x(2) - x(1))^2; 

% Initialize default (uncorrelated) uncertainties
def_m_unc = 0.01;  % 10% 
def_n_unc = 0.005; % 0.5% 
def_S_unc = 0.05; % 5%
def_A_unc = 3 * A_pix_size; 


[~,finder] = min(abs(x - xKP));  
xds = x(finder:length(x));
zds = Z(finder:length(Z)); 
Ads = Area(finder:length(Area)); 
cf = fit(xds',zds','power1'); 
S = abs(differentiate(cf,xds')); 
%p1 = coeffvalues(cf); 
%err1  = confint(cf,0.68); 
%S = abs(p1(1)); 
%S_unc = (err1(2,1) - err1(1,1))/2; 

K = K * 1000; 
K_unc = K_unc * 1000; 
K_unc_p = K_unc / K; 
Cefun = @(K,m,A,S,n) K.*(A.^m).*(S.^(n-1)); 

nsim = 4000; 
[~, maxXI] = max(x); 
leng11 = maxXI - finder; 
num_nodes = leng11 + 1;
Ce_err_out = zeros(nsim,num_nodes); 

%% Calculate Uncorrellated error 
Terr_out = zeros(nsim,1); 
for i = 1:nsim
    % K1 = normrnd(K,K_unc);
    K1 = 1; % changing things so K is out of the integral
    m1 = normrnd(m, def_m_unc); 
    % m1 = m; 
    n1 = normrnd(n,def_n_unc);
    % n1 = n; 
    S1 = abs(S + randn(length(S),1).*std(S).*def_S_unc)';
    A_err = (randn(length(Ads),1) .* def_A_unc)'; 
    A_ds1 = abs(Ads + A_err);
    % S1 = S';
    %S1 = normrnd(S,S_unc); 

%     for j = 1:leng11
%         %nn = finder + j;
%         Ce_err_out(i,j) = Cefun(K1,m1,Ads,j,S1,n1);
%         
%     end
    Ce_err_vec = Cefun(K1,m1,A_ds1,S1,n1);
    Ce_err_out(i,:) = Ce_err_vec * normrnd(K, K_unc);
    Terr_out(i) = trapz(xds,1./Ce_err_vec)*(1/K); 
end
% Ce_mean = nanmean(Ce_err_out);
% Ce_std = nanstd(Ce_err_out);
Ce_err_out = Ce_err_out(:,1); 
if ~isreal(Ce_err_out)
    Ce_err_out(imag(Ce_err_out) ~= 0) = NaN; 
end
if ~isreal(Terr_out)
    Terr_out(imag(Terr_out) ~= 0) = NaN; 
end

%% Generate complete t_err for error correllation
Terr_out2 = zeros(nsim,1); 
for i = 1:nsim
    % K1 = normrnd(K,K_unc);
    K1 = normrnd(1/K, (1/K)*K_unc_p); 
    % K1 = K; 
    m1 = normrnd(m, def_m_unc); 
    % m1 = m; 
    n1 = normrnd(n,def_n_unc);
    % n1 = n; 
    S1 = abs(S + randn(length(S),1).*std(S).*def_S_unc)';
    A_err = (randn(length(Ads),1) .* def_A_unc)'; 
    A_ds1 = abs(Ads + A_err);
    % S1 = S';
    % A_ds1 = Ads; 
    %S1 = normrnd(S,S_unc); 

%     for j = 1:leng11
%         %nn = finder + j;
%         Ce_err_out(i,j) = Cefun(K1,m1,Ads,j,S1,n1);
%         
%     end
    Ce_err_vec = Cefun(1,m1,A_ds1,S1,n1);
    % Ce_err_out(i,:) = Ce_err_vec;
    Terr_out2(i) = trapz(xds,1./Ce_err_vec)*(K1); 
    % Terr_out2(i) = trapz(xds, 1./Ce_err_vec); 
end
if ~isreal(Terr_out2)
    Terr_out2(imag(Terr_out2) ~= 0) = NaN; 
end
Terr_out2(Terr_out2 < 0) = NaN; 
Terr_out2(Terr_out2 > 5000) = NaN; 
% for i = 1:nsim
%     Ce1 = zeros(1,num_nodes);
%     for j = 1:num_nodes
%         Ce1(j) =  normrnd(Ce_mean(j), Ce_std(j));
%     end
%    Terr_out(i) = trapz(xds, 1./Ce1); 
% end
% if ~isreal(Terr_out)
%     Terr_out(imag(Terr_out) ~= 0) = NaN;
% end


%% Output to drain struct
Ce = nanmean(Ce_err_out(:)); 
Ce_unc = nanstd(Ce_err_out(:)); 
tKP = nanmean(Terr_out2); 
tKP_unc = nanstd(Terr_out); 
DrainStruct.DS_nodes = length(S); 
DrainStruct.Ce = Ce; 
DrainStruct.Ce_unc = Ce_unc; 
DrainStruct.A_ds = Ads; 
DrainStruct.tKP = tKP; 
DrainStruct.tKP_uncorr_unc = tKP_unc; 
% DrainStruct.tKP_corr_unc = nanstd(Terr_out2); 
DrainStruct.S = S; 
DrainStruct.n = n; 
% DrainStruct.tKP_tot_unc = DrainStruct.tKP_uncorr_unc + DrainStruct.tKP_corr_unc; 
DrainStruct.tKP_tot_unc = nanstd(Terr_out2); 
% DrainStruct.tKP_corr_unc = tKP * (K_unc / K); 
% DrainStruct.tKP_tot_unc = sqrt((DrainStruct.tKP_corr_unc .^2) + (DrainStruct.tKP_uncorr_unc .^2));
DrainStruct.tKP_corr_unc = sqrt(DrainStruct.tKP_tot_unc^2 - DrainStruct.tKP_uncorr_unc^2); 
DrainStruct.T_all_out = sort(Terr_out2); 