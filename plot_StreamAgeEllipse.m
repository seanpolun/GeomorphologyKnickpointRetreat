function plot_StreamAgeEllipse(DrainStruct1, DrainStruct2)


K = 1.383e-6; % 1/yr
K_unc = 8.622e-7; 


K = K * 1000; 
K_unc = K_unc * 1000; 
Cefun = @(K,m,A,S,n) K.*(A.^m).*(S.^(n-1)); 
nsim = 4000; 
K_rand = randn(nsim, 1).* K_unc + K; 
%% DrainStruct1 
theta_1 = DrainStruct1.Theta; 
n_1 = 1; 
m_1 = theta_1 ./ n_1; 
DrainStruct1.m = m_1; 
m_1_unc = DrainStruct1.sigmaTheta ./ n_1; 
DrainStruct1.n = n_1; 
x_1 = DrainStruct1.MDist; 
Z_1 = DrainStruct1.Elev; 
xKP_1 = DrainStruct1.xKP; 
Area_1 = DrainStruct1.FlowArea; 

% Initialize default (uncorrelated) uncertainties
def_m_unc = 0.09;
def_n_unc = 0.005; 
def_S_unc = 0.05; 
def_A_unc = 2700;

[~,finder_1] = min(abs(x_1 - xKP_1));  
xds_1 = x_1(finder_1:length(x_1));
zds_1 = Z_1(finder_1:length(Z_1)); 
Ads_1 = Area_1(finder_1:length(Area_1)); 
cf_1 = fit(xds_1',zds_1','power1'); 
S_1 = abs(differentiate(cf_1,xds_1')); 
%p1 = coeffvalues(cf); 
%err1  = confint(cf,0.68); 
%S = abs(p1(1)); 
%S_unc = (err1(2,1) - err1(1,1))/2; 

Terr_out1 = zeros(nsim,1);

for i = 1:nsim
    % K_1_1 = normrnd(K,K_unc);
    K_1_1 = K_rand(i); 
    % K1 = K; 
    % K2 = K; 
    m_1_1 = normrnd(m_1, def_m_unc); 
    % m2 = m; 
    % m1 = m; 
    n_1_1 = normrnd(n_1,def_n_unc);
    % n1 = n; 
    S_1_1 = abs(S_1 + randn(length(S_1),1).*std(S_1).*def_S_unc)';
    A_err_1 = (randn(length(Ads_1),1) .* def_A_unc)'; 
    A_ds_1_1 = abs(Ads_1 + A_err_1);
    % S1 = S';
    %S1 = normrnd(S,S_unc); 

%     for j = 1:leng11
%         %nn = finder + j;
%         Ce_err_out(i,j) = Cefun(K1,m1,Ads,j,S1,n1);
%         
%     end
    Ce_err_vec_1 = Cefun(K_1_1,m_1_1,A_ds_1_1,S_1_1,n_1_1);
    Terr_out1(i) = trapz(xds_1,1./Ce_err_vec_1); 

end
%% DrainStruct2 

theta_2 = DrainStruct2.Theta; 
n_2 = 1; 
m_2 = theta_2 ./ n_2; 
DrainStruct2.m = m_2; 
m_2_unc = DrainStruct2.sigmaTheta ./ n_2; 
DrainStruct2.n = n_2; 
x_2 = DrainStruct2.MDist; 
Z_2 = DrainStruct2.Elev; 
xKP_2 = DrainStruct2.xKP; 
Area_2 = DrainStruct2.FlowArea; 

[~,finder_2] = min(abs(x_2 - xKP_2));  
xds_2 = x_2(finder_2:length(x_2));
zds_2 = Z_2(finder_2:length(Z_2)); 
Ads_2 = Area_2(finder_2:length(Area_2)); 
cf_2 = fit(xds_2',zds_2','power1'); 
S_2 = abs(differentiate(cf_2,xds_2')); 

Terr_out2 = zeros(nsim,1);

for i = 1:nsim
    % K_2_1 = normrnd(K,K_unc);
    K_2_1 = K_rand(i); 
    % K1 = K; 
    % K2 = K; 
    m_2_1 = normrnd(m_2, def_m_unc); 
    % m2 = m; 
    % m1 = m; 
    n_2_1 = normrnd(n_2,def_n_unc);
    % n1 = n; 
    S_2_1 = abs(S_2 + randn(length(S_2),1).*std(S_2).*def_S_unc)';
    A_err_2 = (randn(length(Ads_2),1) .* def_A_unc)'; 
    A_ds_2_1 = abs(Ads_2 + A_err_2);
    % S1 = S';
    %S1 = normrnd(S,S_unc); 

%     for j = 1:leng11
%         %nn = finder + j;
%         Ce_err_out(i,j) = Cefun(K1,m1,Ads,j,S1,n1);
%         
%     end
    Ce_err_vec_2 = Cefun(K_2_1,m_2_1,A_ds_2_1,S_2_1,n_2_1);
    Terr_out2(i) = trapz(xds_2,1./Ce_err_vec_2); 

end
%% Plot
figure; 
plot(Terr_out1, Terr_out2, '.k')
