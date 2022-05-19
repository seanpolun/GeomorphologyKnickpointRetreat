function SampleStruct = Calc2StageUplift(SampleStruct)
% SampleStruct = Calc2StageUplift(SampleStruct)
% 
% Performs iterative computation of uplift rate between two end-member 36Cl
% production models. 
% Stage 1: Steady state uplift within the rock column. See Tectonic
% Geomorpholgy, Burbank and Anderson for examples. 
% U = Z*(P0 - N0*lambda) / N0
% Where Z* is Attenuation length / density
% P0 is the production rate
% N0 is the measured concentration of 36Cl
% lambda is the decay constant for 36Cl
% 
% 
% Stage 2: Surface exposure generating 36Cl. This requires that the dated
% surface is elevated above the modern channel. 
% U = h / t, 
% t = -(1/lambda)*ln(1 - (lamda*(N0 - Ne))/P0)
% Where Ne is the produced 36Cl concentration from Stage 1. 
% 
% Adapted from David Horrell's MS thesis. Implemented by Sean Polun, 2019. 

%% Initialize variables
P0 = SampleStruct.P0; 
P0_unc = SampleStruct.P0_unc; 
N_meas = SampleStruct.N_meas; 
N_unc = SampleStruct.N_unc; 
height = SampleStruct.height * 100; % convert m to cm 
height_unc = SampleStruct.height_unc * 100; % convert m to cm
att_len = SampleStruct.att_len; 
att_len_unc = 0.1 * att_len; % 10% unc 

lambda = 2.3028e-6;
density = SampleStruct.density; % g/cc
Z_star = (att_len / density); 
Z_star_unc = sqrt(((SampleStruct.dens_unc ./ density).^2) + ((att_len_unc ./ att_len).^2)).* Z_star; 


%% Calculate Endmember states
St1_Ne = @(P0,lambda,U,Z_star) P0 ./(lambda + (U ./ Z_star)); 
St1_U = @(Ne,P0,lambda,Z_star) (Z_star .* (P0 - (Ne.*lambda)))./Ne;

St2_U = @(Ne,N0,P0,lambda,height) height./((-1./lambda).*log(1 - ((lambda.*(N0 - Ne))./P0))); 

SampleStruct.St1_U_em = St1_U(N_meas,P0,lambda,Z_star); 
SampleStruct.St2_U_em = St2_U(0,N_meas,P0,lambda,height); 
% 
% %% Solve for U via iteration
% tol = 0.00000001; 
% disagree = abs(St1_U_em - St2_U_em); 
% U = St2_U_em; 
% n = 1; 
% while disagree > tol
%     Ne = St1_Ne(P0,lambda,U,Z_star); 
%     U1 = St2_U(Ne,N_meas,P0,lambda,height); 
%     disagree = abs(U1 - U); 
%     U = U1; 
%     n = n + 1; 
% end
% %disp(n)
% %disp('Iterations')

%% Propagate uncertainty via monte carlo simulation
nsim = 4000; 
U_vec = zeros(nsim,1); 
Ne_vec = zeros(nsim,1);
for i = 1:nsim
    P0_1 = normrnd(P0,P0_unc); 
    N_meas_1 = normrnd(N_meas,N_unc); 
    height_1 = normrnd(height,height_unc); 
    Z_star1 = normrnd(Z_star,Z_star_unc); 
    % Ne1 = St1_Ne(P0_1,lambda,U,Z_star1); 
    % U_vec(i,1) = St2_U(Ne1,N_meas_1,P0_1,lambda,height_1); 
    [U_vec(i, 1),Ne_vec(i, 1)] = solveForU(P0_1, N_meas_1, lambda, Z_star1, height_1);
end

U_unc = std(U_vec); 
U_mean = mean(U_vec); 
U = U_mean;
Ne_mean = floor(mean(Ne_vec));
Ne_std = floor(std(Ne_vec));

% U_unc = sqrt((U_unc.^2) + ((U - U_mean).^2)); 

U = U * 10; % Convert cm/yr to m/ka
U_unc = U_unc * 10; % Convert cm/yr to m/ka

SampleStruct.U = U; 
SampleStruct.U_unc = U_unc; 
SampleStruct.Ne_mean = Ne_mean;
SampleStruct.Ne_std = Ne_std;
end

function [U, Ne] = solveForU(P0, N_meas, lambda, Z_star, height)
%% Calculate Endmember states
St1_Ne = @(P0,lambda,U,Z_star) floor(P0 ./(lambda + (U ./ Z_star))); 
St1_U = @(Ne,P0,lambda,Z_star) (Z_star .* (P0 - (Ne.*lambda)))./Ne;

St2_U = @(Ne,N0,P0,lambda,height) height./((-1./lambda).*log(1 - ((lambda.*(N0 - Ne))./P0))); 

U_range = [0:0.00001:0.5];
Ne_range = St1_Ne(P0, lambda, U_range, Z_star);
too_many_atoms = Ne_range < N_meas;
Ne_range = Ne_range(too_many_atoms);
Ne_range = sort(Ne_range);
U1 = St1_U(Ne_range, P0, lambda, Z_star);
U2 = St2_U(Ne_range, N_meas, P0, lambda, height); 
diff = abs(U1 - U2);
[~, I] = min(diff);
Ne = Ne_range(I);
U = St2_U(Ne, N_meas, P0, lambda, height);



% %% Solve for U via iteration
% tol = 0.00000001; 
% n_max = 10000;
% disagree = abs(St1_U_em - St2_U_em); 
% % abs_disagree = disagree; 
% U = St2_U_em; 
% n = 1; 
% while disagree > tol
%     Ne = St1_Ne(P0,lambda,U,Z_star); 
%     U1 = St1_U(Ne,P0,lambda,Z_star); 
%     U2 = St2_U(Ne,N_meas,P0,lambda,height); 
%     disagree = abs(U1 - U2); 
%     U = U2; 
%     n = n + 1;
%     if n > n_max
%         break
%     end
% end
% disp(n)
% U_out = U;
end
