function DrainStruct = thetaFit(DrainStruct, clipping)
% DrainStruct = thetaFit(DrainStruct)
% Calculates the curvature of a selected drainage (Flint's Law), where 
% S = cA^(m/n). The relationship is fit as a power law using MATLAB's curve
% fitting toolbox. 
% In order to fit this relationship with noisy DEM data, the power law
% scaling relationships for Slope - Length and Area - Length are fit
% initially, and the modeled data (with added random noise based on the
% goodness of fit statistics) then has the Flint's law relationship fit. 
% 
% Profiles which are overly dense are interpolated to a lower density
% 
% If the primary knickpoint has not been previously plotted, the program
% will plot the elevation - length profile of the entire trunk stream, and
% the user will identify the primary knickpoint. 
% 
% 

%% Condition Data
AreaNoise = 0.75 ; % 0.25
SlopeNoise = 0.5 ;   % 0.5

nx = 100; % lower limit of points (will interpolate to this if below
ux = 500; % upper limit of points (will downsample to this if above)

Area = DrainStruct.FlowArea;
x = DrainStruct.MDist;
Z = DrainStruct.Elev; 

ASize = length(Area); 
if Area(ASize) < Area(ASize - 1)
    Area(ASize) = []; 
    x(ASize) = []; 
    Z(ASize) = [];
end

if ~isfield(DrainStruct,'xKP')
    plot(x,Z,'ok')
    [KPloc,~] = ginput(1); 
    DrainStruct.xKP = KPloc; 
else 
    KPloc = DrainStruct.xKP; 
end


[~,finder] = min(abs(x - KPloc)); 
finder = finder + 5; 
seg_end = length(Area);
if exist('clipping', 'var')
    if length(clipping) == 1
        seg_end = clipping + finder;
    elseif length(clipping) == 2
        finder = finder - 5;
        seg_end = clipping(2) + finder;
        finder = clipping(1) + finder;
    end
end

A_ds = Area(finder:seg_end);
disp('Measurements Downstream:')
disp(length(A_ds))
A_ds = movmax(A_ds,3); 
Z_ds = Z(finder:seg_end);
x_ds = x(finder:seg_end); 
maxX = max(x_ds); 
if length(x_ds) < nx 
    minX = min(x_ds); 
    xRange = maxX - minX; 
    dx = xRange / nx; 
    xEval = minX:dx:maxX; 
    Z_ds = interp1(x_ds,Z_ds,xEval,'linear');
    A_ds = interp1(x_ds,A_ds,xEval,'linear');
    x_ds = xEval; 
elseif length(x_ds) > ux 
    minX = min(x_ds); 
    xRange = maxX - minX; 
    dx = xRange / ux; 
    xEval = minX:dx:maxX;
    Z_ds = interp1(x_ds,Z_ds,xEval,'spline');
    A_ds = interp1(x_ds,A_ds,xEval,'spline');
    x_ds = xEval; 
else 
    xEval = x_ds; 
end
% xEval = x_ds; 
maxZ = max(Z_ds); 


%% Fit Downstream Z Profile
Z_dsa = Z_ds ./ maxZ; 
x_dsa = x_ds ./ maxX;
%x_dsa = x_dsa - x_dsa(1) + 0.001; 
xEvala = xEval ./ maxX; 
%xEvala = xEvala - xEvala(1) + 0.001; 
% SfitOpt = fitoptions('power2'); 
SfitOpt = fitoptions('power2'); 
SfitOpt.Robust = 'Bisquare';
SfitOpt.Lower = [-Inf -Inf -Inf]; 
SfitOpt.Upper = [Inf 0 Inf]; 
%SfitOpt.Normalize = 'on';

[cf,gof] = fit(x_dsa',Z_dsa','power2',SfitOpt); 
init_rmse = gof.rmse; 
stop = 0; 
zFun = feval(cf,x_dsa); 
last_rmse = init_rmse; 
count = 0; 
while stop == 0
    weight = zeros(size(x_dsa)); 
    weight(Z_dsa < zFun') = 1; 
    weight(Z_dsa >= zFun') = 0.1; 
    %weight = 1./x_dsa; 
    SfitOpt.Weights = weight; 
    [cf1,gof] = fit(x_dsa',Z_dsa','power2',SfitOpt); 
    current_rmse = gof.rmse; 
    diff = current_rmse - last_rmse; 
    disp(diff)
    
    zFun = feval(cf1,x_dsa);
    if abs(diff) < 0.005
    stop = 1; 
    else 
        last_rmse = current_rmse; 
    end
    count = count + 1; 
end
disp('Iterations Elapsed:')
disp(count)
S = abs(differentiate(cf1,xEvala'));
% S = abs(differentiate(cf,x_dsa'));
S = S .* (maxZ./maxX); 
S = abs(S + randn(length(S),1).*gof.rmse.*SlopeNoise); 
%S = S - max(S); 
%S = S ./ max(S); 
DrainStruct.Slope_gof = gof; 
figure;
plot(x_dsa,Z_dsa,'ok')
hold on
plot(cf)
plot(cf1,'--b')


%% Fit A - length profile
%weights = 1./x_ds; 
%weights = (x_ds / max(x_ds)).^2; 
Aoptions = fitoptions('power1'); 
Aoptions.Robust = 'Bisquare'; 
Aoptions.Lower = [-Inf 1.5 ]; 
Aoptions.Upper = [Inf 5 ];
Aoptions.StartPoint = [ 0 1.67]; 
Aoptions.MaxFunEvals = 2400; 
Aoptions.MaxIter = 1200; 
%Aoptions.Normalize = 'on';
%Aoptions.Weights = weights; 



x_ds1 = x_ds - x_ds(1) + 1;
xEval1 = xEval - xEval(1) + 1; 
A_ds1 = abs(A_ds - A_ds(1) + 30^2); 
Amax = max(A_ds1); 
A_ds1 = A_ds1 ./ Amax;
[cA,gof] = fit(x_ds1',A_ds1','power1',Aoptions); 
Aeval = feval(cA,xEval1); 
Afit_unc = gof.rmse*AreaNoise; 
%h = 1/cA.b; 
% Aeval = feval(cA,x_ds1); 

% Aeval = (Aeval + abs(min(Aeval)));
% scF = max(A_ds1) / max(Aeval); 
% Aeval = Aeval.*scF;  
Aeval = abs(Aeval + randn(length(Aeval),1).*Afit_unc);
Aeval = Aeval ./ Amax; 
%weights = 1./(Aeval).^2;
figure; 
plot(x_ds1,A_ds1,'ok')
hold on
plot(cA)
DrainStruct.Area_gof = gof; 

%% Fit Area - Slope profile

ASoptions = fitoptions('power2'); 
ASoptions.Robust = 'Bisquare'; 
ASoptions.Lower = [-Inf 0 -Inf]; 
ASoptions.Upper = [Inf 1 ( max(S) + 0.05)];
ASoptions.StartPoint = [ 0 0.5 0]; 
%ASoptions.Weights = weights; 


[cAS,gof] = fit(Aeval,S,'power2',ASoptions); 
figure; 
plot(Aeval,S,'ok')
hold on
plot(cAS)
%length(Aeval)
coeffAS = coeffvalues(cAS); 
confAS = confint(cAS); 
confAS = confAS(:,2);

theta = abs(coeffAS(2)); 
sigmaTheta = max(confAS) - theta; 

DrainStruct.Flints_gof = gof; 
DrainStruct.Theta = theta; 
DrainStruct.sigmaTheta = sigmaTheta; 
DrainStruct.cA = cA; 
DrainStruct.cAS = cAS; 