list = who; 
outName = 'new_2_25_21.csv'; 
fid = fopen(outName,'w');
header1 = 'Name,Theta,Downstream Nodes,m,Ce,Ce Uncert,Age Kp,Uncert,Uplift Rate,Uplift Uncert'; 
header2 = ' , , , ,[m/ka],[m/ka],[ka],[ka],[m/ka],[m/ka]'; 
fprintf(fid, '%s\n',header1); 
fprintf(fid, '%s\n',header2); 
for i = 1:length(list)
    itemName = list(i);
    itemName = itemName{1};
    
    string = strcat('WorkStruct = ' , itemName , '; ');
    eval(string);
    
    Theta = WorkStruct.Theta;
%     Theta_unc = WorkStruct.Theta_unc;
%     
%     ka = WorkStruct.ka;
%     ka_err = WorkStruct.ka_err;
%     
%     h = WorkStruct.h;
%     h_err = WorkStruct.h_err;
    
    %S = WorkStruct.S;
    DS_Nodes = WorkStruct.DS_nodes; 
    m = WorkStruct.m;
    
    Ce = WorkStruct.Ce;
    Ce_unc = WorkStruct.Ce_unc;
    
    tKP = WorkStruct.tKP;
    tKP_unc = WorkStruct.tKP_unc;
    
    UpliftRate = WorkStruct.UpliftRate;
    UpliftUnc = WorkStruct.UpliftUnc;
    fprintf(fid, '%s,%.2f,%1.2f,%1.2f,%.3f,%.2f,%3.1f,%3.1f,%.3f,%.3f\n', itemName,Theta,DS_Nodes,m,Ce,Ce_unc,tKP,tKP_unc,UpliftRate,UpliftUnc);
end
fclose(fid);
clearvars i string itemName fid outName list Theta DS_Nodes m Ce Ce_unc tKP tKP_unc UpliftRate UpliftUnc header1 header2 WorkStruct