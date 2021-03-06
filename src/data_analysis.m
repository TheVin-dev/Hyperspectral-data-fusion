sample_folder = '3a2\';
maindatafolder = "Data\";
currfolder = pwd;
id = strfind(currfolder, '\');
parentdir = currfolder(1:id(end));
regfolder = append(maindatafolder,"RegistrationOutputs\",sample_folder);
resfolder = append(maindatafolder,'Results\',sample_folder);
regfullfolderout =append(currfolder(1:id(end)),regfolder);
resfullfolder = append(currfolder(1:id(end)),resfolder,"data analysis\");

registrationresult = load(fullfile(regfullfolderout,"fullresult_cube.mat"));
%%
szintcube = registrationresult.size_interpolatedcube;
cube = registrationresult.rawcube; 
tform = registrationresult.reg.Transformation;
heights = registrationresult.heights;
mask = registrationresult.mask;

cutoffwavelength = 650;
[~, cutoffidx] =min(abs(cube.Wavelength-cutoffwavelength));
heights_scaled = heights * 10^6;
n = 1000;
tic
[dtable,coefficients] = perwave_factorized(cube,heights_scaled,cutoffidx,n,mask,tform.T,szintcube);
toc

% ===============================================
% Visualisation only 
% m4 =figure(1); 
% imshowpair(original,distorted,'blend')
% % image(colorize(transformedcube,'Method','rgb'),"CDataMapping","scaled");
% hold on 
% scatter(dtable{10}.Coordinates(:,2),dtable{10}.Coordinates(:,1),'green')
% hold off 
% xlabel("Pixels")
% ylabel("Pixels")
% ==========================================
correlationsfig  = figure(2);
scatter(coefficients(:,1),coefficients(:,2))
actn = length(dtable{1}.Height);
s = sprintf("Each point is based on %i datapoints",actn);
title(s)
xlabel("\lambda [nm]")
ylabel("\rho [-]")
grid on 
% ==================================================

logcorrelationsfig = figure;
scatter(coefficients(:,1),coefficients(:,3))
coef_Wavelength = correlation(coefficients(:,1),coefficients(:,2));
actn = length(dtable{1}.Height);
s = sprintf("Log(I): Each point is based on %i datapoints",actn);
title(s)
xlabel("\lambda [nm]")
ylabel("r [-]")
grid on 

% ==================================================
% Show best correlation wavelength
cutoff_min = 0.1;
cutoff_max = 0.9;
[~,bestwavelengthidx] = max(abs(coefficients(:,3)));
I =  dtable{bestwavelengthidx}.Intensity;
I_scaled = (I - min(min(I))) ./ (max(max(I) - min(min(I))));
I_scaled = I_scaled(I_scaled > cutoff_min & I_scaled<cutoff_max);
filtered_heights = dtable{bestwavelengthidx}.Height(I_scaled > cutoff_min & I_scaled<cutoff_max);
IntensHeightsfig =figure;
wav = str2num(dtable{bestwavelengthidx}.Properties.Description);
s2 = sprintf("At wavelength %.0f nm",cube.Wavelength(wav));
scatter(filtered_heights,I_scaled)
title(s2)
xlabel("h [\mum]")
ylabel("I [-]")
grid on 
% ============================================================
logIntensHeightfig = figure;
wav = str2num(dtable{bestwavelengthidx}.Properties.Description);
linmodel = fitlm(dtable{bestwavelengthidx}.Height,dtable{bestwavelengthidx}.("Natural log Intensity"));
c = corr(dtable{bestwavelengthidx}.Height,dtable{bestwavelengthidx}.("Natural log Intensity"),'rows','complete');
s2 = sprintf("At wavelength %.0f nm \n cor= %.3f",cube.Wavelength(bestwavelengthidx),c);
scatter(dtable{bestwavelengthidx}.Height,dtable{bestwavelengthidx}.("Natural log Intensity"))
hold on 
plot(linmodel)
hold off
title(s2)
xlabel("h [\mum]")
ylabel("log(I) [-]")
grid on 
% =================================================================


%% Saving section 

% saveas(m4,fullfile(resfullfolder,append("Chosen_Points_",num2str(n),".fig")))
% saveas(m4,fullfile(resfullfolder,append("Chosen_Points_",num2str(n),".png")))
% 
% saveas(correlationsfig,fullfile(resfullfolder,append("corr_VS_H_",num2str(n),".fig")))
% saveas(correlationsfig,fullfile(resfullfolder,append("corr_VS_H_",num2str(n),".png")))
% 
% saveas(logcorrelationsfig,fullfile(resfullfolder,append("log_corr_VS_H_",num2str(n),".fig")))
% saveas(logcorrelationsfig,fullfile(resfullfolder,append("log_corr_VS_H_",num2str(n),".png")))
% saveas(IntensHeightsfig,fullfile(resfullfolder,append("I_vs_H_wv",num2str(wav),"_",num2str(n),".fig")))
% saveas(IntensHeightsfig,fullfile(resfullfolder,append("I_vs_H_wv",num2str(wav),"_",num2str(n),".png")))
% saveas(logIntensHeightfig,fullfile(resfullfolder,append("logI_vs_H_wv",num2str(wav),"_",num2str(n),".fig")))
% saveas(logIntensHeightfig,fullfile(resfullfolder,append("logI_vs_H_wv",num2str(wav),"_",num2str(n),".png")))