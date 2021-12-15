close all;
datafolder ="Data\RegistrationInputs\";
regfolder = "Data\RegistrationOutputs\";
currfolder = pwd;
id = strfind(currfolder, '\');
folder = append(currfolder(1:id(end)),datafolder);
regfullfolder =append(currfolder(1:id(end)),regfolder);


fulloriginal = load(fullfile(folder,'lcm_roi.mat'));
original = cell2mat(fulloriginal.lcm_roi(1));
heights = cell2mat(fulloriginal.lcm_roi(2));
fulldistorted = load(fullfile(folder,'hsi_roi_cube.mat'));
cube = fulldistorted.hsi_roi_cube;
distorted = hyperpca(cube,1);
og_distorted = distorted;
distorted = imresize(distorted,size(original));

originalblur =imgaussfilt(original,5);
[fixedfeaturepoints, movingfeaturepoints] = extractPoints(distorted,originalblur);
movingreg = registerImagesSIM(distorted,originalblur,movingfeaturepoints,fixedfeaturepoints);%surfregister(distorted,original);
movingreginv = registerImagesSIM(originalblur,distorted,fixedfeaturepoints,movingfeaturepoints);%surfregister(original,distorted);
close all;

regcube = imresize(cube.DataCube,size(distorted));
regcube = imwarp(regcube, imref2d(size(regcube)), movingreg.Transformation, 'OutputView', imref2d(size(original)), 'SmoothEdges', true);

transformedcube = hypercube(regcube,cube.Wavelength,cube.Metadata);

inlierPtsDistorted = movingfeaturepoints(movingreg.inlierIdx,:);
inlierPtsOriginal  =fixedfeaturepoints(movingreg.inlierIdx,:);
mtransformfigure=figure;


regfig = subplot(2,1,1);
showMatchedFeatures(original,distorted,inlierPtsDistorted,inlierPtsOriginal,'montage');
title('Putatively matched points (inliers only)');
subplot(2,1,2)
imshowpair(original,movingreg.RegisteredImage,'blend')

featfig = figure("Visible",'off'); 
showMatchedFeatures(original,distorted,inlierPtsDistorted,inlierPtsOriginal,'montage');
title('Putatively matched points (inliers only)');
legend('Matched points 1','Matched points 2');



pairfig = figure("Visible",'off');
% imshowpair(original,movingreg.RegisteredImage,'blend')
image(imfuse(original,movingreg.RegisteredImage,'blend'),'CDatamapping','scaled');colormap('gray')
title("Overlay of the original image and the registered image")
xlabel("Pixels")
ylabel("Pixels")

% saveas(featfig,fullfile(regfullfolder,"reg_features_3A2.png"))
% saveas(featfig,fullfile(regfullfolder,"reg_features_3A2.fig"))
% 
% saveas(pairfig,fullfile(regfullfolder,"pair_overlay_3A2.png"))
% saveas(pairfig,fullfile(regfullfolder,"pair_overlay_3A2.fig"))


%%
no_samplepoints = 1000;
[~,~,~] =  EstimateERROR(original,movingreg.RegisteredImage,movingreg.Transformation.T, ... 
    movingreginv.Transformation.T,no_samplepoints,1000,'scale',2.2653e-05*10^6,'show',true,'save',true,'folder',regfullfolder);
%%
tic
data= extractData(transformedcube,cube.DataCube,heights,no_samplepoints,movingreg.Transformation.T,1);
toc
data.data.height = data.data.height * 10^6;

data.data.kernelheight = data.data.kernelheight * 10^6;

% Rescaling intensities
avgI = data.data.avgspectrum;
fullI= data.data.fullspectrum;
rawfull = data.data.('Raw spectrum');
I_scaled = (fullI - min(min(fullI))) ./ (max(max(fullI) - min(min(fullI))));
avgI_scaled = (avgI - min(min(avgI))) ./ (max(max(avgI) - min(min(avgI))));
rawfull = (rawfull - min(min(rawfull))) ./ (max(max(rawfull) - min(min(rawfull))));

data.data.avgspectrum = avgI_scaled;
data.data.fullspectrum = I_scaled;
data.data.('Raw spectrum') = rawfull;
%  Create data plots || No more data manipulation!!! 
datafolder = 'Data\Results\';
folder = pwd;
id = strfind(folder, '\');
folder = append(folder(1:id(end)),datafolder);
plotColors = jet(size(data.data.X,1));

m = figure;

% avgI_scaled = avgI_scaled(avgI_scaled > cutoff_min & avgI_scaled<cutoff_max);
cutoff_min = 0.1;
cutoff_max = 0.9;

for i=1:length(data.data.X)
    avgI_scaled = data.data.avgspectrum(i);

    if (avgI_scaled > cutoff_min && avgI_scaled<cutoff_max)
        scatter(data.data.height(i),data.data.avgspectrum(i));
        
        h(i,:) = round(data.data.height(i),0);
        hold on 
    end

end 
p = [data.data.X,data.data.Y];
s = sprintfc('%d %d',p);
hold off 
title(sprintf('Spectrum vs the depth, used %.0f points',no_samplepoints));
xlabel("Height [\mum]")
ylabel('Intensity [-]')
xlim([min(data.data.height) max(data.data.height)])
ylim([0 1])
grid on 
if length(data.data.X) <= 10
    l = legend(s,'Location','best');
    title(l,"Coordinates")
end
if size(data.data.fullspectrum,1)< 50
m2 = figure;
len_spectrum = length(data.data.fullspectrum(2,:));
for i=1:size(data.data.fullspectrum,1)
    % I_scaled = I_scaled(I_scaled > cutoff_min & I_scaled<cutoff_max);
%     I_scaled = data.data.fullspectrum(i,:);
%     I_scaled = I_scaled(I_scaled > cutoff_min & I_scaled<cutoff_max);

    plot(transformedcube.Wavelength(1:len_spectrum),data.data.('Raw spectrum')(i,:),'Color',plotColors(i,:))
    
    hold on 
end     
hold off 
title('Full spectrums for choosen coordinates');
xlabel("Wavelength [nm]")
ylabel('Intensity [-]')
ylim([0 1])
grid on 


if length(data.data.X) <= 10
    s1 = sprintfc('%d (%d,%d)',[h,p]);
    l2 = legend(s1);
    title(l2,"Height [micrometer] (Coordinates")
end

saveas(m2,fullfile(folder,append("spectrums_diffpoint_",num2str(no_samplepoints),".fig")))
saveas(m2,fullfile(folder,append("spectrums_diffpoint_",num2str(no_samplepoints),".png")))
end

m3 = figure; 
image(data.img,"CDatamapping",'scaled');colormap('gray');
hold on; 
scatter(data.data.X,data.data.Y,'green');
hold off
xlabel('Pixels')
ylabel('Pixels')


m33 = figure; 
image(data.heightimg,"CDatamapping",'scaled');colormap('jet');
hold on; 
scatter(data.data.X,data.data.Y,'black');
hold off
xlabel('Pixels')
ylabel('Pixels')
title("Coordinates on the height map")

% saveas(m,fullfile(folder,append("spec_vs_height_diffpoint_",num2str(no_samplepoints),".fig")))
% saveas(m,fullfile(folder,append("spec_vs_height_diffpoint_",num2str(no_samplepoints),".png")))
% 
% saveas(m3,fullfile(folder,append("points_diffpoint_",num2str(no_samplepoints),".fig")))
% saveas(m3,fullfile(folder,append("points_diffpoint_",num2str(no_samplepoints),".png")))
% saveas(m33,fullfile(folder,append("points_height_",num2str(no_samplepoints),".fig")))
% saveas(m33,fullfile(folder,append("points_height_",num2str(no_samplepoints),".png")))
% 

%% Analysis without taken the average of the spectrum :) 

[~, cutoffidx] =min(abs(transformedcube.Wavelength-650));
heights_scaled = heights * 10^6;
n = 1000;
[dtable,coefficients] = PerWave(cube,heights_scaled,cutoffidx,n,data.mask,movingreg.Transformation.T);

% ===============================================
% Visualisation only 
m4 =figure(1); 
imshowpair(original,distorted,'blend')
% image(colorize(transformedcube,'Method','rgb'),"CDataMapping","scaled");
hold on 
scatter(dtable{10}.Coordinates(:,2),dtable{10}.Coordinates(:,1),'green')
hold off 
xlabel("Pixels")
ylabel("Pixels")
% ==========================================
correlationsfig  = figure(2);
scatter(coefficients(:,1),coefficients(:,2))
actn = length(dtable{1}.Height);
s = sprintf("Each datapoint made with %i points",actn);
title(s)
xlabel("Wavelength [nm]")
ylabel("r [-]")
grid on 
% ==================================================

logcorrelationsfig = figure;
scatter(coefficients(:,1),coefficients(:,3))
coef_Wavelength = correlation(coefficients(:,1),coefficients(:,2));
actn = length(dtable{1}.Height);
s = sprintf("Log Intensity: Each datapoint made with %i points",actn);
title(s)
xlabel("Wavelength [nm]")
ylabel("r [-]")
grid on 

% ==================================================
I =  dtable{1}.Intensity;
I_scaled = (I - min(min(I))) ./ (max(max(I) - min(min(I))));
I_scaled = I_scaled(I_scaled > cutoff_min & I_scaled<cutoff_max);
filtered_heights = dtable{1}.Height(I_scaled > cutoff_min & I_scaled<cutoff_max);
IntensHeightsfig =figure;
wav = str2num(dtable{10}.Properties.Description);
s2 = sprintf("At wavelength %.0f nm",transformedcube.Wavelength(wav));
scatter(filtered_heights,I_scaled)
title(s2)
xlabel("Height [\mum]")
ylabel("Intensity [-]")
grid on 
% ============================================================
logIntensHeightfig = figure;
wav = str2num(dtable{10}.Properties.Description);
s2 = sprintf("Log Intensity: At wavelength %.0f nm",transformedcube.Wavelength(wav));
scatter(dtable{1}.Height,dtable{1}.("Natural log Intensity"))
title(s2)
xlabel("Height [\mum]")
ylabel("log(Intensity) [-]")
grid on 
% =================================================================
%Saving the results in figures and pictures
% saveas(m4,fullfile(folder,append("Chosen_Points_",num2str(n),".fig")))
% saveas(m4,fullfile(folder,append("Chosen_Points_",num2str(n),".png")))
% 
% saveas(correlationsfig,fullfile(folder,append("corr_VS_H_",num2str(n),".fig")))
% saveas(correlationsfig,fullfile(folder,append("corr_VS_H_",num2str(n),".png")))
% 
% saveas(logcorrelationsfig,fullfile(folder,append("log_corr_VS_H_",num2str(n),".fig")))
% saveas(logcorrelationsfig,fullfile(folder,append("log_corr_VS_H_",num2str(n),".png")))
% saveas(IntensHeightsfig,fullfile(folder,append("I_vs_H_wv",num2str(wav),"_",num2str(n),".fig")))
% saveas(IntensHeightsfig,fullfile(folder,append("I_vs_H_wv",num2str(wav),"_",num2str(n),".png")))
% saveas(logIntensHeightfig,fullfile(folder,append("logI_vs_H_wv",num2str(wav),"_",num2str(n),".fig")))
% saveas(logIntensHeightfig,fullfile(folder,append("logI_vs_H_wv",num2str(wav),"_",num2str(n),".png")))
