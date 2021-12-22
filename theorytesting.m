%% Testing the theory 
% I_lambi_i = I_g1,lambda_i * exp(-2*k_i * height)
% 
% 
% 
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
rawcube = registrationresult.rawcube.DataCube;
tform = registrationresult.reg.Transformation.T;
mask = registrationresult.mask;
[yCoordinates, xCoordinates] = find(mask);
invmask = (~mask);
Pmask = [xCoordinates,yCoordinates,zeros(length(xCoordinates),1)] *tform;
Pmask = round(Pmask,0);
[yCoordinatesINV, xCoordinatesINV] = find(invmask);
Pmaskinv = [xCoordinatesINV,yCoordinatesINV,zeros(length(yCoordinatesINV),1)] *tform;
Pmaskinv = round(Pmaskinv,0);

x_t =  round(interp1([1 size(heights,2)], [1 size(rawcube,2)], Pmask(:,1)),0);
y_t =  round(interp1([1 size(heights,1)], [1 size(rawcube,1)], Pmask(:,2)),0);
clear Pmask
Pmask = [x_t,y_t];

x_tinv =  round(interp1([1 size(heights,2)], [1 size(rawcube,2)], Pmaskinv(:,1)),0);
y_tinv =  round(interp1([1 size(heights,1)], [1 size(rawcube,1)], Pmaskinv(:,2)),0);
clear Pmaskinv;
Pmaskinv = [x_tinv,y_tinv];
heights = registrationresult.heights;
avgI_lambds = zeros(84,1);
I_iavg = zeros(84,1);
k_i = zeros(84,1);

image(heights,'CDatamapping','scaled')
colorbar
[x,y] = ginput(1);
x = round(x,0);
y = round(y,0);
close;
%%
for waveband = 1:180 
    lin_indicesmask = sub2ind(size(rawcube),Pmask(:,2),Pmask(:,1),waveband*ones([length(yCoordinates),1]));
    lin_indicesINVmask = sub2ind(size(rawcube),Pmaskinv(:,2),Pmaskinv(:,1),waveband*ones([length(yCoordinatesINV),1]));

    lin_indicesINVmask = lin_indicesINVmask(~isnan(lin_indicesINVmask));
    lin_indicesmask = lin_indicesmask(~isnan(lin_indicesmask));
    I_i = rawcube(lin_indicesINVmask);
    I_iavg(waveband) = average(I_i); 
    L_groundlayer_i  =rawcube(lin_indicesmask);
    
    avgI_lambds(waveband) = average(L_groundlayer_i);
    Iterm = log(I_iavg(waveband) /avgI_lambds(waveband) );

   
    delta = heights(y,x);
    k_i(waveband) = -(1/2*delta) * Iterm;

    if mod(waveband,2) == 0
        fprintf("For waveband %.0f, the average intensity of the ground is: %.3f \n",waveband,avgI_lambds(waveband))
        figure(1)
        hold on 
        plot(L_groundlayer_i)
        figure(2)
        hold on 
        plot(I_i)

        drawnow()

    end


end


diffavgI = avgI_lambds -I_iavg;
figure(3)
subplot(1,3,1)

plot(avgI_lambds)
title("I_{groundlayer_{\lambda_{i}}} average")
grid on 

subplot(1,3,2)

plot(I_iavg)
title("I_{paint}_{\lambda_{i}} average")
grid on 

subplot(1,3,3)
plot(k_i)
title("Calculated k_{\lambda_{i}}")
grid on 