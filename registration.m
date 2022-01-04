
maindatafolder = "Data\";
currfolder = pwd;
id = strfind(currfolder, '\');
parentdir = currfolder(1:id(end));

[lcmfile,lcmpath] = uigetfile(".mat","Select LCM file",append(parentdir,"RegistrationInputs\"));
[hsifile,hsipath] = uigetfile(".mat","Select HSI file",append(parentdir,"RegistrationInputs\"));


sample_folder = '3a2\';
reginfolder = append(maindatafolder,"RegistrationInputs\",sample_folder);
regfolder = append(maindatafolder,"RegistrationOutputs\",sample_folder);
resfolder = append(maindatafolder,'Results\',sample_folder);
regfullfolderin = append(parentdir,reginfolder);
regfullfolderout =append(currfolder(1:id(end)),regfolder);
resfullfolder = append(currfolder(1:id(end)),resfolder);

fulloriginal = load(fullfile(lcmpath,lcmfile));
original = cell2mat(fulloriginal.lcm_roi(1));
original = medfilt2(original,[5,5]);
heights = cell2mat(fulloriginal.lcm_roi(2));
fulldistorted = load(fullfile(hsipath,hsifile));
cube = fulldistorted.hsi_roi_cube;
distorted = hyperpca(cube,1);
og_distorted = distorted;
distorted = imresize(distorted,size(original));

%%
originalblur =imgaussfilt(original,5);
[fixedfeaturepoints, movingfeaturepoints,allFixedfeatures,allMovingfeatures] = extractPoints(distorted,originalblur);

if length(fixedfeaturepoints) <=1 | length(movingfeaturepoints) <=1 
    figure(1000)
    fprintf("There are not enough extracted features for registration \n")
    
end
% registerImages2b1(distorted,originalblur);%
% registerImages2b1inv(originalblur,distorted);%
movingreg = registerImagesAFF(distorted,original);%
movingreginv = registerImagesSIM(original,distorted,movingfeaturepoints,fixedfeaturepoints);%


regcube = imresize(cube.DataCube,size(distorted));
regcube = imwarp(regcube, imref2d(size(regcube)), movingreg.Transformation, 'OutputView', imref2d(size(original)), 'SmoothEdges', true);

transformedcube = hypercube(regcube,cube.Wavelength,cube.Metadata);

% inlierPtsDistorted =movingfeaturepoints(movingreg.inlierIdx,:);
% inlierPtsOriginal  = fixedfeaturepoints(movingreg.inlierIdx,:);


imshowpair(movingreg.RegisteredImage,original,'blend')
%%
grayImagepca = distorted;
level = adaptthresh(grayImagepca);
mask =imbinarize(grayImagepca,level);

fullresult = struct; 
fullresult.reg = movingreg; 
fullresult.reginv = movingreginv; 
fullresult.heights = heights; 
fullresult.rawcube = cube; 
fullresult.size_interpolatedcube = size(transformedcube.DataCube);
fullresult.mask = mask;
name1 = strsplit(hsifile,'.');
name2 = strsplit(name1{1},'_');
namesample = name2{3};
s = sprintf("fullresult_%s",namesample);
save(fullfile(regfullfolderout,s),'-struct',"fullresult")
%%
no_samplepoints = 1000;
iterations = 1000;
registration_accuracyfolder = append(resfullfolder,'registration');
[~,~,~] =  EstimateERROR(original,movingreg.RegisteredImage,movingreg.Transformation.T, ... 
    movingreginv.Transformation.T,no_samplepoints,iterations,'scale',1,'show',true,'save',false,'folder',registration_accuracyfolder);


