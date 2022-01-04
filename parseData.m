
currfolder = pwd;
id = strfind(currfolder, '\');
parentdir = currfolder(1:id(end));
parentdir = append(parentdir, "Data\");
parentdirlcm = append(parentdir, "lcm\");
parentdirhsi = append(parentdir, "hsi\");

[lcmfile,lcmpath] = uigetfile("*.vk4",'Select raw LCM vk4 data',parentdirlcm);
[hsifile,hsipath] = uigetfile("*.hdf5",'Select raw HSI hdf5 data',parentdir);

if all(lcmfile == 0) | all(hsifile ==0)
    fprintf("Did not select files \n")
    return
end
%%
hsidata = hsiContainer(append(hsipath,hsifile));
lcmdata = vk4data(lcmfile,lcmpath);

hsidata = hsidata.PCA(1); 
lcmdata = lcmdata.extractOptical();
lcmdata = lcmdata.extractHeight();

lcmdata.showRGB() 
roi_l = drawrectangle; 
Positionl = round(roi_l.Position);
lcm_roi_img = imcrop(rgb2gray(lcmdata.filtered),Positionl);
heights  =imcrop(lcmdata.h_scaled_corrected,Positionl);
lcm_roi = {lcm_roi_img,heights};%cat(3,lcm_roi_img,heights);
hsidata.showFigure(19,"First PCA component", hsidata.pca);
roi_h = drawrectangle;
Positionh = round(roi_h.Position);
r = Positionh;
hsi_roi_cube = cropData(hsidata.array,r(2):r(2)+r(4),r(1):r(1)+r(3));
close all
%%
s = inputdlg("Enter the pigment name",'pigment name');
reginfolder = 'RegistrationInputs\';
samplefolder= s{1};
folder = append(parentdir,reginfolder,samplefolder);
if not(isempty(s))
    namehsi = fullfile(folder,append('hsi_roi_',s{1},".mat"));
    namelcm = fullfile(folder,append('lcm_roi_',s{1}, ".mat"));
    save(namehsi,"hsi_roi_cube",'-mat')
    save(namelcm,"lcm_roi",'-mat')

end

