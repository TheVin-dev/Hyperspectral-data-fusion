clear all 
filename2 = "3A2_stitched.vk4";
data3A2_full = vk4data(filename2);
filenamehsi = "3A2_2_Normalized.hdf5";
A1 = hsiContainer(filenamehsi);
A1 = A1.PCA(1); 
data3A2_full = data3A2_full.extractOptical();
data3A2_full = data3A2_full.extractHeight();
data3A2_full.showRGB() 
roi_l = drawrectangle; 
Positionl = round(roi_l.Position);
lcm_roi_img = imcrop(rgb2gray(data3A2_full.filtered),Positionl);
heights  =imcrop(data3A2_full.h_scaled,Positionl);
lcm_roi = {lcm_roi_img,heights};%cat(3,lcm_roi_img,heights);
A1.showFigure(19,"First PCA component", A1.pca);
roi_h = drawrectangle;
Positionh = round(roi_h.Position);
r = Positionh;
hsi_roi_cube = cropData(A1.array,r(2):r(2)+r(4),r(1):r(1)+r(3));
close all
s = inputdlg("Enter the pigment name",'pigment name');
datafolder = 'Data\RegistrationResults\';
folder = pwd;
id = strfind(folder, '\');
folder = append(folder(1:id(end)),datafolder);
if not(isempty(s))
    measurement = strsplit("3A2_2_Normalized.hdf5","_");
    measurement = measurement{2};
    namehsi = fullfile(folder,append('hsi_roi_',s{1},"_",measurement,".mat"));
    namelcm = fullfile(folder,append('lcm_roi_',s{1}, ".mat"));
    save(namehsi,"hsi_roi_cube",'-mat')
    save(namelcm,"lcm_roi",'-mat')

end


