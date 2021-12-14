%@Author: Vincent den Ronden
% Own implementation of a simple intensity-based registration algorithm. 
global img1 img2

datafolder = 'Data\Results\';
folder = pwd;
id = strfind(folder, '\');
folder = append(folder(1:id(end)),datafolder);
name = "OwnRegistration_corr_"; 
form = '.fig';


inp = input("Which image for registration? \n Default is checkerboard. \n Choose (C) or (T) ",'s');
switch inp 
    case 'T'
        img =rgb2gray(imread('test.jpg'));
        typ = "Test";
    otherwise 
       
        n=100;
        cols = 3;
        img = checkerboard(n,cols); 
        circleGrayLevel = 0.8; 
        centerx= cols*n; 
        centery = cols*n;
        imageSizeX = size(img,1);
        imageSizeY = size(img,2);
        [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
        radius = n/4;
        circlePixels = (rowsInImage - centery).^2 ...
            + (columnsInImage - centerx).^2 <= radius.^2;
        typ = "Checker";
        img(circlePixels) = circleGrayLevel;
end 
fullsave = fullfile(folder,append(name,typ));
% define transform parameters
tx = 5; 
ty = 10; 
% phi =10;
img2 = imtranslate(img,[tx,ty]);
img1 = imresize(img,size(img2));%reshapeimg(img,img2);

tic
[ximg,param] = Registration(img2,img1);
[~,paraminv] = Registration(img1,img2);
fprintf("Final guess: %2f %2f \n ", param(1), param(2))
T = [1 0 param(1) ; 0 1 param(2); 0 0 1];
invT = [1 0 paraminv(1) ; 0 1 paraminv(2); 0 0 1];
[error,~,~,~] = EstimateERROR(img1,img2,T,invT,10000,1);
toc
errlist = num2cell([[tx,ty]+ param,sqrt(sum(([tx,ty]+ param).^2))]);
[acterrorx,acterrory,actradius] = errlist{:}; 

%%
close all 

subplot(1,3,1)
image(img1,"CDataMapping","scaled");colormap('gray')
daspect([1.5 1.5 1]) 
title('Fixed image')
xlabel('pixels')
ylabel('pixels')
subplot(1,3,2)
image(img2,"CDataMapping","scaled");colormap('gray');
daspect([1.5 1.5 2]) 
title('Moving image')
xlabel('pixels')
ylabel('pixels')
subplot(1,3,3)
image(ximg,"CDataMapping","scaled");colormap('gray');
daspect([1.5 1.5 2]) 
s = sprintf("Correlation: %.2f \n ",error.corr);
serror = sprintf("Actual transform: (tx,ty) = (%.3f,%.3f) \nBest solution: (tx,ty) = (%.3f,%.3f) \nRMSE = %.3f",tx,ty,param(1),param(2),error.rmse);
title('Registered image', s)
xlabel('pixels')
ylabel('pixels')
a = annotation('textbox', [0.269571428571429,0.742857142857143,0.481857142857143,0.201904761904767], 'string', serror);
a.EdgeColor='none';
saveas(gcf,append(fullsave,form))
saveas(gcf,append(fullsave,".png"))



function rimg = reshapeimg(ogimg,targetimg)
    if all(size(ogimg)==size(targetimg))
        rimg = ogimg;
        return 
    end
    t = class(targetimg);
    size_diff = abs(size(targetimg)-size(ogimg));
    diffx = size_diff(1) /2; 
    diffy = size_diff(2)/2;
    rimg = zeros(size(targetimg),t);
    
    rimg(diffx:size(targetimg,1)-diffx-1,diffy:size(targetimg,2)-diffy-1) = ogimg;

end 
function err = measureError(ximg,ogimg,varargin)
    % img1 transfomred
    % img2 og 
    if nargin >2 
        method = varargin{1};   
    else 
        method = "SSD";
    end 
    if ~(all(size(ogimg) == size(ximg)))
        img = imresize(ogimg,size(ximg));%reshapeimg(ogimg,ximg);
    else 
        img = ogimg;
    end 
    switch method 
        case "correlation" 
            err = -1*correlation(ximg,img);
        case "SSD"
            diff = (ximg-img);
            err = sum(diff.^2,'all');
        case "MI" % https://jinjeon.me/post/mutual-info/
            err = MI(ximg,img);
    end
    
end

function err = ErrFunc(param) 
global img2 img1

    % Function to calculate the error after transforming the image
    % according to param 
    tmpimg2 = imtranslate(img2,param(1:2));
    err = measureError(tmpimg2,img1,'correlation');
end 

function [img,param] = Registration(distorted,ogimg) 
    % perform first guess 
    x =linspace(-100,100,1000);
    phi =-30:1:30;
    x_errors = zeros(size(x));
    y_errors = zeros(size(x));
    phi_errors = zeros(size(phi));
    for i=1:length(x)
        xi = x(i); 
        x_errors(i) = measureError(imtranslate(distorted,[xi 0]),ogimg,'correlation');    
        y_errors(i) = measureError(imtranslate(distorted,[0 xi]),ogimg,'correlation');
       
    end 
    for i=1:length(phi)

     phi_errors(i) = measureError(imrotate(distorted,phi(i),'crop'),ogimg,'correlation');

    end
    [~, idX] = min(x_errors);
    [~,idY] = min(y_errors);
%     [~,idP] = min(phi_errors);
    fprintf("Initial guess: %d %d  \n",x(idX),x(idY));
    param = [x(idX),x(idY)];
    
    
    param = fminsearch(@ErrFunc,param);
    
    img = imtranslate(distorted,param(1:2));

end 