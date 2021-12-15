function [data] = extractData(hypercube,rawcube, heightmap, N,tform,autoflag)
% TODO: Add manually point(s) selection 
% Hypercube: contains all data from the hypercube of the sample 
% heigtmap: 3d data containing the recorded height for every pixel in the
% image 
% Transform: The calculated transformation to register lcm with hypercube 
% kernelsize: To average a square of n (must be odd) pixels due to
% inaccuarcy measures. The average isnt of the pixel values but of the
% spectral data and height data. Basically, averaging the square of columns
% centered on (x,y) 

% windowsize = kernelsize ;
% kernel = ones(windowsize);
% kernel = kernel / size(kernel(:));
arguments 
    hypercube;
    rawcube;
    heightmap;
    N; 
    tform;
    
    autoflag (1,1) logical = 1;
    
end 
    datacube = hypercube.DataCube;
    data =struct;
    data.img = colorize(hypercube,"Method","rgb");
    data.heightimg = heightmap;
    sH = size(datacube,[1,2]);
    range = [1,sH(2),1,sH(1)];
    [~, cutoffidx] =min(abs(hypercube.Wavelength-650));
    grayImage = hyperpca(datacube,1);
%     figure;
%     image(grayImage,"CDataMapping","scaled");colormap('gray');
%     daspect([1 1 1])
%     [x,y] = ginput(4);
%     close;
%     filter_points = round([x,y],0);
    invtform = inv(tform);
    % 1 intensity grayImage(xi,yi);
%     for p=1:length(filter_points)
%         filt_int = grayImage(filter_points(p,2),filter_points(p,1));
%     end 
    
%     filt = average(filt_int);
    level = graythresh(grayImage);
    mask =imbinarize(grayImage,level);
    data.mask = mask;

if autoflag

    x = randi([1,sH(2)],N,1);
    y = randi([1,sH(1)],N,1);
    P = [x,y,ones(N,1)]; % original coordinates
    p = P * tform;
    p = round(p(:,1:2)); % Registered coordinate
    
    data.data = getdata(datacube,p,range,cutoffidx,mask,heightmap,tform,rawcube);
    
else 
    figure;
    image(data.img,"CDataMapping","scaled");colormap('gray');
    daspect([1 1 1])
    i=1;
    while i<10
        point = drawpoint;
        points(i,:) = point.Position;
        i= i+1;
    end
    close;
    data.data = getdata(datacube,points,range,cutoffidx,heightmap,invtform);
    
end 




end 

