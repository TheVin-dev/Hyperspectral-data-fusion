
function [datatable,coefficient] =  PerWave(hsicube,heightmap,cutoff,N,mask,tform)
    x = randi([1,size(heightmap,1)],N,1);
    y = randi([1,size(heightmap,2)],N,1);
    P = [x,y];
    sH = size(heightmap,[1,2]);
    range = [50,sH(2)-20,30,sH(1)-20];
    points = round(P,0);
    indices = (points(:,1)>range(3) & points(:,1)<range(4)) & (points(:,2)>range(1) & points(:,2)<range(2));
    
    points = points(indices,:);
    rawcube = hsicube.DataCube(:,:,1:cutoff); 
    indices_mask = sub2ind(size(mask),points(:,1),points(:,2));
    mask_filt = mask(indices_mask);
    points = points(~mask_filt,:);
    coefficient = zeros(size(rawcube,3),2);
    datatable = cell(size(rawcube,3),1);
    for wavelengthidx = 1:size(rawcube,3)
        I = zeros(length(points),1);
        D = zeros(length(points),1);
        for idx=1:length(points(:,1))
            tmpy = points(idx,2);
            tmpx = points(idx,1);
            cubeCoords  = round(tform*[tmpy,tmpx,0]');
            x_new = round(interp1([1 size(heightmap,2)], [1 size(rawcube,2)], cubeCoords(2)),0);
            y_new = round(interp1([1 size(heightmap,1)], [1 size(rawcube,1)], cubeCoords(1)),0);
            
            if ~(xor(isnan(x_new),isnan(y_new)))
            I(idx) = rawcube(y_new,x_new,wavelengthidx);
%             raw_cube_coords(idx) = [x_new,y_new];
%             I(idx) = data(points(idx,1),points(idx,2),wavelengthidx);
            D(idx) = heightmap(points(idx,1),points(idx,2));
            else 
                I(idx) = nan;
                D(idx) = nan;
            end

        end 
%         I = I(~isnan(I));
%         D = D(~isnan(I));
%         points = points(~isnan(I),:);
        
        I = (I - min(min(I))) ./ (max(max(I) - min(min(I))));
            
        log_I = log(abs(I));
        pureI = log_I(xor(~isnan(log_I),isinf(abs(log_I))));
        pureD = D(xor(~isnan(log_I),isinf(abs(log_I))));
        [coefficient(wavelengthidx,2),~] = corr(I,D,'Type','Spearman','rows','complete');
        coefficient(wavelengthidx,3) = corr(pureI,pureD,'rows','complete');
        coefficient(wavelengthidx,1) = hsicube.Wavelength(wavelengthidx);
        datatable{wavelengthidx} = table(points,I,log_I,D,'VariableNames',{'Coordinates','Intensity','Natural log Intensity','Height'});
        
        datatable{wavelengthidx}(datatable{wavelengthidx}.Height == 0, :) = [];
        datatable{wavelengthidx}.Properties.Description = num2str(wavelengthidx);
    end 
        
    
end 