function [datatable,coefficient] =  perwave_factorized(hsicube,heightmap,cutoff,N,mask,tform,sztranscube)
    x = randi([1,size(heightmap,1)],N,1);
    y = randi([1,size(heightmap,2)],N,1);
    P = [x,y];
    sH = size(heightmap,[1,2]);
    range = [50,sH(2)-20,30,sH(1)-20];
    points = round(P,0);
    indices = (points(:,1)>range(3) & points(:,1)<range(4)) & (points(:,2)>range(1) & points(:,2)<range(2));
    
    points = points(indices,:);
    rawcube = hsicube.DataCube(:,:,1:cutoff); 
    indices_mask = sub2ind(sztranscube([1 2]),points(:,1),points(:,2));
    mask_filt = mask(indices_mask);
    points = points(~mask_filt,:);
    coefficient = zeros(size(rawcube,3),2);
    datatable = cell(size(rawcube,3),1);
    P = [points,zeros(length(points),1)];
    cubeCoord2 = P*tform;
    x_newlist =  round(interp1([1 size(heightmap,2)], [1 size(rawcube,2)], cubeCoord2(:,2)),0);
    y_newlist =  round(interp1([1 size(heightmap,1)], [1 size(rawcube,1)], cubeCoord2(:,1)),0);
	P_transformed = [x_newlist,y_newlist];
    nanind = any(isnan(P_transformed),2);
    P_transformed = P_transformed(~nanind,:);
    for wavelengthidx = 1:size(rawcube,3)
        % Refactor effort, not working yet 19-21-2021
        linearnanIndicesC = sub2ind(size(rawcube,[1,2,3]),P_transformed(:,2),P_transformed(:,1),wavelengthidx*ones(length(P_transformed(:,1)),1));
        linearnanIndicesH = sub2ind(size(heightmap),P_transformed(:,2),P_transformed(:,1));
        I_list = rawcube(linearnanIndicesC);
        D_list = heightmap(linearnanIndicesH);
        log_I = log(abs(I_list));
        nanmask = ~xor(isnan(log_I),isinf(abs(log_I)));
        log_I = log_I(nanmask);
        I_list = I_list(nanmask);
        D_list= D_list(nanmask);
        I_list = (I_list - min(min(I_list))) ./ (max(max(I_list) - min(min(I_list))));

        [coefficient(wavelengthidx,2),~] = corr(I_list,D_list,'Type','Spearman','rows','complete');
        coefficient(wavelengthidx,3) = corr(log_I,D_list,'rows','complete');
        coefficient(wavelengthidx,1) = hsicube.Wavelength(wavelengthidx);
        
        
        corrpoints = points(nanmask,:);
        datatable{wavelengthidx} = table(corrpoints,I_list,log_I,D_list,'VariableNames',{'Coordinates','Intensity','Natural log Intensity','Height'});
        
        datatable{wavelengthidx}(datatable{wavelengthidx}.Height == 0, :) = [];
        datatable{wavelengthidx}.Properties.Description = num2str(wavelengthidx);
    end 
        
    
end 