function T = getdata(rawcube,P,rangecoords,cutidx,mask,H,tform)
    
    inParser = inputParser;
 
    addRequired(inParser,'P')
    addRequired(inParser,'rangecoords')
    addRequired(inParser,'cutidx')
    addRequired(inParser,'mask')
    addRequired(inParser,'H')
    addRequired(inParser,'tform');
    addRequired(inParser,'rawcube');

    parse(inParser,rawcube,P,rangecoords,cutidx,mask,H,tform)


    P = round(P,0);
    indices = (P(:,1)>1 & P(:,1)<rangecoords(2)) & (P(:,2)>1 & P(:,2)<rangecoords(4));
    P = P(indices,:);
    indices_mask = sub2ind(size(mask),P(:,2),P(:,1));
    mask_filt = mask(indices_mask);        
    
    P = P(~mask_filt,:);
    avgspectrum = zeros(length(P),1);
    height = zeros(length(P),1);
    for idx=1:length(P)
            tmpx = P(idx,1);
            tmpy = P(idx,2);

            %%Extract data from raw hsi cube 
            % Convert x,y to raw hsi cube coordinate system. 
            % if col or row is nan: skip point
            cubeCoords  = round(tform*[tmpy,tmpx,0]');
            col = round(interp1([1 size(H,2)], [1 size(rawcube,2)], cubeCoords(2)),0);
            row = round(interp1([1 size(H,1)], [1 size(rawcube,1)], cubeCoords(1)),0);
            
            if ~(xor(isnan(col),isnan(row)))
            realfullspectrum(idx,:) = rawcube(row,col,1:cutidx);
            raw_cube_coords(idx,:) = [col,row];
            %%Extracting data from raw hsi cube
            avgspectrum(idx,1) = average(rawcube(row,col,1:cutidx));
            height(idx,1) = averageKernel(tmpy,tmpx,H);
            %height(idx,1) = H(tmpy,tmpx);
            end
            
            
    end

    T = table(P(:,1),P(:,2),raw_cube_coords,height,avgspectrum,realfullspectrum,'VariableNames',{'X','Y','cube_coords','height','avgspectrum','Raw spectrum'});
    
    T(T.height == 0, :) = [];
    T = sortrows(T,'height');
end


function middlePixel = averageKernel(y,x,heights)
    windowidx = 1;
    kernel=[1 1 1; 1 1 1; 1 1 1]/9;
    for i=-1:1:1 
        for j=-1:1:1
            window(windowidx) = heights(y-i,x-j);
            
            windowidx = windowidx+1;
        end
    end 
        window = reshape(window,3,3)';
        P = window*kernel;
        middlePixel = sum(P,'all') ./ length(P(:));

end

function ave = average(x)
    ave = sum(x(:))/numel(x); 
end
