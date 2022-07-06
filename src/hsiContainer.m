
classdef hsiContainer 
    
    properties 
        name; 
        filename;
        mydir;
        array; 
        x_length_per_pixel;
        y_length_per_pixel;
        resized;
        collapsedimg
        length_per_pixel;
        pca;
        d_orig_hsi;
        scaling_len;
        l_orig_h_per_pixel;
        waves  = [  395.47000122   398.45999146   401.45999146   404.45999146   407.45999146...
   410.45999146   413.45999146   416.47000122   419.47000122   422.48001099...
   425.48999023   428.5          431.51000977   434.52999878   437.54000854...
   440.55999756   443.57998657   446.6000061    449.61999512   452.64001465...
   455.67001343   458.69000244   461.72000122   464.75         467.77999878...
   470.80999756   473.8500061    476.88000488   479.92001343   482.95999146...
   486.           489.04000854   492.07998657   495.13000488   498.17999268...
   501.22000122   504.26998901   507.32000732   510.38000488   513.42999268...
   516.48999023   519.53997803   522.59997559   525.65997314   528.7199707...
   531.78997803   534.84997559   537.91998291   540.98999023   544.05999756...
   547.13000488   550.20001221   553.27001953   556.34997559   559.42999268...
   562.5          565.58001709   568.66998291   571.75         574.83001709...
   577.91998291   581.01000977   584.09997559   587.19000244   590.2800293...
   593.38000488   596.4699707    599.57000732   602.66998291   605.77001953...
   608.86999512   611.9699707    615.08001709   618.17999268   621.28997803...
   624.40002441   627.51000977   630.63000488   633.73999023   636.84997559...
   639.9699707    643.09002686   646.21002197   649.33001709   652.46002197...
   655.58001709   658.71002197   661.84002686   664.9699707    668.09997559...
   671.22998047   674.35998535   677.5          680.64001465   683.7800293...
   686.91998291   690.05999756   693.20001221   696.34997559   699.48999023...
   702.64001465   705.78997803   708.94000244   712.09002686   715.25...
   718.40002441   721.55999756   724.7199707    727.88000488   731.03997803...
   734.21002197   737.36999512   740.53997803   743.71002197   746.88000488...
   750.04998779   753.2199707    756.39001465   759.57000732   762.75...
   765.92999268   769.10998535   772.28997803   775.4699707    778.65997314...
   781.84002686   785.0300293    788.2199707    791.40997314   794.60998535...
   797.79998779   801.           804.19000244   807.39001465   810.59002686...
   813.79998779   817.           820.21002197   823.40997314   826.61999512...
   829.83001709   833.03997803   836.25         839.4699707    842.67999268...
   845.90002441   849.11999512   852.34002686   855.55999756   858.78997803...
   862.01000977   865.23999023   868.4699707    871.70001221   874.92999268...
   878.15997314   881.40002441   884.63000488   887.86999512   891.10998535...
   894.34997559   897.59002686   900.84002686   904.08001709   907.33001709...
   910.58001709   913.83001709   917.08001709   920.33001709   923.59002686...
   926.84002686   930.09997559   933.35998535   936.61999512   939.88000488...
   943.15002441   946.40997314   949.67999268   952.95001221   956.2199707...
   959.48999023   962.76000977   966.0300293    969.30999756   972.59002686...
   975.86999512   979.15002441   982.42999268   985.71002197   989. ...
   992.28997803   995.58001709   998.86999512  1002.15997314  1005.45001221...
  1008.75      ]';
    end 
    
    methods 
        function obj = hsiContainer(filename)
            obj.filename = filename;
            obj.mydir = fullfile(filename);
            
            obj = obj.load();
            obj = obj.PCA(1);
            obj = obj.extractCollapsed(obj.array.DataCube);
        end        
        function obj = load(obj)
            data = h5read(obj.mydir,"/Array");
            tmparray = permute(data,[1,3,2]);
            obj.array = hypercube(permute(tmparray,[2,3,1]),obj.waves);
        end
        function obj = extractCollapsed(obj,array)
            res = 0;
            hsidata= array;
            s = size(hsidata);
            
            for i=1:s(1)
                window = hsidata(i,:,:);
                img = window + res;
            end 
            obj.collapsedimg = squeeze(img);
            
        end 
        function obj = showFigure(obj,handle,strtitle,array)
                
                figure(handle);
                image(array,'CDatamapping','scaled'); colormap("gray"); 
                daspect([1 1 1])
                title(strtitle);
                xlabel("pixel")
                ylabel('pixel')
                
        end 
        function showOneWave(wavelength)
                 [~, index] = min(abs(obj.waves-wavelength));
                 img = squeeze(obj.array(index,:,:));
                 f1 = figure("oneWave");
                 image(img,'CDatamapping','scaled'); colormap("gray"); 
                 title(["One wavelength data @ ", num2str(obj.waves(index))]);
        end 
        
        function obj =PCA(obj,num)
            arguments 
                obj ; 
                num single = 1; 
            end 
            obj.pca = double(hyperpca(obj.array,num));
            %obj.pca = uint8(round(tmppca - 1));
        end 
        function showSpectrum(obj)
                % TODO: plot spectrum  average of nxn at point(x,y)
                % for one point:
                obj.showFigure(3,"First PCA component", obj.pca);
                roi_p = drawpoint; 
                x = round(roi_p.Position(1));
                y = round(roi_p.Position(2));
                
                spectraldata = squeeze(obj.array.DataCube(y,x,:));
                figure(2)
                plot(obj.array.Wavelength,spectraldata)
                strtitle = sprintf("Spectrum at point %d,%d",x,y);
                xlabel('Wavelength [nm]')
                ylabel('Spectral intensity')
                title(strtitle)
                grid ON 

        end
        function obj = extractROI(obj,length_lcm,lengthpixellcm)
            roihsi = drawline;
            PositionHSI = round(roihsi.Position)';
            pt1 = PositionHSI(1:2);
            pt2 = PositionHSI(3:4);

            diff_vec = abs(pt2 - pt1);
            l_lcm = length_lcm; %0.0157; % length of pixels in cm 
            l_orig_h_per_pixel = l_lcm ./ diff_vec; 
            fprintf("Length per pixel hsi original x: %.2d,y:%2d m \n",l_orig_h_per_pixel) 
            
            d_orig_hsi = l_orig_h_per_pixel .* size(obj.collapsedimg);
            fprintf("Total size hsi x: %.2d,y:%2d m\n",d_orig_hsi) 
            
            scaling_length = lengthpixellcm./ l_orig_h_per_pixel; 
            fprintf("scaling_length x: %.2d,y:%2d\n",scaling_length) 
           
            obj.resized = imresize(obj.pca,'Scale',scaling_length);
            obj.d_orig_hsi = d_orig_hsi;
            obj.scaling_len = scaling_length;
            obj.l_orig_h_per_pixel = l_orig_h_per_pixel;
        end
    end 
end 