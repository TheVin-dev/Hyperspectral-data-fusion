classdef vk4data
    
    properties (Access= private)
        folder = 'Data\LCM\';
    end 
    properties 
        name;
        mydir;
        filename;
        myoffset;
        A; 
        B;
        %C; not neccessary to be a property. 
        Uncoded;
        Fbase;
        F;
        rgb;
        h_norm;
        h_scaled;
        h_raw;
        LI;
        rgb_size;
        height_size;
        fileID;
        physical_size; 
        measure_conds;
        pixel_area; 
        gray;
        filtered;
    end 
    

    methods 
        function obj = vk4data(filename)
            folder = pwd;
            id = strfind(folder, '\');
            folder = append(folder(1:id(end)),obj.folder);
            obj.mydir = fullfile(folder,filename);
            obj.filename = filename;
            name = split(filename,'.');
            obj.name = name(1);
            obj = obj.load();
            obj = obj.extractMeasureConds();
            
            
        end 
        function obj = load(obj)
            fid = fopen(obj.mydir);
            obj.fileID = fid;
            obj.A = fread(fid,'ubit8');
            obj.B = dec2hex(obj.A);
            obj.myoffset = zeros(1,18);
            for myioffset=1:18
                i1=12+4*(myioffset-1)+1;
                i2=i1+1;
                i3=i1+2;
                i4=i1+3;
                obj.myoffset(myioffset)=16*hex2dec(obj.B(i1,1))+1*hex2dec(obj.B(i1,2))+256*(16*hex2dec(obj.B(i2,1))+1*hex2dec(obj.B(i2,2)))+256*256*(16*hex2dec(obj.B(i3,1))+1*hex2dec(obj.B(i3,2)))+256*256*256*(16*hex2dec(obj.B(i4,1))+1*hex2dec(obj.B(i4,2)));
            end
            

            C=obj.A(6293821:9439548,1)';
            D=dec2hex(C);
            obj.Uncoded= hex2dec(D(4,1))*16^7 +hex2dec(D(4,2))*16^6 +hex2dec(D(3,1))*16^5 + hex2dec(D(3,2))*16^4+ hex2dec(D(2,1))*16^3 + hex2dec(D(2,2))*16^2 + hex2dec(D(1,1))*16^1 + hex2dec(D(1,2))*16^0;
            obj.F=reshape(C,4,1024*768)';
            obj.Fbase=[1 16^2 16^4 16^6]';
            
        end 
        function obj = extractOptical(obj)
            typeoffset=2;
            opticalendoffset=min(obj.myoffset((obj.myoffset-obj.myoffset(typeoffset))>0));
            C=obj.A(obj.myoffset(typeoffset)+21:opticalendoffset)';
            opticalrows=obj.A(obj.myoffset(typeoffset)+1:obj.myoffset(typeoffset)+4)'*obj.Fbase;
            opticalcols=obj.A(obj.myoffset(typeoffset)+5:obj.myoffset(typeoffset)+8)'*obj.Fbase;
%           opticalcodingbase=A(myoffset(typeoffset)+9:myoffset(typeoffset)+12)'*Fbase;
%           largestopticalallowed=A(myoffset(typeoffset)+29:myoffset(typeoffset)+32)'*Fbase;
            obj.rgb_size = [opticalcols,opticalrows];
            OF=reshape(C,3,opticalrows*opticalcols)';
            myoptical=OF;
            myoptical=reshape(myoptical,[opticalrows opticalcols 3]);
            obj.rgb = uint8(permute(myoptical,[2 1 3]));
            [x,y,c] = size(obj.rgb);
            filtered = zeros(x,y,c);
            for i =1:3
                filtered(:,:,i) = medfilt2(obj.rgb(:,:,i));
            end 
            
            obj.filtered = uint8(255 * mat2gray(filtered));
            obj.gray=rgb2gray(obj.rgb);
           
        end 

	function obj = extractHeight(obj)
		typeoffset=7;
		heightendoffset=min(obj.myoffset((obj.myoffset-obj.myoffset(typeoffset))>0));
		C=obj.A(obj.myoffset(typeoffset)+797:heightendoffset)';
		heightrows=obj.A(obj.myoffset(typeoffset)+1:obj.myoffset(typeoffset)+4)'*obj.Fbase; % size in x direction
		heightcols=obj.A(obj.myoffset(typeoffset)+5:obj.myoffset(typeoffset)+8)'*obj.Fbase; % size in y direction
		%heightcodingbase=obj.A(obj.myoffset(typeoffset)+9:obj.myoffset(typeoffset)+12)'*obj.Fbase ; 
		%largestheightallowed=obj.A(obj.myoffset(typeoffset)+29:obj.myoffset(typeoffset)+32)'*obj.Fbase;

		HF=reshape(C,4,heightrows*heightcols)';
		myheight=HF*obj.Fbase;
		myheight=reshape(myheight,[heightrows,heightcols])';
		
		
        obj.h_scaled = obj.measure_conds('z_length_per_digit') * myheight;
        myheight_norm = obj.h_scaled ./max(obj.h_scaled,[],'all');
        obj.h_norm = myheight_norm; 
		obj.h_raw = myheight; 
        obj.height_size = [heightcols,heightrows];
        obj.physical_size = [obj.measure_conds('x_length_per_pixel') * obj.height_size(1), ... 
            obj.height_size(2) * obj.measure_conds('y_length_per_pixel')]; 
        obj.pixel_area = obj.measure_conds('x_length_per_pixel') * obj.measure_conds('y_length_per_pixel');
	end
    function showRGB(obj)
        
       
        figure;
        image(obj.filtered,'CDatamapping','scaled')
        daspect([1 1 1]); 
        title("Optical RGB data")
        
    end 
    function showHeight(obj) 
        figure; 
        image(obj.h_scaled,'CDataMapping','scaled');
        colormap('jet')
        c = colorbar;
        c.Label.String = "Height in meter";
        daspect([1 1 1]);
        title("Height data",'Interpreter','none')
%         figure;
%         image(obj.raw_h,'CDatamapping','scaled')
%     %         colormap('jet')
%     %         colorbar;
%         daspect([1 1 1]);
    end 
    
    function obj = extractMeasureConds(obj)
            offset = 84;
            map = containers.Map();
            id = obj.fileID;
            fseek(obj.fileID, offset,-1);
            map('size')= fread(id,1,'uint32','l');
            map('year') = fread(id,1,'uint32','l');
            map('month') = fread(id,1,'uint32','l');
            map('day')= fread(id,1,'uint32','l');
            map('hour')= fread(id,1,'uint32','l');
            map('minute')= fread(id,1,'uint32','l');
            map('second')= fread(id,1,'uint32','l');
            map('diff_from_UTC')= fread(id,1,'uint32','l');
            map('img_attributes')= fread(id,1,'uint32','l');
            map('user_interface_mode')= fread(id,1,'uint32','l');
            map('color_composite_mode')= fread(id,1,'uint32','l');
            map('img_layer_number')= fread(id,1,'uint32','l');
            map('run_mode')= fread(id,1,'uint32','l');
            map('peak_mode')= fread(id,1,'uint32','l');
            map('sharpening_level')= fread(id,1,'uint32','l');
            map('speed')= fread(id,1,'uint32','l');
            map('distance')= fread(id,1,'uint32','l');
            map('pitch')= fread(id,1,'uint32','l');
            map('optical_zoom')= fread(id,1,'uint32','l');
            map('number_of_lines')= fread(id,1,'uint32','l');
            map('line0_position')= fread(id,1,'uint32','l');
            map('reserved_1') = [fread(id,1,'uint32','l'),fread(id,1,'uint32','l'),fread(id,1,'uint32','l')];
            map('lens_magnification')= fread(id,1,'uint32','l');
            map('PMT_gain_mode')= fread(id,1,'uint32','l');
            map('PMT_gain')= fread(id,1,'uint32','l');
            map('PMT_offset')= fread(id,1,'uint32','l');
            map('ND_filter')= fread(id,1,'uint32','l');
            map('reserved_2')= fread(id,1,'uint32','l');
            map('persist_count')= fread(id,1,'uint32','l');
            map('shutter_speed_mode')= fread(id,1,'uint32','l');
            map('shutter_speed')= fread(id,1,'uint32','l');
            map('white_balance_mode')= fread(id,1,'uint32','l');
            map('white_balance_red')= fread(id,1,'uint32','l');
            map('white_balance_blue')= fread(id,1,'uint32','l');
            map('camera_gain')= fread(id,1,'uint32','l');
            map('plane_compensation')= fread(id,1,'uint32','l');
            map('xy_length_unit')= fread(id,1,'uint32','l');
            map('z_length_unit')= fread(id,1,'uint32','l');
            map('xy_decimal_place')= fread(id,1,'uint32','l');
            map('z_decimal_place')= fread(id,1,'uint32','l');
            %considered in picometers
            map('x_length_per_pixel')= fread(id,1,'uint32','l') * 10^-12;
            map('y_length_per_pixel')= fread(id,1,'uint32','l')*10^-12;
            map('z_length_per_digit')= fread(id,1,'uint32','l')*10^-12;
            map('reserved_3') =[fread(id,1,'uint32','l'),fread(id,1,'uint32','l'),fread(id,1,'uint32','l'), ... 
                fread(id,1,'uint32','l'),fread(id,1,'uint32','l')];
            map('light_filter_type')= fread(id,1,'uint32','l');
            map('reserved_4')= fread(id,1,'uint32','l');
            map('gamma_reverse')= fread(id,1,'uint32','l');
            map('gamma')= fread(id,1,'uint32','l');
            map('gamma_correction_offset')= fread(id,1,'uint32','l');
            map('CCD_BW_offset')= fread(id,1,'uint32','l');
            map('num_aperture')= fread(id,1,'uint32','l');
            map('head_type')= fread(id,1,'uint32','l');
            map('PMT_gain_2')= fread(id,1,'uint32','l');
            map('omit_color_img')= fread(id,1,'uint32','l');
            map('lens_ID')= fread(id,1,'uint32','l');
            map('light_lut_mode')= fread(id,1,'uint32','l');
            map('light_lut_in0')= fread(id,1,'uint32','l');
            map('light_lut_out0')= fread(id,1,'uint32','l');
            map('light_lut_in1')= fread(id,1,'uint32','l');
            map('light_lut_out1')= fread(id,1,'uint32','l');
            map('light_lut_in2')= fread(id,1,'uint32','l');
            map('light_lut_out2')= fread(id,1,'uint32','l');
            map('light_lut_in3')= fread(id,1,'uint32','l');
            map('light_lut_out3')= fread(id,1,'uint32','l');
            map('light_lut_in4')= fread(id,1,'uint32','l');
            map('light_lut_out4')= fread(id,1,'uint32','l');
            % considered in nano meters
            map('upper_position')= fread(id,1,'uint32','l');
            map('lower_position')= fread(id,1,'uint32','l');

            map('light_effective_bit_depth')= fread(id,1,'uint32','l');
            map('height_effective_bit_depth')= fread(id,1,'uint32','l');
            obj.measure_conds = map; 


            end 

    end 
    
    
    
end 