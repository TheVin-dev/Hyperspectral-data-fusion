function [rmse,maxerror,ncorr] = EstimateERROR(FIXED,MOVING,tform,invtform,N,iterations,varargin)
    %% 
    % Calculates the RMSE and MAX error for N random points using the
    % transformation tform and invtform between FIXED and MOVING
    %===========================================================
    % FIXED Reference image 
    % MOVING Transformed image with f 
    % tform transformation matrix for forwards propogation
    % invtform Inverse transformation matrix for backwards propogation
    % N number of points to sample 
    
    % p(x,y) on R corresponds to f(p) on S. 
    % P = f(p) 
    
    % p' = g(P), g is the transformationfrom S to R  
    %g = inv(f); 
    
    % Optional arguments with defaults: 
    %            show logical  = false;
%                scale double = 1;
%      arguments 
%         FIXED;
%         MOVING;
%         tform; 
%         invtform;
%         N double =100;
%         
%         iterations double = 100;
%         varargin;
%      end
    
    defaultshow=false;
    defaultscale = 1;
    defaultsave = false;
    defaultfolder = pwd;
    p = inputParser;
    addRequired(p,'FIXED')
    addRequired(p,'MOVING')
    addRequired(p,'tform')
    addRequired(p,'invtform')
    addRequired(p,'N')
    addRequired(p,'iterations')
    addParameter(p,'show',defaultshow);
    addParameter(p,'scale',defaultscale);
    addParameter(p,'savefigs',defaultsave);
    addParameter(p,'folder',defaultfolder);
    parse(p,FIXED,MOVING,tform,invtform,N,iterations,varargin{:})
    
    scale = p.Results.scale;
    show = p.Results.show;
    savefigs=p.Results.savefigs;
    folder = p.Results.folder;
    rmse = zeros(iterations,3);
    maxerror =zeros(iterations,1);
    coords = zeros(iterations,2);
    
    sS = size(FIXED); 
    sR = size(MOVING); 
    
    ny = min([sS(1),sR(1)]);
    nx = min([sS(2),sR(2)]);   
    % random sample N points inside images
    for iter=1:iterations
        points = cat(2,randi([1,nx],N,1),randi([1,ny],N,1));
        indices_og = sub2ind([ny,nx],points(:,2),points(:,1));
        errors = zeros(length(points),2);
        maxes = zeros(length(points),1);
        
        for i=1:length(points)
            vec = points(i,:);
            pi = [vec,0]; 
            G = pi* (invtform*tform);
            G = G(1:2); % points on the registered image MOVING
            
            if (G(1)> 1 && G(1) < nx) && (G(2)>1 && G(2)<ny)
            dis_points(i,:) = round(G,0);
            errors(i,:) = sqrt((vec - G).^2);
            maxes(i) = sum(abs(vec-G));    
            else
                continue
            end
        end 
        dis_points = dis_points(any(dis_points, 2), :);
        indices_og = indices_og(any(dis_points, 2),:);
        indices_dis = sub2ind([ny,nx],dis_points(:,2),dis_points(:,1));

        tmprmse = (1/ N) * sum(sqrt(errors.^2))* scale;
       
       	tmpmax = max(maxes) * scale;
        totalrms = sqrt(sum(tmprmse.^2,2));
        ncorr(iter) = corr2(FIXED(indices_og),MOVING(indices_dis));
        rmse(iter,:) =[tmprmse,totalrms];
        maxerror(iter) = tmpmax;
%         coords = points;
        
%         error.corr = corr;
%         error.rmse = sqrt(sum(rmse.^2)); 
%         error.maxerror = sqrt(sum(maxerror.^2));
%         error.coords = coords; 


    end
    
    if show
        i = linspace(1,iterations,iterations);
        overallcorrelation = corr2(FIXED,MOVING);
        avgrmse = average(rmse(:,3));
        accuracyfig = figure(102);
        scatter(i,rmse(:,3))
        hold on 
        plot(i,ones(iterations,1)*avgrmse)
        hold off
        s1 = sprintf(" %.0f iterations",iterations);
        title(s1)
        grid ON 
        ylabel('RMSE [\mum]')
        xlabel('iterations [-]')
        legend("",'Average')

        maxfig = figure(100);
        avgmax = average(maxerror);
        scatter(i,maxerror)
        hold on 
        plot(i,ones(iterations,1)*avgmax)
        hold off
        s1 = sprintf("%.0f iterations \ncorr: %.4f",iterations,overallcorrelation);
        title(s1)
        grid on 
        ylabel('MAX [\mum]')
        xlabel('iterations [-]')
        legend("",'Average')
        ncorrfig = figure(101);
        scatter(i,ncorr)
        s1 = sprintf("%.0f iterations",iterations);
        title(s1)
        grid on 
        ylabel('r [-]')
        xlabel('iterations [-]')
        
    end

    if savefigs
        name1 = sprintf("rmse_%0.f.png",iterations);
        name2 = sprintf("max_%.0f.png",iterations);
        name3 = sprintf("corr_%.0f.png",iterations);
        sprintf("Average rmse: %.2f \nAverage max: %.2f", avgrmse,avgmax)
        saveas(accuracyfig,fullfile(folder,name1))
        saveas(maxfig,fullfile(folder,name2))
        saveas(ncorrfig,fullfile(folder,name3))
    end 
    
end 