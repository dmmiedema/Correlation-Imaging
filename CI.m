%% Correlation Imaging (CI) analysis to extract motility of particles moving along straight segments
% Daniël M. Miedema, 2017, Amsterdam, NL

%% Description
% This program calculates the velocity and run length of particles moving
% along a straight segment, as well as the (1D) density of particles on the
% segment. The method is based on spatiotemporal correlations of the
% intensities from the moving objects. It has been developed to measure the
% density, velocity and run length of fluorescently labeled molecular motors moving along
% microtubule segments, and is applicable also at high motor density.

% The main function is divided in 4 analysis parts:
% 1) Segment selection from an image sequence and obtaining the intensity
%   along the segment over time. (uses sub-function 'Im2In').
% 2) Determine the particle density along the segment from the spatially
%   average intensities using fluorescence correlation spectroscopy.
% 3) Calculate the spatiotemporal correlation of intensities. (uses
%   sub-function 'correlation')
% 4) Derive particle velocity and run length from analysis of the
%   evolution of the correlation peak. (uses sub-function 'motility').

%% Input
% directory:        String name of the directory containing the image sequence of interest as separate files (.tiff format).

%% Output
% vel               Velocity (micrometers/second).
% rl                Run length (micrometers).
% dens              Density (particles/micrometer).

%% Sub-functions
% Im2In             Reads a sequence of .tif files from the specified.
%   directory. Allows for the manual selection of a segment and for manual selection of the time window. 
%   Gives the intensity along the segment as output.
% correlation       Calculates the correlation of intensities for a range
%   of temporal and spatial separations.
% motility          Extracts the velocity and run length from the evolution
%   of the correlation peak.
% ----- The sub-functions can be found below the main function -----

%% Main Function
function [dens,vel,rl]=CI(directory)
    % Experimental parameters
    resolution=0.08;                                        % Resolution of images (micrometers/per pixel).
    dt_frame=0.07;                                          % Time between images (seconds).   
    N_pix=512;                                              % Size of (square) image (pixels).
    
    % 1) Read image sequence, select segment and obtain intensity vs. time along the segment.
    [I2,pix_dist]=Im2In(directory,N_pix);
    
    % 2) Calculate density on segment from spatially averaged intensity, using fluorescence correlation spectroscopy.
    seg_length=(size(I2,2)-1)*pix_dist;
    Imean=mean(I2,2);                 
    G2=var(Imean)/(mean(Imean))^2;
    dens=1./(seg_length*G2*resolution);                     % Density (particles/micrometer).
    
    % 3) Calculate correlation of intensity fluctuations
    dx_min=-10;                                             % Select range of dt and dx over which to calculate correlation (in integer values,number of pixels and number of frames, respectively). 
    dx_max=50;    
    dt_min=1;                                                               
    dt_max=100;
    dt_min_fit=2;                                           % Minimum dt value used for Gaussian fitting.
    dt_max_fit=12;                                          % Maximum dt value used for Gaussian fitting: e.g. dt_max_fit = 2 * expected run time of particles /dt_frame. 
    corr0=correlation(I2,dt_min,dt_max,dx_min,dx_max);      % Calculate correlation.
    dt_long_term=dt_max_fit*2;                              % Dt value from which to calculate the long-term correlation.
    corr_bg=mean(corr0(dt_long_term:end,:));                % Calculate long-term correlation.
    corr_bg_ar=repmat(corr_bg,size(corr0,1),1);             
    corr=corr0-corr_bg_ar;                                  % Subtract long-term correlation from correlation at all times to remove background from static objects.
    
    % 4) Calculate velocity and run length from spatiotemporal correlation.
    [vel,rl]=motility(corr,dx_min,dt_min_fit,dt_max_fit,resolution,dt_frame,pix_dist);
end

%% Sub-function Im2In: used to select segment from image sequence and obtain intensity along segment
% -------------- Input ----------------------------  
% directory         Directory of image containing the image sequence.
%   images are expected to be separate and ordered tif files, recognized
%   by .tif extension.
% Npix              Size of (square) image in pixels.
% ------------- Output ---------------------------
% I2                2D Array of intensities along segment (colums) versus
%   time (rows).
% pix_dist          Average distance between pixels along segment in pixel
%   units. To correct for the tilt a segment might have to the horizontal
%   or vertical image axis.
% ------------ Function -------------------------
function [I2,pix_dist]=Im2In(directory,N_pix)
    % Parameters
    stack_min=1;                                            % First image to read in for detection of motion
    stack_movie=50;                                        % Last image to read in for detection of motion
    margin=10;                                              % Margin from selected points included in region of interest containing the segment
    
    % Create movie of first stack_max_pre images to find microtubule
    im_name=strcat(directory,'\img*');
    im=dir(im_name);
    stack_max=length(im);                                   % Number of images in sequence to analyze
    if stack_movie>stack_max                                % Check if total image sequence is longer than the original movie to display.
        stack_movie=stack_max;
    end
    fprintf('Loading initial part of image sequence..\n');
    a1=zeros(N_pix,N_pix,stack_movie-stack_min+1, 'double');% Initialize array of images 
    for i=stack_min:stack_movie                             % Loop through images from stack_min to stack_max
        file=strcat(directory,'\',im(i).name);
        a1(:,:,i-stack_min+1)=double(imread(file));
    end
    fig_handle=figure(1); fig_handle.Units='normalized'; fig_handle.OuterPosition=[0 0.5 1/3 0.5];
    ax=gca; ax.FontSize=12;                     % Displays (first) part of image sequence
    title('Time serie of image sequence'); 
    for i=1:stack_movie-stack_min+1
        imagesc(a1(:,:,i)); 
        ax=gca; ax.XTick=[]; ax.YTick=[]; 
        pause(0.05);
    end
    ax=gca; ax.FontSize=12; ax.XTick=[]; ax.YTick=[]; title('Time serie of image sequence');
    
    % Select segment for analysis from average image: rough determination of begin and end points.
    fig_handle=figure(2); fig_handle.Units='normalized'; fig_handle.OuterPosition=[1/3 0.5 1/3 0.5];
    hold off; imagesc(mean(a1,3));               % Display average image             
    ax=gca; ax.FontSize=12; ax.XTick=[]; ax.YTick=[]; 
    title('Select begin point of segment by clicking on image');
    fprintf('Select begin point of segment by clicking on image\n');  
    [x1,y1]=ginput(1);
    figure(2); hold on; plot(x1,y1,'ok','MarkerSize',5);    % Plot selected begin point of segment
    title('Select end point of segment by clicking on image');
    fprintf('Select end point of segment by clicking on image\n');
    [x2,y2]=ginput(1);
    figure(2); hold on; plot(x2,y2,'ok','MarkerSize',5);    % Plot selected end point of segment
    
    % Display region containing selected segment at higher resolution
    if x1<x2
        x_low=round(max(x1-margin,1));
        x_high=round(min(x2+margin,N_pix));
    else
        x_low=round(max(x2-margin,1));
        x_high=round(min(x1+margin,N_pix));
    end
    if y1<y2
        y_low=round(max(y1-margin,1));
        y_high=round(min(y2+margin,N_pix));
    else
        y_low=round(max(y2-margin,1));
        y_high=round(min(y1+margin,N_pix));
    end
    fprintf('Loading complete image sequence..\n');
    a1=zeros(y_high-y_low+1,x_high-x_low+1,stack_max-stack_min+1,'double');
    for i=stack_min:stack_max                               % Read complete image and store for the region of interest
        file=strcat(directory,'\',im(i).name);
        a0=double(imread(file));
        a1(:,:,i-stack_min+1)=a0(y_low:y_high,x_low:x_high);
    end
    
    % Refine selection of begin and end points of segment
    status='n';       
    while status=='n'                                       % Show selected segment in image and adjust until satisfied.
        figure(2); hold off; imagesc(mean(a1,3)); hold on;
        ax=gca; ax.FontSize=12; ax.XTick=[]; ax.YTick=[]; 
        title('Select refined begin point of segment by clicking on image');
        fprintf('Select refined begin point of segment by clickin on image\n');
        [x1,y1]=ginput(1);
        x1=round(x1); y1=round(y1);
        plot(x1,y1,'ok','MarkerSize',5);
        title('Select refined end point of segment by clicking on image');
        fprintf('Select refined end point of segment by clicking on image\n');
        [x2,y2]=ginput(1);
        x2=round(x2); y2=round(y2);    
        plot(x2,y2,'ok','MarkerSize',5);
        plot([x1 x2],[y1 y2],'k-','LineWidth',2);
        title('Average image with selected segment');
        % Check if ok
        choice = input('Is the segment selection OK? (If not, press n followed by enter. Otherwise press enter)','s');        
        switch choice % Handle response
            case 'n'                                                    
            otherwise    
                status='GOOD';
        end
    end
    
    % Determine the pixel positions along to the segment
    seg_length=sqrt((x2-x1)^2+(y2-y1)^2);                   % Length of segment
    alpha_pre=mod(atan2(y2-y1,x2-x1),pi);                   % Angle of segment
    if alpha_pre>pi/2, alpha=pi-alpha_pre; else alpha=alpha_pre; end
    if alpha<pi/4, pix_dist=1/cos(alpha); else pix_dist=1/sin(alpha); end
    n_pix=round(seg_length/pix_dist)+1;                     % Choose number of pixels along segment to make the discretization of the segment as smooth as possible
    r_pix=zeros(n_pix,2);                                   % Positions of pixels along segment
    for i=0:n_pix
        r_pix(i+1,:)=(1-i/(n_pix-1))*[x1 y1]+(i/(n_pix-1))*[x2 y2];
    end
    r_pix=round(r_pix);
    
    % Calculate intensity along segment over time, store in I2
    n_t=size(a1,3);
    I2=zeros(n_t,n_pix);
    for i1=1:n_t
        for i2=1:n_pix
            I2(i1,i2)=a1(r_pix(i2,2),r_pix(i2,1),i1); 
        end
    end
    
    % Plot time series  of average intensities and check if time range is ok. (Significant drifts in the average intensity hamper a proper the density calculation).
    fig_handle=figure(3); fig_handle.Units='normalized'; fig_handle.OuterPosition=[2/3 0.5 1/3 0.5];
    plot(mean(I2,2)); ax=gca;                    % Plot Intentsity over time
    ax.Title.String='Average intensity on segment versus time';
    ax.XLabel.String='Time (frames)'; ax.YLabel.String='Intensity (a.u.)'; ax.FontSize=12;   
    choice = input('Time limits OK? (If not, press n followed by enter. Otherwise press enter)','s');
    switch choice
        case 'n'
            prompt= {'t_min:', 't_max'};
            dgl_title='Specify timelimits';
            num_lines=1;
            def={num2str(stack_min),num2str(stack_max)};
            answer = inputdlg(prompt,dgl_title,num_lines,def);
            stack_min=str2num(answer{1});
            stack_max=str2num(answer{2});     
    end
    I2=I2(stack_min:stack_max,:);    
    
    % Show movie to find direction of particle movement: check if movement is in direction of arrow 
    fprintf('Determine direction of movement..\n');
    figure(1);  title('Time serie of image sequence'); 
    clf; text(0.1,0.5,'\fontsize{15} \color{black} Determine direction of movement.'); pause(1);
    t_movie=min(stack_max-stack_min+1,100);                 % Number of images in movie
    for i=1:t_movie
        imagesc(a1(:,:,i)); 
        ax=gca; ax.XTick=[]; ax.YTick=[]; 
        axis equal; pause(0.05); hold on;                  
        quiver(x1,y1,x2-x1,y2-y1,'LineWidth',3,'Color','black');% Show arrow 
        hold off;
    end
    ax.FontSize=12;
    title('Was movement in the same direction as the arrow points?')
    choice=input('Movement in direction of arrow? (If not, press n followed by enter. Otherwise press enter)','s'); 
    switch choice
        case 'n'
            I2=I2(:,end:-1:1);                             % If movement was in opposite direction as arrow, adjust the pixel order in I2 to be aligned with the direction of motion 
    end
    
    % Remove background noise from intensities
    histo=histc(reshape(mean(a1,3),[],1),[1:max(max(max(a1)))]);
    [~,Ibg]=max(histo);
    I2=I2-Ibg;
end

%% Sub-function correlation: used to calculate the correlation of intensities as function of temporal separation (dt) and spatial separation (dx) 
% ------------ Input ---------------------------
% I2                2D Array of intensities along segment (colums) versus
%   time (rows).
% dt_min            minimal dt (time separation, in integer image numbers) used to calculate correlation.     
% dt_max            maximal dt used to calculate correlation.     
% dx_min            minimal dx (spatial separation in pixels along segment) used to calculate correlation.     
% dx_max            maximal dx used to calculate correlation.     
% ----------- Output --------------------------
% corr              2D array containing the correlation of intensities
%   separated by dt (rows) and dx (colums).
% ----------- Function ------------------------
function corr=correlation(I2,dt_min,dt_max,dx_min,dx_max)
    Imean=mean(mean(I2));                               % Mean intensity
    Istd=std(reshape(I2,1,[]));                         % Standard deviation of intensities   
    n_pix=size(I2,2);                                   % Number of pixels along segment
    corr=zeros(dt_max-dt_min+1,dx_max-dx_min+1);        % Initialize correlation array       
    
    % Calculate correlation by looping through all combinations of dt and dx
    cc=0;                                               % Array index counter for dt
    for dt=dt_min:dt_max
        cc=cc+1;
        cc2=0;                                          % Array index counter for dx
        for dx=dx_min:dx_max
            cc2=cc2+1;
            cc3=0;                                      % Initialize variable that counts all (dx,dt) combinations, used to normalize
            if dx<0
                for i=-dx+1:n_pix                       % Loop through all possible (dx,dt) combinations in intensity array
                    cc3=cc3+1;
                    corr(cc,cc2)=corr(cc,cc2)+mean((I2(1:end-dt,i)-Imean).*(I2(1+dt:end,i+dx)-Imean));
                end
            else
                for i=1:n_pix-dx
                    cc3=cc3+1;
                    corr(cc,cc2)=corr(cc,cc2)+mean((I2(1:end-dt,i)-Imean).*(I2(1+dt:end,i+dx)-Imean));
                end
            end
            corr(cc,cc2)=corr(cc,cc2)/cc3/Istd^2;       % Normalize correlation by the number of combinations that are summed and the variance of the intensity
        end
    end
end

%% Extract motility parameters from the evolution of the correlation peak
% ----------- Input --------------------------
% corr              2D array containing correlation of intensities
%   separated by dt (rows) and dx (colums).
% dx_min            minimal dx (spatial separation in pixels along segment) used to calculate correlation.     
% dt_min_fit        minimal dt used for extracting velocity and run lengths.
% dt_max_fit        maximal dt used for extracting velocity and run lengths.
% resolution        Resolution of images (micrometers/per pixel).
% dt_frame          Time between images (seconds). 
% pix_dist          Average distance between pixels along segment in pixel
%   units. To correct for the tilt a segment might have to the horizontal
%   or vertical image axis.
% ---------- Output -------------------------
% vel               Velocity (micrometers/second).
% rl                Run length (micrometers).
% --------- Function ------------------------
function [vel,rl]=motility(corr,dx_min,dt_min_fit,dt_max_fit,resolution,dt_frame,pix_dist) 
    corr=corr';
    % Parameters
    dx_fitl=2;                                          % Define which part of the gaussian peak to use for fitting: distance from maximum towards smaller values.
    dx_fith=2;                                          % Define which part of the gaussian peak to use for fitting: distance from maximum towards larger values.
    
    % Initialize arrays to store gaussian fitting results
    dt_range=[dt_min_fit:dt_max_fit];                   % Range of dt over which evolution of correlation peak is analyzed.
    gauss_height=zeros(size(dt_range));                 % Array of peak height from Gaussian fit to correlation peak.    
    gauss_width=zeros(size(dt_range));                  % Array of peak width from Gaussian fit to correlation peak.       
    gauss_dx_mean=zeros(size(dt_range));                % Array of peak dx position from Gaussian fit to correlation peak.    
        
    % Fit Gaussian to correlation function evaluated at each integer dt in dt_range and plot data+fit
    fig_handle=figure(4); fig_handle.Units='normalized'; fig_handle.OuterPosition=[0 0 1/3 0.5];
    hold off; title('Correlation curves with Gaussian fits to peaks'); %
    cc=0;
    for t=dt_range
        cc=cc+1;
        if cc==1                                        % Find dx position of maximum of correlation.
            [~,xpos]=max(corr(-dx_min:-dx_min+10,t));   % Find, for first value of dt used in fitting, dx position of maximum correlation. Range is limited to (-dx:-dx+10) to avoid picking up a noise peak for noisy data. 
            xpos=xpos-dx_min-1;
        else
            [~,xpos_new]=max(corr(xpos-1:xpos+6,t));    % Find, for later values of dt used in fitting, dx position of maximum correlation.
            xpos=xpos_new+xpos-2;
        end
        gauspeakfit=fit([xpos-dx_fitl:xpos+dx_fith]',corr(xpos-dx_fitl:xpos+dx_fith,t),'gauss1'); % Fit Gaussian function to peak of correlation.
        gauss_height(cc)=gauspeakfit.a1;                % Store fit value of peak height
        gauss_dx_mean(cc)=gauspeakfit.b1;               % Store fit value of peak dx position
        gauss_width(cc)=gauspeakfit.c1;                 % Store fit value of peak width
        if mod(t,2)==0, figure(4);plot(corr(:,t),'-'); hold on; plot(gauspeakfit,'--r'); end % Plot the correlation profile of 1 in 2 dt values together with the Gaussian fit.
    end
    figure(4); ax=gca; ax.XLabel.String='dx (pixels)'; ax.YLabel.String='Correlation'; 
    ax.Title.String='Correlation curves with Gaussian fits'; ax.FontSize=12;
    
    % Derive velocity from fit of linear function with zero offset to correlation peak position
    mypropfit=fittype('b*t','independent','t','dependent','y'); % Define linear fitting function.
    fit_vel=fit(dt_range',gauss_dx_mean'+(dx_min-1),mypropfit); % Fit linear function to dx,dt data of peak.
    vel=fit_vel.b*resolution*pix_dist/dt_frame;                 % Calculate velocity in micrometers/second.
    
    % Derive run length from fit of exponential function with zero offset to area below curve.
    fit_rl=fit(dt_range',gauss_width'.*gauss_height','exp1');   % Fit exponential function to decay of peak area versus dt.
    rl=dt_frame./-fit_rl.b*vel;                                 % Calculate run length in micrometers.
    
    % Plot dx vs dt of peak and linear fit for velocity calculation
    fig_handle=figure(5); fig_handle.Units='normalized'; fig_handle.OuterPosition=[1/3 0 1/3 0.5];
    hold off; title('Peak position with linear fit to obtain velocity');
    plot(dt_range,gauss_dx_mean+(dx_min-1),'o'); hold on;
    plot(fit_vel,'r');
    ax=gca; ax.XLabel.String='dt (frames)'; ax.YLabel.String='Peak position (pixels)'; 
    ax.FontSize=12;
    ax.Title.String=strcat('Velocity is:',num2str(vel),' \mum');
    
    % Plot area vs dt of peak and exponential fit for velocity calculation
    fig_handle=figure(6); fig_handle.Units='normalized'; fig_handle.OuterPosition=[2/3 0 1/3 0.5];
    hold off; title('Area of peak with exponential fit to obtain run length');
    semilogy(dt_range,gauss_height.*gauss_width,'o'); hold on;
    plot(fit_rl,'r');
    ax=gca; ax.XLabel.String='dt (frames)'; ax.YLabel.String='Area below peak'; 
    ax.FontSize=12; 
    ax.Title.String=strcat('Run length is:',num2str(rl),' \mum'); 
end