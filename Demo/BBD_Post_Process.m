%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is a demonstration using for the phase reconstruction using BBD method in paper of:
% Binbin Wang and David W. McComb, Ultramicroscopy 2023, "Phase imaging in scanning transmission electron microscopy using bright-field balanced divergency method"
%
%The TIE calculation is based on:
% Main file for Experiment 3 (multiple cells) for the TIE tutorial.
% Related Reference:
% "Transport-of-intensity equation��A tutoral" Chao Zuo, et. al.
% 
% Last modification: Binbin Wang, 12/18/2022
% wang.10255@osu.edu or binbin1.wang@intel.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear everything existing.
clc
close all
clearvars
addpath('./Functions');
fpath = fullfile('.','MoS2_monolayer_Def4nm');
addpath(fpath)

%% Setting up system parameters

I0 = double(imread('1.tif'))/255;
Iz = double(imread('2.tif'))/255;
%% Dose balance
I0 = (mean2(Iz)/mean2(I0))*I0;

% Valid domain within the apecture
Apecture = ones(size(I0));
Apecture(I0<max(max(I0))/10)=0;


eV=input('What is the electron energy in keV?                 : ')*1000;
c1=9.7846113e-7;
c2=12.263868;
ke=sqrt(eV+c1*eV^2)/c2;

Pixelsize = 11.9e-9;    % pixel size (m)
lambda = 1/ke;        % wavelength (m)
k = 2*pi/lambda;        % wave number
dz=input('What is the defocus (nm)?                 : ')*1e-9;
% dz = -4e-9;            % defocus distance (m)

% Axial intensity derivative
dIdz = (Iz-I0)/(dz);

% Check if satisfying Energy Conservation
disp(['mean2(Iz-I0)) = ' num2str(mean2(Iz-I0))]);
if abs(mean2(Iz-I0))>max(max(I0))/255    % Need to check the threshold
    warndlg('Conservation of Energy may not be satisfied!');
end

%% Show the images.
figure
iptsetpref('ImshowBorder','tight');
subplot(3,1,1)
% imshow(I0,[0,1]);
imshow(I0,[]);
title('Defocused intensity');
subplot(3,1,2)
imshow(Iz,[]);
title('Overfocus intensity');
subplot(3,1,3)
imshow(dIdz,[]);
title('Intensity Derivative');


%% Solve TIE with iteration parameters.
%% Solve TIE with iteration parameters.
RegPara=input('Low pass filter (1 is none), e.g. 1e-2                 : ');
%RegPara = 1e-2;
IntThr = 1e-2;
JudgeFlag_DCT = 1;
JudgeFlag_US = 1;
MaxIterNum = 10;
Method = 'Fresnel';  %'TIE','Angular Spectrum','Fresnel'

%% Solve TIE with FFT-TIE.
Phi_FFT = TIE_FFT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr);

%% Solve TIE with DCT-TIE.
Phi_DCT = TIE_DCT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr);

%% Solve TIE with US-TIE.
Phi_US = Universal_Solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_US,Method);

%% Show Reconstruction Result
figure;
iptsetpref('ImshowBorder','tight');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',scrsz);

subplot(2,2,1)
imshow(Phi_FFT.*Apecture,[]);
colormap(gca,gray);colorbar
title('Reconstructed Phase (FFT-TIE)');
R = georasterref('RasterSize',size(Phi_FFT.*Apecture));
geotiffwrite(fullfile(fpath,'Reconstructed Phase (FFT-TIE).tif'),Phi_FFT.*Apecture,R)
subplot(2,2,2)
imshow(Phi_DCT.*Apecture,[]);
colormap(gca,gray);colorbar
R = georasterref('RasterSize',size(Phi_FFT.*Apecture));
geotiffwrite(fullfile(fpath,'Reconstructed Phase (DCT-TIE).tif'),Phi_DCT.*Apecture,R)
title('Reconstructed Phase (DCT-TIE)');
subplot(2,2,3)
imshow(Phi_US.*Apecture,[]);
colormap(gca,gray);colorbar
R = georasterref('RasterSize',size(Phi_FFT.*Apecture));
geotiffwrite(fullfile(fpath,'Reconstructed Phase (US-TIE).tif'),Phi_US.*Apecture,R)
title('Reconstructed Phase (US-TIE)');

FieldShow(Phi_FFT.*Apecture,fpath,'FFT_TIE')
FieldShow(Phi_US.*Apecture,fpath,'US_TIE')

% FieldShow(Phi_US.*Apecture)

%% make diff to potential maps make 1st derivetive to electric field and 2ed deri to charge 
%% Binbin Wang
function FieldShow(Im,fpath,filename)
    mkdir(fullfile(fpath,filename))
    [numRows,numCols] = size(Im);
    u=numRows-1;
    v=numCols-1;
    
    dify=diff(Im,1,1);
    difx=diff(Im,1,2);
    elecfl_inty=sqrt(dify(1:end,1:end-1).^2+difx(1:end-1,1:end).^2);
    theta=atan2(dify(1:end,1:end-1),difx(1:end-1,1:end));
    theta=theta.*(theta>=0)+(theta+2*pi).*(theta<0);
    color_index=int16((theta./(2*pi)).*256);
    threshold=max(max(elecfl_inty));
    mask=elecfl_inty./threshold;
    c=hsv(256);

    electric_vect=zeros(u,v,3);
    for i=1:u
        for j=1:v
            if color_index(i,j)==0
                color_index(i,j)=1;
            end
        electric_vect(i,j,:)=c(color_index(i,j),:);  % give electric field mapping
        end
    end
    figure;
    subplot(2,2,1)
    imagesc(Im);axis image; colormap(gray);
    title('Reconstructed Phase');
    subplot(2,2,2)
    imagesc(elecfl_inty);axis image; colormap(gray);
    title('Electric field amplitude');
    R = georasterref('RasterSize',size(elecfl_inty));
    geotiffwrite(fullfile(fpath,filename,'Electric field amplitude.tif'),elecfl_inty,R)
    subplot(2,2,3)
    electric_vect=electric_vect.*mask;
    colormap(hsv);imagesc(electric_vect); axis image; 
    title('Electric field mapping');
    R = georasterref('RasterSize',size(electric_vect));
    geotiffwrite(fullfile(fpath,filename,'Electric field vector.tif'),electric_vect,R)
    subplot(2,2,4)
    difyy=diff(Im,2,1);
    difxx=diff(Im,2,2);
    charg_inty=sqrt(difyy(1:end,1:end-2).^2+difxx(1:end-2,1:end).^2);
    imagesc(charg_inty);axis image; colormap(gray);
    R = georasterref('RasterSize',size(charg_inty));
    geotiffwrite(fullfile(fpath,filename,'Charge density.tif'),charg_inty,R)
    title('Charge density');
%     figure(5);imagesc(Im*Im);axis image; colormap(gray);
end 



