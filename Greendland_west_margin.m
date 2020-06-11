%% Supplementary material from
% Determining the evolution of an alpine glacier drainage system by solving inverse problems
% Inigo Irarrazaval 2020.

% I provide an overlook of introducing a gaussian random feld on Shreve's
% equation to add variability and produce a large variety of subglacial
% channels.

% To run the example:

% First. Download/clone Topotoolbox2
% Schwanghart, W., & Scherler, D. (2014). Short Communication:
% TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences.
% Earth Surf. Dynam., 2(1), 1-7. doi:10.5194/esurf-680 2-1-2014
% ==> https://topotoolbox.wordpress.com/download/

% Second. Download glacier bed and thickness from:
% Lindbäck, K., Pettersson, R., Doyle, S. H., Helanow, C., Jansson, P., Kristensen, S. S., Stenseng, L., Forsberg, R., and Hubbard, A. L.:
% High-resolution ice thickness and bed topography of a land-terminating section of the Greenland Ice Sheet,
% Earth Syst. Sci. Data, 6, 331–338, https://doi.org/10.5194/essd-6-331-2014, 2014.
% ==> Geotiff download links here
% http://store.pangaea.de/Publications/Lindbaeck_et_al_2014/GAP_DEM_bed_V1_1_geotiff.zip
% http://store.pangaea.de/Publications/Lindbaeck_et_al_2014/GAP_DEM_thickness_V1_geotiff.zip
%%
close all
clear
clc
set(groot,'defaultFigureColor','w')
%% load Russel Glacier
B=GRIDobj('C:\Users\iirarraz\Documents\01_UNIL\50_SSDS_GORNER_2\Russel_glacier\GAP_DEM_bed_V1_1.tif');
H=GRIDobj('C:\Users\iirarraz\Documents\01_UNIL\50_SSDS_GORNER_2\Russel_glacier\GAP_DEM_thickness_V1.tif');
H.refmat=B.refmat; % Note! there is a 0.3725e-8 misalignement. So we just regularize it.
[ny,nx]=size(B.Z);

%% Define some variables
rng(2)        % fix random seed for repetibility.
rho_i=900;    % define constants
rho_w=1000;
g_grav=9.8;
f=1;          % flotation factor in Shreves equation (case 1 and 2)
a=3;          % Parameters a and b define channel radii r(ui)=a*exp(ui*b)
b=3;          % In this exmple are used for plotting visualization only.
c=200;        % c densification factor (flow accumulation area threshold for channels)

%% CASE 1. Inspect channels without adding a gaussian random field
% multidirectional flow accumulation algorithm

phi=rho_i*g_grav*H+f*rho_w*g_grav*B;
phi_2=fillsinks(phi);          % pre-process
FD = FLOWobj(phi_2,'Multi');   % multi direction flow algorithm
Fpp=flowacc(FD);

fig=figure;
fig.Position=[241  393   1079 481];
s1=subplot(121);
imagesc(log10(Fpp))
axis equal
caxis([0 5])
title('Shreves equation (Multi)')

s2=subplot(122);

try
    St=STREAMobj(FD,Fpp>=c);
    plot(St)
    linkaxes([s1,s2])
catch
    disp('Warning! some versions of TopoToolbox does not support this feature anymore')
    disp('Please refer to Case1 figure')
%    open('Case1.fig')
end

% Zoom to the channels to see closely.
% One of the challenges is that multdirection algorithm produces many
% channels that divide downstream. This is not likely in a
% subglacial drainage system as channels tend to run in low pressure areas
% capturing water form the surrouding. Several channels near to
% each others is unestable.

%% CASE 2. Inspect channels without adding a gaussian random field
% single direction flow accumulation algorithm

phi=rho_i*g_grav*H+f*rho_w*g_grav*B;
phi_2=fillsinks(phi);
FD = FLOWobj(phi_2);
Fpp=flowacc(FD);

fig=figure;
fig.Position=[241  393   1079 481];
s1=subplot(121);
imagesc(log10(Fpp))
axis equal
caxis([0 5])
title('Shreves equation (D8)')

s2=subplot(122);
linkaxes([s1,s2])

St=STREAMobj(FD,Fpp>=c); % stream object save for later
ustr=streamorder_shreve_acc(St,Fpp);
ui=ustr/max(ustr);
r=a*exp(ui*b);
sc=scatter(St.x,St.y,r,'filled');
sc.MarkerFaceAlpha = 0.1;  %
axis equal

% One of the problems with single flow direction is that no uncertainity on
% the data is considered. Therefore only one network is generated. Previous
% studies have modified the flotation fator to increase variability, but is
% homogeneously distributed.
%% CASE 3. Inspect channels when adding gaussian random field
% radomly sampling values of flotation (f), shift (s), integral scales
% (lxy) and gaussian random field variance (varGRF).

WN=rand(2000,size(H.Z,2)); % generate a white noise
n=50; % number of channel networks to generate
ns=4; % number of channel networks to save and plot

F=zeros(n,1); S=zeros(n,1);LXY=zeros(n,1); VAR=zeros(n,1); % pre-allocate
FACC=zeros(size(H.Z));
STT=cell(n,1); FPP=cell(n,1); cc=round(linspace(1,n,ns));

% for jj=1:n
parfor jj=1:n

    % generate random channel network model parameters
    varGRF=10+(10-rand*20);    % variance of gaussian random field (fixed in Irarrazaval et al., 2020)
    lxy=15+(10-rand*20);        % Integral scales (in pixels)
    f=0.9+rand*.2;              % flotation factor
    shift=rand*(size(WN,1)-ny-1)+1; % shift
    
    phi_R = interp_GRF(shift,lxy,varGRF,nx,ny,WN); % compute GRF
    phi=rho_i*g_grav*H+f*rho_w*g_grav*B;
    phi_2=phi+phi_R*1000*9.8; % equivalent to phi* perturbed shreves equation
    
    phi_2=fillsinks(phi_2);
    FD = FLOWobj(phi_2);
    Fpp=flowacc(FD);
    
    FACC=FACC+Fpp.Z; % save realizations for later
    F(jj)=f;
    S(jj)=shift;
    LXY(jj)=lxy;
    VAR(jj)=varGRF;
    
    if ismember(jj,cc) % save only if explicited
        STT{jj}=STREAMobj(FD,Fpp>=c); % stream object save for later
        FPP{jj}=Fpp; % flow object save for later
    end
    disp(num2str(jj))
end

fig=figure;
fig.Position=[241  393   1079 481];
s1=subplot(121);
FACC2=Fpp; FACC2.Z=FACC;
imagesc(log10(FACC2/n))
axis equal
title('Shreves equation* (Eq.3)')
caxis([0 5])

s2=subplot(122);
linkaxes([s1,s2])
for jj=cc
    St=STT{jj}; Fpp=FPP{jj};
    ustr=streamorder_shreve_acc(St,Fpp);
    ui=ustr/max(ustr);
    r=a*exp(ui*b);
    sc=scatter(St.x,St.y,r,'filled');
    sc.MarkerFaceAlpha = 0.1;  %
    hold on
    axis equal
end
% many different channels network where generated.

fig=figure;
fig.Position=[587   381   810   390];
subplot(141)
histogram(F,10,'Normalization','probability')
xlabel('flotation f')
subplot(142)
histogram(S,10,'Normalization','probability')
xlabel('shift s')
subplot(143)
histogram(LXY,10,'Normalization','probability')
xlabel('Integral scale l_x_y')
subplot(144)
histogram(VAR,10,'Normalization','probability')
xlabel('GRF variance')
suptitle('Parameter distributions')
%% Until now, you have only produced some "stochastic" subglacial channel networks.
% Following steps are:
% Step 2: Define water recharge, boundary conditions and use water flow
% equations to test each channel network. The output is a water pressure field and
% transit speeds in subglacial channels.

% Step 3: The output pressure and transit speeds can be
% compared with field observations through a Likelihood function (misfit function).
% the best fitting models will have a high likelihood.
% As many channel networks would perform poorly, it is more efficitent to
% find an algorithm to navigate the model space. In Irarrazaval et
% al., 2019 a Markov-chain Monte Carlo is prefered.