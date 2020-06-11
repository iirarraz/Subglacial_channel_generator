function RR = interp_GRF(shift,lxy,varGRF,nx,ny,SEED)
% RR = interp_GRF(t_seed,lxy,varGRF,nx,ny,SEED)
% function patch to interpolate between t_seed 

% E.g.
% t_seed=200.1;
% lxy=30;
% varGRF=60;
% nx=199; ny=102;

t_dis1=floor(shift);
t_dis2=ceil(shift);
if t_dis1==t_dis2
seed1=SEED(t_dis1:t_dis1+ny-1,:);
RR=MGSimulFFTseed(nx,ny,1,0,varGRF,2,lxy,lxy,1,seed1);
else
seed1=SEED(t_dis1:t_dis1+ny-1,:);
seed2=SEED(t_dis2:t_dis2+ny-1,:);
RR1=MGSimulFFTseed(nx,ny,1,0,varGRF,2,lxy,lxy,1,seed1); 
RR2=MGSimulFFTseed(nx,ny,1,0,varGRF,2,lxy,lxy,1,seed2);
RR=RR2*(shift-t_dis1)+RR1*(t_dis2-shift);
end
end