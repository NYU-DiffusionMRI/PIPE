%% Generating big library of Standard Model basis functions
clc,clear,close all

root = '/Users/coelhs01/Documents/SantiagoCoelho/Git/PIPE';

% % Large experiment
M=5e4; N_b = 1000; N_beta=1;

% Compute synthetic dwis according to SM
lb=[0.01, 1, 1, 0.1,   0,  50,  50, 0.05]; % Since data is fixed TE, compartmental T2s are irrelevant
ub=[0.99, 3, 3, 1.0, 0.2, 150, 120, 0.95];

Lmax=6;
[f,Da,Depar,Deperp,f_w,T2a,T2e,~,~] = SMI.Get_uniformly_distributed_SM_prior(M,lb,ub,Lmax);
f_extra=1-f-f_w;
kernel=[f, Da, Depar, Deperp, 1-f-f_extra, T2a,T2e];

bmax=10; nametag='25_06_20_largesvd_LTE_b10k_GitHub_v0.mat'; % broad range for library

b = PIPE.get_Chebyshev_nodes(N_b,[0 bmax]);
NB=N_b;

b=b(:)';beta=0*b+1;

D_FW=3;
tic
K0 = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(0,b,beta,b*0,kernel,D_FW)';
K2 = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(2,b,beta,b*0,kernel,D_FW)';
K4 = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(4,b,beta,b*0,kernel,D_FW)';
K6 = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(6,b,beta,b*0,kernel,D_FW)';
t=toc; fprintf('K_ell library generated in %f s\n',t)

Lmax=6;
Nfibers=2;%round(rand(M,1))+1;
CS_phase=1;
tic
[plm,pl] = rand_REALplm_powerLaw_DNnorm_3Ea_Nfibers(M,Lmax,[],[],Nfibers);
plm=plm';
p2=pl(:,1);
if Lmax>2, p4=pl(:,2); end
if Lmax>4, p6=pl(:,3); end
t=toc; fprintf('random fODFs generated in %f s\n',t)

SNR=50;
NdwiTypical=150;

sigma_0=1/(SNR*sqrt(NdwiTypical));
sigma_2=1/(SNR*sqrt(5*NdwiTypical));
sigma_4=1/(SNR*sqrt(9*NdwiTypical));
sigma_6=1/(SNR*sqrt(13*NdwiTypical));

tic
[U0,S0,V0]=svd(K0,'econ');
[~, ~, n0_MP] = MP(K0+randn(NB,M)*sigma_0, 25);
clearvars K0
[U2,S2,V2]=svd(K2,'econ');
[~, ~, n2_MP] = MP(K2+randn(NB,M)*sigma_2, 25);
clearvars K2
[U4,S4,V4]=svd(K4,'econ');
[~, ~, n4_MP] = MP(K4+randn(NB,M)*sigma_4, 25);
clearvars K4
[U6,S6,V6]=svd(K6,'econ');
[~, ~, n6_MP] = MP(K6+randn(NB,M)*sigma_6, 25);
clearvars K6
t=toc; fprintf('ALL K_ell svds computed in %f s\n',t)

% Get target ul basis functions
n0=n0_MP;n2=n2_MP;n4=n4_MP;n6=n6_MP;
disp([n0 n2 n4 n6])

% Building full system for inversion
U0=U0(:,1:n0);
U2=U2(:,1:n2);
U4=U4(:,1:n4);
U6=U6(:,1:n6);
S0V0t=S0(1:n0,1:n0)*V0(:,1:n0)';
S2V2t=S2(1:n2,1:n2)*V2(:,1:n2)';
S4V4t=S4(1:n4,1:n4)*V4(:,1:n4)';
S6V6t=S6(1:n6,1:n6)*V6(:,1:n6)';
clearvars V0 V2 V4 V6

save(fullfile(root,nametag))



