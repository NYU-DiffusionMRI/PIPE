%%
clc,clear,close all
root = '/Users/coelhs01/Documents/SantiagoCoelho/Git/PIPE';
nametag='25_06_20_largesvd_LTE_b10k_GitHub_v0.mat'; % broad range for library
LIB = load(fullfile(root,nametag));


flag_random_b = 1; % Random b-values
% flag_random_b = 0; % Shelled protocol


% Noise propagation SMI test (not for library, only for SMI regression test)
M = 4e4;
lb=[0.2 1 1 0.15 0 50 30 0.1];
ub=[0.8 2.5 2.5 0.8 0.2 150 100 0.9];
Lmax=6;
[f,Da,Depar,Deperp,f_w,T2a,T2e,~,~] = SMI.Get_uniformly_distributed_SM_prior(M,lb,ub,Lmax);
f_extra=1-f-f_w;
kernel=[f, Da, Depar, Deperp, 1-f-f_extra, T2a,T2e];
% Generate fODF
Lmax=6;
Nfibers=2;%round(rand(M,1))+1;
[plm,pl] = rand_REALplm_powerLaw_DNnorm_3Ea_Nfibers(M,Lmax,[],[],Nfibers);
plm=plm';

% Setting up target protocol
NB_target = 200; bmax = 8;

if flag_random_b
    % Random b-values
    b_target=[0 rand(1,NB_target-1)*bmax];
else
    % Shelled protocol
    b_target=[0 ones(1,30) 2*ones(1,30) 3*ones(1,34) 4*ones(1,45) 8*ones(1,60)];
end

beta_target=ones(1,NB_target);
dirs_target=rand_sph([],NB_target);
D_FW = LIB.D_FW;


[~,id_b_sorted]=sort(b_target);
tprime=linspace(0,1,length(b_target));
t=tprime.^(0.75);
cmap=[1 0 0]'*t+[0 0 1]'*(1-t);
figure('Position',[480 541 1111 472]), subplot(121), hold on, subplot(122), hold on
start=1;
for ii=1:length(b_target)
    subplot(121), plot(ii,b_target(id_b_sorted(ii)),'.','LineWidth',2,'MarkerSize',10,'Color',cmap(:,ii)')
    subplot(122), plot3(dirs_target(id_b_sorted(ii),1),dirs_target(id_b_sorted(ii),2),dirs_target(id_b_sorted(ii),3),'.','LineWidth',2,'MarkerSize',25,'Color',cmap(:,ii)')
end
subplot(121), set(gca,'FontSize',25), ylabel('b-value [$\mathrm{ms}/\mu\mathrm{m}^2$]','interpreter','latex'), xlabel('Measurements sorted','interpreter','latex'), set(gca,'FontSize',25),
subplot(122), set(gca,'FontSize',25), axis equal, axis off, campos([4.1104  -16.7049    1.7235])


tic
S_target = SMI.SM_wFW_b_beta_TE_RealSphHarm_quadInt(f,Da,Depar,Deperp,f_extra,T2a,T2e,plm,b_target,dirs_target,beta_target,0*b_target,1202,LIB.CS_phase);
t=toc; fprintf('Generating forward SM signals in %f s\n',t)

% Interpolating U for new b,beta values
tic
U0_target = Chebyshev_interpolation_U_b(LIB.U0(:,1:LIB.n0),b_target,LIB.bmax);
U2_target = Chebyshev_interpolation_U_b(LIB.U2(:,1:LIB.n2),b_target,LIB.bmax);
U4_target = Chebyshev_interpolation_U_b(LIB.U4(:,1:LIB.n4),b_target,LIB.bmax);
U6_target = Chebyshev_interpolation_U_b(LIB.U6(:,1:LIB.n6),b_target,LIB.bmax);
t=toc; fprintf('Interpolating u(b) for TARGET b in %f s\n',t)

tic
Lmax=6;
Ylm=SMI.get_even_SH(dirs_target,Lmax,LIB.CS_phase);
Y00_target=Ylm(:,1);
Y2m_target=Ylm(:,2:6);
Y4m_target=Ylm(:,7:15);
Y6m_target=Ylm(:,16:28);
ALPHA_target=zeros(NB_target,LIB.n0+LIB.n2*5+LIB.n4*9+LIB.n6*13);
for ii=1:NB_target
    UY0=kron(U0_target(ii,:),Y00_target(ii,:));
    UY2=kron(U2_target(ii,:),Y2m_target(ii,:));
    UY4=kron(U4_target(ii,:),Y4m_target(ii,:));
    UY6=kron(U6_target(ii,:),Y6m_target(ii,:));
    ALPHA_target(ii,:)=[UY0 UY2 UY4 UY6];
end
t=toc; fprintf('Computing ALPHA(b,g) in %f s\n',t)


%% Run polynomial regression on linearly estimated gamma_nlm
clc,close all
Mtest = 5e3; Mtrain = M-Mtest; SNR_dwi = 50;

id_train = 1:(M-Mtest);
id_test = (Mtrain+1):M;

S_target_noisy=S_target+randn(NB_target,M)/SNR_dwi;
GAMMA_hat=ALPHA_target\(S_target_noisy(:,id_train)/(4*pi));
GAMMA_gt=ALPHA_target\(S_target(:,id_train)/(4*pi));

GAMMA_hat_test=ALPHA_target\(S_target_noisy(:,id_test)/(4*pi));
GAMMA_gt_test = ALPHA_target\(S_target(:,id_test)/(4*pi));

New_Nell = [LIB.n0 LIB.n2 LIB.n4 LIB.n6];

% regression_n = [5 4 0 0]; % Choosing how many basis functions we use for the regression
regression_n = [4 3 0 0]; % Choosing how many basis functions we use for the regression
regression_ids = [1:regression_n(1) LIB.n0+(1:regression_n(2)*5)  LIB.n0+5*LIB.n2+(1:regression_n(3)*9) LIB.n0+5*LIB.n2+9*LIB.n4+(1:regression_n(4)*13) ];

options.training_data = GAMMA_hat(regression_ids,:);   
options.test_data = GAMMA_hat_test(regression_ids,:);
tag_train = ['$N_\ell$=[',num2str(regression_n),']'];

% PR fitting
options.training_objective = kernel(1:Mtrain,1:5)';
% options.Degree = 1; reg_tag = 'Linear PR';
% options.Degree = 2; reg_tag = 'Quadratic PR';
options.Degree = 3; reg_tag = 'Cubic PR';
options.strategy = 'PR';
tic
out = PIPE.DataDriven_regression(options);  
t=toc; fprintf('PR fit done in %f s\n',t)
kernel_hat = out.fits;

kernel_gt = kernel(id_test,:);
param_lims=[0.2 0.8;1 2.5;1 2.5;0 1;0 0.2];
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$'};
figure('Position',[72 742 2137 401])
for ii=1:5
    kernel_target_fit(:,ii) = kernel_hat(:,ii);

    id_mins=kernel_target_fit(:,ii)<param_lims(ii,1);
    id_maxs=kernel_target_fit(:,ii)>param_lims(ii,2);
    kernel_target_fit(id_mins,ii)=param_lims(ii,1);
    kernel_target_fit(id_maxs,ii)=param_lims(ii,2);

    RMSE=sqrt(mean((kernel_gt(:,ii)-kernel_target_fit(:,ii)).^2));
    RMSEtag=['RMSE=',num2str(RMSE,3)];
    subplot(1,5,ii), 
    h = scatter_kde(kernel_gt(:,ii),kernel_target_fit(:,ii),'filled'); hold on, plot(param_lims(ii,:),param_lims(ii,:),'r-.'), axis([param_lims(ii,:) param_lims(ii,:)])
    title([paramNames{ii},' - ',RMSEtag],'interpreter','latex')
    xlabel('Ground truth','interpreter','latex'), ylabel('NN output','interpreter','latex'), set(gca,'FontSize',20)
    fprintf('================== Done with param %d/5\n',ii)
end
sgt = sgtitle([reg_tag,' $\gamma_{n\ell m}$ (SNR = ',num2str(SNR_dwi),') - ',tag_train]); sgt.FontSize = 30; sgt.Interpreter='latex';


%% Run polynomial regression on denoised ratios of <gamma_nlm/gamma_1lm> (rotational invariants of gamma)
clc,close all

Mtest = 5e3; Mtrain = M-Mtest; SNR_dwi = 100;

id_train = 1:(M-Mtest);
id_test = (Mtrain+1):M;

S_target_noisy=S_target+randn(NB_target,M)/SNR_dwi;
GAMMA_hat=ALPHA_target\(S_target_noisy(:,id_train)/(4*pi));
GAMMA_gt=ALPHA_target\(S_target(:,id_train)/(4*pi));

GAMMA_hat_test=ALPHA_target\(S_target_noisy(:,id_test)/(4*pi));
GAMMA_gt_test = ALPHA_target\(S_target(:,id_test)/(4*pi));

New_Nell = [LIB.n0 LIB.n2 LIB.n4 LIB.n6];
N2_svd = 4; N4_svd = 3; N6_svd = 2;

id_2m = (New_Nell(1)+1):(New_Nell(1)+N2_svd*5);
id_4m = (New_Nell(1)+5*New_Nell(2))+(1:N4_svd*9);
id_6m = (New_Nell(1)+5*New_Nell(2)+9*New_Nell(3))+(1:N6_svd*13);
GAMMA_hat_dn = PIPE.rank1denoisegamma(GAMMA_hat,id_2m,2);
GAMMA_hat_dn = PIPE.rank1denoisegamma(GAMMA_hat_dn,id_4m,4);
GAMMA_hat_dn = PIPE.rank1denoisegamma(GAMMA_hat_dn,id_6m,6);
GAMMA_hat_test_dn = PIPE.rank1denoisegamma(GAMMA_hat_test,id_2m,2);
GAMMA_hat_test_dn = PIPE.rank1denoisegamma(GAMMA_hat_test_dn,id_4m,4);
GAMMA_hat_test_dn = PIPE.rank1denoisegamma(GAMMA_hat_test_dn,id_6m,6);

% Taking ratios
[GAMMA_hat_RI,~] =         PIPE.ComputeRatioRotationallyInvariantGamma(GAMMA_hat,New_Nell);
[GAMMA_hat_dn_RI,~] =      PIPE.ComputeRatioRotationallyInvariantGamma(GAMMA_hat_dn,New_Nell);
[GAMMA_hat_test_RI,~] =    PIPE.ComputeRatioRotationallyInvariantGamma(GAMMA_hat_test,New_Nell);
[GAMMA_hat_test_dn_RI,~] = PIPE.ComputeRatioRotationallyInvariantGamma(GAMMA_hat_test_dn,New_Nell);
[GAMMA_gt_RI,~] =          PIPE.ComputeRatioRotationallyInvariantGamma(GAMMA_gt,New_Nell);
[GAMMA_gt_test_RI,~] =     PIPE.ComputeRatioRotationallyInvariantGamma(GAMMA_gt_test,New_Nell);

% regression_n = [5 4 0 0]; % Choosing how many basis functions we use for the regression
regression_n = [4 3 0 0]; % Choosing how many basis functions we use for the regression
regression_ids = [1:regression_n(1) LIB.n0+(1:regression_n(2))  LIB.n0+LIB.n2+(1:regression_n(3)) LIB.n0+LIB.n2+LIB.n4+(1:regression_n(4)) ];


options.training_data = GAMMA_hat_dn_RI(regression_ids,:);   
options.test_data = GAMMA_hat_test_dn_RI(regression_ids,:);
tag_train = ['$N_\ell$=[',num2str(regression_n),']'];

% PR fitting
options.training_objective = kernel(1:Mtrain,1:5)';
% options.Degree = 1; reg_tag = 'Linear PR';
% options.Degree = 2; reg_tag = 'Quadratic PR';
options.Degree = 3; reg_tag = 'Cubic PR';
options.strategy = 'PR';
tic
out = PIPE.DataDriven_regression(options);  
t=toc; fprintf('PR fit done in %f s\n',t)
kernel_hat = out.fits;


kernel_gt = kernel(id_test,:);
param_lims=[0.2 0.8;1 2.5;1 2.5;0 1;0 0.2];
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$p_2$'};
% Gamma gt
figure('Position',[72 742 2137 401])
for ii=1:5
    kernel_target_fit(:,ii) = kernel_hat(:,ii);
    
    id_mins=kernel_target_fit(:,ii)<param_lims(ii,1);
    id_maxs=kernel_target_fit(:,ii)>param_lims(ii,2);
    kernel_target_fit(id_mins,ii)=param_lims(ii,1);
    kernel_target_fit(id_maxs,ii)=param_lims(ii,2);
    
    RMSE=sqrt(mean((kernel_gt(:,ii)-kernel_target_fit(:,ii)).^2));
    RMSEtag=['RMSE=',num2str(RMSE,3)];
    subplot(1,5,ii), 
    %     hold on, histogram(kernel_target(:,ii)), histogram(kernel_target_fit(:,ii))
    h = scatter_kde(kernel_gt(:,ii),kernel_target_fit(:,ii),'filled'); hold on, plot(param_lims(ii,:),param_lims(ii,:),'r-.'), axis([param_lims(ii,:) param_lims(ii,:)])
    title([paramNames{ii},' - ',RMSEtag],'interpreter','latex')
    xlabel('Ground truth','interpreter','latex'), ylabel('NN output','interpreter','latex'), set(gca,'FontSize',20)
    fprintf('================== Done with param %d/5\n',ii)
end
sgt = sgtitle([reg_tag,' $\langle\gamma_{n\ell m}/\gamma_{1\ell m}\rangle_m$ (SNR = ',num2str(SNR_dwi),') - ',tag_train]); sgt.FontSize = 30; sgt.Interpreter='latex';



