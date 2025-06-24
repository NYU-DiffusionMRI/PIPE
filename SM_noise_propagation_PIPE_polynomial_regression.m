%% Run SM noise propagation using basis functions from library
clc,clear,close all
root = '/Users/coelhs01/Documents/SantiagoCoelho/Git/PIPE';
% nametag='25_06_20_largesvd_LTE_b10k_GitHub_v0.mat'; % broad range for library
nametag='25_06_20_largesvd_LTE_b15k_GitHub_v0.mat'; % broad range for library
library_path = fullfile(root,nametag);
LIB = load(library_path);

flag_random_b = 1; % Random b-values
% flag_random_b = 0; % Shelled protocol

% Noise propagation SMI test (not for library, only for SMI regression test)
M = 5e3;
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
dirs_target = PIPE.rand_dirs_S2(NB_target);
D_FW = LIB.D_FW;

% Display protocol
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
% S_target = SMI.SM_wFW_b_beta_TE_RealSphHarm_quadInt(f,Da,Depar,Deperp,f_extra,T2a,T2e,plm,b_target,dirs_target,beta_target,0*b_target,1202,LIB.CS_phase);
S_target = SMI.SM_wFW_b_beta_TE_RealSphHarm(f,Da,Depar,Deperp,f_extra,T2a,T2e,plm,b_target,dirs_target,beta_target,0*b_target,LIB.CS_phase,LIB.D_FW);
t=toc; fprintf('Generating forward SM signals in %f s\n',t)

Nl_gamma_fit = [5 4 2 0]; mask = [];
% Add random noise
SNR_dwi = 100;
dwi_synthetic = permute(reshape(S_target+randn(NB_target,M)/SNR_dwi, [NB_target 10 10 M/100]),[2 3 4 1]);

tic
[~, gamma_hat] = PIPE.compute_alpha_gamma(dwi_synthetic,b_target,dirs_target,library_path,Nl_gamma_fit,mask);
t=toc; fprintf('Computing alpha(b,g) and gamma(kernel,plm) in %f s\n',t)

% Estimate SM kernel parameters
Nl_kernel = Nl_gamma_fit;
sigma_dwi_norm = 1/SNR_dwi; PRdegree = 3;
% Use isocenter protocol for gamma noise propagation
tic
kernel_fit = PIPE.SM_LTE_gamma2kernel(gamma_hat,b_target,dirs_target,library_path,Nl_gamma_fit,Nl_kernel,mask,sigma_dwi_norm,PRdegree);
t=toc; fprintf('Computing Standard Model parameters from gamma_nlm in %f s\n',t)

kernel_hat = reshape(kernel_fit, [M 5]);
kernel_gt = kernel(:,1:5);
param_lims=[0.2 0.8;1 2.5;1 2.5;0 1;0 0.2];
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$'};
figure('Position',[72 742 2137 401])
for ii=1:5
    kernel_target_fit = kernel_hat(:,ii);
    id_mins=kernel_target_fit<param_lims(ii,1);
    id_maxs=kernel_target_fit>param_lims(ii,2);
    kernel_target_fit(id_mins)=param_lims(ii,1);
    kernel_target_fit(id_maxs)=param_lims(ii,2);

    RMSE=sqrt(mean((kernel_gt(:,ii)-kernel_target_fit).^2));
    RMSEtag=['RMSE=',num2str(RMSE,3)];
    subplot(1,5,ii), 
    h = scatter_kde(kernel_gt(:,ii),kernel_target_fit,'filled'); hold on, plot(param_lims(ii,:),param_lims(ii,:),'r-.'), axis([param_lims(ii,:) param_lims(ii,:)])
    title([paramNames{ii},' - ',RMSEtag],'interpreter','latex')
    xlabel('Ground truth','interpreter','latex'), ylabel('NN output','interpreter','latex'), set(gca,'FontSize',20)
    fprintf('================== Done with param %d/5\n',ii)
end
