%% PIPE 
clc,clear,close all
root_data = '/Volumes/labspace/Santiago/PIPE_GE/MAGNUS_data';
parpool

tic
nii_dwi = load_untouch_nii(fullfile(root_data,'dwi.nii'));
nii_b = load_untouch_nii(fullfile(root_data,'bval_actual.nii'));
nii_g = load_untouch_nii(fullfile(root_data,'bvec_actual.nii'));
nominal_protocol = load(fullfile(root_data,'nominal_protocol.mat'));
nii_mask = load_untouch_nii(fullfile(root_data,'brain_mask.nii'));
t=toc; fprintf('================== Done loading data in %f seconds\n',t)
mask = logical(nii_mask.img);

root_git = '/Users/coelhs01/Documents/SantiagoCoelho/Git/PIPE';
nametag='25_06_20_largesvd_LTE_b15k_GitHub_v0.mat'; % broad range for library
library_path = fullfile(root_git,nametag);
Nl_gamma_fit=[4 3 0 0];
tic
[alpha_nlm_flat, gamma_hat_flat] = PIPE.compute_alpha_gamma(nii_dwi.img,nominal_protocol.b,nominal_protocol.g,library_path,Nl_gamma_fit,mask);
t=toc; fprintf('Computing alpha(b,g) and gamma(kernel,plm) in %f s (flat L field - isocenter)\n',t)

tic
[alpha_nlm, gamma_hat] = PIPE.compute_alpha_gamma(nii_dwi.img,nii_b.img,nii_g.img,library_path,Nl_gamma_fit,mask);
t=toc; fprintf('Computing alpha(b,g) and gamma(kernel,plm) in %f s (spatially varying L)\n',t)

% Computing signal rotational invariants through projections on many dirs
new_b = [1 2 4 6];
tic
[S0_b,S2_b] = PIPE.project_gamma_onto_RotInvs(gamma_hat,mask,new_b,library_path,Nl_gamma_fit);
t=toc; fprintf('Computing signal rotational invariants in %f s\n',t)

% Estimate SM kernel parameters
Nl_kernel = Nl_gamma_fit;
sigma_dwi_norm = 1/50; PRdegree = 3;
% Use isocenter protocol for gamma noise propagation
tic
kernel_fit = PIPE.SM_LTE_gamma2kernel(gamma_hat,nominal_protocol.b,nominal_protocol.g,library_path,Nl_gamma_fit,Nl_kernel,mask,sigma_dwi_norm,PRdegree);
t=toc; fprintf('Computing Standard Model parameters from gamma_nlm in %f s\n',t)




%% Display results
clc,close all
% plot gamma00
allMAPS = gamma_hat(:,:,:,1:4);
slice=34;
clims = [-3.5   -2.5 ; -0.75    0.5 ;  -0.25    0.25 ; -0.1    0.1];
allMAPS=flip(permute(allMAPS,[2 1 3 4]),1);
figure('Position',[2606 799 1772 374]), colormap gray
nametags={'$\hat{\gamma}_{n=1,\ell=0,m=0}$','$\hat{\gamma}_{n=2,\ell=0,m=0}$','$\hat{\gamma}_{n=3,\ell=0,m=0}$','$\hat{\gamma}_{n=4,\ell=0,m=0}$'};
PIPE.plotSlices(allMAPS, slice,clims,nametags,1,[],1,0),

% plot gamma2m
allMAPS = gamma_hat(:,:,:,[ 5:9 10:14 15:19]);
slice=34; clims = repelem([ -0.4    0.4 ; -0.06    0.06 ; -0.02    0.02 ],5,1);
allMAPS=flip(permute(allMAPS,[2 1 3 4]),1);
figure('Position',[93 32 2151 1305]), colormap gray
nametags = {'$\hat{\gamma}_{n=1,\ell=2,m=-2}$','$\hat{\gamma}_{n=1,\ell=2,m=-1}$','$\hat{\gamma}_{n=1,\ell=2,m=0}$','$\hat{\gamma}_{n=1,\ell=2,m=1}$','$\hat{\gamma}_{n=1,\ell=2,m=2}$',...
            '$\hat{\gamma}_{n=2,\ell=2,m=-2}$','$\hat{\gamma}_{n=2,\ell=2,m=-1}$','$\hat{\gamma}_{n=2,\ell=2,m=0}$','$\hat{\gamma}_{n=2,\ell=2,m=1}$','$\hat{\gamma}_{n=2,\ell=2,m=2}$',...
            '$\hat{\gamma}_{n=3,\ell=2,m=-2}$','$\hat{\gamma}_{n=3,\ell=2,m=-1}$','$\hat{\gamma}_{n=3,\ell=2,m=0}$','$\hat{\gamma}_{n=3,\ell=2,m=1}$','$\hat{\gamma}_{n=3,\ell=2,m=2}$'};
PIPE.plotSlices(allMAPS, slice,clims,nametags,3,[],1,1),
% colormap(cbrewer2('iRdBu', 256));

% Plot signal rotational invariants
allMAPS=cat(4,S0_b,S2_b);
clims = [ repmat([0 5e-2],length(new_b),1) ; repmat([0 1e-2],length(new_b),1) ];
allMAPS=flip(permute(allMAPS,[2 1 3 4]),1);
slice=32;
figure('Position',[2678 677 1603 607]), colormap gray
nametags = {'$S_{0}(b=1)$','$S_{0}(b=2)$','$S_{0}(b=4)$','$S_{0}(b=6)$','$S_{2}(b=1)$','$S_{2}(b=2)$','$S_{2}(b=4)$','$S_{2}(b=6)$'};
PIPE.plotSlices(allMAPS, slice,clims,nametags,2,[],1,0),

% Plot SM kernel parameters
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$'};
slice=32;
allMAPS=flip(permute(kernel_fit,[2 1 3 4]),1);
c_lims=[0 1;1 3;1 3;0 1.5;0 1];
figure('Position',[152 803 2200 321]), colormap gray
PIPE.plotSlices(allMAPS, slice,c_lims,paramNames,1,[],1),

