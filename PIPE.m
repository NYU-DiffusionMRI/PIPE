classdef PIPE
    % PIPE: Protocol Independent Parameter Estimation class
    %
    % PIPE is an implementation of the Standard Model framework for
    % modeling white matter microstructure with diffusion MRI (supports
    % diffusion relaxometry data too)
    %
    % =====================================================================
    %
    %  Authors: Santiago Coelho (santiago.coelho@nyulangone.org), Els Fieremans, Dmitry Novikov
    %  Copyright (c) 2025 New York University
    %              
    %   Permission is hereby granted, free of charge, to any non-commercial entity ('Recipient') obtaining a 
    %   copy of this software and associated documentation files (the 'Software'), to the Software solely for
    %   non-commercial research, including the rights to use, copy and modify the Software, subject to the 
    %   following conditions: 
    % 
    %     1. The above copyright notice and this permission notice shall be included by Recipient in all copies
    %     or substantial portions of the Software. 
    % 
    %     2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
    %     NOT LIMITED TO THE WARRANTIESOF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    %     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    %     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE
    %     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
    % 
    %     3. In no event shall NYU be liable for direct, indirect, special, incidental or consequential damages
    %     in connection with the Software. Recipient will defend, indemnify and hold NYU harmless from any 
    %     claims or liability resulting from the use of the Software by recipient. 
    % 
    %     4. Neither anything contained herein nor the delivery of the Software to recipient shall be deemed to
    %     grant the Recipient any right or licenses under any patents or patent application owned by NYU. 
    % 
    %     5. The Software may only be used for non-commercial research and may not be used for clinical care. 
    % 
    %     6. Any publication by Recipient of research involving the Software shall cite the references listed
    %     below.
    %
    %  REFERENCES:
    %

    methods     ( Static = true )
        % =================================================================
        function x = get_Chebyshev_nodes(n,bounds)
        % x = get_Chebyshev_nodes(n,bounds)
        %
        % Get n Chebyshev nodes between bounds(1) and bounds(2)
        % 
        k=n:(-1):1;
        x=cos(pi*(2*k-1)./(2*n));
        if exist('bounds', 'var') && ~isempty(bounds)
            x=1/2*(x+1)*(bounds(2)-bounds(1))+bounds(1);
        end
        end
        % =================================================================
        function dirs = rand_dirs_S2(Ndirs)
        % dirs = rand_dirs_S2(Ndirs)
        %
        % returns randomly sampled points on the unit sphere in 3D
        dirs = randn(Ndirs,3);
        dirs = dirs./repmat(sqrt(sum(dirs.^2,2)),1,3);
        end
        % =================================================================
        function [alpha_nlm, gamma_nlm] = compute_alpha_gamma(dwi,target_b,target_dirs,library_path,Nl_fit,mask)
        % [alpha_nlm, gamma_nlm] = compute_alpha_gamma(dwi,target_b,target_dirs,library_path,Nl_fit,mask)
        
        LIB = load(library_path);
        if max(target_b(:))>100
            target_b = target_b/1000;
        end
        
        if ~exist('Nl_fit', 'var') || isempty(Nl_fit)
            Nl_fit = [ LIB.n0 LIB.n2 LIB.n4 LIB.n6 ];
        end
        total_nlm_lib = LIB.n0+LIB.n2*5+LIB.n4*9+LIB.n6*13;
        Nl_lib = [ LIB.n0 LIB.n2 LIB.n4 LIB.n6 ];
        nlm_fit = PIPE.get_nlm_from_Nl(Nl_fit,Nl_lib);

        if isvector(target_b)
            b0_ids = target_b<1e-3;
            NB_target = length(target_b);
        else
            b_isocenter = squeeze(target_b(round(end/2),round(end/2),round(end/2),:));
            b0_ids = b_isocenter<1e-3;
            [~,~,~,NB_target] = size(target_b);
        end

        sz = size(dwi);
        if ~exist('mask', 'var') || isempty(mask)
            mask = true(sz(1:3));
        else
            mask = logical(mask);
        end
        b0 = mean(dwi(:,:,:,b0_ids),4);
        dwi = RICEtools.vectorize(dwi./b0,mask);

        threshold_in_b = LIB.bmax; % Do not interpolate basis functions outside library range
        CS_phase = LIB.CS_phase;
        % Interpolating U for new b,beta values
        if isvector(target_b)
            alpha_nlm_all = zeros(NB_target,total_nlm_lib);
            U0_target = PIPE.Chebyshev_interpolation_U_b(LIB.U0(:,1:LIB.n0),target_b,LIB.bmax);
            U2_target = PIPE.Chebyshev_interpolation_U_b(LIB.U2(:,1:LIB.n2),target_b,LIB.bmax);
            U4_target = PIPE.Chebyshev_interpolation_U_b(LIB.U4(:,1:LIB.n4),target_b,LIB.bmax);
            U6_target = PIPE.Chebyshev_interpolation_U_b(LIB.U6(:,1:LIB.n6),target_b,LIB.bmax);
            Lmax=6;
            Ylm=SMI.get_even_SH(target_dirs,Lmax,CS_phase);
            Y00_target=Ylm(:,1);
            Y2m_target=Ylm(:,2:6);
            Y4m_target=Ylm(:,7:15);
            Y6m_target=Ylm(:,16:28);
            for ii=1:NB_target
                U0Y0=kron(U0_target(ii,:),Y00_target(ii,:));
                U2Y2=kron(U2_target(ii,:),Y2m_target(ii,:));
                U4Y4=kron(U4_target(ii,:),Y4m_target(ii,:));
                U6Y6=kron(U6_target(ii,:),Y6m_target(ii,:));
                alpha_nlm_all(ii,:)=[U0Y0 U2Y2 U4Y4 U6Y6];
            end
            alpha_nlm = alpha_nlm_all(:,nlm_fit);
            % gamma_nlm = alpha_nlm\(dwi/(4*pi));
            gamma_nlm = alpha_nlm\dwi;
        else
            bval_masked = PIPE.vectorize(target_b,mask);
            gx = PIPE.vectorize(target_dirs(:,:,:,:,1),mask);
            gy = PIPE.vectorize(target_dirs(:,:,:,:,2),mask);
            gz = PIPE.vectorize(target_dirs(:,:,:,:,3),mask);
            bvec_masked = permute(cat(3,gx,gy,gz),[1 3 2]);
            Nvox_mask = size(bval_masked,2);
            alpha_nlm = zeros(NB_target,length(nlm_fit),Nvox_mask);
            gamma_nlm = zeros(length(nlm_fit),Nvox_mask);
            Lmax=6;
            parfor ii=1:Nvox_mask  %parfor loop takes 4.36s for 1000 vox, for loop takes 19.76
                current_b = bval_masked(:,ii)';
                if any(current_b>threshold_in_b)
                    continue
                end
                current_dirs = bvec_masked(:,:,ii);
                Ylm = SMI.get_even_SH(current_dirs,Lmax,CS_phase);
                Y00 = Ylm(:,1);
                Y2m = Ylm(:,2:6);
                Y4m = Ylm(:,7:15);
                Y6m = Ylm(:,16:28);
                U0_target = PIPE.Chebyshev_interpolation_U_b(LIB.U0(:,1:LIB.n0),current_b,LIB.bmax);
                U2_target = PIPE.Chebyshev_interpolation_U_b(LIB.U2(:,1:LIB.n2),current_b,LIB.bmax);
                U4_target = PIPE.Chebyshev_interpolation_U_b(LIB.U4(:,1:LIB.n4),current_b,LIB.bmax);
                U6_target = PIPE.Chebyshev_interpolation_U_b(LIB.U6(:,1:LIB.n6),current_b,LIB.bmax);
                alpha_nlm_local = zeros(NB_target,total_nlm_lib);
                for ll=1:NB_target
                    U0Y0 = kron(U0_target(ll,:),Y00(ll,:));
                    U2Y2 = kron(U2_target(ll,:),Y2m(ll,:));
                    U4Y4 = kron(U4_target(ll,:),Y4m(ll,:));
                    U6Y6 = kron(U6_target(ll,:),Y6m(ll,:));
                    alpha_nlm_local(ll,:) = [U0Y0 U2Y2 U4Y4 U6Y6];
                end
                alpha_nlm(:,:,ii) = alpha_nlm_local(:,nlm_fit);
                % gamma_nlm(:,ii) = alpha_nlm_local(:,nlm_fit)\(dwi(:,ii)/(4*pi));
                gamma_nlm(:,ii) = alpha_nlm_local(:,nlm_fit)\dwi(:,ii);
            end
        end
        gamma_nlm = PIPE.vectorize(gamma_nlm,mask);

        end
        % =================================================================
        function [alpha_nlm] = compute_alpha(target_b,target_dirs,library_path,Nl_fit,mask)
        % [alpha_nlm] = compute_alpha(target_b,target_dirs,library_path,Nl_fit,mask)
        
        LIB = load(library_path);
        if max(target_b(:))>100
            target_b = target_b/1000;
        end


        if ~exist('Nl_fit', 'var') || isempty(Nl_fit)
            Nl_fit = [ LIB.n0 LIB.n2 LIB.n4 LIB.n6 ];
        end
        total_nlm_lib = LIB.n0+LIB.n2*5+LIB.n4*9+LIB.n6*13;
        Nl_lib = [ LIB.n0 LIB.n2 LIB.n4 LIB.n6 ];
        nlm_fit = PIPE.get_nlm_from_Nl(Nl_fit,Nl_lib);


        if isvector(target_b)
            b0_ids = target_b<0.1;
            NB_target = length(target_b);
        else
            b_isocenter = squeeze(target_b(end/2,end/2,end/2,:));
            b0_ids = b_isocenter<0.1;
            [~,~,~,NB_target] = size(target_b);
            sz = size(target_b);
            if ~exist('mask', 'var') || isempty(mask)
                mask = true(sz(1:3));
            else
                mask = logical(mask);
            end
        end

        threshold_in_b = LIB.bmax; % Do not interpolate basis functions outside library range
        CS_phase = LIB.CS_phase;
        % Interpolating U for new b,beta values
        if isvector(target_b)
            alpha_nlm_all = zeros(NB_target,total_nlm_lib);
            U0_target = PIPE.Chebyshev_interpolation_U_b(LIB.U0(:,1:LIB.n0),target_b,LIB.bmax);
            U2_target = PIPE.Chebyshev_interpolation_U_b(LIB.U2(:,1:LIB.n2),target_b,LIB.bmax);
            U4_target = PIPE.Chebyshev_interpolation_U_b(LIB.U4(:,1:LIB.n4),target_b,LIB.bmax);
            U6_target = PIPE.Chebyshev_interpolation_U_b(LIB.U6(:,1:LIB.n6),target_b,LIB.bmax);
            Lmax=6;
            Ylm=SMI.get_even_SH(target_dirs,Lmax,CS_phase);
            Y00_target=Ylm(:,1);
            Y2m_target=Ylm(:,2:6);
            Y4m_target=Ylm(:,7:15);
            Y6m_target=Ylm(:,16:28);
            for ii=1:NB_target
                U0Y0=kron(U0_target(ii,:),Y00_target(ii,:));
                U2Y2=kron(U2_target(ii,:),Y2m_target(ii,:));
                U4Y4=kron(U4_target(ii,:),Y4m_target(ii,:));
                U6Y6=kron(U6_target(ii,:),Y6m_target(ii,:));
                alpha_nlm_all(ii,:)=[U0Y0 U2Y2 U4Y4 U6Y6];
            end
            alpha_nlm = alpha_nlm_all(:,nlm_fit);
        else
            bval_masked = PIPE.vectorize(target_b,mask);
            gx = PIPE.vectorize(target_dirs(:,:,:,:,1),mask);
            gy = PIPE.vectorize(target_dirs(:,:,:,:,2),mask);
            gz = PIPE.vectorize(target_dirs(:,:,:,:,3),mask);
            bvec_masked = permute(cat(3,gx,gy,gz),[1 3 2]);
            Nvox_mask = size(bval_masked,2);
            alpha_nlm = zeros(NB_target,length(nlm_fit),Nvox_mask);
            Lmax=6;
            parfor ii=1:Nvox_mask  %parfor loop takes 4.36s for 1000 vox, for loop takes 19.76
                current_b = bval_masked(:,ii)';
                if any(current_b>threshold_in_b)
                    continue
                end
                current_dirs = bvec_masked(:,:,ii);
                Ylm = SMI.get_even_SH(current_dirs,Lmax,CS_phase);
                Y00 = Ylm(:,1);
                Y2m = Ylm(:,2:6);
                Y4m = Ylm(:,7:15);
                Y6m = Ylm(:,16:28);
                U0_target = PIPE.Chebyshev_interpolation_U_b(LIB.U0(:,1:LIB.n0),current_b,LIB.bmax);
                U2_target = PIPE.Chebyshev_interpolation_U_b(LIB.U2(:,1:LIB.n2),current_b,LIB.bmax);
                U4_target = PIPE.Chebyshev_interpolation_U_b(LIB.U4(:,1:LIB.n4),current_b,LIB.bmax);
                U6_target = PIPE.Chebyshev_interpolation_U_b(LIB.U6(:,1:LIB.n6),current_b,LIB.bmax);
                alpha_nlm_local = zeros(NB_target,total_nlm_lib);
                for ll=1:NB_target
                    U0Y0 = kron(U0_target(ll,:),Y00(ll,:));
                    U2Y2 = kron(U2_target(ll,:),Y2m(ll,:));
                    U4Y4 = kron(U4_target(ll,:),Y4m(ll,:));
                    U6Y6 = kron(U6_target(ll,:),Y6m(ll,:));
                    alpha_nlm_local(ll,:) = [U0Y0 U2Y2 U4Y4 U6Y6];
                end
                alpha_nlm(:,:,ii) = alpha_nlm_local(:,nlm_fit);
            end
        end
        end
        % =================================================================
        function kernel = SM_LTE_gamma2kernel(gamma_nlm,target_b,target_dirs,library_path,Nl_gamma,Nl_kernel,mask,sigma,MLTraining)
        % kernel = SM_LTE_gamma2kernel(gamma_nlm,target_b,target_dirs,library_path,Nl_gamma,Nl_kernel,mask,sigma,MLTraining)
        LIB = load(library_path);
        sz = size(gamma_nlm);
        if ~exist('mask', 'var') || isempty(mask)
            mask = true(sz(1:3));
        else
            mask = logical(mask);
        end
        if ~exist('sigma', 'var') || isempty(sigma)
            sigma = 1/100;
            fprintf('default SNR = 100\n')
        end

        if ~isfield(MLTraining,'Degree')
            Degree = 2; fprintf('default quadratic regression\n')
        else
            Degree = MLTraining.Degree;
        end
        if ~isfield(MLTraining,'Mtrain')
            Mtrain = 40000; fprintf('default training size = 40000\n')
        else
            Mtrain = MLTraining.Mtrain;
        end
        if ~isfield(MLTraining,'KernelBounds')
            lb=[0.2 1 1 0.15 0 50 30 0.1]; ub=[0.8 2.5 2.5 0.8 0.2 150 100 0.9];
            lb=[0.2 1 1 0.15 0 50 30 0.1]; ub=[0.8 2.5 2.5 0.8 0.2 150 100 0.9];
        else
            lb = MLTraining.KernelBounds(1,:);
            ub = MLTraining.KernelBounds(2,:);
        end

        alpha_target = PIPE.compute_alpha(target_b,target_dirs,library_path,Nl_gamma,mask);
        nlm_regression = PIPE.get_nlm_from_Nl(Nl_kernel,Nl_gamma);

        % Noise propagation SMI test (not for library, only for SMI regression test)
        Mtest = 5e3; M = Mtest + Mtrain;
        Lmax = 6;
        [f,Da,Depar,Deperp,f_w,T2a,T2e,~,~] = SMI.Get_uniformly_distributed_SM_prior(M,lb,ub,Lmax);
        f_extra=1-f-f_w;
        kernel_gt=[f, Da, Depar, Deperp, f_w, T2a,T2e];
        % Generate fODF
        Lmax = 6;
        Nfibers=2;%round(rand(M,1))+1;
        [plm,pl] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_Nfibers(M,Lmax,[],[],Nfibers);
        plm=plm';
        
        % Setting up target protocol
        NB_target = length(target_b);
        target_beta = ones(1,NB_target);  % The current library only supports beta = 1
        
        % S_target = SMI.SM_wFW_b_beta_TE_RealSphHarm_quadInt(f,Da,Depar,Deperp,f_extra,T2a,T2e,plm,target_b,target_dirs,target_beta,0*target_b,1202,LIB.CS_phase,0,LIB.D_FW);
        S_target = SMI.SM_wFW_b_beta_TE_RealSphHarm(f,Da,Depar,Deperp,f_extra,T2a,T2e,plm,target_b,target_dirs,target_beta,0*target_b,LIB.CS_phase,LIB.D_FW);
        mask_train = true([10 10 M/100]);
        S_target_noisy = PIPE.vectorize(S_target+randn(NB_target,M)*sigma,mask_train);
        
        id_train = 1:Mtrain;
        [alpha_target, GAMMA_hat] = PIPE.compute_alpha_gamma(S_target_noisy,target_b,target_dirs,library_path,Nl_gamma,mask_train);
        GAMMA_hat = PIPE.vectorize(GAMMA_hat,mask_train);
        options.training_data = GAMMA_hat(nlm_regression,id_train);   
        gamma_nlm = PIPE.vectorize(gamma_nlm,mask);
        options.test_data = gamma_nlm(nlm_regression,:);

        % PR fitting
        options.training_objective = kernel_gt(id_train,1:5)';
        options.Degree = Degree;
        options.strategy = 'PR';
        out = PIPE.DataDriven_regression(options);  
        kernel = PIPE.vectorize(out.fits',mask);

        % % % % Testing local function
        % % % id_test = (Mtrain+1):M;
        % % % options.test_data = GAMMA_hat(nlm_regression,id_test);
        % % % out_test = PIPE.DataDriven_regression(options);  
        % % % kernel_test_hat = out_test.fits;
        % % % kernel_test_gt = kernel_gt(id_test,:);
        % % % param_lims=[0.2 0.8;1 2.5;1 2.5;0 1;0 0.2];
        % % % paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$'};
        % % % figure('Position',[72 742 2137 401])
        % % % for ii=1:5
        % % %     kernel_target_fit(:,ii) = kernel_test_hat(:,ii);
        % % %     id_mins=kernel_target_fit(:,ii)<param_lims(ii,1);
        % % %     id_maxs=kernel_target_fit(:,ii)>param_lims(ii,2);
        % % %     kernel_target_fit(id_mins,ii)=param_lims(ii,1);
        % % %     kernel_target_fit(id_maxs,ii)=param_lims(ii,2);
        % % %     RMSE=sqrt(mean((kernel_test_gt(:,ii)-kernel_target_fit(:,ii)).^2));
        % % %     RMSEtag=['RMSE=',num2str(RMSE,3)];
        % % %     subplot(1,5,ii), 
        % % %     h = scatter_kde(kernel_test_gt(:,ii),kernel_target_fit(:,ii),'filled'); hold on, plot(param_lims(ii,:),param_lims(ii,:),'r-.'), axis([param_lims(ii,:) param_lims(ii,:)])
        % % %     title([paramNames{ii},' - ',RMSEtag],'interpreter','latex')
        % % %     xlabel('Ground truth','interpreter','latex'), ylabel('NN output','interpreter','latex'), set(gca,'FontSize',20)
        % % %     fprintf('================== Done with param %d/5\n',ii)
        % % % end
        end
        % =================================================================
        function [S0_b,S2_b] = project_gamma_onto_RotInvs(gamma_nlm,mask,target_b,library_path,Nl_fit)
        % 
        LIB = load(library_path);
        if ~exist('Nl_fit', 'var') || isempty(Nl_fit)
            Nl_fit = [ LIB.n0 LIB.n2 LIB.n4 LIB.n6 ];
        end
        Nl_lib = [ LIB.n0 LIB.n2 LIB.n4 LIB.n6 ];
        nlm_fit = PIPE.get_nlm_from_Nl(Nl_fit,Nl_lib);

        sz = size(gamma_nlm);
        if ~exist('mask', 'var') || isempty(mask)
            mask = true(sz(1:3));
        else
            mask = logical(mask);
        end
        Lmax = 6;
        dirs = [-0.130068563916747	0.0932848110280071	0.301874114147997	-0.160588543132534	-0.963953191993350	0.978837928613104	-0.549522980273078	0.592328122242095	-0.729164145170031	0.238633693018993	0.434651955799002	0.652635346397925	-0.734000052801543	-0.403599416152604	0.853698627179037	0.179812898840564	-0.511230028936703	-0.347227091988192	-0.151596735206056	0.719922327180808	-0.822392486806688	0.704707871446281	0.829703699851419	-0.862550086857427	-0.335321484736507	-0.314225258918558	0.194088504613566	0.327542356254105	0.505947381273477	-0.834296964584067	0.309476571411950	-0.117340713967168	-0.363227137830263	0.521389418163999	0.945240558338114	-0.524417749263304	-0.763218546997936	-0.145525976291921	-0.0122722916686933	0.685740007952837	-0.956922291701663	0.484894981408397	-0.223076366317911	-0.0912265656525496	0.939489850665938	-0.548677930291736	0.0899830129379917	-0.529095236280458	0.857923466019652	0.0313910318301205	-0.503828318019078	0.488811191935837	0.294732743377875	-0.956802130789944	0.300000411642112	0.00897408722442505	-0.709091302055868	0.650561378894920	-0.711236989891623	0.789508073435427 ; 0.282325584363997	-0.332500776697028	0.898485458170594	-0.968673668717596	0.261165366605327	-0.200111124452651	0.264272708506881	-0.285537557395265	-0.355775017966924	0.854283276620088	-0.848130577636177	0.307985185967729	-0.537419604056184	0.910143183919012	0.515134135617358	-0.483539314532076	0.478397264082337	-0.766740371788774	0.531267969479908	-0.684110754074832	-0.564282076688453	0.157821133257434	-0.133437282194303	0.257480038001757	-0.543060152023286	0.874112698528370	0.230046094506407	0.475172375772998	-0.530090908461732	0.190037107797886	-0.938443721190627	0.870791399069994	-0.419677486106270	0.848661761051783	0.130543357216993	-0.753119539980642	0.630507380213516	-0.179058930954968	-0.966700744188533	0.632133714249398	-0.162507387356043	-0.696844732289984	0.0147999107837224	0.687075735781306	0.218072671767995	-0.00591585331993446	-0.729496857445582	0.657056695425554	-0.330180891156839	0.998046624894397	-0.857385809697511	-0.101631079322322	0.0246017271241047	-0.166433335585371	0.635666172830116	-0.816956746605835	0.630753698155324	0.578870483786831	-0.154577580948680	-0.533415663960950 ; 0.950460116519448	-0.938478117766921	-0.318741118563859	0.189427145224130	0.0508615270291595	-0.0428000861997289	-0.792580866341241	0.753402746821313	0.584588561288273	0.461794373984095	-0.302906257113693	-0.692251565441700	-0.415239800191922	-0.0935248418712461	0.0763902891366770	0.856654570239605	0.713988736067953	-0.539946801878329	-0.833530308074974	0.117065447446015	0.0725005901498665	0.691721985930035	-0.542020536671319	-0.435558696036320	0.769834510241357	0.370390978466020	-0.953629092875851	0.816656119894876	-0.680456373435304	0.517527267441655	0.153452383184467	-0.477444757172385	-0.831827418398286	0.0890297138886364	0.299129936247056	0.397236683548934	-0.141272407124139	0.973015975958025	-0.255615066149348	-0.360787484264772	-0.240605645609469	0.528473533949100	-0.974688615626135	0.720836074998898	-0.264202820430936	0.836012877586388	-0.678039373753192	-0.536968090245386	0.393634228145500	0.0540142356934065	-0.105102803442199	-0.866449503637717	0.955262982117209	0.238389654393286	-0.711286348620636	0.576629118180567	0.315181372565219	0.491608436960495	-0.685746101467561	-0.303553506686020 ];
        Ndirs = 60;
        Ylm=SMI.get_even_SH(dirs,Lmax,LIB.CS_phase);
        Y00=Ylm(:,1);
        Y2m=Ylm(:,2:6);
        Y4m=Ylm(:,7:15);
        Y6m=Ylm(:,16:28);
        gamma_nlm = PIPE.vectorize(gamma_nlm,mask);
        S0_b = []; S2_b = []; %S4_b = [];
        for kk = 1:length(target_b)
            U0S0_new_b_flat = PIPE.Chebyshev_interpolation_U_b(LIB.U0(:,1:LIB.n0),target_b(kk)*ones(1,Ndirs),LIB.bmax);
            U2S2_new_b_flat = PIPE.Chebyshev_interpolation_U_b(LIB.U2(:,1:LIB.n2),target_b(kk)*ones(1,Ndirs),LIB.bmax);
            U4S4_new_b_flat = PIPE.Chebyshev_interpolation_U_b(LIB.U4(:,1:LIB.n4),target_b(kk)*ones(1,Ndirs),LIB.bmax);
            U6S6_new_b_flat = PIPE.Chebyshev_interpolation_U_b(LIB.U6(:,1:LIB.n6),target_b(kk)*ones(1,Ndirs),LIB.bmax);
            alpha_nlm_flat_newb = [];
            for ll=1:Ndirs
                US0Y0=kron(U0S0_new_b_flat(ll,:),Y00(ll,:));
                US2Y2=kron(U2S2_new_b_flat(ll,:),Y2m(ll,:));
                US4Y4=kron(U4S4_new_b_flat(ll,:),Y4m(ll,:));
                US6Y6=kron(U6S6_new_b_flat(ll,:),Y6m(ll,:));
                alpha_nlm_flat_newb(ll,:)=[US0Y0 US2Y2 US4Y4 US6Y6];
            end
            dwi_newb = alpha_nlm_flat_newb(:,nlm_fit) * gamma_nlm ; 
        
            % Spherical harmonics fit
            dwi_newb = PIPE.vectorize(dwi_newb,mask);
            [~,Sl,~,table_4D_sorted] = SMI.Fit2D4D_LLS_RealSphHarm_wSorting_norm_varL(dwi_newb,[],ones(1,Ndirs),dirs,[],[],4,[]);
            currentS0 = squeeze(Sl(:,:,:,1,:));
            currentS2 = squeeze(Sl(:,:,:,2,:));
        %     S4=squeeze(Sl(:,:,:,3,:));
            S0_b = cat(4,S0_b,currentS0);
            S2_b = cat(4,S2_b,currentS2);
        end
        
        end
        % =================================================================
        function U_target = Chebyshev_interpolation_U_b(U,b_target,bmax)
        % U_target = Chebyshev_interpolation_U_b(U,b_target,bmax)
        Nkeep=size(U,2);
        NB_target=length(b_target);
        U_target=zeros(NB_target,Nkeep);
        for ii=1:Nkeep
            y=U(:,ii)';
            U_interp=chebyshevInterpolate(y,b_target,[0 bmax]);
            U_target(:,ii)=U_interp';
        end
        end
        % =================================================================
        function [gamma_ratios_average,gamma_ratios_all] = ComputeRatioRotationallyInvariantGamma(gamma,n_ell)
        % [gamma_ratios_average,gamma_ratios_all] = ComputeRatioRotationallyInvariantGamma(gamma,n_ell)
        
        Lmax = (length(n_ell)-1)*2;
        M_target = size(gamma,2);
        gamma_ratios_average = zeros(sum(n_ell),M_target);
        gamma_ratios_all = zeros(size(gamma));
        L_index = repelem(0:2:Lmax,n_ell(:)');
        n_index = L_index * 0;
        unique_vals = unique(L_index, 'stable'); % Preserve original order
        for ii = 1:length(unique_vals)
            idx = find(L_index == unique_vals(ii));
            n_index(idx) = 1:length(idx);
        end
        % disp([L_index ; n_index])
        past=1;
        for jj=1:sum(n_ell)
            if L_index(jj)==0
                ids=past;
            elseif L_index(jj)==2
                ids=past:(past+4);
            elseif L_index(jj)==4
                ids=past:(past+8);
            elseif L_index(jj)==6
                ids=past:(past+12);
            elseif L_index(jj)==8
                ids=past:(past+16);
            end
        
            if L_index(jj) == 0
                gamma_ratios_aux = gamma(ids,:);
            else
                if n_index(jj) == 1
                    gamma_1lm = gamma(ids,:);
                    gamma_ratios_aux = sqrt(sum(gamma_1lm.^2,1));
                    current_ri = gamma_ratios_aux;
                else
                    gamma_ratios_aux = (gamma(ids,:) ./ gamma_1lm ) .* current_ri;
                end
            end
            gamma_ratios_average(jj,:)=mean(gamma_ratios_aux,1);
            gamma_ratios_all(ids,:)=repmat(mean(gamma_ratios_aux,1),length(ids),1);
        
            past=ids(end)+1;
            % display(ids)
        end

        end
        % =================================================================
        function gamma_nlm_dn = rank1denoisegamma(gamma_nlm_noisy,id_l,ell)
        % gamma_nlm_dn = rank1denoisegamma(gamma_nlm_noisy,id_l,ell)
        Nelems = size(gamma_nlm_noisy,2);
        Nell = length(id_l)/(2*ell+1);
        gamma_nellm_noisy = gamma_nlm_noisy(id_l,:);
        gamma_nellm_noisy = reshape(gamma_nellm_noisy,[2*ell+1 Nell Nelems]);
        
        gamma_nlm_dn = gamma_nlm_noisy;
        
        gamma_nellm_dn = 0 * gamma_nellm_noisy;
        for ii = 1:Nelems
            gamma_nellm_dn(:,:,ii) = PIPE.LowRankDenoising(gamma_nellm_noisy(:,:,ii),1);
        end
        gamma_nellm_dn = reshape(gamma_nellm_dn,[(2*ell+1)*Nell Nelems]);
        gamma_nlm_dn(id_l,:) = gamma_nellm_dn;
        end
        % =================================================================
        function Xdn = LowRankDenoising(X,p)
            % Xdn = LowRankDenoising(X,p)
            % 
            % Fixed p low-rank denoising, does svd and keeps p largest singular values
            [U,S,V] = svd(X);
            Xdn= U(:,1:p)*S(1:p,1:p)*V(:,1:p)';
        end
        % =================================================================
        function out = DataDriven_regression(options)
        % [NET_TRAINED,TR] = Train_NN_for_function_approximation(training_data,training_objective,Nneurons,NhiddenLayers,Nrounds)
        %
        % training_data is [DIM x TrainingSize]
        % training_objective is [1 x TrainingSize] or [DIM_output x TrainingSize]
        %
        % options.strategy = 'NN' or 'PR'
        %
        %
        %
        % Nneurons is the # of neurons per layer (NhiddenLayers is # of layers) 
        % Nrounds is the amount of times the NN is trained (random initializations)
        %
        % By: Santiago Coelho (26/02/2020)
        training_data=options.training_data;
        training_objective=options.training_objective;
        test_data=options.test_data;
        DIM=size(training_data,1);
        if strcmp(options.strategy,'NN')||strcmp(options.strategy,'nn')
            if isfield(options, 'Nneurons')
                Nneurons=options.Nneurons;
            else
                Nneurons=ceil(1.5*DIM);
            end
            if isfield(options, 'NhiddenLayers')
                NhiddenLayers=options.NhiddenLayers;
            else
                NhiddenLayers=3;
            end
            if isfield(options, 'Nrounds')
                Nrounds=options.Nrounds;
            else
                Nrounds=5;
            end
            if isfield(options, 'useGPU')
                useGPU=options.useGPU;
            else
                useGPU='no';
            end
            trainFcn = 'trainscg';  % Levenberg-Marquardt backpropagation.
            % Create a Fitting Network
            hiddenLayerSize = Nneurons*ones(1,NhiddenLayers);
            net = fitnet(hiddenLayerSize,trainFcn);
            % Setup Division of Data for Training, Validation, Testing
            net.divideParam.trainRatio = 70/100;
            net.divideParam.valRatio = 15/100;
            net.divideParam.testRatio = 15/100;
            net.trainParam.showWindow = false;
            for round_id=1:Nrounds
                % =========================================================================
                % Train the Network
                [NET_trained{round_id},tr{round_id}] = train(net,training_data,training_objective,'useGPU',useGPU,'showResources','no');
                RMSE_best_perf(round_id)=tr{round_id}.best_tperf;
                % =========================================================================
            end
            [~,id_best_net]=min(RMSE_best_perf);
            out.NET_TRAINED=NET_trained{id_best_net};
            out.TR=tr{id_best_net};
            out.fits = out.NET_TRAINED(test_data,'useGPU',useGPU,'showResources','no');
            out.RMSE_training=RMSE_best_perf(id_best_net);
        elseif strcmp(options.strategy,'PR')||strcmp(options.strategy,'pr')
            if isfield(options, 'Degree')
                Degree=options.Degree;
            else
                Degree=2;
            end
            if isfield(options, 'svdFLAG')
                svdFLAG=options.svdFLAG;
                keep_sv=options.keep_sv;
            else
                svdFLAG=0;
            end
            if svdFLAG&&(DIM>1)
                % =========================================================================
                % Applying PCA for dimensionality reduction
                [coeff, score] = pca(training_data');
                reducedDimension = coeff(:,1:keep_sv);
                data_training = (training_data' * reducedDimension)';
                % =========================================================================
                % Applying training transformation to test dataset (optimal)
                data_testing = (test_data' * reducedDimension)';
            else
                data_training=training_data;
                data_testing=test_data;
            end
            % =========================================================================
            % Polynomial Regression FITTING
            % =========================================================================
            X_moments = PIPE.Compute_extended_moments(data_training',Degree);
            X_moments_test = PIPE.Compute_extended_moments(data_testing',Degree);        
            pinvX=pinv(X_moments);
            coeffs=pinvX*(training_objective');   
            training_obj=X_moments*coeffs;
            clearvars pinvX
            % Applying PR
            out.fits=X_moments_test*coeffs;
            out.RMSE_training=sqrt(mean((training_obj-training_objective').^2));
        end
        out.strategy=options.strategy;
        
        end
        % =================================================================
        function plotSlices(ARRAY_4D, slice,clims,names,Nrows,positions,nanTransparent,colorbar_flag)
            %
            sz=size(ARRAY_4D);
            if length(sz)==4
                Nplots=sz(4);
            elseif length(sz)==3
                Nplots=sz(3);
            end

            if ~exist('positions', 'var') || isempty(positions)
                positions=1:Nplots;
            end
            if ~exist('nanTransparent', 'var') || isempty(nanTransparent)
                nanTransparent=0;
            end
            if ~exist('colorbar_flag', 'var') || isempty(colorbar_flag)
                colorbar_flag=1;
            end

            if isvector(clims)
                clims=repmat(clims,Nplots,1);
            end

            for ii=1:Nplots
                subplot(Nrows,ceil(Nplots/Nrows),positions(ii))
                if ~nanTransparent
                    if length(sz)==4
                        imagesc(ARRAY_4D(:,:,slice,ii),clims(ii,:))
                    elseif length(sz)==3
                        imagesc(ARRAY_4D(:,:,ii),clims(ii,:))
                    end
                else
                    if length(sz)==4
                        h=pcolor(ARRAY_4D(:,:,slice,ii)); caxis(clims(ii,:)),set(h, 'EdgeColor', 'none');
                    elseif length(sz)==3
                        h=pcolor(ARRAY_4D(:,:,ii)); caxis(clims(ii,:)),set(h, 'EdgeColor', 'none');
                    end
                end

                if isempty(names)
                    title(['case ',num2str(ii)],'interpreter','latex')
                else
                    title(names{ii},'interpreter','latex')
                end
                set(gca,'FontSize',30), axis off, grid off

                if colorbar_flag
                cb=colorbar('south'); 
            %     cb_pos =  cb.Position;cb_pos(3)=0.03;
            %     set(cb,'position',cb_pos)
                cb.Ticks=clims(ii,:);    
                end
            end
        end        
        % =================================================================
        function ids_nlm = get_nlm_from_Nl(Nl_fit,Nl_lib)
        % ids_nlm = get_nlm_from_Nl(Nl_fit,Nl_lib)
        %
        ell = 0:2:6;
        elems_per_sector = cumsum((2*ell+1).*Nl_lib);
        ids_nlm = [ 1:Nl_fit(1),...
                    (elems_per_sector(1) + 1):(elems_per_sector(1) + Nl_fit(2)*5 ),...
                    (elems_per_sector(2) + 1):(elems_per_sector(2) + Nl_fit(3)*9 ),...
                    (elems_per_sector(3) + 1):(elems_per_sector(3) + Nl_fit(4)*13 )];
        % elems_per_sector = (2*ell+1).*Nl_lib;
        % ids_nlm = [ 1:Nl_fit(1),...
        %             (elems_per_sector(1) + 1):(elems_per_sector(1) + Nl_fit(2)*5 ),...
        %             (elems_per_sector(2) + 1):(elems_per_sector(2) + Nl_fit(3)*9 ),...
        %                                 (elems_per_sector(3) + 1):(elems_per_sector(3) + Nl_fit(4)*13 )];
        end
        % =================================================================
        function [s, mask] = vectorize(S, mask)
            if nargin == 1
               mask = ~isnan(S(:,:,:,1));
            end
            if ismatrix(S)
                n = size(S, 1);
                [x, y, z] = size(mask);
                s = NaN([x, y, z, n], 'like', S);
                for i = 1:n
                    tmp = NaN(x, y, z, 'like', S);
                    tmp(mask(:)) = S(i, :);
                    s(:,:,:,i) = tmp;
                end
            else
                for i = 1:size(S, 4)
                   Si = S(:,:,:,i);
                   s(i, :) = Si(mask(:));
                end
            end
        end
        % =================================================================
        function [plm,pl] = rand_REALplm_powerLaw_DNnorm_3Ea_Nfibers(M,Lmax,bounds_p2,bounds_lambda,N)
        %
        % [plm,pl] = rand_REALplm_powerLaw_DNnorm_3Ea_Nfibers(M,Lmax,bounds_p2,bounds_lambda,N)
        %
        % Output: plm is [Nplm x M] where M is the number of samples (input)
        %         and Nplm = Lmax*(Lmax+3)/2 is the number of plm coefficients
        %
        % By: Santiago Coelho (2021/10/06)
        if ~exist('Lmax', 'var') || isempty(Lmax)
            Lmax=2;
        end
        if ~exist('bounds_p2', 'var') || isempty(bounds_p2)
            bounds_p2=[0.02 0.9];
        end
        if ~exist('bounds_lambda', 'var') || isempty(bounds_lambda)
            bounds_lambda=[0.5 0.9];
        end
        if ~exist('N', 'var') || isempty(N)
            N=1;
        end
        if N>3, error('N can only be 1, 2, or 3'), end
        if N==1
            [plm,pl] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda);
        elseif N==2
            f=rand(1,M)*0.7+0.15;
            [plm1,~] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda);
            [plm2,~] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda);
            plm=plm1.*f+plm2.*(1-f);
        elseif N==3
            f1_sq=rand(1,M)*0.7+0.15;
            f2_sq=rand(1,M)*0.7+0.15;
            [f1,f2] = PIPE.From_unit_square_to_unit_triangle(f1_sq,f2_sq);
            [plm1,~] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda);
            [plm2,~] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda);
            [plm3,~] = PIPE.rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda);
            plm=plm1.*f1+plm2.*f2+plm3.*(1-f1-f2);
        end
        p2=sqrt(sum(plm(1:5,:).^2,1));
        pl=p2';
        if Lmax>=4
            p4=sqrt(sum(plm(6:14,:).^2,1));
            pl=[p2;p4]';
        end
        if Lmax>=6
            p6=sqrt(sum(plm(15:27,:).^2,1));
            pl=[p2;p4;p6]';
        end
        if Lmax>=8
            p8=sqrt(sum(plm(28:44,:).^2,1));
            pl=[p2;p4;p6;p8]';
        end
        if Lmax>=10
            p10=sqrt(sum(plm(45:65,:).^2,1));
            pl=[p2;p4;p6;p8;p10]';
        end
        end
        % =================================================================
        function [plm,pl] = rand_REALplm_powerLaw_DNnorm_3Ea_helper(M,Lmax,bounds_p2,bounds_lambda)
        %
        % plm = rand_REALplm_powerLaw_DNnorm_3Ea(M,Lmax,bounds_p2,bounds_lambda)
        %
        % Output: plm is [Nplm x M] where M is the number of samples (input)
        %         and Nplm = Lmax*(Lmax+3)/2 is the number of plm coefficients
        %
        % By: Santiago Coelho (2021/04/27)
        if ~exist('Lmax', 'var') || isempty(Lmax)
            Lmax=2;
        end
        if ~exist('bounds_p2', 'var') || isempty(bounds_p2)
            bounds_p2=[0.02 0.9];
        end
        if ~exist('bounds_lambda', 'var') || isempty(bounds_lambda)
            bounds_lambda=[0.5 0.9];
        end
        l=4:2:10;
        p2=rand(M,1)*(bounds_p2(2)-bounds_p2(1))+bounds_p2(1);
        if any(~isfinite(bounds_lambda))
            % This is too much
        %     p4=p2.*10.^(rand(M,1)*-1);
        %     p6=p4.*10.^(rand(M,1)*-1);
        %     p8=p6.*10.^(rand(M,1)*-1);
        %     p10=p8.*10.^(rand(M,1)*-1);
        %     pl=[ones(M,1) p2 p4 p6 p8 p10];
            
            % This perturbation from the exponential is better (plus combined with N>1)
            bounds_lambda=[0.5 0.9];
            lambda=rand(M,1)*(bounds_lambda(2)-bounds_lambda(1))+bounds_lambda(1);
            C=p2./lambda.^2;
            pl=[ones(M,1) p2 C.*lambda.^l];
            pert=[ones(M,1) rand(M,5)*0.2+0.9];
            pl=pert.*pl;
        else
            lambda=rand(M,1)*(bounds_lambda(2)-bounds_lambda(1))+bounds_lambda(1);
            C=p2./lambda.^2;
            pl=[ones(M,1) p2 C.*lambda.^l];
        end
        p2m=[zeros(M,2) pl(:,2) zeros(M,2)]';
        p4m=[zeros(M,4) pl(:,3) zeros(M,4)]';
        p6m=[zeros(M,6) pl(:,4) zeros(M,6)]';
        p8m=[zeros(M,8) pl(:,5) zeros(M,8)]';
        p10m=[zeros(M,10) pl(:,6) zeros(M,10)]';
        
        p2m_r=p2m*0;
        p4m_r=p4m*0;
        p6m_r=p6m*0;
        p8m_r=p8m*0;
        p10m_r=p10m*0;
        
        s = PIPE.rand_sph(M)';
        [phi,theta,~]=cart2sph(s(1,:),s(2,:),s(3,:)); theta=pi/2-theta;
        gamma=rand(M,1)*2*pi;
        r = eul2rotm([phi', theta', gamma],'ZYZ'); 
        for ii=1:M
        %     R0=rotz(gamma(ii)*180/pi);
        %     R1=roty(theta(ii)*180/pi);
        %     R2=rotz(phi(ii)*180/pi);
        %     r=R0*R2*R1;
        %     r=eul2rotm([phi(ii), theta(ii), gamma(ii)],'ZYZ');
            [R] = SHRotate(r(:,:,ii), Lmax);
            p2m_r(:,ii)=R{3}*p2m(:,ii);
            if Lmax>=4
                p4m_r(:,ii)=R{5}*p4m(:,ii);
            end
            if Lmax>=6
                p6m_r(:,ii)=R{7}*p6m(:,ii);
            end
            if Lmax>=8
                p8m_r(:,ii)=R{9}*p8m(:,ii);
            end
            if Lmax>=10
                p10m_r(:,ii)=R{11}*p10m(:,ii);
            end
        end
        if Lmax==2
            plm=p2m_r;
        elseif Lmax==4
            plm=[p2m_r;p4m_r];
        elseif Lmax==6
            plm=[p2m_r;p4m_r;p6m_r];
        elseif Lmax==8
            plm=[p2m_r;p4m_r;p6m_r;p8m_r];
        elseif Lmax==10
            plm=[p2m_r;p4m_r;p6m_r;p8m_r;p10m_r];
        end
        % if nargout==2
            pl=pl(:,2:(Lmax/2+1));
        % end
        end
        % =================================================================
        function [x_triangle,y_triangle] = From_unit_square_to_unit_triangle(x_square,y_square)
        % [x_triangle,y_triangle] = From_unit_square_to_unit_triangle(x_square,y_square)
        %
        % From x,y \in [0,1] to xp,yp whose sum is in [0,1]
        %
        % Inverse mapping given by From_unit_triangle_to_unit_square
        %
        % By: Santiago Coelho (22/09/2019)
        sz_x=size(x_square);
        sz_y=size(y_square);
        x=x_square(:);
        y=y_square(:);
        flag1=x>y;
        flag2=x<y;
        xp=0*x;
        yp=0*y;
        xp(flag1)=x(flag1)-y(flag1)/2;
        yp(flag1)=y(flag1)/2;
        xp(flag2)=x(flag2)/2;
        yp(flag2)=y(flag2)-x(flag2)/2;
        x_triangle=reshape(xp,sz_x);
        y_triangle=reshape(yp,sz_y);
        end
        % =================================================================
        function Moments_Extended = Compute_extended_moments(Moments,degree)
        % Moments_Extended = Compute_extended_moments(Moments,degree)
        %
        % This function computes the extended matrix of moments for polynomial
        % interpolation in multiple dimensions.
        %
        % Moments is [NxD] where N is the number of samples and D the dimension of
        % the problem (the number of different moments)
        %
        % e.g. in 1D this extended matrix is [1 X X^2 ... X^degree]
        %
        %
        % By: Santiago Coelho
        N=size(Moments,1);
        N_moments=size(Moments,2);
        % Preparing data points for polynomial fitting
        x_0=ones(N,1);
        if degree>=1
            x_1=Moments;
            X=[x_0 x_1];
        else
                X=x_0;
        end
        if degree>=2
            for current_degree=2:degree
                combinations=nchoosek(1:(N_moments+current_degree-1),current_degree);
                NtotalCombs=size(combinations,1);
                flags=zeros(NtotalCombs,N_moments+current_degree-1);
                sums=zeros(NtotalCombs,N_moments+current_degree-1);
                which=zeros(NtotalCombs,current_degree);
                for ii=1:NtotalCombs
                    flags(ii,combinations(ii,:))=1;
                    sums(ii,:)=cumsum(~flags(ii,:))+1;
                    which(ii,:)=sums(ii,combinations(ii,:));
                end
                Currentx=ones(N,NtotalCombs);
                for ii=1:NtotalCombs
                    for jj=1:current_degree
                        current_id=which(ii,jj);
                        Currentx(:,ii)=Currentx(:,ii).*Moments(:,current_id);
                    end
                end
                X=[X Currentx];
            end
        end
        Moments_Extended=X;
        end
        % =================================================================
        function s = rand_sph(nsamples)
        % s = rand_sph(nsamples)
        % 
        % This function draws randomly from the unit sphere in 3D
        if isempty(nsamples)
            nsamples = 1;
        end
        N = randn(nsamples,3);
        den = sqrt(sum(N.^2,2));
        s = N./repmat(den,1,3);        
        end
        % =================================================================
    end
end

