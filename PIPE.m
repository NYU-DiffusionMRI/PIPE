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
        % =================================================================
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
                else
                    gamma_ratios_aux = gamma(ids,:)./gamma_1lm;
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
            X_moments = Compute_extended_moments(data_training',Degree);
            X_moments_test = Compute_extended_moments(data_testing',Degree);        
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
    end
end

