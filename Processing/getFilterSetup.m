function  filter = getFilterSetup(params,data,ftype)

switch ftype

    case 'objf'
        ref_window_power = params.reg_window_power;
        
        % normalization factor
        reg_scale = 0.5 * data.obj.base_target_sz/data.setup.featureRatio;
        
        % construct grid
        wrg = -(data.obj.use_sz(1)-1)/2:(data.obj.use_sz(1)-1)/2;
        wcg = -(data.obj.use_sz(2)-1)/2:(data.obj.use_sz(2)-1)/2;
        [wrs, wcs] = ndgrid(wrg, wcg);
        
        % construct the regukarization window
        reg_window = (params.reg_window_edge - params.reg_window_min) * (abs(wrs/reg_scale(1)).^ref_window_power + abs(wcs/reg_scale(2)).^ref_window_power) + params.reg_window_min;
        filter.reg_window_SR = reg_window;
        
        if params.enableSSAR ==1
            reg_window = Get_SSAR(reg_window, params, data);
        end
        filter.reg_window = reg_window;
        
        % compute the DFT and enforce sparsity
        reg_window_dft = fft2(reg_window) / prod(data.obj.use_sz);
        reg_window_dft_sep = cat(3, real(reg_window_dft), imag(reg_window_dft));
        reg_window_dft_sep(abs(reg_window_dft_sep) < params.reg_sparsity_threshold * max(abs(reg_window_dft_sep(:)))) = 0;
        reg_window_dft = reg_window_dft_sep(:,:,1) + 1i*reg_window_dft_sep(:,:,2);
        
        % do the inverse transform, correct window minimum
        reg_window_sparse = real(ifft2(reg_window_dft));
        reg_window_dft(1,1) = reg_window_dft(1,1) - data.obj.support_sz * min(reg_window_sparse(:)) + params.reg_window_min;
        
        % construct the regularizsation matrix
        regW = cconvmtx2(reg_window_dft);
        
        regW_dfs = real(data.setup.dfs_matrix * regW * data.setup.dfs_matrix');
        
        filter.WW_block = regW_dfs' * regW_dfs;
        
        % If the filter size is small enough, remove small values in WW_block.
        % It takes too long time otherwise.
        if data.obj.support_sz <= 120^2
            filter.WW_block(0<abs(filter.WW_block) & abs(filter.WW_block)<0.00001) = 0;
        end
        
    case 'conf'
        % create weight window
        ref_window_power = params.reg_window_power;
        
        % normalization factor
        reg_scale = 0.5 * data.obj.base_target_sz/data.setup.featureRatio;
        
        % construct grid
        wrg = -(data.obj.use_sz(1)-1)/2:(data.obj.use_sz(1)-1)/2;
        wcg = -(data.obj.use_sz(2)-1)/2:(data.obj.use_sz(2)-1)/2;
        [wrs, wcs] = ndgrid(wrg, wcg);
        
        % construct the regukarization window
        reg_window = (params.reg_window_edge - params.reg_window_min) * (abs(wrs/reg_scale(1)).^ref_window_power + abs(wcs/reg_scale(2)).^ref_window_power) + params.reg_window_min;
        
        maxmax = max(max(reg_window));
        minmin = min(min(reg_window));
        mmdiff = maxmax + minmin;
        reg_window(:,:) = mmdiff - reg_window(:,:);

        % compute the DFT and enforce sparsity
        reg_window_dft = fft2(reg_window) / prod(data.obj.use_sz);
        reg_window_dft_sep = cat(3, real(reg_window_dft), imag(reg_window_dft));
        reg_window_dft_sep(abs(reg_window_dft_sep) < params.reg_sparsity_threshold * max(abs(reg_window_dft_sep(:)))) = 0;
        reg_window_dft = reg_window_dft_sep(:,:,1) + 1i*reg_window_dft_sep(:,:,2);
        
        % do the inverse transform, correct window minimum
        reg_window_sparse = real(ifft2(reg_window_dft));
        reg_window_dft(1,1) = reg_window_dft(1,1) - data.obj.support_sz * min(reg_window_sparse(:)) + params.reg_window_min;
        
        % construct the regularizsation matrix
        regW = cconvmtx2(reg_window_dft);
        
        regW_dfs = real(data.setup.dfs_matrix * regW * data.setup.dfs_matrix');
        
        filter.WW_block = regW_dfs' * regW_dfs;
        
        % If the filter size is small enough, remove small values in WW_block.
        % It takes too long time otherwise.
        if data.obj.support_sz <= 120^2
            filter.WW_block(0<abs(filter.WW_block) & abs(filter.WW_block)<0.00001) = 0;
        end
        
    case 'mixf'
        % else use a scaled identity matrix
        filter.WW_block = params.lambda * speye(data.obj.support_sz);
        params.reg_window_min = sqrt(params.lambda);
end

% create block diagonal regularization matrix
WW = eval(['blkdiag(filter.WW_block' repmat(',filter.WW_block', 1, data.setup.feature_dim-1) ');']);

WW_L = tril(WW);
filter.WW_U = triu(WW, 1);

% create indexing vectors
% first create the indices for the symmetrix (real) part of the spectrum
index_i_sym = zeros(2*data.setup.feature_dim, length(data.setup.dft_sym_ind), data.setup.feature_dim);
index_j_sym = zeros(size(index_i_sym));

index_i_sym_re = repmat(bsxfun(@plus, data.obj.support_sz*(0:data.setup.feature_dim-1)', 1:length(data.setup.dft_sym_ind)), [1 1 data.setup.feature_dim]); %index for the Real-Real part
index_i_sym(1:2:end, :, :) = index_i_sym_re;
index_i_sym(2:2:end, :, :) = NaN; % these will be zero

index_j_sym_re = permute(index_i_sym_re, [3 2 1]);
index_j_sym(1:2:end, :, :) = index_j_sym_re;
index_j_sym(2:2:end, :, :) = NaN; % these will be zero

% create the indices for the remaining part
index_i = zeros(2*data.setup.feature_dim, 2*length(data.setup.dft_pos_ind), data.setup.feature_dim);
index_j = zeros(size(index_i));

index_i_re = repmat(bsxfun(@plus, data.obj.support_sz*(0:data.setup.feature_dim-1)', (length(data.setup.dft_sym_ind)+1:2:data.obj.support_sz)), [1 1 data.setup.feature_dim]); %index for the Real-Real part
index_i(1:2:end, 1:2:end, :) = index_i_re;
index_i(2:2:end, 1:2:end, :) = index_i_re + 1;
index_i(1:2:end, 2:2:end, :) = index_i_re;
index_i(2:2:end, 2:2:end, :) = index_i_re + 1;

index_j_re = permute(index_i_re, [3 2 1]);
index_j(1:2:end, 1:2:end, :) = index_j_re;
index_j(2:2:end, 1:2:end, :) = index_j_re;
index_j(1:2:end, 2:2:end, :) = index_j_re + 1;
index_j(2:2:end, 2:2:end, :) = index_j_re + 1;

% concatenate the results
index_i = cat(2, index_i_sym, index_i);
index_j = cat(2, index_j_sym, index_j);

index_i = index_i(:);
index_j = index_j(:);

% the imaginary part of the autocorrelations (along the diagonal) will be zero
zero_ind = (index_i == index_j-1) | (index_i == index_j+1);
index_i(zero_ind) = NaN;
index_j(zero_ind) = NaN;

% indexing masks for upper and lower triangular part
filter.data_L_mask = index_i >= index_j;
filter.data_U_mask = index_i < index_j;

data_L_i = index_i(filter.data_L_mask);
data_L_j = index_j(filter.data_L_mask);
data_U_i = index_i(filter.data_U_mask);
data_U_j = index_j(filter.data_U_mask);

% extract the linear indeces from the data matrix and regularization matrix
WW_L_ind = find(WW_L);
data_L_ind = sub2ind(size(WW_L), data_L_i, data_L_j);

% compute the linear indeces of the non-zeros in the matrix
[L_ind, ~, data_WW_in_L_index] = unique([data_L_ind; WW_L_ind]);

% compute the corresponding indices for the values in the data and reg
% matrix
filter.data_in_L_index = uint32(data_WW_in_L_index(1:length(data_L_ind)));
WW_in_L_index = data_WW_in_L_index(length(data_L_ind)+1:end);

% create the arrays of values in the regularization matrix
nnz_L = length(L_ind);
WW_L_vec = zeros(nnz_L, 1, 'single');
WW_L_vec(WW_in_L_index) = full(WW_L(WW_L_ind));

% precompute the data part of the regularization matrix
filter.WW_L_vec_data = WW_L_vec(filter.data_in_L_index);

% initialize the content vectors for the matrices
filter.L_vec = WW_L_vec;

% preallocate the matrices
mat_size = data.setup.feature_dim* data.obj.support_sz;
[L_i, L_j] = ind2sub(size(WW_L), L_ind);
filter.AL = sparse(double(L_i), double(L_j), ones(double(nnz_L),1), double(mat_size), double(mat_size));
filter.AU_data = sparse(double(data_U_i), double(data_U_j), ones(length(double(data_U_i)),1), double(mat_size), double(mat_size));

filter.isinitial=0;
end