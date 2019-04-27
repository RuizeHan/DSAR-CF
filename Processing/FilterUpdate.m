function [params, data] = FilterUpdate(params, data)

%this is the training code used to update/initialize the tracker

%Compute coefficients for the tranlsation filter
pixels = get_pixels(data.seq.im,data.obj.pos,round(data.obj.sz*data.obj.currentScaleFactor),data.obj.sz);
xl = bsxfun(@times,get_features(pixels,params.t_features,params.t_global),data.setup.cos_window);

xlf = fft2(xl);
xlf_reshaped = reshape(xlf, [data.obj.support_sz, data.setup.feature_dim]);

% new rhs sample
xyf_corr = bsxfun(@times, data.obj.yf_vec, conj(xlf_reshaped));
xy_dfs = real(data.setup.dfs_matrix * double(xyf_corr));
new_hf_rhs = xy_dfs(:);

xlf_reshaped_sym = xlf_reshaped(data.setup.dft_sym_ind, :);    % extract the symmetric part of the spectrum x_0
xlf_reshaped_pos = xlf_reshaped(data.setup.dft_pos_ind, :);    % extract the positive part of the spectrum x_+

% compute autocorrelation
xxlf_sym = bsxfun(@times, conj(permute(xlf_reshaped_sym, [2 1])), permute(xlf_reshaped_sym, [3 1 2]));
xxlf_pos = bsxfun(@times, conj(permute(xlf_reshaped_pos, [2 1])), permute(xlf_reshaped_pos, [3 1 2]));
xxlf_pos_real = real(xxlf_pos);

% partition the real and imaginary parts
data.setup.xxlf_sep(1:2:end, 1:data.setup.num_sym_coef, :) = real(xxlf_sym);
data.setup.xxlf_sep(1:2:end, data.setup.num_sym_coef+1:2:end, :) = xxlf_pos_real;
data.setup.xxlf_sep(2:2:end, data.setup.num_sym_coef+1:2:end, :) = imag(xxlf_pos);
data.setup.xxlf_sep(1:2:end, data.setup.num_sym_coef+2:2:end, :) = -imag(xxlf_pos);
data.setup.xxlf_sep(2:2:end, data.setup.num_sym_coef+2:2:end, :) = xxlf_pos_real;

switch params.selector   
    case -1
        ftype = 'conf';
    case 0
        ftype = 'mixf';        
    case 1
        ftype = 'objf';        
end
       
filter = getfield(data,ftype);

if data.seq.frame == 1 || ~filter.isinitial
    filter.hf_rhs = new_hf_rhs;
    filter.hf_autocorr = data.setup.xxlf_sep(:);
    hf_init_autocorr = double(sum(xlf_reshaped .* conj(xlf_reshaped), 2));
    switch params.init_strategy
        case 'const_reg'       % exact solution for constant regularization
            hf_init = bsxfun(@rdivide, xyf_corr, hf_init_autocorr + params.reg_window_min^2);
            hf_init = real(data.setup.dfs_matrix * hf_init);
            filter.hf_vec = hf_init(:);
        case 'indep'           % independent filters for each feature
            A_init = real(data.setup.dfs_matrix * spdiags(hf_init_autocorr, 0, double(data.obj.support_sz), double(data.obj.support_sz)) * data.setup.dfs_matrix') + data.setup.feature_dim * filter.WW_block;
            b_init = reshape(filter.hf_rhs, data.obj.support_sz, data.setup.feature_dim);
            hf_init = A_init \ b_init;
            filter.hf_vec = hf_init(:);
    end
    filter.isinitial = 1;
else
    filter.hf_rhs = (1 - params.learning_rate) * filter.hf_rhs + params.learning_rate * new_hf_rhs;
    filter.hf_autocorr = (1 - params.learning_rate) * filter.hf_autocorr + params.learning_rate * data.setup.xxlf_sep(:);
end

% add the autocorrelation to the matrix vectors with the regularization
filter.L_vec(filter.data_in_L_index) = filter.hf_autocorr(filter.data_L_mask) + filter.WW_L_vec_data;

% update the matrices with the new non-zeros
filter.AL = setnonzeros(filter.AL, double(filter.L_vec));
filter.AU_data = setnonzeros(filter.AU_data, double(filter.hf_autocorr(filter.data_U_mask)));

time = toc;
%----------------------------Update the weight map------------------------
if params.selector == 1 && params.enableDSAR ==1 && rem(data.seq.frame+params.updateW-1,params.updateW) == 0 
   [params, data, filter] = update_w(params, data, filter);
end
%----------------------------Update the weight map------------------------
data.time.timeW = data.time.timeW + toc - time;

% do Gausss-Seidel
for iter = 1:params.num_GS_iter
    filter.hf_vec = filter.AL \ (filter.hf_rhs - filter.AU_data * filter.hf_vec - filter.WW_U * filter.hf_vec);
end

% reconstruct the filter
filter.hf = reshape(single(data.setup.dfs_matrix' * reshape(filter.hf_vec, [data.obj.support_sz, data.setup.feature_dim])), [data.obj.use_sz, data.setup.feature_dim]);

data = setfield(data,ftype,filter);

%% 
% if params.enableCSR && data.seq.frame == 1
%     data = RidgeRrgression(data, params);
% end

end