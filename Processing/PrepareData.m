function [params, data] = PrepareData(params)

data.obj = [];
data.objf = [];
data.conf = [];
data.mixf = [];
data.seq = [];
data.setup =[];
data.PSR_avg = 0;
data.max_res_avg = 0;
data.APCE_avg = 0;
data.time.timeS = 0;
data.time.timeW = 0;

data.obj.pos = floor(params.init_pos);
data.obj.target_sz = floor(params.wsize);
data.seq.num_frames = numel(params.s_frames);
init_target_sz = data.obj.target_sz;
data.setup.featureRatio = params.t_global.cell_size;

search_area = prod(init_target_sz / data.setup.featureRatio * params.search_area_scale);

if isfield(params.t_global, 'cell_selection_thresh')
    if search_area < params.t_global.cell_selection_thresh * params.translation_model_max_area
        params.t_global.cell_size = min(data.setup.featureRatio, max(1, ceil(sqrt(prod(init_target_sz * params.search_area_scale)/(params.t_global.cell_selection_thresh * params.translation_model_max_area)))));
        
        data.setup.featureRatio = params.t_global.cell_size;
        search_area = prod(init_target_sz / data.setup.featureRatio * params.search_area_scale);
    end
end

if search_area > params.translation_model_max_area
    data.obj.currentScaleFactor = sqrt(search_area / params.translation_model_max_area);
else
    data.obj.currentScaleFactor = 1.0;
end

% target size at the initial scale
data.obj.base_target_sz = data.obj.target_sz / data.obj.currentScaleFactor;

%window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        data.obj.sz = floor( data.obj.base_target_sz * params.search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        data.obj.sz = repmat(sqrt(prod(data.obj.base_target_sz * params.search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        data.obj.sz = data.obj.base_target_sz + sqrt(prod(data.obj.base_target_sz * params.search_area_scale) + (data.obj.base_target_sz(1) - data.obj.base_target_sz(2))/4) - sum(data.obj.base_target_sz)/2; % const padding
    case 'custom'
        data.obj.sz = [data.obj.base_target_sz(1)*2 data.obj.base_target_sz(2)*4]; % for testing
end

% set the size to exactly match the cell size
data.obj.sz = round(data.obj.sz / data.setup.featureRatio) * data.setup.featureRatio;
data.obj.use_sz = floor(data.obj.sz/data.setup.featureRatio);

% construct the label function
output_sigma = sqrt(prod(floor(data.obj.base_target_sz/data.setup.featureRatio))) * params.output_sigma_factor;
rg = circshift(-floor((data.obj.use_sz(1)-1)/2):ceil((data.obj.use_sz(1)-1)/2), [0 -floor((data.obj.use_sz(1)-1)/2)]);
cg = circshift(-floor((data.obj.use_sz(2)-1)/2):ceil((data.obj.use_sz(2)-1)/2), [0 -floor((data.obj.use_sz(2)-1)/2)]);
[rs, cs] = ndgrid( rg,cg);
data.obj.y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf = fft2(data.obj.y);

if params.interpolate_response == 1
    data.obj.interp_sz = data.obj.use_sz * data.setup.featureRatio;
else
    data.obj.interp_sz = data.obj.use_sz;
end

% construct cosine window
data.setup.cos_window = single(hann(data.obj.use_sz(1))*hann(data.obj.use_sz(2))');

% the search area size
data.obj.support_sz = prod(data.obj.use_sz);

% Calculate feature dimension
data.seq.im = imread(params.s_frames{1});
if size(data.seq.im,3) == 3
    if all(all(data.seq.im(:,:,1) == data.seq.im(:,:,2)))
        data.seq.colorImage = false;
    else
        data.seq.colorImage = true;
    end
else
    data.seq.colorImage = false;
end

data.setup.feature_dim = 0;
for n = 1:length(params.t_features)
    
    if ~isfield(params.t_features{n}.fparams,'useForColor')
        params.t_features{n}.fparams.useForColor = true;
    end;
    
    if ~isfield(params.t_features{n}.fparams,'useForGray')
        params.t_features{n}.fparams.useForGray = true;
    end;
    
    if (params.t_features{n}.fparams.useForColor && data.seq.colorImage) || (params.t_features{n}.fparams.useForGray && ~data.seq.colorImage)
        data.setup.feature_dim= data.setup.feature_dim+ params.t_features{n}.fparams.nDim;
    end;
end;

if size(data.seq.im,3) > 1 && data.seq.colorImage == false
    data.seq.im = data.seq.im(:,:,1);
end

% compute the indices for the real, positive and negative parts of the
% spectrum
[data.setup.dft_sym_ind, data.setup.dft_pos_ind, data.setup.dft_neg_ind] = partition_spectrum2(data.obj.use_sz);

% the discrete fourier series output indices
dfs_sym_ind = (1:length(data.setup.dft_sym_ind))';
dfs_real_ind = dfs_sym_ind(end) - 1 + 2 * (1:length(data.setup.dft_pos_ind))';
dfs_imag_ind = dfs_sym_ind(end) + 2 * (1:length(data.setup.dft_pos_ind))';

% construct the transformation matrix from dft to dfs (the real fourier
% series)
data.setup.dfs_matrix = dft2dfs_matrix(double(data.setup.dft_sym_ind), double(data.setup.dft_pos_ind), data.setup.dft_neg_ind, dfs_sym_ind, dfs_real_ind, dfs_imag_ind);

% create vectorized desired correlation output
data.obj.yf_vec = single(yf(:));
data.setup.num_sym_coef = length(data.setup.dft_sym_ind);

%% generate filters setup for objf conf
data.objf = getFilterSetup(params,data,'objf');
data.conf = getFilterSetup(params,data,'conf');
data.mixf = getFilterSetup(params,data,'mixf');

if params.number_of_scales > 0
    scale_exp = (-floor((params.number_of_scales-1)/2):ceil((params.number_of_scales-1)/2));
    data.setup.scaleFactors  = params.scale_step .^ scale_exp;
    %force reasonable scale changes
    data.setup.min_scale_factor = params.scale_step ^ ceil(log(max(5 ./ data.obj.sz)) / log(params.scale_step));
    data.setup.max_scale_factor = params.scale_step ^ floor(log(min([size(data.seq.im,1) size(data.seq.im,2)] ./ data.obj.base_target_sz)) / log(params.scale_step));
end

if params.interpolate_response >= 3
    % Pre-computes the grid that is used for gradient ascent.
    data.obj.ky = circshift(-floor((data.obj.use_sz(1) - 1)/2) : ceil((data.obj.use_sz(1) - 1)/2), [1, -floor((data.obj.use_sz(1) - 1)/2)]);
    data.obj.kx = circshift(-floor((data.obj.use_sz(2) - 1)/2) : ceil((data.obj.use_sz(2) - 1)/2), [1, -floor((data.obj.use_sz(2) - 1)/2)])';
end

% initialize the projection matrix
data.obj.rects = zeros(data.seq.num_frames, 4);

data.setup.xxlf_sep = zeros(2*data.setup.feature_dim, length(data.setup.dft_sym_ind) + 2 * length(data.setup.dft_pos_ind), data.setup.feature_dim, 'single');
data.setup.multires_pixel_template = zeros(data.obj.sz(1), data.obj.sz(2), size(data.seq.im,3), params.number_of_scales, 'uint8');
end