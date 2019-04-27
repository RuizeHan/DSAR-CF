function results=default(seq, res_path, bSaveImage, parameters)

close all

s_frames = seq.s_frames;

%translation filter feature parameters
hog_params.nDim = 31;

grayscale_params.colorspace='gray';
grayscale_params.nDim = 1;

cn_params.tablename = 'CNnorm';
cn_params.nDim = 10;
cn_params.useForGray = false;

params.t_features = {
    struct('getFeature',@get_colorspace, 'fparams',grayscale_params),...
    struct('getFeature',@get_fhog,'fparams',hog_params),...
    struct('getFeature',@get_table_feature, 'fparams',cn_params),...
};
params.t_global.cell_size = 4; %global feature-cell size
params.t_global.cell_selection_thresh = 0.75^2;

% Translation filter parameters
params.search_area_shape = 'square';
params.search_area_scale = 3.5;         % the scaling of the target size to get the search area
params.output_sigma_factor = 1/16;		% spatial bandwidth (proportional to target)
params.lambda = 1e-2;					% the minimum value of the regularization window
params.learning_rate = 0.025;			% learning rate
params.refinement_iterations = 1;       % number of iterations used to refine the resulting position in a frame
params.translation_model_max_area = 50^2;   % the maximal search area in number of feature grid cells
params.interpolate_response = 4;        % correlation score interpolation strategy
params.gradient_ascent_iterations = 5;
params.gradient_step_size = 0.00;
params.init_strategy = 'indep';         % strategy for initializing the filter, 'const_reg' or 'indep'
params.num_GS_iter = 4;

% Regularization window parameters
params.use_reg_window = 1;              % wather to use windowed regularization or not
params.reg_window_min = 0.1;					% the minimum value of the regularization window
params.reg_window_edge = 12.0;
params.reg_window_power = 2;            % the degree of the polynomial to use (e.g. 2 is a quadratic window)
params.reg_sparsity_threshold = 0.05;   % a relative threshold of which DFT coefficients that should be set to zero

% Scale parameters
params.number_of_scales = 7;
params.scale_step = 1.01;

% Other parameters
params.visualization = 0;
params.debug = 0;


params.wsize = [seq.init_rect(1,4), seq.init_rect(1,3)];
params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);
params.s_frames = s_frames;

results = tracker(params);
