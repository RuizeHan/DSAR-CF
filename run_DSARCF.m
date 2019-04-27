function results = run_DSARCF(seq, ~, ~, ~)

% By RuizeHan
% Dynamic Saliency-Aware Regularization for Correlation Filter based Object Tracking
% In TIP 2019

close all
addpath(genpath('Processing/'));

seq.format = 'otb';
% set the parameters of the method 
params = SetParams(seq);

% prepare the data of the method 
[params, data] = PrepareData(params);

time = 0;
timeF = 0;
timeD = 0;

for frame = 1:data.seq.num_frames
    data.seq.frame = frame;
    data.seq.im = imread(params.s_frames{data.seq.frame});
    if size(data.seq.im,3) > 1 && data.seq.colorImage == false
        data.seq.im = data.seq.im(:,:,1);
    end

    tic();
    % Detection process of CF tracking
    [params, data] = Detection(params, data);
    timeD = timeD + toc();
    
    time = toc;
    
    if rem(frame+params.update-1,params.update) == 0
    % initial the correlation filter
    if params.selector == -1 && params.selector_last == 1
        data.conf = getFilterSetup(params,data,'conf');
    end
    
    % update the correlation filter
    [params, data] = FilterUpdate(params, data);
    end
    timeF = timeF + toc() -time;
    timeW = data.time.timeW;
    timeS = data.time.timeS;
    time = timeD + timeF;
    
    Visualization(params.visualization, params.selector, data.seq.frame, data.seq.im, data.obj.pos, data.obj.target_sz);
    
end

% get the running fps of DSARCF
fps = numel(params.s_frames) / time;

% compute the running time of each component
% disp(['fps: ' num2str(fps)])
results.time.timeF = timeF;
results.time.timeD = timeD;
results.time.time = time;
results.time.timeW = timeW;
results.time.timeS = timeS;
results.type = 'rect';
results.res = data.obj.rects; % each row is a rectangle
results.fps = fps;

end