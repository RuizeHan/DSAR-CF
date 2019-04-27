function [features, support_sz] = get_cnn( im, fparam, gparam )
% extract neural net features from the image
persistent current_net;
persistent net_;
persistent average_;
persistent output_layer;
persistent net_info;

if isempty(net_)
    net = load(['networks/' fparam.nn_name]);
    current_net = fparam.nn_name;
    net_.layers = net.layers(1:fparam.output_layer);
    average_ = single(net.normalization.averageImage);
    output_layer = fparam.output_layer;
    net_info = vl_simplenn_display(net_);
    if isfield(fparam, 'scale') && fparam.scale ~= 1
        average_ = imresize(average_, fparam.scale);
    end
else
    if ~strcmp(current_net,fparam.nn_name) || fparam.output_layer ~= output_layer
        net = load(['networks/' fparam.nn_name]);
        current_net = fparam.nn_name;
        net_.layers = net.layers(1:fparam.output_layer);
        average_ = single(net.normalization.averageImage);
        output_layer = fparam.output_layer;
        net_info = vl_simplenn_display(net_);
        if isfield(fparam, 'scale') && fparam.scale ~= 1
            average_ = imresize(average_, fparam.scale);
        end
    end;
end;

if size(im,3) == 1
    im = repmat(im, [1 1 3]);
end

%preprocess the image
im_ = imresize(single(im),[size(average_,1),size(average_,2)]) - average_;
% im_ = single(mexResize(im,[size(average_,1),size(average_,2)],'auto')) - average_;

result = vl_simplenn(net_,im_);

features = result(fparam.output_layer + 1).x;

if fparam.output_layer > 0
    stride = net_info.receptiveFieldStride(:,end)';
    support_sz = [size(features,1) size(features,2)] .* stride .* [size(im,1) size(im,2)] ./ [size(im_,1) size(im_,2)];
else
    support_sz = [size(im,1) size(im,2)];
end

