function data = RidgeRrgression(data, params)

% Get multi-resolution image
for scale_ind = 1:params.number_of_scales
    data.setup.multires_pixel_template(:,:,:,scale_ind) = ...
        get_pixels(data.seq.im, data.obj.pos, round(data.obj.sz*data.obj.currentScaleFactor*data.setup.scaleFactors(scale_ind)), data.obj.sz);
end

xt = bsxfun(@times,get_features(data.setup.multires_pixel_template,params.t_features,params.t_global),data.setup.cos_window);

xtf = fft2(xt);

responsef = permute(sum(bsxfun(@times, data.objf.hf, xtf), 3), [1 2 4 3]);

if params.interpolate_response == 2
    % use dynamic interp size
    data.interp_sz = floor(size(data.y) * data.setup.featureRatio * data.obj.currentScaleFactor);
end
responsef_padded = resizeDFT2(responsef, data.obj.interp_sz);

% response
response = ifft2(responsef_padded, 'symmetric');

resp = fftshift(response);
max_response = max(resp(:));
[~,~,k] = size(response);
max_layer = 1;
for layer = 1 : k
    if(max(resp(:,:,layer)) == max_response)
        max_layer = layer;
    end
end

data.reg = jud_overlap([data.obj.pos([2,1]), data.obj.target_sz([2,1])],resp(:,:,max_layer));

end