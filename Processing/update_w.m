function [params,data, filter]=update_w(params,data,filter) %responsemap,enbel_planewin,reg_window0,reg_window,
     
     srcImg0 = data.seq.im;
     sz_0=size(srcImg0);  
     get_scal = params.sal_reg;
     
    % obtain the saliency detection region that is k-times of the bounding
    % box
    xs = floor(data.obj.pos(2)) + (1:get_scal*data.obj.target_sz(2)) - floor(get_scal*data.obj.target_sz(2)/2);
    ys = floor(data.obj.pos(1)) + (1:get_scal*data.obj.target_sz(1)) - floor(get_scal*data.obj.target_sz(1)/2);
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(data.seq.im,2)) = size(data.seq.im,2);
    ys(ys > size(data.seq.im,1)) = size(data.seq.im,1);
    
    srcImg = get_pixels(srcImg0,data.obj.pos,round(get_scal*data.obj.target_sz),get_scal*data.obj.target_sz); 
    
    %enlarge the saliencyregion if it is too small
    flag_lager=0;
    if min(size(srcImg))<128
        if(data.obj.target_sz(1) <= data.obj.target_sz(2))
        srcImg=imresize(srcImg,[128,128*double(get_scal*data.obj.target_sz(2))/double(get_scal*data.obj.target_sz(1))]);
        elseif(data.obj.target_sz(1)>data.obj.target_sz(2))
        srcImg=imresize(srcImg,[128*double(get_scal*data.obj.target_sz(1))/double(get_scal*data.obj.target_sz(2)),128]);
        end
        flag_lager=1;
    end
    
    % the saliency detection needs 3-channel image
    img_channels=ndims(srcImg);  
    if img_channels==2
        I =srcImg(:,:);
        clear I3;
        I3(:,:,1)=I;
        I3(:,:,2)=I;
        I3(:,:,3)=I;   
        srcImg=I3;
    end

    time = toc;
    saliency_img=get_saliency_SCA(srcImg,200);
    % other methods to get the saliency map
    %saliency_img=get_saliency_wCtr(srcImg);
    %saliency_img=get_saliency_RA(srcImg);
    data.time.timeS = data.time.timeS + toc - time;
   
    if flag_lager ==1
        saliency_img=imresize(saliency_img,floor(get_scal*data.obj.target_sz));
    end
   
    mask_0=zeros(sz_0(1),sz_0(2));        
    mask_0(ys,xs)=double(saliency_img(:,:)); 
   
    tar_xs = floor(data.obj.pos(2)) + (1:data.obj.target_sz(2)) - floor(data.obj.target_sz(2)/2);
    tar_ys = floor(data.obj.pos(1)) + (1:data.obj.target_sz(1)) - floor(data.obj.target_sz(1)/2);  
    tar_xs(tar_xs < 1) = 1;
    tar_ys(tar_ys < 1) = 1;
    tar_xs(tar_xs > size(data.seq.im,2)) = size(data.seq.im,2);
    tar_ys(tar_ys > size(data.seq.im,1)) = size(data.seq.im,1); 
    if params.SaliencyBS == 1
        target_slcimg = double(get_pixels(mask_0,data.obj.pos,round(data.obj.target_sz),data.obj.target_sz)); 
        mask_0=zeros(sz_0(1),sz_0(2));   
        % slc_thr=sum(sum(target_slcimg(:,:)))/(length(tar_ys)*length(tar_xs));    
        mask_0(tar_ys,tar_xs)=double(target_slcimg(:,:));    
    end
    maskpixels = double(get_pixels2(mask_0,data.obj.pos,round(data.obj.sz*data.obj.currentScaleFactor),data.obj.sz)); 
    
    %   savepath = (strcat('D:\tracking\OTB\tracker_benchmark_v1.0\trackers\DSARCF\case_CLE/saliency/singer2/sal_',num2str(data.seq.frame),'_BS.jpg'));
    %   imwrite(maskpixels,savepath);
    %   figure(2) 
    %   imshow(maskpixels)
    %   hold on;
    
    % compute the ratio of boundingbox in W
    ratio_tar(1) = length(tar_xs)/(data.obj.currentScaleFactor*data.obj.sz(1));
    ratio_tar(2) = length(tar_ys)/(data.obj.currentScaleFactor*data.obj.sz(2));
  
    % normlize the values of mask_w
    maskpixels=(max(maskpixels(:))-maskpixels)/(max(maskpixels(:))-min(maskpixels(:)));
            
    saliency=imresize(maskpixels,size(filter.reg_window_SR));
    filter.reg_window = energy_w(params,filter,saliency,ratio_tar);%reg_window,,responsemap,reg_window0,,enbel_planewin,);

    % compute the DFT and enforce sparsity
    reg_window_dft = fft2(filter.reg_window) / prod(data.obj.use_sz);
    reg_window_dft_sep = cat(3, real(reg_window_dft), imag(reg_window_dft));
    reg_window_dft_sep(abs(reg_window_dft_sep) < params.reg_sparsity_threshold * max(abs(reg_window_dft_sep(:)))) = 0;
    reg_window_dft = reg_window_dft_sep(:,:,1) + 1i*reg_window_dft_sep(:,:,2);

    % do the inverse transform, correct window minimum
    reg_window_sparse = real(ifft2(reg_window_dft));
    reg_window_dft(1,1) = reg_window_dft(1,1) - data.obj.support_sz * min(reg_window_sparse(:)) + params.reg_window_min;

    % construct the regularizsation matrix
    regW = cconvmtx2(reg_window_dft);

    regW_dfs = real(data.setup.dfs_matrix * regW * data.setup.dfs_matrix');

    WW_block = regW_dfs' * regW_dfs;

    % If the filter size is small enough, remove small values in WW_block.
    % It takes too long time otherwise.
    if data.obj.support_sz <= 120^2
        WW_block(0<abs(WW_block) & abs(WW_block)<0.00001) = 0;
    end

    % create block diagonal regularization matrix
    WW = eval(['blkdiag(WW_block' repmat(',WW_block', 1, data.setup.feature_dim-1) ');']);

    % upper and lower triangular parts of the regularization matrix
    filter.WW_U = triu(WW, 1); 

end