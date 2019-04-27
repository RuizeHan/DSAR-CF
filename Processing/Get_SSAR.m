function  reg_window = Get_SSAR(reg_window, params, data)
%**************************************   saliency   ****************S***********************

    srcImg0 = data.seq.im;
    sz_0=size(srcImg0);
    get_scal = params.sal_reg;
    xs = floor(data.obj.pos(2)) + floor(1:get_scal*data.obj.target_sz(2)) - floor(get_scal*data.obj.target_sz(2)/2);
    ys = floor(data.obj.pos(1)) + floor(1:get_scal*data.obj.target_sz(1)) - floor(get_scal*data.obj.target_sz(1)/2);  
       
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(data.seq.im,2)) = size(data.seq.im,2);
    ys(ys > size(data.seq.im,1)) = size(data.seq.im,1);    
    srcImg = get_pixels(srcImg0,data.obj.pos,round(get_scal*data.obj.target_sz),get_scal*data.obj.target_sz); 
    
    %enlarge the saliencyregion if it is too small
    flag_lager=0;
    if min(size(srcImg))<128
        if(data.obj.target_sz(1)<data.obj.target_sz(2))
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
    I3(:,:,1)=I;
    I3(:,:,2)=I;
    I3(:,:,3)=I;   
    srcImg=I3;
    end
           
    %   path_img='D:\tracking\OTB\tracker_benchmark_v1.0\trackers\SRDCF_saliency\saliency_resSCA';
    %   srcName=fullfile(path_img,strcat(params.seqname,'src.jpg'));
    %   imwrite(srcImg,srcName);
    %   figure(2);imshow(srcImg);
    tic
    saliency_img=get_saliency_SCA(srcImg,200);
    data.time.timeS = data.time.timeS + toc;    
    %   saliency_img=get_saliency_wCtr(srcImg);
    %   saliency_img=get_saliency_RA(srcImg);
    if flag_lager ==1
        saliency_img=imresize(saliency_img,floor(get_scal*data.obj.target_sz));
    end
        
    mask_0=zeros(sz_0(1),sz_0(2));   
    % figure(2);imshow(saliency_img);
    mask_0(ys,xs)=double(saliency_img(:,:)); 
    maskpixels = double(get_pixels(mask_0,data.obj.pos,round(data.obj.sz*data.obj.currentScaleFactor),data.obj.sz));     
    % figure(1) 
    % imshow(maskpixels);
    if params.SaliencyBS == 1
        target_slcimg = double(get_pixels(mask_0,data.obj.pos,round(data.obj.target_sz),data.obj.target_sz)); 
        tar_xs = floor(data.obj.pos(2)) + (1:data.obj.target_sz(2)) - floor(data.obj.target_sz(2)/2);
        tar_ys = floor(data.obj.pos(1)) + (1:data.obj.target_sz(1)) - floor(data.obj.target_sz(1)/2);  
        mask_BS=zeros(sz_0(1),sz_0(2)); 
        mask_BS(tar_ys,tar_xs)=double(target_slcimg(:,:));
        maskpixels= double(get_pixels(mask_BS,data.obj.pos,round(data.obj.sz*data.obj.currentScaleFactor),data.obj.sz));     
    % figure(2) 
    % imshow(maskpixels);
    end

    maskpixels((maskpixels(:)<0))=0;
    maskpixels=(max(maskpixels(:))-maskpixels)/(max(maskpixels(:))-min(maskpixels(:)));
% **************************************   saliency   **************E***********************

    mask_w=imresize(maskpixels,size(reg_window));
    reg_window= reg_window.*double(mask_w); 
    % reg_window_S = double(mask_w);
    % reg_window_Snorm = normalizing(double(mask_w),min(reg_window(:)),max(reg_window(:)));  
    % reg_window = reg_window_S;
     
end