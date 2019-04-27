 
function res=get_saliency_SCA(img,pixels)

%     imgpath=fullfile(path_img,strcat('img (',num2str(i),').jpg'));   
%     img=imread(imgpath);
    im_size=size(img);
    mask=zeros(im_size(1),im_size(2));
    box_l=round(im_size(1)/4);
    box_t=round(im_size(2)/4);
    mask(box_l:box_l*3,box_t:box_t*3)=1;
    sal=mask;
    [superpixel, ~] = slicmex(img,pixels,5);
    superpixel=superpixel+1;
    res=SCAfun(img, superpixel, sal);
%     savepath=fullfile('results',strcat('res(',num2str(i),').jpg')); 
%     imwrite(res,savepath);
%     level = graythresh(res);  
%     bw_slcimg=im2bw(res,level);  
%     savepath1=fullfile('results',strcat('bw(',num2str(i),').jpg')); 
%     imwrite(bw_slcimg,savepath1);

