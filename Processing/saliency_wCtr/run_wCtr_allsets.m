clear all;
close all;
clc;

addpath('F:\Saliency_Dataset_Code\Running_code\HS\hsaliencyexe');

setname={'DUT_OMRON_image'};%;'ECSSD';'MSRA1000';'MSRA5000';'PASCAL_S'
savepath1='D:\tracking\saliency\wCtr\Saliency_Dataset_Code\Running_result\wCtr\';
savepath1_SF='D:\tracking\saliency\wCtr\Saliency_Dataset_Code\Running_result\SF\';
savepath1_GS='D:\tracking\saliency\wCtr\Saliency_Dataset_Code\Running_result\GS\';
imagepath1='D:\tracking\saliency\wCtr\Saliency_Dataset_Code\Running_dataset\';
imagetype='*.jpg';

for i=1:length(setname)
    imagepath=[imagepath1 setname{i} '\' 'images'];
    images=dir([imagepath imagetype]);
    number_images=length(images);
    savepath=[savepath1 setname{i} ];
    if(~exist(savepath,'file'))
        mkdir(savepath);
    end
    savepath_SF=[savepath1_SF setname{i}];
    savepath_GS=[savepath1_GS setname{i}];
    if(~exist(savepath_SF,'file'))
        mkdir(savepath_SF);
    end
    if(~exist(savepath_GS,'file'))
        mkdir(savepath_GS);
    end
    %imagepath='F:\Saliency_Dataset&Code\Running_code\HS\hsaliencyexe\src\';
    param_str1=[imagepath];
    param_str2=[savepath];
    demo(param_str1,param_str2,savepath_SF,savepath_GS);
%     dos(['HSaliency.exe' ' "' param_str1 '" "' param_str2 '\"']);
%     for j=1:number_images
%         param_str1=[imagepath images(j).name];
%         param_str2=[savepath '\'];
%         dos(['HSaliency.exe ' param_str1 ' ' param_str2]);
%         %eval(['!HSaliency ' imagepath images(j).name ' ' savepath '\']);
%     end
end
