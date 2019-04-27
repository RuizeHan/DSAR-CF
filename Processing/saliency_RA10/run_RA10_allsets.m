clear all;
close all;
clc;

%addpath('F:\Saliency_Dataset_Code\Running_code\HS\hsaliencyexe');

setname={'ECSSD'};
savepath1='C:\Users\scm\Desktop\tietu\runningtime\res\';
imagepath1='C:\Users\scm\Desktop\tietu\runningtime\';
imagetype='*.jpg';
inputimagetype='.jpg';

%tempfile='F:\Saliency_Dataset_Code\Running_code_2\DSR\temp_bmp\';
tic;
for k = 1:10
for i=1:length(setname)
    imagepath=[imagepath1 setname{i} '\' 'images\'];
    images=dir([imagepath imagetype]);
    number_images=length(images);
    savepath=[savepath1 setname{i} '\'];
    if(~exist(savepath,'file'))
        mkdir(savepath);
    end
    
    %imagepath='F:\Saliency_Dataset&Code\Running_code\HS\hsaliencyexe\src\';
    
    for j=1:number_images
        image_name=[imagepath images(j).name];
        savename=[savepath images(j).name(1:length(images(j).name)-4) '_RA10.png'];
        
        if(~exist(savename,'file'))
            try
                img=imread(image_name);
                [L,a,b]=RGB2Lab(img);
                [salMat,salMatInd]=saliencyMeasure({L,a,b});
                imwrite(salMat,savename,'png');
                disp(image_name);
            catch
                errorfile=[savepath 'errorfile\'];
                if(~exist(errorfile,'file'))
                    mkdir(errorfile);
                end
                copyfile(image_name,[errorfile images(j).name]);
            end
        end
    end
%     param_str1=[imagepath imagetype];
%     param_str2=[savepath];
%     dos(['HSaliency.exe' ' "' param_str1 '" "' param_str2 '\"']);
%     for j=1:number_images
%         param_str1=[imagepath images(j).name];
%         param_str2=[savepath '\'];
%         dos(['HSaliency.exe ' param_str1 ' ' param_str2]);
%         %eval(['!HSaliency ' imagepath images(j).name ' ' savepath '\']);
%     end
end
end
toc;
