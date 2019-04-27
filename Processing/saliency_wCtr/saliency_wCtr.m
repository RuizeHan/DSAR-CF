clear all;
close all;
clc;


imagepath1='D:\tracking\Saliency_Detection\wCtr\Saliency_Dataset_Code\Running_dataset\DUT_OMRON_image\images\0001.jpg';
srcImg = imread(imagepath1);
saliency_img=get_saliency_wCtr(srcImg);
imshow(saliency_img);