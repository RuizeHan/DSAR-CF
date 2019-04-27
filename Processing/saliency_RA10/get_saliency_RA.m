
function saliency_img = get_saliency_RA(img)

[L,a,b]=RGB2Lab(img);
[salMat,salMatInd]=saliencyMeasure({L,a,b});
saliency_img = salMat;
% imwrite(salMat,'saliency.png');

end

          