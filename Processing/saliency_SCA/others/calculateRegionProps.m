function regions = calculateRegionProps(sup_num,sulabel_im)
% function regions = calculateRegionProps(sup_num,sulabel_im) can  
% calculate the region pixel index for each superpixel.
%
% Input: sup_num - the number of superpixels
%        sublabel_im - the superpixel labels for all pixels
%
% Output: regions - the region pixel index for each superpixel

for r = 1:sup_num
	indxy = find(sulabel_im==r);
	[indx indy] = find(sulabel_im==r);
	regions{r}.pixelInd = indxy;
    regions{r}.pixelIndxy = [indx indy];
end