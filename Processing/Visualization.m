function Visualization(bSaveImage, slector, frame, im, pos, target_sz)
% figure(1)
if bSaveImage
    rect_position_vis = [pos([2,1]) - (target_sz([2,1]) - 1)/2, target_sz([2,1])];
    im_to_show = double(im)/255;
    if size(im_to_show,3) == 1
        im_to_show = repmat(im_to_show, [1 1 3]);
    end
    if frame == 1
        imagesc(im_to_show);
        hold on;
        rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
        text(10, 10, int2str(frame), 'color', [0 1 1]);
        hold off;
        axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
    else
        imagesc(im_to_show);
        hold on;
        if slector == 1
            edgecolor = 'g';           
        elseif slector == 0
            edgecolor = 'b';
        else
            edgecolor = 'r';
        end
        rectangle('Position',rect_position_vis, 'EdgeColor',edgecolor, 'LineWidth',2);
        text(10, 10, int2str(frame), 'color', [0 1 1]);
        hold off;
        axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
    end
    drawnow
end

end