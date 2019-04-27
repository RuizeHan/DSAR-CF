function s=getSelector(params, data, ss)

    lamda = params.s_lamda;

    ref_window_power = params.reg_window_power;
        
    % normalization factor
    reg_scale = 0.5 * data.obj.base_target_sz/data.setup.featureRatio;

    % construct grid
    wrg = -(data.obj.use_sz(1)-1)/2:(data.obj.use_sz(1)-1)/2;
    wcg = -(data.obj.use_sz(2)-1)/2:(data.obj.use_sz(2)-1)/2;
    [wrs, wcs] = ndgrid(wrg, wcg);

    % construct the regukarization window
    reg_window = (params.reg_window_edge - params.reg_window_min) * (abs(wrs/reg_scale(1)).^ref_window_power + abs(wcs/reg_scale(2)).^ref_window_power) + params.reg_window_min;
    
    W_o = reg_window;
    maxmax = max(max(reg_window));
    minmin = min(min(reg_window));
    mmdiff = maxmax + minmin;
    W_c(:,:) = mmdiff - reg_window(:,:);
    
    switch params.selector   
    case -1 % context filter
          F = data.conf.hf;
    case 0 % mix filter
          F = data.mixf.hf;
    case 1 % object filter
          F = data.objf.hf;
    end
    
    f = real(ifft2(F));
    P = sum(sum(sum(W_c.*(W_o - W_c) .* f)));
    Q = sum(sum(sum((W_o - W_c).^2 .* f)));
    
    PP = P/(P+Q);
    QQ = Q/(P+Q);
    
    s= (lamda * ss - PP)./ (QQ + lamda) ;
    
end
