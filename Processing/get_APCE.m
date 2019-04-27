function [max_res,APCE] = get_APCE(C)

C_max=max(max(C));
C_min=min(min(C));
max_res = C_max;
APCE=(C_max-C_min)^2/mean((C(:)-C_min).^2);

end