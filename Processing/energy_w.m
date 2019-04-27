function W = energy_w(params,filter,saliency,ratio_tar)

W_SR = filter.reg_window_SR;
W = filter.reg_window;
S = saliency;

% learning rate parameters setting
params.deta1=0.02;
params.deta2=0.1;
lambda_1=1;
lambda_2=1;
lambda_3 = 1;
zeta =params.zeta;

% initialization for the level-set function
W_size=size(W);
bdbox_size(1) = round(W_size(1) * ratio_tar(1));
bdbox_size(2) = round(W_size(2) * ratio_tar(2));
xs = floor(W_size(1)/2) + (1:bdbox_size(1)) - floor(bdbox_size(1)/2);
ys = floor(W_size(2)/2) + (1:bdbox_size(2)) - floor(bdbox_size(2)/2);
WW=W;WW(ys,xs,:)=-1;
omega = find(WW==-1);
% no_omega = find(WW~=-1);
% zeta = mean(W(omega));
omega_in = omega((W(omega)<zeta));  
omega_out = omega((W(omega)>zeta)); 
mu_in = mean(saliency(omega_in)); % interior mean
mu_out = mean(saliency(omega_out)); % exterior mean

WEF = W_EneryFun(W,mu_in,mu_out);
max_its = 10;
its=1;
NewW = W;

% update the level-set function by iterative method
while ((its < max_its) )  %&& ~stop             
   %-- Minimized energy from the BSpline coefficients
   %-- Compute new feature parameters   
    new_mu_in = sum(saliency(omega) .* heavyside(NewW(omega)) ) / sum( heavyside(NewW(omega)) );
    new_mu_out = sum(saliency(omega) .* ( 1 - heavyside(NewW(omega))) ) / sum( 1 - heavyside(NewW(omega)) );
    dW_dt1 = dirac(NewW(omega)).*(-lambda_1*(S(omega)-mu_in).^2+lambda_2*(S(omega)-mu_out).^2);
    %dW_dt1_rect=reshape(dW_dt1,[bdbox_size(2),bdbox_size(1)]);  
    %dW_dt2=-2*lambda_3*(NewW(no_omega)-W_S(no_omega))-2*lambda_4*(NewW(no_omega)-W_C(no_omega));    
    dW_dt2 =-2*lambda_3*(NewW-W_SR);
    NewW(omega) = NewW(omega) + params.deta1 * dW_dt1;
    NewW= NewW + params.deta2 * dW_dt2;
    New_WEF = W_EneryFun(NewW,new_mu_in,new_mu_out);   

    %-- Update feature parameters
    if ( New_WEF < WEF )
    WEF = New_WEF; 
    W = NewW;  
    end 
    its=its+1;
end

%-- Compute the regularized heaviside function
function H = heavyside(w)
    epsilon = params.epsilon;
    H = 1-0.5 * ( 1 + (2/pi) * atan((w-zeta)/epsilon) );
end

%-- Compute the regularized dirac function
function d = dirac(w)
    epsilon = params.epsilon;
    d = -(1/(pi*epsilon)) ./ ( 1 + ((w-zeta)/epsilon).^2 );
end

function WEF = W_EneryFun(W,mu_in,mu_out) 
    W1=sum((S(:)-mu_in).^2 .* heavyside(W(:))); 
    W2=sum((S(:)-mu_out).^2 .*(1- heavyside(W(:)))); 
    W3=sum(sum((W-W_SR).^2)); 
    WEF = lambda_1* W1 + lambda_2* W2 + lambda_3* W3 ;
end
   
end
