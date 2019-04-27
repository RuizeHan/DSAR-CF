function res=SCAfun(image, superpixel, sal)

%%--------------------------set parameters---------------------------%%
% SLIC algorithm parameters
spnumber=200;  %the number of superpixels
compactness=20;
% Single-layer Cellular Automata parameters
theta = 10;  % control the strength of similarity between neighbors
a=0.6;
b=0.2;     % a and b control the strength of coherence

%%-------------------------Read an input image-----------------------%%
input_im=image;
input_imlab = RGB2Lab(image);
[m,n,r]=size(input_imlab);

%%------------------------generate superpixels----------------------%%
sulabel=superpixel; % superpixel label matrix
supNum = max(sulabel(:));% the actual superpixel number

% read saliency maps of other algorithms as prior maps
input_im_prior=double(sal);
input_im_prior=mat2gray(input_im_prior);
sal_prior=input_im_prior; % remove the image frame

% compute the saliency of each superpixel in prior maps
S_prior=zeros(1,supNum);
for i=1:supNum
    S_prior(i)=mean(sal_prior(find(sulabel==i)));
end

% label the boundary superpixels
boundsp = extract_bg_sp(sulabel,m,n);

% compute the feature (mean color in lab color space)
% for each superpixel
input_vals=reshape(input_im, m*n, r);
rgb_vals=zeros(supNum,1,3);
inds=cell(supNum,1);
for i=1:supNum
    inds{i}=find(sulabel==i);
    rgb_vals(i,1,:)=mean(input_vals(inds{i},:),1);
end
lab_vals = colorspace('Lab<-', rgb_vals);
seg_vals=reshape(lab_vals,supNum,3);

% compute Impact Factor Matrix
impfactor=AdjcProcloop(sulabel,supNum);
edges=[];
for i=1:supNum
    indext=[];
    ind=find(impfactor(i,:)==1);
    for j=1:length(ind)
        indj=find(impfactor(ind(j),:)==1);
        indext=[indext,indj];
    end
    indext=[indext,ind];
    indext=indext((indext>i));
    indext=unique(indext);
    if(~isempty(indext))
        ed=ones(length(indext),2);
        ed(:,2)=double(i)*ed(:,2);
        ed(:,1)=indext;
        edges=[edges;ed];
    end
end
weights = makeweights(edges,seg_vals,10);
F = adjacency(double(edges), double(weights),double(supNum));

% calculate a row-normalized impact factor matrix
D_sam=sum(F,2);
D=diag(D_sam);
F_normal=inv(D)*F;   % the row-normalized impact factor matrix

% compute Coherence Matrix
C=a*normalization(1./max(F'),0)+b;
C_normal=diag(C);

%%-----------------Single-layer Cellular Automata---------------%%
S_prior=normalization(S_prior,0);
S_N1=S_prior';
diff = setdiff(1:supNum, boundsp);

% step1: decrease the saliency value of boundary superpixels
for lap=1:5
    S_N1(boundsp) = S_N1(boundsp) - 0.6;
    neg_Ind = find(S_N1 < 0);
    if numel(neg_Ind) > 0
        S_N1(neg_Ind) = 0.001;
    end
    S_N1=C_normal*S_N1+(1-C_normal).*diag(ones(1,supNum))*F_normal*S_N1;
    S_N1(diff)=normalization(S_N1(diff),0);
end

% step2: control the ratio of foreground larger than a threshold
for lap = 1:5
    S_N1(boundsp) = S_N1(boundsp) - 0.6;
    neg_Ind = find(S_N1 < 0);
    if numel(neg_Ind) > 0
        S_N1(neg_Ind) = 0.001;
    end
    most_sal_sup = find(S_N1 > 0.93);
    if numel(most_sal_sup) < 0.02*supNum
        sal_diff = setdiff(1:supNum, most_sal_sup);
        S_N1(sal_diff) = normalization(S_N1(sal_diff), 0);
    end
    S_N1=C_normal*S_N1+(1-C_normal).*diag(ones(1,supNum))*F_normal*S_N1;
    S_N1(diff)=normalization(S_N1(diff),0);
end

% step3: simply update the saliency map according to rules
for lap = 1:10
    S_N1 = C_normal*S_N1+(1-C_normal).*diag(ones(1,supNum))*F_normal*S_N1;
    S_N1 = normalization(S_N1, 0);
end

image_sam=zeros(m,n);
image_sam(:)=S_N1(sulabel(:));
res=image_sam;
