

clear all; close all; 

% path = 'Images\3063.jpg';
% im = im2double(imread(path));
% alpha = 0.7;
% k = 2;

path = 'Images\100007.jpg';
im = im2double(imread(path));
%E = evalclusters(ToVector(im2double(im)),'kmeans','CalinskiHarabasz','klist',[2:50]);
%im = imnoise(im,'salt & pepper');
alpha = 0.95;
k = 44; %E.OptimalK;
sz = size(im);
Hc=ones(sz(1:2));
Vc=Hc;
i_ground = 0; % rank of the bakground for plotting, 0: the darkest; 
%k-1 the brightest; 99: nowhere

diff=10000;
an_energy=999999999;
iter=0;
iter_v=0;
energy_global_min=99999999;

distance = 'sqEuclidean'; % Feature space distance

% Initialization: cluster the data into k regions
tic,
disp('Start kmeans');
data = ToVector(im);
[idx c] = kmeans(data, k, 'distance', distance,'EmptyAction','drop','maxiter',100);
c=c(find(isfinite(c(:,1))),:);
k=size(c,1);
k_max=k;
kmean_time=toc,
k_means_result = reshape(idx, size(im(:,:,1)));

Dc = zeros([sz(1:2) k],'single');   
c,


tic
while iter < 5
    iter=iter+1;
    clear Dc
    clear K
    c;
    for ci=1:k
        K=kernel_RBF(im,c(ci,:));
        Dc(:,:,ci)=1-K;
    end   
    clear Sc
    clear K
    %% The smoothness term
    Sc = alpha*(ones(k) - eye(k)); 
    gch = GraphCut('open', Dc, Sc, Vc, Hc);
    [gch L] = GraphCut('swap',gch);
    [gch se de] = GraphCut('energy', gch);
    nv_energy=se+de;
    gch = GraphCut('close', gch);
 
    if (nv_energy<=energy_global_min)
        diff=abs(energy_global_min-nv_energy)/energy_global_min;
        energy_global_min=nv_energy;
        L_global_min=L;
        k_max=k;
        
        nv_energy;
        iter_v=0;
        % Calculate region Pl of label l
        if size(im, 3)==3 % Color image
        for l=0:k-1
            Pl=find(L==l);
            card=length(Pl);
            K1=kernel_RBF(im(Pl),c(l+1,1));K2=kernel_RBF(im(Pl),c(l+1,2));K3=kernel_RBF(im(Pl),c(l+1,3));
            smKI(1)=sum(im(Pl).*K1); smKI(2)=sum(im(Pl+prod(sz(1:2))).*K2); smKI(3)=sum(im(Pl+2*prod(sz(1:2))).*K3);
            smK1=sum(K1);smK2=sum(K2);smK3=sum(K3);
            if (card~=0)
                c(l+1,1)=smKI(1)/smK1;c(l+1,2)=smKI(2)/smK2;c(l+1,3)=smKI(3)/smK3;
            else
                c(l+1,1)=999999999;c(l+1,2)=999999999;c(l+1,3)=999999999;
            end
        end
        end
        
        if size(im, 1)==1 % Gray-level image
        for l=0:k-1
            Pl=find(L==l);
            card=length(Pl);
            K=kernel_RBF(im(Pl),c(l+1,1));
            smKI=sum(im(Pl).*K);
            smK=sum(K);
            if (card~=0)
                c(l+1,1)=smKI/smK;
            else
                c(l+1,1)=999999999;
            end
        end
        end
        
        
        c=c(find(c(:,1)~=999999999),:);
        c_global_min=c;
        k_global=length(c(:,1));
        k=k_global;

    else
        iter_v=iter_v+1;
        
        if size(im, 3)==3 % Color image 
        for l=0:k-1           
            Pl=find(L==l);
            card=length(Pl);
            K1=kernel_RBF(im(Pl),c(l+1,1));K2=kernel_RBF(im(Pl),c(l+1,2));K3=kernel_RBF(im(Pl),c(l+1,3));
            smKI(1)=sum(im(Pl).*K1); smKI(2)=sum(im(Pl+prod(sz(1:2))).*K2); smKI(3)=sum(im(Pl+2*prod(sz(1:2))).*K3);
            smK1=sum(K1);smK2=sum(K2);smK3=sum(K3);
            % Calculate contour Cl of region Pl
            if (card~=0)
                c(l+1,1)=smKI(1)/smK1;c(l+1,2)=smKI(2)/smK2;c(l+1,3)=smKI(3)/smK3;
            else
                c(l+1,1)=999999999;c(l+1,2)=999999999;c(l+1,3)=999999999;
                area(l+1)=999999999;
            end
        end
        end
        
        if size(im, 3)== 1 % Gray-level image 
        for l=0:k-1           
            Pl=find(L==l);
            card=length(Pl);
            K=kernel_RBF(im(Pl),c(l+1,1));
            smKI=sum(im(Pl).*K);
            smK=sum(K);
            % Calculate contour Cl of region Pl
            if (card~=0)
                c(l+1,1)=smKI/smK;
            else
                c(l+1,1)=999999999;
                area(l+1)=999999999;
            end
        end
        end
              
        c=c(find(c(:,1)~=999999999),:);
        k=length(c(:,1));
    end
end

L=L_global_min;
energy_global_min;
c=c_global_min;

size(c,1)
iter;



if size(im, 3)==3 % Color image 
img=zeros(sz(1),sz(2),3);
j=1;
imagesc(im); axis off; hold on; colormap winter

for i=0:k_max-1
    LL=(L_global_min==i);
    is_zero=sum(sum(LL));
    if is_zero
         img(:,:,1)=img(:,:,1)+LL*c(j,1);
         img(:,:,2)=img(:,:,2)+LL*c(j,2);
         img(:,:,3)=img(:,:,3)+LL*c(j,3);
         j=j+1;
    end
    if i~=i_ground
        color=[rand rand rand];
        contour(LL,[1 1],'LineWidth',2.5,'Color',color); hold on;
    end
end
figure(2);
imagesc(img); axis off; colormap winter
end

if size(im, 3)==1 % Gray-level image 
img=zeros(sz(1),sz(2));
j=1;
imagesc(im); axis off; hold on; colormap winter; 

for i=0:k_max-1
    LL=(L_global_min==i);
    is_zero=sum(sum(LL));
    if is_zero
         img(:,:,1)=img(:,:,1)+LL*c(j,1);
         j=j+1;
    end
    if i~=i_ground
        color=[rand rand rand];
        contour(LL,[1 1],'LineWidth',2.5,'Color',color); hold on;
    end
end
figure(2);

imagesc(img); axis off;
end
toc

