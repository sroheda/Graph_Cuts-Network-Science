clc
clear all
close all

%% initialize
Im = imread('len_top.jpg');
[nRow,nCol,dim] = size(Im);              
N  = nRow*nCol;                          
Vx = reshape(Im,N,dim);                  
% figure(); imshow(Im);
%% Texture features
m    = 1.5;                  
ST   = 10; %10               
dim2 = 4;                    
In = im2double(rgb2gray(Im));
T = zeros(nRow,nCol,dim2);   
disp('nRow:'); disp(nCol)
for i = 1:nRow               
    for j = 1:nCol   
Vi  = (i-floor(m)):(i+floor(m));  
Vj  = ((j-floor(m)):(j+floor(m)));
Vi  = Vi(Vi>=1 & Vi<=nRow);           
Vj  = Vj(Vj>=1 & Vj<=nCol);
blk = In(Vi,Vj);                      
T(i,j,1) = mean(blk(:));              
T(i,j,2) = var(blk(:));               
T(i,j,3) = skewness(blk(:));          
T(i,j,4) = kurtosis(blk(:));          
    end
end
T(:,:,1) = (T(:,:,1)-min(min(T(:,:,1))))/(max(max(T(:,:,1)))-min(min(T(:,:,1))));
T(:,:,2) = (T(:,:,2)-min(min(T(:,:,2))))/(max(max(T(:,:,2)))-min(min(T(:,:,2))));
T(:,:,3) = (T(:,:,3)-min(min(T(:,:,3))))/(max(max(T(:,:,3)))-min(min(T(:,:,3))));
T(:,:,4) = (T(:,:,4)-min(min(T(:,:,4))))/(max(max(T(:,:,4)))-min(min(T(:,:,4))));
T = uint8(255*T);                        % normalization
% figure(); imshow(T)
%% Compute weight matrix W
r  = 1.5;                                
SI = 5; %5                               
SX = 6; %6                               
sNcut = 2.5; %0.22                       
sArea = 600; %30                         
W = sparse(N,N);                         
F = reshape(Im,N,1,dim);                 
G = reshape(T,N,1,dim2);                 
X = cat(2,repmat((1:nRow)',1,nCol),repmat((1:nCol),nRow,1));
X = reshape(X,N,1,2);                    
for Ic = 1:nCol                          
    if ((mod(Ic,100) == 0) || Ic == 1)
        disp(Ic)
    end
    for Ir = 1:nRow    
Vc = (Ic-floor(r)):(Ic+floor(r));        
Vr = ((Ir-floor(r)):(Ir+floor(r)))';
Vc = Vc(Vc>=1 & Vc<=nCol);               
Vr = Vr(Vr>=1 & Vr<=nRow);
VN = length(Vc)*length(Vr);
I = Ir +(Ic-1)*nRow;                    
V = repmat(Vr,1,length(Vc))+repmat((Vc-1)*nRow,length(Vr),1);
V = reshape(V,length(Vc)*length(Vr),1); 
XV = X(V,1,:);                          
XI = repmat(X(I,1,:),length(V),1);
DX = XI-XV;                             
DX = sum(DX.*DX,3);                     
constraint = find(sqrt(DX)<=r);         
V = V(constraint);
DX = DX(constraint);
FV = F(V,1,:);                          
FI = repmat(F(I,1,:),length(V),1);      
DF = FI-FV;                             
DF = sum(DF.*DF,3);                     
GV = G(V,1,:);                          
GI = repmat(G(I,1,:),length(V),1);      
DG = GI-GV;                             
DG = sum(DG.*DG,3);                     
W(I, V) = exp(-DF/(SI*SI)) .* exp(-DX/(SX*SX)) .* exp(-DG/(ST*ST));
    end
end
tic
%% NcutPartition
Seg = (1:N)';         
id = 'ROOT';          
d = sum(W, 2);        
D = spdiags(d, 0, N, N);      
disp('Computed D')
warning off;                  
[U,S] = eigs(D-W, D, 2, 'sm');  
U2 = U(:, 2);                    
t = mean(U2);                    
t = fminsearch('NcutValue', t, [], U2, W, D);
disp('Found Optimal Threshold')
A = find(U2 > t);
B = find(U2 <= t);

x = (U2 > t);
x = (2 * x) - 1;
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
ncut = (y' * (D - W) * y) / ( y' * D * y );
disp('ncut done')
%% recursive 2-way repartition
if (length(A) < sArea || length(B) < sArea) 
    Seg{1}   = Seg;
    %Id{1}   = id;   
    %Ncut{1} = ncut; 
    return;
end
[SegA IdA NcutA] = NcutPartition(Seg(A), W(A, A), sNcut, sArea, [id '-A']); 
[SegB IdB NcutB] = NcutPartition(Seg(B), W(B, B), sNcut, sArea, [id '-B']); 
Seg  = [SegA SegB];                      
Id   = [IdA IdB];
Ncut = [NcutA NcutB];
%% show
Io   = zeros(size(Im),'uint8');
for k=1:length(Seg)
 [r, c] = ind2sub(size(Im),Seg{k});
 for i=1:length(r)
 Io(r(i),c(i),1:dim) = uint8(round(mean(Vx(Seg{k}, :))));
 end
end
figure()
subplot(121); imshow(Im); title('Original'); 
subplot(122); imshow(Io);  title(['ncut',': ',num2str(length(Seg))]);
toc
