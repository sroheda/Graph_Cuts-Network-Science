function [Seg Id Ncut] = NcutPartition(I, W, sNcut, sArea, id)
%% NcutPartition - Partitioning
N = length(W);
d = sum(W, 2);
D = spdiags(d, 0, N, N);    % diagonal matrix
warning off;                
[U,S] = eigs(D-W, D, 2, 'sm');
U2 = U(:, 2);

t = mean(U2);
t = fminsearch('NcutValue', t, [], U2, W, D);
A = find(U2 > t);
B = find(U2 <= t);

ncut = NcutValue(t, U2, W, D);
if (length(A) < sArea || length(B) < sArea) || ncut > sNcut
    Seg{1}   = I;
%    Id{1}   = id;          
%    Ncut{1} = ncut;        
    return;
end

[SegA IdA NcutA] = NcutPartition(I(A), W(A, A), sNcut, sArea, [id '-A']);

[SegB IdB NcutB] = NcutPartition(I(B), W(B, B), sNcut, sArea, [id '-B']);

Seg   = [SegA SegB];
Id   = [IdA IdB];
Ncut = [NcutA NcutB];
end
