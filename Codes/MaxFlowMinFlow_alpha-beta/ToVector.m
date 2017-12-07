function v = ToVector(im)
sz = size(im);
if size(im,3)==3
    v = reshape(im, [prod(sz(1:2)) 3]);
elseif size(im,3)== 2 
    v = reshape(im, [prod(sz(1:2)) 2]);
else
    v = reshape(im, [prod(sz(1:2)) 1]);
end