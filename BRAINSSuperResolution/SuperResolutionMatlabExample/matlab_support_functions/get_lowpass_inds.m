function ind = get_lowpass_inds(k,siz)
%GET_LOWPASS_INDS return center k-space indices in rectangule with
%dimensions siz = height x width
%input:    k = 2xprod(res) linear array of kspace indicies
%               (as obtained from get_kspace_inds)
%        siz = dimensions of center k-space indicies to return
%output: ind = linear index set cooresponding to center kspace rectangle
if mod(siz(2),2)
    kxL = (siz(2)-1)/2;
else
    kxL = siz(2)/2;
end
if mod(siz(1),2)
    kyL = (siz(1)-1)/2;
else
    kyL = siz(1)/2;
end
if mod(siz(3),2)
    kzL = (siz(3)-1)/2;
else
    kzL = siz(3)/2;
end
    kxR = floor((siz(2)-1)/2);
    kyR = floor((siz(1)-1)/2);
    kzR = floor((siz(3)-1)/2);

ind = find((-kxL <= k(1,:)) & (k(1,:) <= kxR) & (-kyL <= k(2,:)) & (k(2,:) <= kyR) & (-kzL <= k(3,:)) & (k(3,:) <= kzR));
end
