function ind = get_lowpass_inds(k,siz)
%GET_LOWPASS_INDS return center k-space indices in rectangule with
%dimensions siz = height x width
%input:    k = 2xprod(res) linear array of kspace indicies
%               (as obtained from get_kspace_inds)
%        siz = dimensions of center k-space indicies to return
%output: ind = linear index set cooresponding to center kspace rectangle

% NOTE: siz(2) = Row, siz(1) = Column, siz(3) = Slice 
%       (i.e. fortran indexing)

%                                         U  L
% if res(2) = 7, cf(2) = 3 ->  [ 0  1  2  3 -3 -2 -1 ]
% if res(2) = 6, cf(2) = 3 ->  [ 0  1  2    -3 -2 -1 ]
%                                      U     L
cf = floor(siz/2); %CenterFrequency

kxL = -floor(cf(2));
kyL = -floor(cf(1));
kzL = -floor(cf(3));

kxR = floor((cf(2)-mod(siz(2)+1,2)));
kyR = floor((cf(1)-mod(siz(1)+1,2)));
kzR = floor((cf(3)-mod(siz(3)+1,2)));

% NOTE: linear indexed memory is organized in C order Column, Row, Slice
ind = find((kxL <= k(1,:)) & (k(1,:) <= kxR) & (kyL <= k(2,:)) & (k(2,:) <= kyR) & (kzL <= k(3,:)) & (k(3,:) <= kzR));
ind = int64(ind); % Index locations should be integers
end
