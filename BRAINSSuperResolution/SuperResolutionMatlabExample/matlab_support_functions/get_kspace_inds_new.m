function k = get_kspace_inds_new( res )
%GET_KSPACE_INDS define k-space index sets based on fft indexing
%input:  res = fourier domain resoultion
%outoput:   k = 2xprod(res) linear array of kspace indicies


% if res(2) = 7, cf(2) = 3 ->  [ 0  1  2  3 -3 -2 -1 ]
% if res(2) = 6, cf(2) = 3 ->  [ 0  1  2    -3 -2 -1 ]
%                                         ^ remove one positive element                                        
cf =   int64( floor(res/2)  ); %CenterFrequency
indx = int64( [ 0:(cf(2) - mod(res(2)+1,2)), (-cf(2)):-1] );
indy = int64( [ 0:(cf(1) - mod(res(1)+1,2)), (-cf(1)):-1] );
indz = int64( [ 0:(cf(3) - mod(res(3)+1,2)), (-cf(3)):-1] );

[kx,ky,kz] = meshgrid(indx,indy,indz);
k = [kx(:)'; ky(:)'; kz(:)'];
k = int64(k);
end