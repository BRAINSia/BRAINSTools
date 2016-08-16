function k = get_kspace_inds_new( res )
%GET_KSPACE_INDS define k-space index sets based on fft indexing
%input:  res = fourier domain resoultion
%outoput:   k = 2xprod(res) linear array of kspace indicies

%NOTE: This uses a convention that is compatible with DFT's that
%      only store 1/2 of the complex numbers in a r2c_dft scheme
%NOTE: This scheme is shifted by one form the ifftshift functions

% if res(2) = 7, cf(2) = 3 ->  [ 0  1  2  3 -3 -2 -1 ]
% if res(2) = 6, cf(2) = 3 ->  [ 0  1  2  3    -2 -1 ]
%                                            ^ remove one negative element                                        
cf =   int64( floor(res/2)  ); %CenterFrequency
indx = int64( [ 0:(cf(2) ), (-cf(2) + mod(res(2)+1,2)):-1] );
indy = int64( [ 0:(cf(1) ), (-cf(1) + mod(res(1)+1,2)):-1] );
indz = int64( [ 0:(cf(3) ), (-cf(3) + mod(res(3)+1,2)):-1] );

[kx,ky,kz] = meshgrid(indx,indy,indz);
k = [kx(:)'; ky(:)'; kz(:)'];
k = int64(k);
end
