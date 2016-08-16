function  [] = WriteComplex(cmplx, template, prefix)
   ss=size(template.data);
   if 0 %half hermatian
        x_size = floor((ss(1)/2+1));
        reshapecmplx = reshape(cmplx, ss);
        halfcmplx= reshapecmplx(1:x_size,:,:);
        WriteFile(single(real(halfcmplx)),template,[prefix '_real.nii']);
        WriteFile(single(imag(halfcmplx)),template,[prefix '_imag.nii']);
   else
        reshapecmplx = reshape(cmplx, ss);
        WriteFile(single(real(reshapecmplx)),template,[prefix '_real.nii']);
        WriteFile(single(imag(reshapecmplx)),template,[prefix '_imag.nii']);
   end

end
