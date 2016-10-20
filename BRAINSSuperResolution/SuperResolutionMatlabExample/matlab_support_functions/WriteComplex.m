function  [] = WriteComplex(cmplx, template, prefix)
   ss=size(template.data);
   if 0 %half hermatian
        x_size = floor((ss(1)/2+1));
        reshapecmplx = reshape(cmplx, ss);
        halfcmplx= reshapecmplx(1:x_size,:,:);
        WriteITKImageFromMatlabStructure(single(real(halfcmplx)),template,[prefix '_real.nii']);
        WriteITKImageFromMatlabStructure(single(imag(halfcmplx)),template,[prefix '_imag.nii']);
   else
        reshapecmplx = reshape(cmplx, ss);
        WriteITKImageFromMatlabStructure(single(real(reshapecmplx)),template,[prefix '_real.nii']);
        WriteITKImageFromMatlabStructure(single(imag(reshapecmplx)),template,[prefix '_imag.nii']);
   end

end
