function  [] = WriteITKImageFromMatlabStructure(data, template, filename)
   %return
   template.data = single(real(data));
   itkSaveWithMetadata(filename,template); 
end
