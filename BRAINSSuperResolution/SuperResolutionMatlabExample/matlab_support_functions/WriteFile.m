function  [] = WriteFile(data, template, filename)
   %return
   template.data = single(real(data));
   itkSaveWithMetadata(filename,template); 
end
