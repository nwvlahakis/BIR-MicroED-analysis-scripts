function [stack]=read_tvips(filename1,header_size,image_size,new_size,scalefactor)
%The header for .tvips files are 180 bytes.
%Image sizes are typically 2048 x 2048.
%New size represents the size of the region of interest would like to sample from
%the raw data - centered at the raw image center.
%scalefactor ranges from 0 to 1 - resizes images to be smaller.

fclose('all');



if exist('header_size','var')==0
    header_size=180;
end
if exist('image_size','var')==0
    image_size=2048;
end
if exist('new_size','var')==0
    new_size=image_size;
end
if exist('scalefactor','var')==0 || scalefactor>1
    scalefactor=1;
end

Y=image_size;
YC=floor((Y.*scalefactor+1)/2);
X=image_size;
XC=floor((X.*scalefactor+1)/2);
hnew_size=floor(new_size/2);

%tvips files from the XF416 detector are sometimes split into multiple parts if size exceeds a certain amount
% This chunk reads the second tvips file if it exists in the same directory
% as the first

[filepath,name,ext] = fileparts(filename1);

if isfile(strcat(filepath,'/',name(1:end-1),num2str(1),ext))
   filename2 = strcat(filepath,'/',name(1:end-1),num2str(1),ext);
   filename = strcat(filepath,'/',name(1:end-3),'cat',ext);
   command = ['cat ' filename1 ' ' filename2,' > ', filename];
   status = system(command);
else 
    filename = filename1;
end

fid=fopen(filename,'r');
ct=0;

while ~feof(fid) && ct<2000
    if ct==0 
       fseek(fid,256,'bof');
    end
    fseek(fid,header_size,'cof');
    imagein=single(fread(fid,[Y X],'uint16','ieee-le'));
    
    if size(imagein,1)== image_size && size(imagein,2)== image_size
        
        imagein=imresize(imagein,scalefactor);
        imagein=imagein(YC-hnew_size+1:YC+hnew_size,XC-hnew_size+1:XC+hnew_size);
        
        if ct==0
            stack= imagein; 
        else
   
            stack=cat(3,stack,imagein);
        end
        
        if mod(ct,10)==1
        	figure(11);
        	subplot(1,2,1), imagesc(imagein), axis image, colormap gray, caxis([0 100]), title(strcat('frame no. ',num2str(ct)));
        	subplot(1,2,2), histogram(imagein(imagein<60000),10), title(strcat('histogram of image ',num2str(ct))); drawnow();
            pause(0.1);
        end
        
        ct=ct+1;
    end
end


fclose(fid);

end