function [] = diffmontage2(pathout,rootname,imgarray,sampleframe,diffarray,montagesize)

%produces imaging and diffraction montage for paired image/diffraction
%image stacks collected in alternating fashion on the same target (such as
%a crystal exhibiting beam-induced bend contour propagation).


imgstack = zeros(size(imgarray,1),size(imgarray,2),montagesize);

for aa = 1:montagesize
    imgstack(:,:,aa) = imfuse(imgarray(:,:,(aa-1)*10+1),sum(imgarray,3),'blend');
end


diffstack = zeros(size(diffarray,1),size(diffarray,2),montagesize);

for bb = 1:montagesize
    diffstack(:,:,bb) = diffarray(:,:,bb);
end


f9 = figure(009);
ax(1) = subplot(2,1,1);
montage(imgstack,'Size',[1 5])
clim(ax(1),'auto');
ax(2) = subplot(2,1,2);
montage(diffstack,'Size',[1 5])
clim(ax(2),[0 2]);


exportgraphics(f9,strcat(pathout,'/',rootname,'_kfiltmontage.pdf'),'ContentType','vector');
saveas(f9, strcat(pathout,'/',rootname,'_kfiltmontage.png'));
end