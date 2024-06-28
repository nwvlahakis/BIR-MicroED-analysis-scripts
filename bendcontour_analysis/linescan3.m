function[stackscans,xi,yi]=linescan3(stackin)

[Y,X,Z]=size(stackin);

stacksum=squeeze(sum(stackin,3));
figure(005); 
subplot(1,2,1), imagesc(stacksum), axis image; 
[xi,yi]=getpts();

c = improfile(stacksum,xi,yi); %size(c),
stackscans=zeros(max(size(c)),Z,'single');

for kk=1:Z
    c=improfile(stackin(:,:,kk),xi,yi);
    stackscans(:,kk)=c;
end

subplot(1,2,2), imagesc(stackscans);