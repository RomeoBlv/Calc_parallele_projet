clear all
close all
clc

%%%%%%%%%%%%%Dimensions 180x150%%%%%%%%%%%%%%%%%

I=imread('barbara.gif');
fileID = fopen('image_barbara.bin','w');
fwrite(fileID,I','double');
fclose(fileID);


fid=fopen('image_barbara.bin','r','l');
mydata=fread(fid,'double');
for(i=1:512)
    for(j=1:512)
        A(i,j)=mydata((i-1)*512+j);
    end
end
figure;imshow(A,[]);



