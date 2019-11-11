function [ a ] = loadfield3d(fn)

% fn:filename; n:Ncell; l:depth

fid=fopen(fn,'r');
%fid=fopen(fn,'r','b'); % big endian
p1=fread(fid,'real*4'); % real*4
fclose(fid);
n=round(numel(p1)^(1/3))
a=reshape(p1,[n,n,n]);
%disp(n)
end