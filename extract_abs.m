function [flag,x,y,z] = extract_abs(fn)

fid = fopen(fn, 'r');

rows = textscan(fid,'%d%*[^\n]',1);
cols = textscan(fid,'%d%*[^\n]',1);

% rows = rows{1};
% cols = cols{1};
rows = 480;
cols = 640;

len = rows * cols;

textscan(fid,'%*[^\n]',1);

flag = textscan(fid,'%u8',len);
flag = reshape(flag{:},cols,rows)';
x = textscan(fid,'%f',len);
x = reshape(x{:},cols,rows)'; 
y = textscan(fid,'%f',len);
y = reshape(y{:},cols,rows)';
z = textscan(fid,'%f',len);
z = reshape(z{:},cols,rows)';
% why transpose? Because the data in file represents the face image from
% 1st to last line, but in matlab 'reshape' will shape the matrix from
% 1st col to last col.

fclose(fid);