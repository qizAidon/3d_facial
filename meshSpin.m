function [spinImage, vertices3d, dsc_lmk, shape_lmk] = meshSpin( pts_file, abs_file, L, resolution, E, index )
%MESHSPIN ------calculating spin image for a mesh
% Input:
%   pts_file --- read landmark 2-D coordinates, if empty, input '' instead
%   abs_file --- read 3-D coordinates and flag of points on a mesh
%          L --- number of landmarks to extract from pts_file
% resolution --- parameter for calculate spin image descriptor
%          E --- parameter for calculate spin image descriptor, it is the
%          number of edges taken into account when calculting bin size
% Output:
%  spinImage --- spin image descriptor for each valid point on the mesh
% vertices3d --- 3-D coordinates for each valid point on the mesh
%    dsc_lmk --- returns the descriptor of landmarks when pts_file ~empty
%  shape_lmk --- returnd the 3-D coordinates of landmarks when pts_file ~empty
%      index --- indicate the No. of the current file

[flag,x,y,z]=extract_abs(abs_file);
% downsampling
flag = flag(1:2:end,1:2:end);
x = x(1:2:end,1:2:end);
y = y(1:2:end,1:2:end);
z = z(1:2:end,1:2:end);

[rows,cols,val] = find(flag==1);
cnt = length(rows);
vertices2d = zeros(cnt,2);% 2-D coordinate, row-wise
vertices3d = zeros(cnt,3);% 3-D coordinate, row-wise
for i = 1:cnt
    vertices3d(i,:) = [x(rows(i),cols(i)),y(rows(i),cols(i)),z(rows(i),cols(i))];
    vertices2d(i,:) = [rows(i),cols(i)];
end
mesh.vertices = vertices3d;
mesh.vertexNormals = surfnorm(vertices3d);% surface normals
side_len = 0;
for i = 1:(E-1)
    for j = (i+1):E
        delta = vertices3d(i,:) - vertices3d(j,:);
        side_len = side_len + sqrt(delta*delta');
    end
end

vertices3d = vertices3d'; % column-wise

side_len = side_len/(E*(E-1)/2);% bin_size
temp_idx = find(abs_file=='.');
file_name = abs_file((find(abs_file=='U')+15):temp_idx(2)-1);
file_name = strcat(file_name,'.mat');
if exist(strcat('.\spin_imgs\',file_name),'file') == 0 
    spinImage = calcSpinImages(mesh,side_len,resolution, index);
    save(strcat('.\spin_imgs\',file_name),'spinImage');
else
    load(strcat('.\spin_imgs\',file_name));
end

dsc_lmk = [];
shape_lmk = [];
if ~isempty(pts_file) && L>0 
    shape = read_shape(pts_file,L);
    dsc_lmk = zeros(resolution^2,L);
    shape_lmk = zeros(3,L);
    delta_radius = 1;

    shape = round(shape/2);% downsampling
    for i = 1:L
        idx = find((abs(vertices2d(:,1)-shape(i,2))<=delta_radius).*(abs(vertices2d(:,2)-shape(i,1))<=delta_radius)==1);
        if ~isempty(idx)
            if ~isempty(spinImage)
                dsc_lmk(:,i) = spinImage(idx(1)).spinIm(:);
            end
            shape_lmk(:,i) = vertices3d(:,idx(1));
        else
            disp('Landmark spin image absent...')
            dsc_lmk = [];
            shape_lmk = [];
            break
        end
    end
end

end

