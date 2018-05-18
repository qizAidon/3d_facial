function spinImages = calcSpinImages(mesh, binSize, spinImageSize, index)

% spinImages = calcSpinImages(mesh, binSize, spinImageSize)
%
% This is a quick implementation of the spin images. A. Johnson and M. Hebert, 
% "Using Spin Images for Efficient Object Recognition in Cluttered 3D
% Scenes", IEEE TRANS. ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE.
%
% calculates spin images of a mesh which are used for mesh registration and
% 3D object recognition. 
%
%
% Input Argument: mesh - data structure .. must have vertexNormals
%                      - mesh.vertices, 顶点
%                      - mesh.vertexNormals, 顶点处的法向量
%              binSize - resolution, 通常取所有边的平均值
%        spinImageSize - size of spin image, 通常取10*10或者20*20
%
% Output Argument: spinImages - data structure of spinIm
%                             - spinImages(k).spinIm, 即顶点k的spin image
%
% Copyright : This code is written by Ajmal Saeed Mian {ajmal.mian@uwa.edu.au}
%              Computer Science, The University of Western Australia. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion this copyright information.
%
% Disclaimer : This code is provided as is without any warrantly.

vertices3d = mesh.vertices;
vertices_norm = mesh.vertexNormals;
N = length(mesh.vertices);
for ii = 1 : N
    sprintf('No. %d picture, Total %d points, it is No. %d vertex now...\n',index,N,ii)
    spinIm = zeros(spinImageSize);
    Pt = vertices3d(ii,:);
    Normal = vertices_norm(ii,:);
    if isnan(Normal(1))
        continue;
    end
    for jj = 1 : length(vertices3d)
        thisPt = vertices3d(jj,:);
        thisNormal = vertices_norm(jj,:);
        if acos(Normal*thisNormal')>pi/3
            continue;
        end
        alpha = norm(cross((thisPt - Pt),Normal));
        beta = Normal(1)*(thisPt(1)-Pt(1)) + Normal(2)*(thisPt(2)-Pt(2)) + Normal(3)*(thisPt(3)-Pt(3));
        %beta = beta/sqrt(Normal(1)^2 + Normal(2)^2 + Normal(3)^2);
        % shifting the center up
         
        beta = beta - spinImageSize*binSize/2;
        if beta > 0 || beta < -spinImageSize*binSize
            continue;
        end        
        b = floor(-beta/binSize) + 1;
        
        a = floor(alpha/binSize) + 1;
        if a > spinImageSize
            continue;
        end
        spinIm(a,b) = spinIm(a,b) + 1;        
    end
    spinImages(ii).spinIm = spinIm;
end  

end