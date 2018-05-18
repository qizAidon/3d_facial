clear;
clc;
addpath(genpath('.\files'));

L = 5;
num_sample = 49;
E = 100; % 选取100个vertices,用于计算边的平均长度
resolution = 10;
dsc_lmk_sample = zeros(L*resolution^2,num_sample);% matrix storing descriptors of lamdmarks in samples
shape_lmk_sample = zeros(L*3,num_sample);
%% ----------- ① calculate descriptors of all landmarks in samples ----------- %%
pts_path = '.\files\UND-G-45lr\';
abs_path = '.\files\UND-HID-data-G\';
pts_list = dir(strcat(pts_path,'*.pts'));
abs_list = dir(strcat(abs_path,'*.abs'));
cnt = 1;
if ~isempty(pts_list) && ~isempty(abs_list)
    % 以第一个样本为参照样本
    pts_file = strcat(pts_path,pts_list(1).name);
    disp(pts_file)
    file_name = pts_list(1).name(1:(find(pts_list(1).name=='.')-1));
    abs_file = strcat(abs_path,file_name,'.abs');

    [~,~,dsc_lmk,shape_lmk] = meshSpin(pts_file,abs_file,L,resolution,E, 1);
    dsc_lmk_sample(:,1) = dsc_lmk(:);
    shape_lmk_sample(:,1) = shape_lmk(:);

    ref_shape = shape_lmk_sample(:,1);
    ref_shape = reshape(ref_shape,3,5);
    ref_shape = ref_shape';
    for k = 2:num_sample
        pts_file = strcat(pts_path,pts_list(k).name);
        disp(pts_file)
        file_name = pts_list(k).name(1:(find(pts_list(k).name=='.')-1));
        abs_file = strcat(abs_path,file_name,'.abs');

        [~,~,dsc_lmk,shape_lmk] = meshSpin(pts_file,abs_file,L,resolution,E, k);
        shape_lmk = shape_lmk';
        [d,Z,tform] = procrustes(ref_shape,shape_lmk,'scaling',0,'reflection',0);
        shape_lmk = tform.b * shape_lmk * tform.T + repmat(tform.c(1,:),size(shape_lmk,1),1); 
        shape_lmk = shape_lmk';
        dsc_lmk_sample(:,k) = dsc_lmk(:);
        shape_lmk_sample(:,k) = shape_lmk(:);
    end
    
end
dsc_template = mean(dsc_lmk_sample,2);
dsc_template = reshape(dsc_template,resolution^2,L);
shape_template = mean(shape_lmk_sample,2);
shape_template = reshape(shape_template,3,L);

%% ----------------------- ②calculate phi and lambda ------------------------ %%
%{
shape_lmk_sample = shape_lmk_sample - repmat(mean(shape_lmk_sample,2),1,size(shape_lmk_sample,2));
% for eig(), each row is an observation
[V,D,W] = eig(cov(shape_lmk_sample'));% V is right eigenvectors, W is left eigenvectors; 
eigenvalues = diag(D);
eigenvalues = eigenvalues(end:-1:1);
V = V(:,end:-1:1); V = V';
pc = V*shape_lmk_sample;
%}
mean_sample = mean(shape_lmk_sample,2);
[Phi,pc,Lambda] = princomp(shape_lmk_sample'); % Phi is row-wise
% Phi = Phi';
beta = 4;
M = 9;% the number of retained principle components

%% ③输入一个mesh,希望找到所有L个landmarks.对所有的vertices,获取3D坐标,计算spin image
% 取第50个样本做测试
pts_file = strcat(pts_path,pts_list(50).name);
disp(pts_file)
file_name = pts_list(50).name(1:(find(pts_list(50).name=='.')-1));
abs_file = strcat(abs_path,file_name,'.abs');

[spin_image,vertices3d,dsc_lmk,shape_lmk] = meshSpin(pts_file,abs_file,L,resolution,E, k);
% 坐标对齐（这里是作为测试，实际中用算法1中预测得到的5个关键点作为测试样本与参考样本对齐的参照点）
% #########################################################################
shape_lmk = shape_lmk';
[d,Z,tform] = procrustes(ref_shape,shape_lmk,'scaling',0,'reflection',0);
shape_lmk = tform.b * shape_lmk * tform.T + repmat(tform.c(1,:),size(shape_lmk,1),1); 
vertices3d = vertices3d';
vertices3d = tform.b * vertices3d * tform.T + repmat(tform.c(1,:),size(vertices3d,1),1);
vertices3d = vertices3d';

%% ④设定acceptance radius = 10mm, ek = 50，找到每个landmark template的candidate set
accept_r = 10;
ek = 50;
num_pt = size(vertices3d,2);
spin_img_mat = zeros(resolution^2,num_pt);
candidate_flag = zeros(L,num_pt);
% score_mat = zeros(L,num_pt); % score = -dist_dsc
% rd_mat = zeros(L,num_pt); % RD
% highest_idx = zeros(L,2);
for i = 1:num_pt
    spin_img_mat(:,i) = spin_image(i).spinIm(:);
end
for i = 1:L
    dsc_temp = dsc_template(:,i);
    shape_temp = shape_template(:,i);
    
    dist_shape = vertices3d - repmat(shape_temp,1,num_pt);
    dist_shape = sqrt(sum(dist_shape.^2));
    idx1 = dist_shape<=accept_r;
%    candidate_flag(i,idx1) = 1;
    
    dist_dsc = spin_img_mat - repmat(dsc_temp,1,num_pt);
    dist_dsc = sqrt(sum(dist_dsc.^2));
    MAX = max(dist_dsc) + 1;
    dist_dsc(~idx1) = MAX;
%    score_mat(i,:) = -dist_dsc;
    
    rd_vec = zeros(1,num_pt);
    for j = 1:num_pt
        rd_vec(j) = sum(dist_dsc<=dist_dsc(j));
    end
    idx2 = rd_vec<=ek;
    idx = logical(idx1.*idx2);
    candidate_flag(i,idx) = 1;
    
    lmk_candidate(i).set = find(idx>0);% set里面存储第i个landmark的包含了哪些点作为候选点
end

%% ⑤利用candidate set确定输入mesh对应的 L 个landmark,即Algorithm 1
fix_init = 4;
all_comb_lmk = nchoosek(1:L,fix_init);
all_lmk_types = size(all_comb_lmk,1);
index_all = 1:3*L;
lmk_all = 1:L;
cnt_comb_valid = 1;
for i = 1:all_lmk_types
    fix_lmk = all_comb_lmk(i,:);% fix_lmk指示哪四个landmark作为fixed
    fix1 = lmk_candidate(fix_lmk(1)).set;
    fix2 = lmk_candidate(fix_lmk(2)).set;
    fix3 = lmk_candidate(fix_lmk(3)).set;
    fix4 = lmk_candidate(fix_lmk(4)).set;
    
    for x1 = 1:length(fix1)
        comb_all(cnt_comb_valid).lmk = fix1(x1);
        
        lmk_fix(1).idx = fix_lmk(1);% idx指示该点来自于哪一个landmark的candidate set
        lmk_fix(1).pt = vertices3d(:,fix1(x1));% pt为该点的三维坐标
        
        for x2 = 1:length(fix2)
            comb_all(cnt_comb_valid).lmk = [comb_all(cnt_comb_valid).lmk,fix2(x2)];
            
            lmk_fix(2).idx = fix_lmk(2);
            lmk_fix(2).pt = vertices3d(:,fix2(x2));
            
            for x3 = 1:length(fix3)
                comb_all(cnt_comb_valid).lmk = [comb_all(cnt_comb_valid).lmk,fix3(x3)];
                
                lmk_fix(3).idx = fix_lmk(3);
                lmk_fix(3).pt = vertices3d(:,fix3(x3));
                
                for x4 = 1:length(fix4)
                    sprintf('[x1, x2, x3, x4] = [%d, %d, %d, %d]\n',x1,x2,x3,x4)
                    comb_all(cnt_comb_valid).lmk = [comb_all(cnt_comb_valid).lmk,fix4(x4)];
                    
                    lmk_fix(4).idx = fix_lmk(4);
                    lmk_fix(4).pt = vertices3d(:,fix4(x4));
                    
                    % line 7
                    yf = [lmk_fix(1).pt;lmk_fix(2).pt;lmk_fix(3).pt;lmk_fix(4).pt];
                    index_fix = [(fix_lmk(1)-1)*3+1:fix_lmk(1)*3,(fix_lmk(2)-1)*3+1:fix_lmk(2)*3,(fix_lmk(3)-1)*3+1:fix_lmk(3)*3,(fix_lmk(4)-1)*3+1:fix_lmk(4)*3];
                    yf = yf - mean_sample(index_fix);
                    Phif = Phi(index_fix,:);
                    Phig = Phi;
                    Phig(index_fix,:) = [];
                    
                    % 改变lambda的顺序,lambda_guess,lambda_fix
                    lambda_guess = Lambda;
                    lambda_guess(index_fix) = [];
                    lambda_fix = Lambda(index_fix);
                    lambda_temp = [lambda_guess;lambda_fix];
                    
                    yg = -pinv(Phig*pinv(diag(lambda_temp))*Phig')*(Phig*pinv(diag(lambda_temp))*Phif')*yf; % line 8
                    y = [yg;yf];
                    b = [Phig;Phif]*y;
         
                    [lambda_val,lambda_idx] = sort(lambda_temp,'descend');
                    cond = sum(b(lambda_idx(1:M)).^2./lambda_val(1:M));% condition after 4-tuple inference
                    sprintf('Condition is %d...\n',cond)
                    
                    gama_k = [];
                    idx_k = [];
                    gama_comb = [];
                    
                    guess_lmk = 1:L;
                    guess_lmk(fix_lmk) = [];% guess_lmk指示剩下的landmark作为guess
                    while cond < beta^2 && ~isempty(guess_lmk) % line 9
                       for j = 1:length(guess_lmk)
                           candi_temp = lmk_candidate(guess_lmk(j)).set;
                           cost_ck = zeros(1,length(candi_temp));
                           index_k = (guess_lmk(j)-1)*3+1:guess_lmk(j)*3;
                           
                           index_fix_new = [index_fix,index_k]; % line 12
                           lambda_guess_new = Lambda;
                           lambda_guess_new(index_fix_new) = [];
                           lambda_fix_new = Lambda(index_fix_new);
                           lambda_temp_new = [lambda_guess_new;lambda_fix_new];

                           Phig_test = Phi;
                           Phig_test(index_fix_new,:) = [];
                           Phif_test = Phi(index_fix_new,:);

                           fix_lmk_new = [fix_lmk,guess_lmk(j)]; % 每次加入一个guess landmark
                           guess_lmk_new = 1:L;
                           guess_lmk_new(fix_lmk_new) = [];
                           
                           mean_sample_g = mean_sample;
                           mean_sample_g(index_fix_new) = [];
                           mean_sample_f = mean_sample(index_fix_new);
                           mean_sample_new = [mean_sample_g;mean_sample_f];
                           
                           for k = 1:length(candi_temp) % line 11
                               ck = vertices3d(:,candi_temp(k));
                               ck = ck - mean_sample(index_k);
                               yf_test = [yf;ck]; % line 12
                               
                               delta_x_test_in = zeros(1,length(fix_lmk_new));
                               if ~isempty(guess_lmk_new)
                                   yg_test =  -pinv(Phig_test*pinv(diag(lambda_temp_new))*Phig_test')*(Phig_test*pinv(diag(lambda_temp_new))*Phif_test')*yf_test; % line 13
                                   y_test = [yg_test;yf_test];% line 13
                                   x_test = y_test + mean_sample_new; 
                                   
                                   delta_x_test_nin = zeros(1,length(guess_lmk_new)); 
                                   
                                   for p = 1:length(guess_lmk_new)
                                       x_test_lk = x_test((guess_lmk_new(p)-1)*3+1:guess_lmk_new(p)*3);
                                       candidate_lk = lmk_candidate(guess_lmk_new(p)).set;% lk 的candidates对应的点的标号
                                       candidate_lk = vertices3d(:,candidate_lk); % lk 的candidates对应的3-D坐标
                                       dist = candidate_lk - repmat(x_test_lk,1,size(candidate_lk,2));
                                       dist = dist.^2;
                                       dist = sum(dist);
                                       dist = min(dist);
                                       delta_x_test_nin(p) = dist;
                                   end
                                   
                                   delta_x_test = [delta_x_test_in,delta_x_test_nin];
                               else 
                                   y_test = yf_test;
                                   x_test = y_test + mean_sample_new;
                                   delta_x_test = delta_x_test_in;
                               end
                               cost_ck(k) = median(delta_x_test);% line 15
                           end % line 16
                           cost_ck_min = min(cost_ck);
                           idx_min = find(cost_ck==min(cost_ck));
                           gama_k = [gama_k,cost_ck_min];
                           idx_k = [idx_k,idx_min(1)]; 
                       end
                       [~,idx_add] = find(gama_k==min(gama_k));
                       gama_comb = [gama_comb,min(gama_k)];
                       ck_add = vertices3d(:,candi_temp(idx_k(idx_add(1))));
                       yf = [yf;ck_add-mean_sample((guess_lmk(idx_add(1))-1)*3+1:guess_lmk(idx_add(1))*3)];
                       comb_all(cnt_comb_valid).lmk = [comb_all(cnt_comb_valid).lmk,candi_temp(idx_k(idx_add(1)))];
                       
                       % calculate yg_test using formula (17), get y_test 
                       index_k = (guess_lmk(idx_add(1))-1)*3+1:guess_lmk(idx_add(1))*3;
                       index_fix = [index_fix,index_k];
                       lambda_guess_new = Lambda;
                       lambda_guess_new(index_fix) = [];
                       lambda_fix_new = Lambda(index_fix);
                       lambda_new = [lambda_guess_new;lambda_fix_new];
                       
                       mean_sample_g = mean_sample;
                       mean_sample_g(index_fix) = [];
                       mean_sample_f = mean_sample(index_fix);
                       mean_sample_new = [mean_sample_g;mean_sample_f];
                       
                       Phig = Phi;
                       Phig(index_fix,:) = [];
                       Phif = Phi(index_fix,:);
                       yg =  -pinv(Phig*pinv(diag(lambda_new))*Phig')*(Phig*pinv(diag(lambda_new))*Phif')*yf;
                       y = [yg;yf];
                       
                       b = [Phig;Phif]*y;
                   
                       % calculate condition
                       [lambda_val,lambda_idx] = sort(lambda_new,'descend');
                       cond = sum(b(lambda_idx(1:M)).^2./lambda_val(1:M));
                       if cond < beta^2
                          fix_lmk_new;
                          guess_lmk = 1:L;
                          guess_lmk(fix_lmk_new) = [];
                       end
                    end
                    
                    if ~isempty(gama_comb)
                        score(cnt_comb_valid) = length(yf) + exp(-gama_comb(end)); % line 22
                        
                        lmk_idx = comb_all(cnt_comb_valid).lmk;
                        lmk_find = [vertices3d(:,lmk_idx(1))';vertices3d(:,lmk_idx(2))';...
                            vertices3d(:,lmk_idx(3))';vertices3d(:,lmk_idx(4))';vertices3d(:,lmk_idx(5))'];

                        lmk_dist = lmk_find - shape_lmk;
                        lmk_dist = lmk_dist.^2;
                        lmk_dist = sum(lmk_dist,2);
                        lmk_dist = sqrt(lmk_dist);
                        lmk_dist = mean(lmk_dist);
                        dist2temp(cnt_comb_valid) = lmk_dist;        
                    else
                        score(cnt_comb_valid) = length(yf);
                        dist2temp(cnt_comb_valid) = 0;
                    end
                    score_cond(cnt_comb_valid) = cond;
                    sprintf('No.%d score is %d...',cnt_comb_valid,score(cnt_comb_valid))

                    cnt_comb_valid = cnt_comb_valid + 1;
                    comb_all(cnt_comb_valid).lmk = comb_all(cnt_comb_valid-1).lmk(1:3);

                end
            end
        end
    end
end
[score_max,idx_max] = max(score);
set_optim = comb_all(idx_max(1)).lmk;