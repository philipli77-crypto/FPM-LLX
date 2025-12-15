clc;
clear all;
load("dataFPM_HEpathology_v4.mat");
F = @(x) fftshift(fft2(ifftshift(x)));
IF = @(x) fftshift(ifft2(ifftshift(x)));

rng(42);
datacell_raw=cell(1,size(dataFPM,3));
for i=1:size(dataFPM,3)
    datacell_raw{i}=dataFPM(:,:,i);
    
end
figure('Color','w');
tiledlayout(4,5,'TileSpacing','compact','Padding','compact'); % 替代 subplot，更整洁

k_0_temp = 2 * pi / lambda;
real_kx_temp = k_0_temp .* kx;
real_ky_temp = k_0_temp .* ky;
k_NA_limit_temp = k_0_temp * NA;


filter_radius = 0.96 * k_NA_limit_temp; 
valid_indices = find(sqrt(real_kx_temp.^2 + real_ky_temp.^2) < filter_radius);

fprintf('原始图片数量: %d\n', length(kx));
fprintf('剔除边缘后保留数量: %d\n', length(valid_indices));


% fprintf('原始图片数量: %d\n', length(kx));
% 
% valid_indices = setdiff(1:length(kx), [16 18 20]);
% 
% 
% fprintf('剔除第18张后保留数量: %d\n', length(valid_indices));

% 执行剔除
dataFPM = dataFPM(:, :, valid_indices);
kx = kx(valid_indices);
ky = ky(valid_indices);

% % 计算所有图像的全局最小值和最大值
% allMin = inf;
% allMax = -inf;
% for i = 1:length(valid_indices)
%     currentMin = min(datacell_raw{i}(:));
%     currentMax = max(datacell_raw{i}(:));
%     allMin = min(allMin, currentMin);
%     allMax = max(allMax, currentMax);
% end
% 
% figure('Color','w');
% tiledlayout(4,5,'TileSpacing','compact','Padding','compact');
% 
% for i = 1:length(valid_indices)
%     nexttile;
%     imagesc(datacell_raw{i});
%     axis image off;
%     colormap(gca,'gray');       % 灰度色图，可改为其他
%     clim([allMin allMax]);    
%     title(sprintf('Image %d', i), 'FontSize', 10);% 统一颜色范围
% end
% 
% % 如果需要统一的颜色条，可以加上：
% cb = colorbar('Location','eastoutside');
% cb.Layout.Tile = 'east';
%%
%物domain \delta_x
ob_dx=dpix_c/ mag;

%波数
wavelength=lambda;
k_0=2*pi/wavelength;
%EPRY参数
a_EPRY=0.5;
b_EPRY=0.5;
%把波数重定义成三角函数
alpha=kx;
beta=ky;

% real_kx=kx;
% real_ky=ky;
real_kx=k_0.*alpha;
real_ky=k_0.*beta;
%用以比较kx kychange  作为copy
real_kx_1=real_kx;
real_ky_1=real_ky;
%N.A.=a*sin\theta = n*sqrt(1-(k_0/kz)^2) 我们得到 在频域fx fy中的半径 令n=1
R_f=NA/wavelength;% 截止频率 单位1/m

N=size(dataFPM,1);
N_LED=length(kx);
%空域长度
L_x=N*ob_dx;
L_y=N*ob_dx;

%spatial domain 坐标
x_obj=(-N/2+1:N/2)*ob_dx;
y_obj=(-N/2+1:N/2)*ob_dx;
[X,Y]=meshgrid(x_obj,y_obj);

%频域
d_fx=1/L_x;
d_fy=1/L_y;
B_fx=1/ob_dx;
B_fy=1/ob_dx;

%频域网格
fx_f=(-N/2+1:N/2)*d_fx;
fy_f=(-N/2+1:N/2)*d_fy;
[FX,FY]=meshgrid(fx_f,fy_f);
KX=2*pi*FX;
KY=2*pi*FY;

PixelX=FX./d_fx;
PixelY=FY./d_fy;

%%
%calibration -1 bright field 
%找到所有的bright field indices ？适当缩小避免边界位置？
bf_indices = find(sqrt(real_kx.^2+ real_ky.^2) < 1*R_f*2*pi);
%对所有Iquta 归一化,I_quta对应原文的自相关
Sum_I=zeros(N,N);
for i=1:length(bf_indices)
    bf_led_idx=bf_indices(i);
    I_r=dataFPM(:,:,bf_led_idx);
    Sum_I=Sum_I+abs(F(I_r));
end
Sum_I_mean=Sum_I./length(bf_indices);

Iquta=zeros(512,512,length(valid_indices));
I_qutaf=zeros(512,512,length(valid_indices));
for i=1:length(valid_indices)
    Iquta(:,:,i)=F(double(dataFPM(:,:,i)));
    epsilon = 1e-6; 
    I_qutaf(:,:,i)= abs(Iquta(:,:,i)) ./ (Sum_I_mean + epsilon);
end

%%
% % 直接求和（不太行）
% % change=zeros(length(valid_indices),2);
% % 
% % %通过梯度法进行准直
% % %直接把校准后的kx ky返回real_kx real_ky
% % for i=1:length(bf_indices)
% %     I_raw=dataFPM(:,:,i);
% %     %直接得到normalize 后的频域两个圆
% %     Spec=abs(F(I_raw))./Sum_I_mean;
% % 
% %     %高斯模糊
% %     Spec=imgaussfilt(Spec,2);
% %     Spec=Spec.*pupil_OTF;
% %     %卷积圆环找到最大值
% %     FIND_center=abs(IF(F(Spec).*F(pupil_ring)));
% %     % 找到前 k 个最大值及其线性索引
% %     k=5;
% %     [val, idx] = maxk(FIND_center(:), k);
% % 
% %     % 转换为二维坐标 (行,列)
% %     [row, col] = ind2sub(size(FIND_center), idx); update_ky = (row(1)-256)*d_fy*2*pi; update_kx = (col(1)-256)*d_fx*2*pi;
% %     %更新kx ky
% % 
% %     if update_kx*real_kx(i)>0
% %         real_kx(i)=update_kx;
% %     else
% %         real_kx(i)=-update_kx;
% %     end
% %     if update_ky*real_ky(i)>0
% %         real_ky(i)=update_ky;
% %     else
% %         real_ky(i)=-update_ky;
% %     end  
% %     change(i,:) = [real_kx(i)- real_kx_1(i), real_ky(i)- real_ky_1(i)];
% % 
% %     disp([real_kx(i),real_ky(i),change(i,:)]);
% % 
% % end
% %粗略绘图显示两个圆,通过绝度值寻找边缘%设置OTFpupil过滤掉圆外杂信息
% OTF=2*R_f;
% pupil_OTF=double(sqrt(FX.^2+FY.^2)<=OTF);
% %设置圆环pupil
% pupil_ring=double(sqrt(FX.^2 + FY.^2) <= (R_f/d_fx+sqrt(2))*d_fx)-double(sqrt(FX.^2 + FY.^2) <= (R_f/d_fx-sqrt(2))*d_fx);
% %平均值mean(|Iquta|) 
% I_raw = double(dataFPM(:,:,3)); 
% %归一化
% I_no_DC = I_raw ; 
% Spec = F(I_no_DC);
% % 添加一个极小值，防止分母太小导致边缘数值飞出天际
% epsilon = 1e-5 * max(Sum_I_mean(:)); 
% I_normal = abs(Spec) ./ (Sum_I_mean + epsilon);
% 
% % level = graythresh(I_normal); 
% % 
% % % 4. 生成二值 Mask
% % %    Mask 为 1 的地方是信号，为 0 的地方是背景
% % Otsu_Mask = imbinarize(I_normal, level);
% % 
% % % 5. 膨胀 Mask (安全措施)
% % %    往外扩 4 个像素，防止边缘被误杀
% % Otsu_Mask = imdilate(Otsu_Mask, strel('disk', 4));
% % 
%
%R' calibration
% R_update_list=zeros(length(valid_indices),1);
% R_update=0;
% Radius_pix = R_f / d_fx;
% pixel_x=real_kx_1/2/pi/d_fx;  pixel_y=real_ky_1/2/pi/d_fy;  
% for ID=1:length(valid_indices)
%     %统一到pixel 
%     [Gx,Gy]=gradient(I_qutaf(:,:,ID));
% 
% 
%     shift_pixel_x=pixel_x(ID); shift_pixel_y=pixel_y(ID);
%     %对应第ID图的 r_map  kx-f-pixel
%     r_map = sqrt((PixelX-shift_pixel_x).^2+(PixelY-shift_pixel_y).^2);
%     range=5;
%     thick_ring_mask = double(abs(r_map - Radius_pix) <= range);
%     %构造R向量
%     r_norm=r_map;
%     r_norm(r_norm==0) = 1e-6; 
%     Vec_x =(PixelX-shift_pixel_x) ./ r_norm.*thick_ring_mask; % 指向圆心的 X 分量的厚圆环
%     Vec_y =(PixelY-shift_pixel_y) ./ r_norm.*thick_ring_mask; % 指向圆心的 Y 分量的厚圆环
%     %环状梯度图
%     gradient_map=(Vec_x.*Gx+Vec_y.*Gy).*thick_ring_mask;
%     Gradsum=0;
%     gradient_map_iteration=gradient_map;
%     for R_prime=Radius_pix-range:Radius_pix+range
%         Inner_pupil=double(r_map<R_prime);
%         Grad_save=sum((Inner_pupil.*gradient_map_iteration),'all');
% 
%         if Grad_save > Gradsum
%             Gradsum =Grad_save;
%             R_update=R_prime;
%         end
% 
%         gradient_map_iteration=gradient_map_iteration.*(1-Inner_pupil);
%         R_update_list(ID)=R_update;
%     end
% end
% R_update=mean(R_update_list);
% R_f=R_update*d_fx;

%%
%%找到梯度最大的圆心

% I_normal=I_qutaf(:,:,1);
% OTF = 2 * R_f*1.1; %适当放大OTF避免边界位置的问题
% pupil_OTF = double(sqrt(FX.^2+FY.^2) <= OTF);
% I_normal = I_normal .* pupil_OTF;
% I_normal = imgaussfilt(I_normal, 2);
% 
% 
% %梯度矢量 
% [Gx,Gy]=gradient(I_normal);
% %构造矢量卷积核
% %构造圆环
% radius_pix = R_update;
% r_map = sqrt(FX.^2 + FY.^2) / d_fx;
% ring_mask = double(abs(r_map - radius_pix) <= 0.75);
% %构造R向量
% r_norm = sqrt(FX.^2 + FY.^2);
% r_norm(r_norm==0) = 1e-6; % 防除零
% vec_x = FX ./ r_norm; % 指向圆心的 X 分量
% vec_y = FY ./ r_norm; % 指向圆心的 Y 分量
% %X Y卷积核
% Kernel_X = vec_x .* ring_mask;
% Kernel_Y = vec_y .* ring_mask;
% 
% K_sum = sum(ring_mask(:));
% Kernel_X = Kernel_X / K_sum;
% Kernel_Y = Kernel_Y / K_sum;
% 
% Conv_X = real(IF(F(Gx) .* F(Kernel_X)));
% Conv_Y = real(IF(F(Gy) .* F(Kernel_Y)));
% 
% FIND_center_grad=Conv_Y+Conv_X;
% figure('Name', 'Gradient Vector Search');
% imagesc(-N/2+1 : N/2, -N/2+1 : N/2, FIND_center_grad); 
% axis image; axis xy; colormap jet; colorbar;
% title('Gradient search center');
% xlabel('k_y shift (pixels)'); ylabel('k_x shift (pixels)');
% 
% % 找到前 k 个最大值
% k = 10;                 
% [val, idx] = maxk(FIND_center_grad(:), k);
% [row,col]=ind2sub(size(FIND_center_grad),idx);
% 
% row_y=row-256; col_x= col-256;
% row_ky=row_y*d_fy*2*pi; col_kx=col_x*2*pi*d_fx;
% %标记圆心
% hold on;
% plot( col_x(1:2), row_y(1:2),'w+', 'MarkerSize', 15, 'LineWidth', 2);
% col_kx=col_kx(1,:); row_ky=row_ky(1,:);
% 
% % 显示结果
% disp(table(val, row_y, col_x));
% disp(['col_kx is '  num2str(col_kx) ', and row_ky is ' num2str(row_ky)]);
% figure; 
% imagesc( I_normal); 
% axis image; colormap jet; colorbar;
% title('normalized I-quta');
% %viscircles([real_kx_1(3)/d_fx/2/pi, real_ky_1(3)/d_fy/2/pi], R_f/d_fx, 'Color','b','LineWidth',1);
% 
% viscircles([ col(1:2),row(1:2)], R_f/d_fx, 'Color','k','LineWidth',1);
% viscircles([ col(1:2),row(1:2)], R_update, 'Color','b','LineWidth',1);
% %viscircles([ 256+shift_pixel_y,256+shift_pixel_x], R_f/d_fx, 'Color','r','LineWidth',1);
% %这里的x y 由于习惯问题存在颠倒

%%
changes=zeros(length(valid_indices),2);
col_list=zeros(length(valid_indices),1);
row_list=zeros(length(valid_indices),1);
%梯度卷积搜寻法
for i=1:length(bf_indices)
    ii_bf=bf_indices(i);

 
    I_normal = I_qutaf(:,:,ii_bf);
    OTF = 2 * R_f*1.1; %适当放大OTF避免边界位置的问题
    pupil_OTF = double(sqrt(FX.^2+FY.^2) <= 1.2*OTF);
    I_normal = I_normal .* pupil_OTF;
    I_normal = imgaussfilt(I_normal, 2);
    
    

    %梯度矢量 
    [Gx,Gy]=gradient(I_normal);
    %构造矢量卷积核
    %构造圆环
    radius_pix = R_f/d_fx ;
    r_map = sqrt(FX.^2 + FY.^2) / d_fx;
    ring_mask = double(abs(r_map - radius_pix) <= 0.75);
    %构造R向量
    r_norm = sqrt(FX.^2 + FY.^2);   
    r_norm(r_norm==0) = 1e-6; % 防除零
    vec_x = FX ./ r_norm; % 指向圆心的 X 分量
    vec_y = FY ./ r_norm; % 指向圆心的 Y 分量
    %X Y卷积核
    Kernel_X = vec_x .* ring_mask;
    Kernel_Y = vec_y .* ring_mask;
    %归一化 避免圈住的像素数浮动
    K_sum = sum(ring_mask(:));
    Kernel_X = Kernel_X / K_sum;
    Kernel_Y = Kernel_Y / K_sum;

    Conv_X = real(IF(F(Gx) .* F(Kernel_X)));
    Conv_Y = real(IF(F(Gy) .* F(Kernel_Y)));
    FIND_center_grad=Conv_Y+Conv_X;
    % 找到前 k 个最大值
    k = 10;                 
    [val, idx] = maxk(FIND_center_grad(:), k);
    [row,col]=ind2sub(size(FIND_center_grad),idx);
    [y,x]=deal(col,row);
    
    row_y=y-256; col_x=x-256;
    row_ky=row_y*d_fy*2*pi; col_kx=col_x*2*pi*d_fx;

    col_kx=col_kx(1,:); row_ky=row_ky(1,:);
    %由于对称性 调整正确的正负号
    if col_kx*real_kx_1(ii_bf)<0
        col_kx=-col_kx;
        col_idx=-col_x(1);
    else 
        col_idx=col_x(1);
    end
    if row_ky*real_ky_1(ii_bf)<0
        row_ky=-row_ky;
        row_idx=-row_y(1);
    else 
        row_idx=row_y(1);
    end
    real_kx(ii_bf)=col_kx;
    real_ky(ii_bf)=row_ky;
    col_list(ii_bf)=col_idx;
    row_list(ii_bf)=row_idx;
    %计算change shift
    changes(ii_bf,:)=[real_kx(ii_bf)-real_kx_1(ii_bf), real_ky(ii_bf)-real_ky_1(ii_bf)];

end


%%
%高分辨率图像网格
max_kshift=max(sqrt(real_kx.^2+real_ky.^2));
max_f=max_kshift/(2*pi)+R_f;
N_hr=ceil(N*double(max_f/R_f));
if mod(N_hr,2)~=0
    N_hr=N_hr+1;

end

%ob_dx_high_resolution=ob_dx/up_sample_frate;
%初始化高分辨率 频域网格 全零
frequency_domain_hr_result=zeros(N_hr,N_hr);

%高分辨率不改变spatial domain L大小和d_fx d_fy 大小,因此，设置frequency domain B
B_f_h_r=N_hr*d_fx;

%初始aperture
Pupil_LR_initial = double(sqrt(FX.^2 + FY.^2) <= R_f);
%aperture会不断变化-EPRY-用在循环中
Pupil=Pupil_LR_initial;

%频域平移
dkx_LR = 2 * pi * d_fx;
dky_LR = 2 * pi * d_fy;
k_idx_x = round(real_kx / dkx_LR);
k_idx_y = round(real_ky / dky_LR);
%循环迭代次数
N_iter=50;
%total_amplitude=sum(sqrt(dataFPM,3));
%%
% figure('Name','Try sum(sqrt)');

% imagesc(sum(sqrt(dataFPM),3));
% axis image; colormap hsv; colorbar
% 
 figure('Name','Try sqrt(sum)')
  imagesc(sqrt(sum(dataFPM,3)))
  axis image; colormap gray; colorbar
% 

%%
%初始化
[~,center_idx]=min(real_kx.^2 + real_ky.^2);
initial_amplitude = sum(sqrt(dataFPM),3);


patchFT_center = F(initial_amplitude);
% 获取实际的偏移量
center_shift_x = k_idx_x(center_idx);
center_shift_y = k_idx_y(center_idx);
idx_hr_y = (N_hr/2 - N/2 + 1 ) : (N_hr/2 + N/2);
idx_hr_x = (N_hr/2 - N/2 + 1 ) : (N_hr/2 + N/2);

% 将频谱填入正确的位置
frequency_domain_hr_result(idx_hr_x, idx_hr_y) = patchFT_center .* Pupil_LR_initial;
%%%%% measurement 全部加起来 作为initialization
% %
% 看camera
% figure('Name', 'total ');
% subplot(1, 2, 1);
% imagesc(abs(total_amplitude));
% axis image; 
% title('1');
% subplot(1, 2, 2);
% imagesc(center_amplitude);
% axis image; 
% title('2)')
% colormap hsv; 
% 
% h=colorbar('southoutside');
% h.Position = [0.25 0.05 0.5 0.03];
%%

%FPM  phase retrieval
%R_f=R_update*d_fx;
for iter_time=1:N_iter
    fprintf('Iteration：%d /%d\n',iter_time,N_iter);
    idx_rand = randperm(N_LED);
    
    for p_idx =1:N_LED
        led_idx = idx_rand(p_idx);
        
        measured_amplitude = sqrt(dataFPM(:,:,led_idx));
        
        %波矢->频率 
        %fx_led=real_kx(led_idx)/(2*pi);
        %fy_led=real_ky(led_idx)/(2*pi);

        %构建aperture，包含在当前角度下位移了多少
        fx_shift=k_idx_x(led_idx);
        fy_shift=k_idx_y(led_idx); 
        fy_range = (N_hr/2 - N/2 + 1 + fy_shift) : (N_hr/2 + N/2 + fy_shift);
        fx_range = (N_hr/2 - N/2 + 1 + fx_shift) : (N_hr/2 + N/2 + fx_shift);
        
        patch_FT_old=frequency_domain_hr_result(fx_range,fy_range);

        patch_FT_pupil=patch_FT_old.*Pupil;
        %逆变换回spatial domain
        g_xy_est=IF(patch_FT_pupil);
        %3.projection约束
        g_xy_new=measured_amplitude.*exp(1j*angle(g_xy_est)) ;
        patch_FT_updated=F(g_xy_new);

        %EPRY-1 迭代高分辨率频域Sn，实际上也就是在迭代更新frequency_domain_hr_result
        %定义pupil normalize项
        P_normal=conj(Pupil)/max(max(abs(Pupil.^2)));
        %定义diff_FT 作为\phi之差项
        diff_FT=patch_FT_updated-patch_FT_pupil;
        frequency_domain_hr_result(fx_range,fy_range)=frequency_domain_hr_result(fx_range,fy_range)+a_EPRY*P_normal.*diff_FT;


        
        %EPRY-2 迭代pupil，Sn对应patch_FT_old中的内容, \phi_n对应 patch_FT_updated
        %EPRY 中Sn都是直接存储在频域平面的量，phi_n是经过了pupil调制的
        %这里定义S_normal 作为Sn对应项
        S_normal=conj(patch_FT_old)/max(max(abs(patch_FT_old).^2));
        
        Pupil_update=b_EPRY*S_normal.*diff_FT;
        Pupil=Pupil+Pupil_update;
        %把pupil约束在圆内
        Pupil=Pupil.*Pupil_LR_initial;
        
        %SC grid search
        if iter_time>10
            search_range = [-1, 0, 1];
            min_error = inf;
            best_nx = 0;
            best_ny = 0;
            
            % 重新提取更新后的物体频谱（作为优化参考）
            current_obj_spectrum = frequency_domain_hr_result;
            
            for nx = search_range
                for ny = search_range
                    % 候选角度
                    cand_fx = fx_shift + nx;
                    cand_fy = fy_shift + ny;
                    
                    % 提取候选位置的频谱块
                    cand_fy_range = (N_hr/2 - N/2 + 1 + cand_fy) : (N_hr/2 + N/2 + cand_fy);
                    cand_fx_range = (N_hr/2 - N/2 + 1 + cand_fx) : (N_hr/2 + N/2 + cand_fx);
                    
                    % 确保范围在有效范围内
                    if min(cand_fy_range) < 1 || max(cand_fy_range) > N_hr || ...
                       min(cand_fx_range) < 1 || max(cand_fx_range) > N_hr
                        continue;
                    end
                    
                    % 提取候选频谱块
                    patch_FT_cand = current_obj_spectrum(cand_fx_range, cand_fy_range);
                    
                    % 用当前瞳孔函数调制
                    patch_FT_pupil_cand = patch_FT_cand .* Pupil;
                    
                    % 计算预测图像
                    g_xy_cand = IF(patch_FT_pupil_cand);
                    
                    % 计算误差（公式7）
                    error = norm(abs(g_xy_cand) - measured_amplitude, 'fro')^2;
                    
                    if error < min_error
                        min_error = error;
                        best_nx = nx;
                        best_ny = ny;
                    end
                end
            end
            
            % 3.2 更新角度（如果找到更优角度）
            if best_nx ~= 0 || best_ny ~= 0
                k_idx_x(led_idx) = k_idx_x(led_idx) + best_nx;
                k_idx_y(led_idx) = k_idx_y(led_idx) + best_ny;
            end
        end

    end

    %返回每一LED的最终error
    %fprintf('LED %d, final error: %.6f\n', led_idx, error);

end
%%
object_highresolution_recovery=IF(frequency_domain_hr_result);
amplitude=abs(object_highresolution_recovery);

figure('Name', 'FPM ');
subplot(1, 2, 1);
imagesc(amplitude);
axis image; 

title('Amplitude');

subplot(1, 2, 2);
imagesc(angle(object_highresolution_recovery));
axis image; 
title('Phase');

colormap gray; 
h = colorbar('southoutside'); % 在整个 figure 下方放一个 colorbar
h.Position = [0.25 0.05 0.5 0.03]; % 手动调整位置
