
load("FPMsimulatonData.mat");
F = @(x) fftshift(fft2(ifftshift(x)));
IF = @(x) fftshift(ifft2(ifftshift(x)));


datacell_simulation=cell(1,size(measurements,3));
for i=1:size(measurements,3)
    datacell_simulation{i}=measurements(:,:,i);
end
allMin = inf;
allMax = -inf;
for i = 1:20
    currentMin = min(datacell_simulation{i}(:));
    currentMax = max(datacell_simulation{i}(:));
    allMin = min(allMin, currentMin);
    allMax = max(allMax, currentMax);
end

figure('Color','w');
tiledlayout(4,5,'TileSpacing','compact','Padding','compact');

for i = 1:20
    nexttile;
    imagesc(datacell_simulation{i});
    axis image off;
    colormap(gca,'gray');       % 灰度色图，可改为其他
    clim([allMin allMax]);     % 统一颜色范围
end

% 如果需要统一的颜色条，可以加上：
cb = colorbar('Location','eastoutside');
cb.Layout.Tile = 'east';
%%
% 预处理归一化
% 假设 A 是 M×N×K
% 找到每一层的最大值 (1×1×K)
layerMax = squeeze(max(max(measurements,[],1),[],2));

% 扩展成 M×N×K，方便逐层除法
L_measurements =  measurements./ reshape(layerMax, [1 1 size(measurements,3)]);



%%
%物domain \delta_x
ob_dx=pixelSize/ magnification;

%波数
k_0=2*pi/wavelength;
%把波数重定义成三角函数
alpha=kx;
beta=ky;
% real_kx=kx;
% real_ky=ky;
real_kx=k_0.*alpha;
real_ky=k_0.*beta;

%N.A.=a*sin\theta = n*sqrt(1-(k_0/kz)^2) 我们得到 在频域fx fy中的半径 令n=1
R_f=NA/wavelength;

N=size(measurements,1);
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

%高分辨率图像网格
max_kshift=max(sqrt(real_kx.^2+real_ky.^2));
max_f=max_kshift/(2*pi)+R_f;
up_sample_frate=ceil(max_f/(N/2*d_fx));
aaa=up_sample_frate;
up_sample_frate=max(4,up_sample_frate);

N_high_resolution=N*up_sample_frate;
ob_dx_high_resolution=ob_dx/up_sample_frate;
%初始化高分辨率 频域网格 全零
frequency_domain_hr_result=zeros(N_high_resolution,N_high_resolution);

%高分辨率不改变spatial domain L大小和d_fx d_fy 大小,因此，设置frequency domain B
B_f_h_r=N_high_resolution*d_fx;
%初始aperture
Pupil_LR = double(sqrt(FX.^2 + FY.^2) <= R_f);
%频域平移
dkx_LR = 2 * pi * d_fx;
dky_LR = 2 * pi * d_fy;
k_idx_x = round(real_kx / dkx_LR);
k_idx_y = round(real_ky / dky_LR);
%循环迭代次数
N_iter=100;

%初始化
[~,center_idx]=min(real_kx.^2 + real_ky.^2);
center_amplitude = sqrt(measurements(:,:,center_idx));
%total_amplitude=sum(sqrt(measurements,3));
%%
% figure('Name','Try sum(sqrt)');
% imagesc(sum(sqrt(measurements),3));
% axis image; colormap hsv; colorbar
% %%
% figure('Name','Try sqrt(sum)')
% imagesc(sqrt(sum(measurements,3)))
% axis image; colormap hsv; colorbar
% 
%%
patchFT_center = F(center_amplitude);
%patchFT_center=F(total_amplitude);
idx_hr_y = (N_high_resolution/2 - N/2 + 1) : (N_high_resolution/2 + N/2);
idx_hr_x = (N_high_resolution/2 - N/2 + 1) : (N_high_resolution/2 + N/2);
frequency_domain_hr_result(idx_hr_y, idx_hr_x) = patchFT_center .* Pupil_LR;

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
for iter_time=1:N_iter
    fprintf('迭代：%d /%d\n',iter_time,N_iter);
    idx_rand = randperm(N_LED);

    for p_idx =1:N_LED
        led_idx = idx_rand(p_idx);
        measured_amplitude = sqrt(measurements(:,:,led_idx));
        
        %波矢->频率 
        %fx_led=real_kx(led_idx)/(2*pi);
        %fy_led=real_ky(led_idx)/(2*pi);

        %构建aperture
        fx_shift=k_idx_x(led_idx);
        fy_shift=k_idx_y(led_idx); 
        fy_range = (N_high_resolution/2 - N/2 + 1 + fy_shift) : (N_high_resolution/2 + N/2 + fy_shift);
        fx_range = (N_high_resolution/2 - N/2 + 1 + fx_shift) : (N_high_resolution/2 + N/2 + fx_shift);
        patch_FT_old=frequency_domain_hr_result(fx_range,fy_range);
        patch_FT_pupil=patch_FT_old.*Pupil_LR;
        %逆变换回spatial domain
        g_xy_est=IF(patch_FT_pupil);
            
        %为观察收敛情况，定义error，并计算
        %g_xy_estamplitude=abs(g_xy_est);
        %error=sqrt(sum(sum((g_xy_estamplitude-measured_amplitude).^2)));
        %fprintf('LED %d, 迭代 %d, 误差: %.6f\n', led_idx, iter_time, error);

        %3.projection
        g_xy_new=measured_amplitude.*exp(1j*angle(g_xy_est)) ;
        
        
        %4. 下一次循环  高分辨率频谱叠加
        patch_FT_updated=F(g_xy_new);
        frequency_domain_hr_result(fx_range,fy_range)=patch_FT_old.*(1-Pupil_LR)+patch_FT_updated.*Pupil_LR;
        
    end

    %返回每一LED的最终error
    %fprintf('LED %d, final error: %.6f\n', led_idx, error);

end
%%
object_highresolution_recovery=IF(frequency_domain_hr_result);
figure('Name', 'FPM ');
subplot(1, 2, 1);
imagesc(abs(object_highresolution_recovery));
axis image; 
title('Amplitude');

subplot(1, 2, 2);
imagesc(angle(object_highresolution_recovery));
axis image; 
title('Phase');

colormap gray; 
h = colorbar('southoutside'); % 在整个 figure 下方放一个 colorbar
h.Position = [0.25 0.05 0.5 0.03]; % 手动调整位置
