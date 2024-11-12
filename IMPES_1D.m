%% 此程序基于<IMPES>方法，求解一维两相水驱油过程压力与含水饱和度随岩心长度的变化
% 作者：ZhengHe-upc
% 参考：《油藏数值模拟方法与应用》 第六章 第一节 一维油水两相水驱油的有限差分数值模拟方法
% 时间：2024年11月7日
% 1)模型假设：
%         1.油水两相互不混溶
%         2.等温渗流，符合达西渗流定律
%         3.一维流动，不考虑重力
%         4.流体和岩石不可压缩，不考虑毛管力pcow=0
%         5.油藏全区孔渗保持恒定，各向同性
% 2)初始条件：
%         岩心中饱和油和束缚水，在左端以恒定流量注水,注水为稳定驱替，注入、产出量均为Qv
%         p(x,0)=pi, Sw(x,0)=Swc
%         qv(x=0)=Qv, qv(x=L)=qvw+qvo=Qv

%% 初始化
clear; clc;

%% 输入参数，单位采取<cm,D,atm,s>单位制，参照<渗流力学>
% 时空参数
N    = 40;      % 一维油藏，一共有N个网格，[/]
L    = 100;     % 岩心长度，[cm]
nx   = L/N;     % 单个网格长度，[cm]
Tmax = 1100;    % 模拟最大时长，[s]
nt   = 10;      % 时间步长，[s]
% 岩心参数
A    = 10;      % 岩心横截面积，[cm2]
fai  = 0.3;     % 孔隙度，[/]
k    = 1;       % 渗透率，[um2]
Swc  = 0.2;     % 束缚水含水饱和度，[/]
% 流体参数
uo   = 2;       % 油相流体粘度，[mPa·s]
uw   = 1;       % 水相流体粘度，[mPa·s]
% 注入参数
Pi   = 1;       % 岩心的初始压力值，[10-1Mpa atm]
Qv   = 0.1;     % 注入流量值，[cm3/s]
% 油水相渗数据
Sat  = [0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80];
krw0 = [0.00,0.03,0.06,0.10,0.14,0.17,0.27,0.35,0.42,0.52,0.65,0.79,0.90];
kro0 = [0.85,0.75,0.62,0.49,0.31,0.19,0.14,0.10,0.07,0.05,0.03,0.01,0.00];

% 绘制油水相对渗透率曲线
figure;
plot(Sat, krw0, 'b-', 'LineWidth', 2);
hold on;
plot(Sat, kro0, 'r-', 'LineWidth', 2);
xlabel('Water Saturation');
ylabel('Relative Permeability');
title('Oil - Water Relative Permeability Curves');
legend('Original Water', 'Original Oil')
axis([0.2 0.8 0 1]);
grid on;
 
%% 初始化变量，提前分配内存空间
krw = 0;                    % 水相相对渗透率
kro = 0;                    % 油相相对渗透率
fc  = zeros(1, N-1);        % 油和水流动系数之和
fcw = zeros(1, N-1);        % 水相流动系数
fco = zeros(1, N-1);        % 油相流动系数
Sw  = Swc*ones(1, N);       % 储存含水饱和度
P   = zeros(1, N);          % 储存压力
a   = zeros(1, N);          % 压力方程系数矩阵的下三角
b   = zeros(1, N);          % 压力方程系数矩阵的主对角线
c   = zeros(1, N);          % 压力方程系数矩阵的上三角
d   = zeros(1, N);          % 压力方程的右边项

%% 主循环
for t = 0: nt: Tmax  
    
    % 计算不同Sw下的水相、油相以及总的流动系数fcw、fco、fc；使用MATLAB自带的插值函数<interp1>求取渗透率值
    for i = 1: N-1
        krw    = interp1(Sat, krw0, Sw(i), 'linear');
        kro    = interp1(Sat, kro0, Sw(i), 'linear');
        fcw(i) = k * krw / uw;  
        fco(i) = k * kro / uo;
        fc(i)  = fcw(i) + fco(i);
    end
    
    % 构造压力系数矩阵及右边项
    % 构造下三角
    a = [0, fc(1:N-2), 1];
    % 构造主对角
    b = [1, -(fc(1:N-2) + fc(2:N-1)), -1];
    % 构造上对角
    c = [-1, fc(2:N-1), 0];
    % 构造右边项
    d = zeros(1, N);
    d(1) = nx * Qv / A / fc(1);
    d(N) = nx * Qv / A / fc(N-1);
    
    % 使用追赶法求解压力
    P = solve_tridiagonal_matrix(a, b, c, d);
    
    % 显式求饱和度
    % 第1个网格处饱和度
    Sw(1) = Sw(1) + nt * (Qv / A - fcw(1) * (P(1) - P(2)) / nx) / (fai * nx);
    % 第2：N-1处网格饱和度
    for i = 2:N-1
        Sw(i) = Sw(i) + nt * (fcw(i) * (P(i+1) - P(i)) - fcw(i-1) * (P(i) - P(i-1))) / (fai * nx * nx);
    end
    % 第N个网格处饱和度
    Sw(N) = Sw(N) + nt * (fcw(N-1) * (P(N-1) - P(N)) / nx - Qv * fcw(N-1) / (A * (fco(N-1) + fcw(N-1)))) / (fai * nx);
    
    % 每隔100步保存一次压力P、饱和度数据Sw，以.mat形式保存
    if mod(t, 100) == 0 
        save(['data_t_', num2str(t), '.mat'], 'P', 'Sw');
    end
        
    % 求解水的突破时间，需要增加最大迭代步长
    if (Sw(N) ~= Swc)
        disp(t)
    end
end

%% 结果可视化
t_values = 0: 100: Tmax;                     % 创建一个包含所有要绘制的时间步值的向量
num_timesteps = length(t_values);            % 获取时间步的数量、用于控制循环
P_all = cell(num_timesteps, 1);              % 创建元胞数组P_all，用于存储不同时间步下的压力数据
Sw_all = cell(num_timesteps, 1);             % 创建元胞数组Sw_all，用于存储不同时间步下的饱和度数据
xx = linspace(1.25, 98.75, N);               % 生成横轴坐标，N是网格数，注意xx要在每一个网格的中心位置，符合物理场景

for i = 1:num_timesteps
    load(['data_t_', num2str(t_values(i)), '.mat']);        % 加载对应时间步的数据文件
    P_all{i} = P;
    Sw_all{i} = Sw;
end

% 绘制压力变化图
figure;
for i = 1:num_timesteps
    plot(xx, P_all{i} ,'linewidth', 1.5);
    hold on;
end
xlabel('Length x [cm]');
ylabel('Pressure P [atm]');
title('Pressure along length x at different time steps');
legend(arrayfun(@(t) ['t = ', num2str(t)], t_values, 'UniformOutput', false), 'FontSize', 12, 'FontName', 'Times New Roman');

% 绘制含水饱和度变化图
figure;
for i = 1:num_timesteps
    plot(xx, Sw_all{i},'linewidth', 1.5);
    hold on;
end
ylim([0.18 0.80])
xlabel('Length x [cm]');
ylabel('Water Saturation [/]');
title('Water Saturation along length x at different time steps');
legend(arrayfun(@(t) ['t = ', num2str(t)], t_values, 'UniformOutput', false), 'FontSize', 12, 'FontName', 'Times New Roman');
