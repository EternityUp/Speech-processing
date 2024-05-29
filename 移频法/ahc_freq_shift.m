clear all
close all
clc

%[x,fs] = wavread('man.wav');
[x,fs] = audioread('man.wav');
%[x,fs] = audioread('timit.wav');
x = x(1:fs*4);   % notice!
%plot(x)
g = load('path.txt');

% 希尔伯特变换所用滤波器，方法：先得到低通，然后移频
h = fir2(200,[0,0.48,0.5,1],[1,1,0,0]); 
h = h(:); 
figure
freqz(h,1)
h = h.*exp(2*pi*5i*(1:length(h))'/4);
figure
freqz(h,1)

h_dummy = zeros(size(h)); h_dummy((end+1)/2) = 1;

K = 0.12;                                % 增益

g = g(:);                               % 反馈声学路径g
c = [0,0,0,0,1]';                       % 扩音系统内部传递路径c

xs1 = zeros(size(c));
xs2 = zeros(size(g));
xs3 = zeros(size(h_dummy));

y1 = zeros(size(x));                    % 先分配y1和y2空间，避免运行中临时分配空间占用大量的运算量
y2 = zeros(size(x));
temp = 0;
yfb = zeros(size(x));

for i = 1:length(x)                     
    % 卷积形成反馈回路
    % 构建用于与前馈路径卷积的输入信号缓冲区
    % 将最新的sample值存储于缓冲区的最前端
    % 这种更新方式可以将卷积中的错位乘积变成直接相乘
    xs1 = [x(i)+temp; xs1(1:end-1)];    % 等待与c卷积的信号缓存
    
    % 前向路径
    y1(i) = K*(xs1'*c);                 % 馈给扬声器的信号
    
    % 通过一个只有时延的滤波器，为了将y1与y2的群时延条件相同
    xs3 = [y1(i); xs3(1:end-1)];        
    y1(i) = xs3' * h_dummy;
    
    y1(i) = min(1,y1(i));               % 幅度约束，啸叫则出现截止
    y1(i) = max(-1,y1(i));
    
    % 计算反馈信号
    xs2 = [y1(i); xs2(1:end-1)];        % 等待与g卷积的信号缓存
    % temp作为单样点缓存，待下一采样点处理时与输入信号混合
    temp = xs2'*g;  
    yfb(i) = temp;
end

audiowrite('man_howling.wav',y1,fs);
%audiowrite('timit_howling.wav',y1,fs);

figure
subplot 221
plot(x)
axis tight
legend('original')
subplot 222
plot(y1)
axis tight
legend('yout')
subplot 223
plot(yfb)
axis tight
legend('yfb')
subplot 224
plot(g)
axis tight
legend('g')


xs1 = zeros(size(c));
xs2 = zeros(size(g));
xs3 = zeros(size(h));
temp = 0;
f_shift = 3;                            % 移频频率为3Hz

for i = 1:length(x)
    xs1 = [x(i)+temp; xs1(1:end-1)];
    y2(i) = K*(xs1'*c);
    
    xs3 = [y2(i); xs3(1:end-1)];
    y2(i) = xs3' * h;                   % 通过滤波器得到信号频谱的正半轴部分
    y2(i) = y2(i)*exp(2*pi*1i*i/fs*f_shift);    % 频移f_shift
    y2(i) = real(y2(i));                % 取实部，恢复出频谱在负半轴部分的信号

    y2(i) = min(1,y2(i));
    y2(i) = max(-1,y2(i));
    xs2 = [y2(i); xs2(1:end-1)];
    temp = xs2'*g;
end

audiowrite('man_howling_suppression.wav',y2,fs);
%audiowrite('timit_howling_suppression.wav',y2,fs);
figure
subplot(211)
plot(y1)           % 未处理信号
axis tight
legend('yout')
subplot(212)
plot(y2)           % 啸叫抑制信号
axis tight
legend('freq\_shift')
