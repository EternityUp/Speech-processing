function [yn,en,yfb,S] = HANLMSadapt(un,S)

% HANLMSadapt       NLMS Algorithm for Hearing Aids
%
%                   Perform over the entire length of input sequence for hearing
%                   aids. The history of output, square error and coefficients of 
%                   FIR filters are passed out to extenal. Refer to Fig. 1.9. 
%                   
% Parameters:
% un                Input signal
% S                 Adptive filter parameters as defined in HANLMSinit.m
% yn                History of output signal
% en                History of error signal
% yfb               History of feedback path signal
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

mu = S.step;                      % Step size of NLMS algorithm
w = S.coeffs;                     % Weight vector of FIR filter
M = length(S.coeffs);             % Length of FIR filter
fftap = zeros(1,M);               % Feedforward taps 
M_fb = length(S.fb);              % Length of feedback filter
fbtap = zeros(1,M_fb);            % Feedback taps
delay = S.delay;

ITER = length(un);                % Length of input sequence
yn = zeros(1,ITER);               % Initialize output sequence to zero
en = zeros(1,ITER);               % Initialize error sequence to zero
en(1) = un(1);                    %   and let the 1st sample as 1st input sample
yfb = zeros(1,ITER);              % Initialize output of feedback path to zero
u = zeros(1,ITER);                % Initialize input of adaptive filter to zero
u1 = zeros(1,ITER);               % Initialize input of feedforward path to zero
alpha = S.alpha;


for n = 1:ITER
    %% 将前向输出延迟一段时间赋值给u
    if n > (delay+1)
      u(n) = S.ff(end)*en(n-delay-1);
    else
      u(n) = 0;
    end
   
    %% 将最新数据置于fbtap的最前端，同时更新fbtap
    % 上述排布方式使得应用反馈路径响应时可以直接用点乘代替
    % 而不是写成显式的循环错位相乘
    for k = 1:1:M_fb-1            % Updates of tapped-delay line for feedback plant
       fbtap(M_fb-k+1) = fbtap(M_fb-k);
   end
   fbtap(1) = un(n);
   %% 应用反馈路径响应函数
   yfb(n) = S.fb'*fbtap';
   
   %if n >delay                   % Delayed input before adaptive filtering
   %    u(n) = u1(n-delay);
   %else 
   %    u(n) = 0;
   %end
   
   %% 类似地，更新前向路径缓冲区，将最新数据置于最前端
   for k = 1:1:M-1               % Updates of tapped-delay line of adaptive filter
       fftap(M-k+1) = fftap(M-k);
   end
   fftap(1) = un(n);
   %% 滤波，滤波器系数是作用在前向路径上
   yn(n) = w'*fftap';             % Compute filter output by inner product
   %% 计算滤波器输出的误差信号
   en(n) = un(n)+yfb(n)-yn(n);    % Compute error signal
   %% 滤波器权重更新
   w = w + ((mu*en(n))/(fftap*fftap'+alpha))*fftap'; % NLMS algorithm 
   S.iter = S.iter + 1;
end


S.coeffs = w;                    % Coefficient values at final iteration 
