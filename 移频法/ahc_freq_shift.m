clear all
close all
clc

%[x,fs] = wavread('man.wav');
[x,fs] = audioread('man.wav');
%[x,fs] = audioread('timit.wav');
x = x(1:fs*4);   % notice!
%plot(x)
g = load('path.txt');

% ϣ�����ر任�����˲������������ȵõ���ͨ��Ȼ����Ƶ
h = fir2(200,[0,0.48,0.5,1],[1,1,0,0]); 
h = h(:); 
figure
freqz(h,1)
h = h.*exp(2*pi*5i*(1:length(h))'/4);
figure
freqz(h,1)

h_dummy = zeros(size(h)); h_dummy((end+1)/2) = 1;

K = 0.12;                                % ����

g = g(:);                               % ������ѧ·��g
c = [0,0,0,0,1]';                       % ����ϵͳ�ڲ�����·��c

xs1 = zeros(size(c));
xs2 = zeros(size(g));
xs3 = zeros(size(h_dummy));

y1 = zeros(size(x));                    % �ȷ���y1��y2�ռ䣬������������ʱ����ռ�ռ�ô�����������
y2 = zeros(size(x));
temp = 0;
yfb = zeros(size(x));

for i = 1:length(x)                     
    % ����γɷ�����·
    % ����������ǰ��·������������źŻ�����
    % �����µ�sampleֵ�洢�ڻ���������ǰ��
    % ���ָ��·�ʽ���Խ�����еĴ�λ�˻����ֱ�����
    xs1 = [x(i)+temp; xs1(1:end-1)];    % �ȴ���c������źŻ���
    
    % ǰ��·��
    y1(i) = K*(xs1'*c);                 % �������������ź�
    
    % ͨ��һ��ֻ��ʱ�ӵ��˲�����Ϊ�˽�y1��y2��Ⱥʱ��������ͬ
    xs3 = [y1(i); xs3(1:end-1)];        
    y1(i) = xs3' * h_dummy;
    
    y1(i) = min(1,y1(i));               % ����Լ����Х������ֽ�ֹ
    y1(i) = max(-1,y1(i));
    
    % ���㷴���ź�
    xs2 = [y1(i); xs2(1:end-1)];        % �ȴ���g������źŻ���
    % temp��Ϊ�����㻺�棬����һ�����㴦��ʱ�������źŻ��
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
f_shift = 3;                            % ��ƵƵ��Ϊ3Hz

for i = 1:length(x)
    xs1 = [x(i)+temp; xs1(1:end-1)];
    y2(i) = K*(xs1'*c);
    
    xs3 = [y2(i); xs3(1:end-1)];
    y2(i) = xs3' * h;                   % ͨ���˲����õ��ź�Ƶ�׵������Ჿ��
    y2(i) = y2(i)*exp(2*pi*1i*i/fs*f_shift);    % Ƶ��f_shift
    y2(i) = real(y2(i));                % ȡʵ�����ָ���Ƶ���ڸ����Ჿ�ֵ��ź�

    y2(i) = min(1,y2(i));
    y2(i) = max(-1,y2(i));
    xs2 = [y2(i); xs2(1:end-1)];
    temp = xs2'*g;
end

audiowrite('man_howling_suppression.wav',y2,fs);
%audiowrite('timit_howling_suppression.wav',y2,fs);
figure
subplot(211)
plot(y1)           % δ�����ź�
axis tight
legend('yout')
subplot(212)
plot(y2)           % Х�������ź�
axis tight
legend('freq\_shift')
