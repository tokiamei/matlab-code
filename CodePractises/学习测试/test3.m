%%%%%%%%%%%%求函数极小值
% x1=-50;x2=5;  
% yx=@(x)(sin(x).^2.*exp(-0.1*x)-0.5*sin(x).*(x+0.1));    %定义函数句柄
[xc0,fc0,exitflag]=fminbnd(yx,x1,x2) %找到的其中一个极小值 

ezplot(yx,[-50,5]);                     %ezplot同样可用于函数句柄，只需指定x的范围即可
xlabel('x'),grid on 

xx=[-23,-20,-18];                       %观察最小值疑似存在的两个区段
fc=fc0;xc=xc0;                          %暂时设立最小值点和最小值为初始搜索的结果
for k=1:2
 [xw,fw]=fminbnd(yx,xx(k),xx(k+1));  %分别计算两个区段的极小值
 if fw<fc, xc=xw;fc=fw;end           %若有更小的极小值则更换最小值
end
fprintf('函数最小值%6.5f发生在x=%6.5f处',fc,xc)