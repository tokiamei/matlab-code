%主函数
function [a,b] = myFirstCode(x,y)   %此处的a,b就是输出形参；x,y就是输入形参
a = mySub1(x,y);                    %myMain就是函数名称，一般也是文件名
b = mySub2(x,y);
end
%子函数
function z = mySub1(x,y)
z = x - y;
end
%又一个子函数
function z = mySub2(x,y)
z =x + y;
end  