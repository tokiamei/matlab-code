function p = draw(flag)
t = (0:50)/50*2*pi;
x = sin(t);
y = cos(t);
p = @cirline;
feval(p,flag,x,y,t);
end

function cirline(wd,x,y,t)
switch wd
    case 'line'
        y=x+y;
    case 'circle'
        plot(x,y);
    otherwise
        error('输入的只能是“line”或“circle”');
end
end