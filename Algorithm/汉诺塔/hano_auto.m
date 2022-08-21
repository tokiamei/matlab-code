function hano_auto(level)
if nargin < 1
   level=6;
end
hold on,axis equal
axis(0.5+[0,120,0,60])
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
set(gca,'color','k')
ground=[120.5,1.5;120.5,0.5;0.5,0.5;0.5,1.5];
column_1=[21,1.5;21,44;19,44;19,1.5];
column_2=[61,1.5;61,44;59,44;59,1.5];
column_3=[101,1.5;101,44;99,44;99,1.5];
step=0;
arrow=[1];k=1.5;control=1;
A=[];B=[];C=[];
A=level:-1:1;
block_1=@(l,h)[20+2*l,3*h-1.5;20+2*l,3*h+1.5;20-2*l,3*h+1.5;20-2*l,3*h-1.5];
block_2=@(l,h)[60+2*l,3*h-1.5;60+2*l,3*h+1.5;60-2*l,3*h+1.5;60-2*l,3*h-1.5];
block_3=@(l,h)[100+2*l,3*h-1.5;100+2*l,3*h+1.5;100-2*l,3*h+1.5;100-2*l,3*h-1.5];
block__=@(l,p)[p+2*l,45-1.5;p+2*l,45+1.5;p-2*l,45+1.5;p-2*l,45-1.5];
set(gcf,'tag','co','CloseRequestFcn',@clo);
    function clo(~,~)
        control=0;
        delete(findobj('tag','co'));
        clf
        close
    end 


    function draw(~,~)
        delete(findobj(gcf,'type','text'))
        delete(findobj(gcf,'type','patch'))
        text(5,-4,'step=')
        text(20,-4,num2str(step))
        text(5,-10,'level=')
        text(20,-10,num2str(level))
        text(34,-4,'use the key(\leftarrow \rightarrow) to move the claw')
        text(34,-10,'use the key(\uparrow \downarrow) to pick or land the block')
        fill(ground(:,1),ground(:,2),'w','EdgeColor','none')
        fill(column_1(:,1),column_1(:,2),'w','EdgeColor','none')
        fill(column_2(:,1),column_2(:,2),'w','EdgeColor','none')
        fill(column_3(:,1),column_3(:,2),'w','EdgeColor','none')
        text(40*arrow(1)-20-2,55,'\downarrow','Color','w','Fontsize',25)
        if length(arrow)==2
            p=40*arrow(1)-20;
            block=block__(arrow(2),p);
            fill(block(:,1),block(:,2),'w','EdgeColor','none')  
        end
        for i=1:length(A)
            block=block_1(A(i),k*i);
            fill(block(:,1),block(:,2),'w','EdgeColor','none')   
        end
        for i=1:length(B)
            block=block_2(B(i),k*i);
            fill(block(:,1),block(:,2),'w','EdgeColor','none')   
        end
        for i=1:length(C)
            block=block_3(C(i),k*i);
            fill(block(:,1),block(:,2),'w','EdgeColor','none')   
        end
        if length(C)==level
            pause(0.2)
            buttonName1=questdlg('Congratulations! You win','you win','close','restart','next level','close');
            if isempty(buttonName1)
                buttonName1='end';
            end
            if strcmp(buttonName1,'restart')
                hano_auto(level)
            end
            if strcmp(buttonName1,'close')
                close;
            end
            if strcmp(buttonName1,'next level')
                if level<10
                    hano_auto(level+1)
                end
                if level==10
                    msgbox('it is already the highest level')
                end
            end
        end
    end

   function key(event)
   switch event
       case 'leftarrow'
           if arrow(1)-1~=0
               arrow(1)=arrow(1)-1;step=step+1;
           end
       case 'rightarrow'
           if arrow(1)-3~=0
               arrow(1)=arrow(1)+1;step=step+1;
           end
       case 'uparrow'
           if length(arrow)==1
               if (arrow(1)==1)&&(~isempty(A))
                   arrow=[arrow(1),A(end)];
                   A(end)=[];
               end
               if (arrow(1)==2)&&(~isempty(B))
                   arrow=[arrow(1),B(end)];
                   B(end)=[];
               end
               if (arrow(1)==3)&&(~isempty(C))
                   arrow=[arrow(1),C(end)];
                   C(end)=[];
               end
           end
       case 'downarrow'
           if length(arrow)==2
               if (arrow(1)==1)
                   if isempty(A)
                       A=[A,arrow(2)];
                       arrow(2)=[];
                   else
                       if A(end)>arrow(2)
                           A=[A,arrow(2)];
                           arrow(2)=[];
                       end
                   end
               end
               if (arrow(1)==2)
                   if isempty(B)
                       B=[B,arrow(2)];
                       arrow(2)=[];
                   else
                       if B(end)>arrow(2)
                           B=[B,arrow(2)];
                           arrow(2)=[];
                       end
                   end           
               end
               if (arrow(1)==3)
                   if isempty(C)
                       C=[C,arrow(2)];
                       arrow(2)=[];
                   else
                       if C(end)>arrow(2)
                           C=[C,arrow(2)];
                           arrow(2)=[];
                       end
                   end
               end
           end
   end
   pause(0.2)
   if control==1
       draw()
   else
       pause
   end
   end
draw()
matrix=[1 2 3];
digui(level,matrix)
    function digui(n,structdigui)
            if n==1
                move(structdigui(1),structdigui(3));
            else
                digui(n-1,[structdigui(1),structdigui(3),structdigui(2)])
                move(structdigui(1),structdigui(3))
                digui(n-1,[structdigui(2),structdigui(1),structdigui(3)])
            end
    end
    function move(a,b)
        while arrow(1)~=a
            if arrow(1)>a
                key('leftarrow')
            end
            if arrow(1)<a
                key('rightarrow')
            end
        end
        key('uparrow')
        while arrow(1)~=b
            if arrow(1)>b
                key('leftarrow')
            end
            if arrow(1)<b
                key('rightarrow')
            end
        end
        key('downarrow')
    end
        

end
