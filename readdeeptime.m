function readdeeptime(mul1,mul2,phase1,phase2,fn,deglabel,timelabel,hor,ver,label,interv)
fileid=fopen(fn);
A=textscan(fileid,'%s %f %f %s %f %f %f %s');
%phase='PcP';
fclose(fileid);
nphase=length(A{1});
deg=A{1,2};
time=A{1,5};
phasename=A{1,4};
p=A{1,6};
i=1;
while(i<nphase)
    j=i+1;
    while((deg(j)<=(deg(i)+0.001))&&(strcmp(char(phasename{i,1}),char(phasename{j,1}))))
        phasename{j,1}='q';
        j=j+1;
        if(j>nphase)
            break;
        end
    end
    i=j;
end
n1=0;
n2=0;
for i=1:1:nphase
     if(strcmp(phase1,char(phasename{i,1})))
         %plot(deg(i),time(i)/60,'k.','MarkerSize',0.5);
         %hold on
         n1=n1+1;
         deg1(n1)=deg(i);
         time1(n1)=time(i);
         p1(n1)=p(i);
     end
     if(strcmp(phase2,char(phasename{i,1})))
         %plot(deg(i),time(i)/60,'k.','MarkerSize',0.5);
         %hold on
         n2=n2+1;
         deg2(n2)=deg(i);
         time2(n2)=time(i);
         p2(n2)=p(i);
     end
end

for j=1:interv:n2
    np1=0;
    np2=0;
    pmin=min(p1);
    pmax=max(p1);
    for i=1:1:n1
        if ((p1(i)<=p2(j))&&(p1(i)>pmin))
            np1=i;
            pmin=p1(i);
        end
        if ((p1(i)>=p2(j))&&(p1(i)<pmax))
            np2=i;
            pmax=p1(i);
        end
    end
    if((np1>0)&&(np2>0))
        deg0=deg1(np1);
        time0=time1(np1);
        if(abs(p1(np1)-p1(np2))>1e-3)
            deg0=deg0+(deg1(np1)-deg1(np2))*(p2(j)-p1(np2))/(p1(np1)-p1(np2));
            time0=time0+(time1(np1)-time1(np2))*(p2(j)-p1(np2))/(p1(np1)-p1(np2));
        end
        degdraw=deg0*mul1+deg2(j)*mul2;
        if(degdraw>180)
            degdraw=360-degdraw;
        end
        timedraw=time0*mul1+time2(j)*mul2;
        timedraw=timedraw/60;
        plot(degdraw,timedraw,'r.','MarkerSize',8);
        hold on
%         if(strcmp(sign,'+'))
%             plot(deg0+deg2(j),(time0+time2(j))/60,'r.','MarkerSize',4);
%             hold on
%         else
%             plot(deg0-deg2(j),(time0-time2(j))/60,'r.','MarkerSize',4);
%             hold on
%         end
    end
end
text('Interpreter','latex','HorizontalAlignment',hor,'VerticalAlignment',...
        ver,'String',label,'Position',[deglabel,timelabel],'FontSize',8,'Color','g');
hold on
%xlim([0 180]);
%ylim([0 60]);