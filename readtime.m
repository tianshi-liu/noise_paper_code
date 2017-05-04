%clear all
%clc;
function readtime(phase,fn,mul,deglabel,timelabel,hor,ver,label,interv)
fileid=fopen(fn);
A=textscan(fileid,'%s %f %f %s %f %f %f %s');
%phase='PcP';
fclose(fileid);
nphase=length(A{1});
deg=A{1,2};
time=A{1,5};
phasename=A{1,4};
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
% for i=1:1:nphase
%     if(strcmp(phase,char(phasename{i,1})))
%         plot(deg(i),time(i),'.','MarkerSize',0.5);
%         hold on
%     end
% end
% hold off
if(strcmp(phase,'all'))
    plot(deg,time/60,'.','MarkerSize',0.01);
    hold on
elseif(strcmp(phase,'deep'))
    plot(deg*2,time*2/60,'.','MarkerSize',0.01);
    hold on
else
     ncount=0;
     for i=1:1:nphase
         if(strcmp(phase,char(phasename{i,1})))
             ncount=ncount+1;
             degout(ncount)=deg(i)*mul;
             if(degout(ncount)>180)
                 degout(ncount)=360-degout(ncount);
             end
             timeout(ncount)=time(i)*mul/60;
             %plot(deg(i),time(i)/60,'b.','MarkerSize',0.5);
             %hold on
         end
     end
     plot(degout(1:interv:length(degout)),timeout(1:interv:length(timeout)),'r.','MarkerSize',8);
     hold on
end
text('Interpreter','latex','HorizontalAlignment',hor,'VerticalAlignment',...
     ver,'String',label,'Position',[deglabel,timelabel],'FontSize',8,'Color','g');
hold on
%xlim([0 180]);
%ylim([0 60]);