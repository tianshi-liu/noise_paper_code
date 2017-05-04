clear all
clc;
A=load('noisestack_full_deeptt');
N=1780;
nps=3600;
dt=1;
domega=pi/(nps*dt);
dsigma=domega;
omegamax1=0.045*2*pi;
omegamax2=0.055*2*pi;
omegamin1=0.002*2*pi;
omegamin2=0.004*2*pi;
Fc1=0.01;
Fc2=0.04;
n1=ceil(omegamax1/domega);
n2=floor(omegamax2/domega);
n1m=ceil(omegamin1/domega);
n2m=floor(omegamin2/domega);
for i=1:1:n2-n1
    omega=domega*(n1+i);
    win(i)=erf(-(omegamax2-omega)/(omegamax1-omega));
end
for i=1:1:n2m-n1m
    omega=domega*(n1m+i);
    winm(i)=erf(-(omegamin1-omega)/(omegamin2-omega));
end
tt=0:2*nps-1;
tt=tt*dt;
d0=1;
dd=0.1;
degree=d0+dd:dd:d0+N*dd;
Hd=seismo_filter(Fc1,Fc2);
for j=1:1:N
    fp(1:nps)=A(nps+1:2*nps,j)+1i*A(3*nps+1:4*nps,j);
    fp(nps+1:2*nps)=A(1:nps,j)+1i*A(2*nps+1:3*nps,j);
    fp(n1+1:n2)=fp(n1+1:n2).*win;
    fp(n2+1:nps)=0;
%     fp(n1m+1:n2m)=fp(n1m+1:n2m).*winm;
%     fp(1:n1m)=0;
%     fp(nps*2-n2m+1:nps*2-n1m)=fp(nps*2-n2m+1:nps*2-n1m).*winm(n2m-n1m:-1:1);
%     fp(nps*2-n1m+1:nps*2)=0;
    yp=ifft(fp);
    yp=real(yp);
    for i=1:1.5*nps
        yp(i)=yp(i)*exp(tt(i)*dsigma);
    end
    yp=yp(1.5*nps:-1:1);
    yp=filter(Hd,yp);
    yp=yp(1.5*nps:-1:1);
    yp=filter(Hd,yp);
    B(1:nps,j)=yp(1:nps);
end
bmax=max(max(B));
bmin=min(min(B));
if bmax+bmin>0
    bmin=-bmax;
else
    bmax=-bmin;
end
% reg=4;
% for i=1:1:nps
%     for j=1:1:N
%         if(B(i,j)>0)
%             B(i,j)=B(i,j)/bmax;
%             B(i,j)=real(B(i,j)^(1/reg));
%         else
%             B(i,j)=B(i,j)/bmin;
%             B(i,j)=-real(B(i,j)^(1/reg));
%         end
%     end
% end
rsat=500;
%rsat=700;
%rsat=300;
%rsat=3000;
for i=1:1:nps
    for j=1:1:N
        if(B(i,j)>bmax/rsat)
            B(i,j)=bmax/rsat;
        end
        if(B(i,j)<bmin/rsat)
            B(i,j)=bmin/rsat;
        end
    end
end
f=figure(2);
set(gcf,'PaperPositionMode','Manual');
set(gcf,'PaperUnits','inches');
set(gcf,'Papersize',[6,6]);
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition',[0.05,0.05,0.9,0.9]);
h=pcolor(degree,tt(1:nps)/60,B);
colormap(gray(1024));
set(h,'edgecolor','None');
hold on
% readtime('S','S_traveltime',1,94,24.5,'center','middle','a',6);
% readtime('ScS','S_traveltime',1,10,16,'center','top','c',12);
% readtime('ScS','S_traveltime',2,20,32,'center','top','d',6);
% readtime('S','S_traveltime',2,95,30.5,'center','middle','b',3);
% readtime('S','S_traveltime',1,95,24.5,'center','middle','a');
% readtime('ScS','S_traveltime',1,10,16,'center','top','c');
% readtime('ScS','S_traveltime',2,20,32,'center','top','d');
% readtime('S','S_traveltime',2,95,30.5,'center','middle','b');
% readtime('S','S_traveltime',1,58,17,'right','bottom','S(a)');
% readtime('ScS','S_traveltime',1,10,16,'center','top','ScS(c)');
% readtime('ScS','S_traveltime',2,20,32,'center','top','ScSScS(d)');
% readtime('SS','S_traveltime',1,114,34,'right','bottom','SS(b)');
% readtime('P','P_traveltime',1,80,10,'center','bottom','P');
% readtime('S','S_traveltime',1,58,16,'right','bottom','S');
% readtime('PP','P_traveltime',1,100,15,'right','bottom','PP');
% readtime('SS','S_traveltime',1,114,34,'right','bottom','SS');
% readtime('ScS','S_traveltime',1,10,16,'center','bottom','ScS');
% readtime('PcP','P_traveltime',1,10,8,'center','bottom','PcP');
% readtime('PKP','P_traveltime',1,160,22,'center','bottom','PKP');
% readtime('PKIKP','P_traveltime',1,170,18,'center','bottom','PKIKP');
%%%%%%%%%%%%%%%%%%%%%%%200km%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readtime('S','S_traveltime',1,95,24.5,'center','middle','a',6);
readtime('ScS','S_traveltime',1,10,16,'center','top','c',12);
readtime('ScS','S_traveltime',2,20,32,'center','top','d',6);
readtime('S','S_traveltime',2,95,30.5,'center','middle','b',3);
readdeeptime(1,1,'S','S','deep_traveltime',95,29,'center','middle','b${}^\prime$',3);
readdeeptime(4,2,'s','S','deep_traveltime',95,32.5,'center','middle','b${}^{\prime\prime}$',3);
readdeeptime(1,1,'ScS','ScS','deep_traveltime',20,29.5,'center','top','d${}^\prime$',6);
readdeeptime(1,1,'sScS','sScS','deep_traveltime',20,34,'center','top','d${}^{\prime\prime}$',6);
readdeeptime(3,1,'s','S','deep_traveltime',95,26.5,'center','middle','a${}^{\prime\prime}$',6);
readdeeptime(-1,1,'s','S','deep_traveltime',95,22.5,'center','middle','a${}^\prime$',6);
readdeeptime(1,1,'sScS','s','deep_traveltime',10,18.5,'center','top','c${}^{\prime\prime}$',1);
readdeeptime(1,-1,'ScS','s','deep_traveltime',10,14,'center','top','c${}^\prime$',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%300km%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readtime('S','S_traveltime',1,94,24.5,'center','middle','a',6);
% readtime('ScS','S_traveltime',1,10,16,'center','top','c',12);
% readtime('ScS','S_traveltime',2,20,32,'center','top','d',6);
% readtime('S','S_traveltime',2,95,30.5,'center','middle','b',3);
% readdeeptime(1,1,'S','S','ddeep_traveltime',95,28,'center','middle','b${}^\prime$',2);
% readdeeptime(4,2,'s','S','ddeep_traveltime',95,33,'center','middle','b${}^{\prime\prime}$',2);
% readdeeptime(1,1,'ScS','ScS','ddeep_traveltime',20,28.5,'center','top','d${}^\prime$',6);
% readdeeptime(1,1,'sScS','sScS','ddeep_traveltime',20,35,'center','top','d${}^{\prime\prime}$',6);
% readdeeptime(3,1,'s','S','ddeep_traveltime',94,27,'center','middle','a${}^{\prime\prime}$',6);
% readdeeptime(-1,1,'s','S','ddeep_traveltime',94,21,'center','middle','a${}^\prime$',6);
% readdeeptime(1,1,'sScS','s','ddeep_traveltime',10,19.5,'center','top','c${}^{\prime\prime}$',2);
% readdeeptime(1,-1,'ScS','s','ddeep_traveltime',10,13,'center','top','c${}^\prime$',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  readdeeptime('+','P','P','deep_traveltime');
%  readdeeptime('+','pP','pP','deep_traveltime');
% readdeeptime('+','ScS','ScS','deep_traveltime');
% readdeeptime('+','sScS','sScS','deep_traveltime');
%  readdeeptime('+','pP','p','deep_traveltime');
%  readdeeptime('-','P','p','deep_traveltime');
%  readdeeptime('+','pPcP','p','deep_traveltime');
%  readdeeptime('-','PcP','p','deep_traveltime');
% text('HorizontalAlignment','center','VerticalAlignment',...
 %            'bottom','String','(\uppercase\expandafter{\romannumeral4}) noise sources in whole Earth','Interpreter','latex','Position',[75,61],'FontSize',12);
% text('HorizontalAlignment','center','VerticalAlignment',...
%           'bottom','String','(\uppercase\expandafter{\romannumeral1}) noise sources on free surface','Interpreter','latex','Position',[75,61],'FontSize',12);
text('HorizontalAlignment','center','VerticalAlignment',...
           'bottom','String','(\uppercase\expandafter{\romannumeral2}) noise sources at 200 $$\mathrm{km}\;{}$$ deep','Interpreter','latex','Position',[75,61],'FontSize',12);
% text('HorizontalAlignment','center','VerticalAlignment',...
%           'bottom','String','(\uppercase\expandafter{\romannumeral3}) noise sources at 300 $$\mathrm{km}\;{}$$ deep','Interpreter','latex','Position',[75,61],'FontSize',12);
xlim([0 180]);
ylim([0 60]);
xlabel('$$\Delta (deg)$$','Interpreter','latex','FontSize',14);
ylabel('time(min)','Interpreter','latex','FontSize',14);
print(f,'lowq_noisestack_full_deeptt.eps','-depsc','-r200');
print(f,'lowq_noisestack_full_deeptt.png','-dpng','-r200');