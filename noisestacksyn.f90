program main
    integer i,j,N,nps
    complex(8), allocatable :: uf(:,:),y(:)
    real(8) degree
    character(100) uffn,form
    N=1780
    !N=2
    nps=3600
    allocate(uf(1:2*nps,1:N))
    do i=1,N
        degree=1.0+0.1*i
        allocate(y(1:2*nps))
        call syn(degree,y,nps)
        do j=1,2*nps
            uf(j,i)=y(j)
        enddo
        deallocate(y)
        write(*,*)degree
    enddo
    uffn='noisestack_full_isott'
    open(unit=1,file=uffn,form='formatted')
    write(form,*)N
    form='('//trim(form)//'E15.4E5)'
    do j=1,2*nps
        write(unit=1,fmt=form)(dreal(uf(j,i)),i=1,N)
    enddo
    do j=1,2*nps
        write(unit=1,fmt=form)(dimag(uf(j,i)),i=1,N)
    enddo
    close(unit=1)
    end
    
subroutine syn(degree1,ufp,nps)
    implicit none
    character(50) modelfn,eigenfn_S,eigenfn_T,glfn,ititle
    integer ifanis,tref,ifdeck,n,nic,noc,ngrid,jcom,normin,normax,lmin,lmax,nmode,length,i,j,lm,l,nps,k,ncurr,n1,n2
    integer nn,ll,m,lcurr,ni,ns
    pointer (pnn, nn)
    pointer (pll, ll)
    real(8), allocatable :: gl(:,:,:),glt(:,:) !time*modeorder(1~6)*direction(1~3)
    real(4), allocatable :: rec(:)
    real(8) degp,domega,rhon,rn,c,degree1,pi,omegamax,omegamin,dt,depth
    real(8) r(200),rho(200),vpv(200),vsv(200),qkappa(200),qshear(200),vph(200),vsh(200),eta(200)
    real(4) omega,q,gc,ur,vr,wr,us,vs,ws
    real(8), allocatable :: omeganl(:),qnl(:),unl(:),vnl(:),wnl(:),usnl(:),vsnl(:),wsnl(:)
    complex(8) ufp(1:2*nps)
    complex(8) ampsum
    !real(8) momtens(1:3,1:3)
    !!open
    
    rhon=5515
    rn=6371000
    pi=3.1415926535
    depth=300000
    c=(1e25)/(rhon*(rn**3))
    c=sqrt(c)
    omegamax=0.4*2*pi
    omegamin=0.001*2*pi
    normin=0
    normax=510
    modelfn='prem_noocean.txt'
    open(unit=1,file=modelfn,form='formatted')
    !eigenfn_S='enoocean_resort_S'
    !eigenfn_S='eprem_fund_S'
    eigenfn_S='epremfull_resort_S'
    open(unit=2,file=eigenfn_S,form='unformatted')
    !eigenfn_T='enoocean_resort_T'
    !eigenfn_T='eprem_fund_T'
    eigenfn_T='epremfull_resort_T'
    open(unit=3,file=eigenfn_T,form='unformatted')
    !!open done
    !!input
    dt=1
    !nps=6000
    domega=2*pi/(dt*nps*2.0)
    read(1,*) ititle
    read(1,*) ifanis,tref,ifdeck
    read(1,*) n,nic,noc
    read(1,*) (r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)
    ns=n
    do while(r(ns)>(r(n)-depth))
        ns=ns-1
    enddo
    if(depth*2<(2*r(n)-r(ns)-r(ns+1)))then
        ns=ns+1
    endif
    do i=1,2*nps
        ufp(i)=0      
    enddo
    allocate(omeganl(normin:normax))
    allocate(qnl(normin:normax))
    allocate(unl(normin:normax))
    allocate(vnl(normin:normax))
    allocate(wnl(normin:normax))
    allocate(usnl(normin:normax))
    allocate(vsnl(normin:normax))
    allocate(wsnl(normin:normax))
    
    
        

    !!input legendre
    degp=degree1/180*pi
    lmin=0
    lmax=3000
    !normin=0
    !normax=300
    allocate(gl(lmin:lmax,-1:1,-1:1))
    if(dabs(degp)<1e-5)then
        do l=lmin,lmax,1
            do m=-1,1,1
                do k=-1,1,1
                    gl(l,m,k)=0
                enddo
            enddo
            gl(l,0,0)=1
            gl(l,1,1)=1
            gl(l,-1,-1)=1
        enddo
            
    else if(dabs(degp-pi)<1e-5)then
        do l=lmin,lmax,1
            do m=-1,1,1
                do k=-1,1,1
                    gl(l,m,k)=0
                enddo
            enddo
            if(mod(l,2)==0)then
                gl(l,0,0)=-1
                gl(l,-1,1)=1
                gl(l,1,-1)=1
            else
                gl(l,0,0)=1
                gl(l,-1,1)=-1
                gl(l,1,-1)=-1
            endif
        enddo
        
            
    else
        
      do l=lmin,lmax,1
        allocate(glt(-l:l,-1:1))
        call GL_generator_MR(l,glt,degp,-1,1)
!        write(*,*)l
        do m=-1,1,1
            do k=-1,1,1
                if(abs(m)>abs(l))then
                    gl(l,m,k)=0
                else
                    
                   gl(l,m,k)=glt(m,k)
                endif
                
            enddo
        enddo
        deallocate(glt)
      enddo
    endif
    
    !!
    lcurr=-1
    ncurr=0
    length=n*6+5
    allocate(rec(1:length))
    do while(.not.eof(2))    
        read(unit=2)(rec(i),i=1,length)
        pnn=loc(rec(1))
        pll=loc(rec(2))
        omega=rec(3)
        q=rec(4)
        ur=rec(5+n)*c
        us=rec(5+ns)*c
        vr=rec(5+n*3)*c
        vs=rec(5+n*2+ns)*c
        wr=0
        ws=0
        if(ll>lcurr)then
            do n1=1,ncurr
                ampsum=0
                do n2=1,ncurr
                    call ampsumcal(lcurr,omeganl(n1),qnl(n1),unl(n1),vnl(n1),wnl(n1),usnl(n1),vsnl(n1),wsnl(n1),omeganl(n2),qnl(n2),unl(n2),vnl(n2),wnl(n2),usnl(n2),vsnl(n2),wsnl(n2),gl,lmax,ampsum)   
                enddo
                call add(lcurr,ampsum,omeganl(n1),qnl(n1),domega,nps,ufp)
            enddo            
            lcurr=ll
            ncurr=0
        endif
        if((omega<omegamax).and.(omega>omegamin))then
            ncurr=ncurr+1
            omeganl(ncurr)=omega
            qnl(ncurr)=q
            unl(ncurr)=ur
            usnl(ncurr)=us
            vnl(ncurr)=vr
            vsnl(ncurr)=vs
            wnl(ncurr)=wr
            wsnl(ncurr)=ws
        endif
    enddo
    deallocate(rec)
    
    !!
    lcurr=-1
    ncurr=0
    length=n*2+5
    allocate(rec(1:length))
    do while(.not.eof(3))      
        read(unit=3)(rec(i),i=1,length)
        pnn=loc(rec(1))
        pll=loc(rec(2))
        omega=rec(3)
        q=rec(4)
        ur=0
        us=0
        vr=0
        vs=0
        wr=rec(5+n)*c
        ws=rec(5+ns)*c
        if(ll>lcurr)then
            do n1=1,ncurr
                ampsum=0
                do n2=1,ncurr
                    call ampsumcal(lcurr,omeganl(n1),qnl(n1),unl(n1),vnl(n1),wnl(n1),usnl(n1),vsnl(n1),wsnl(n1),omeganl(n2),qnl(n2),unl(n2),vnl(n2),wnl(n2),usnl(n2),vsnl(n2),wsnl(n2),gl,lmax,ampsum)                 
                enddo
                call add(lcurr,ampsum,omeganl(n1),qnl(n1),domega,nps,ufp)
            enddo            
            lcurr=ll
            ncurr=0
        endif
        if((omega<omegamax).and.(omega>omegamin))then
            ncurr=ncurr+1
            omeganl(ncurr)=omega
            qnl(ncurr)=q
            unl(ncurr)=ur
            usnl(ncurr)=us
            vnl(ncurr)=vr
            vsnl(ncurr)=vs
            wnl(ncurr)=wr
            wsnl(ncurr)=ws
        endif 
    enddo
    deallocate(rec)
    deallocate(gl)
    deallocate(omeganl)
    deallocate(qnl)
    deallocate(unl)
    deallocate(vnl)
    deallocate(wnl)
    deallocate(usnl)
    deallocate(vsnl)
    deallocate(wsnl)
    !close(unit=4)
    !write(*,*)nmode
    !!
    close(unit=1)
    close(unit=2)
    close(unit=3)
    !print *, 'input name of uz file'
    !read(*,*)uzfn
    !print *, 'input name of ur file'
    !read(*,*)urfn
    !print *, 'input name of ut file'
    !read(*,*)utfn
    !print *, 'input name of fundamental mode file'
    !read(*,*)uffn
    !open(unit=1,file=uffn,form='formatted')
    
    !do j=1,3
    !    do k=1,3
    !        write(unit=1)(uf(i,j,k),i=1,nps)
    !    enddo
    !enddo
    !do i=1,2*nps
    !    write(unit=1,fmt='(9E15.4E5)')((dreal(ufp(i,j,k)),k=1,3),j=1,3)
    !enddo
    !do i=1,2*nps
    !    write(unit=1,fmt='(9E15.4E5)')((dimag(ufp(i,j,k)),k=1,3),j=1,3)
    !enddo
    !do i=1,2*nps
    !    write(unit=1,fmt='(9E15.4E5)')((dreal(ufn(i,j,k)),k=1,3),j=1,3)
    !enddo
    !do i=1,2*nps
    !    write(unit=1,fmt='(9E15.4E5)')((dimag(ufn(i,j,k)),k=1,3),j=1,3)
    !enddo
        
    !close(unit=1)   
    !deallocate(ufp)
    !deallocate(ufn)
    
    end
    
!subroutine add(ll,omega1,q1,u1,v1,w1,us1,vs1,ws1,omega2,q2,u2,v2,w2,us2,vs2,ws2,gl,domega,nps,lmax,ufp)
!    integer ll,nps,lmax,i,j,k,nspec,iomegac1,iomegac2
!    real(8) domega,pi,omega1,omega2,q1,q2,u1,u2,v1,v2,w1,w2,us1,us2,vs1,vs2,ws1,ws2,sigma1,sigma2,sigma0,amp,tw
!    real(8) gl(0:lmax,-1:1,-1:1)
!    real(8), allocatable :: comp(:,:) 
!    complex(8) ufp(-nps:nps-1)
!    complex(8) c1,c2,tp,tn,const
!    real(8) omegan2,rhon,g
!    !complex(8) s(-2:2),a(1:3)
!    !complex(8) e1,e2
!    rhon=5515
!    g=6.6723e-11
!    pi=3.1415926535
!    nspec=500
!    tw=3600*4
!    iomegac1=int(omega1/domega)
!    iomegac2=-int(omega2/domega)
!    omegan2=rhon*pi*g
!    sigma1=0.5*omega1/q1
!    sigma2=0.5*omega2/q2
!    sigma0=domega
!    allocate(comp(1:3,1:3))
!    amp=us1*us2+vs1*vs2+ws1*ws2
!    !if (abs(omega1-omega2)<1e-10)then
!    !    amp=1/(omega1*omega2)*omegan2
!    !else
!    !    amp=0
!    !endif
!    comp(1,1)=u1*u2*gl(ll,0,0)
!    comp(1,2)=u1*v2*gl(ll,-1,0)
!    comp(1,3)=-u1*w2*gl(ll,-1,0)
!    comp(2,1)=v1*u2*gl(ll,-1,0)
!    comp(2,2)=v1*v2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5-w1*w2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5
!    comp(2,3)=-v1*w2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5-w1*v2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5
!    comp(3,1)=-w1*u2*gl(ll,-1,0)
!    comp(3,2)=-w1*v2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5-v1*w2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5
!    comp(3,3)=-v1*v2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5+w1*w2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5
!    !const=((sigma1+sigma2)*(sigma1+sigma2)+(omega1-omega2)*(omega1-omega2))*((sigma1+sigma2)*(sigma1+sigma2)+(omega1+omega2)*(omega1+omega2))
!    !c1=dcmplx((sigma1+sigma2)*2*omega1,-(sigma1+sigma2)*(sigma1+sigma2)+omega1*omega1-omega2*omega2)/const*omega1*omega2*omega2
!    !c2=dcmplx((sigma1+sigma2)*2*omega2,-(sigma1+sigma2)*(sigma1+sigma2)+omega2*omega2-omega1*omega1)/const*omega2*omega1*omega1
!    !c1=dcmplx((sigma1+sigma2)*2*omega1,-(sigma1+sigma2)*(sigma1+sigma2)+omega1*omega1-omega2*omega2)/const*omega1*omega1*omega1*omega2*omega2
!    !c2=dcmplx((sigma1+sigma2)*2*omega2,-(sigma1+sigma2)*(sigma1+sigma2)+omega2*omega2-omega1*omega1)/const*omega2*omega2*omega2*omega1*omega1
!    const=0.5*((1-dexp(-(sigma1+sigma2)*tw)*dcmplx(dcos((omega1-omega2)*tw),-dsin((omega1-omega2)*tw)))/dcmplx(sigma1+sigma2,omega1-omega2)-(1-dexp(-(sigma1+sigma2)*tw)*dcmplx(dcos((omega1+omega2)*tw),-dsin((omega1+omega2)*tw)))/dcmplx(sigma1+sigma2,omega1+omega2))
!    !const=0.5*((1-dexp(-(sigma1+sigma2)*tw)*dcmplx(dcos((omega1-omega2)*tw),dsin((omega1-omega2)*tw)))/dcmplx(sigma1+sigma2,omega2-omega1)+(1-dexp(-(sigma1+sigma2)*tw)*dcmplx(dcos((omega1+omega2)*tw),dsin((omega1+omega2)*tw)))/dcmplx(sigma1+sigma2,-omega1-omega2))
!    c1=const*omega1*omega2
!    !c1=const*omega1*omega2*omega1*omega2
!    if (iomegac1+nspec<nps-1)then
!        do i=-nspec,nspec
!            omega=domega*(i+iomegac1)
!            tp=c1/dcmplx((sigma1+sigma0),omega-omega1)
!            !do j=1,3
!            !    do k=1,3
!                    ufp(i+iomegac1)=ufp(i+iomegac1)+(2.0*ll+1.0)*tp*amp*comp(1,1)/(4*pi*omegan2*omegan2)
!            !    enddo
!            !enddo
!        enddo
!    endif
!    !if (iomegac2-nspec>-nps)then
!    !    do i=-nspec,nspec
!    !        omega=domega*(i+iomegac2)
!    !        tn=c2/dcmplx((sigma2+sigma0),-omega-omega2)
!    !        do j=1,3
!    !            do k=1,3
!    !                ufn(i+iomegac2,j,k)=ufn(i+iomegac2,j,k)+(2.0*ll+1.0)*tn*amp*comp(j,k)/(4*pi*omegan2*omegan2)
!    !            enddo
!    !        enddo
!    !    enddo
!    !endif
!    deallocate(comp)
!    return
!    end
subroutine ampsumcal(ll,omega1,q1,u1,v1,w1,us1,vs1,ws1,omega2,q2,u2,v2,w2,us2,vs2,ws2,gl,lmax,ampsum)
    integer ll,lmax,i,j,k,itype
    real(8) omega1,omega2,q1,q2,u1,u2,v1,v2,w1,w2,us1,us2,vs1,vs2,ws1,ws2,sigma1,sigma2,const,amp,tw
    real(8) gl(0:lmax,-1:1,-1:1)
    real(8), allocatable :: comp(:,:) 
    real(8) c1,c2,rhon,g,pi,omegan2
    complex(8) consttemp,tp,ampsum
    itype=3
    tw=3600*4
    rhon=5515
    g=6.6723e-11
    pi=3.1415926535
    omegan2=rhon*pi*g
    sigma1=0.5*omega1/q1
    sigma2=0.5*omega2/q2
    allocate(comp(1:3,1:3))
    !amp=us1*us2+vs1*vs2+ws1*ws2
    amp=1/(omega1*omega2)*omegan2
    if (abs(omega1-omega2)<1e-10)then
        amp=1/(omega1*omega2)*omegan2
    else
        amp=0
    endif
    comp(1,1)=u1*u2*gl(ll,0,0)
    comp(1,2)=u1*v2*gl(ll,-1,0)
    comp(1,3)=-u1*w2*gl(ll,-1,0)
    comp(2,1)=v1*u2*gl(ll,-1,0)
    comp(2,2)=v1*v2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5-w1*w2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5
    comp(2,3)=-v1*w2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5-w1*v2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5
    comp(3,1)=-w1*u2*gl(ll,-1,0)
    comp(3,2)=-w1*v2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5-v1*w2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5
    comp(3,3)=-v1*v2*(gl(ll,1,-1)+gl(ll,-1,-1))*0.5+w1*w2*(gl(ll,1,-1)-gl(ll,-1,-1))*0.5
    if(itype==3)then
        const=((sigma1+sigma2)*(sigma1+sigma2)+(omega1-omega2)*(omega1-omega2))*((sigma1+sigma2)*(sigma1+sigma2)+(omega1+omega2)*(omega1+omega2))
        c1=((sigma1+sigma2)*2*omega1)/const*omega1*omega2*omega2
        c2=((sigma1+sigma2)*(sigma1+sigma2)+omega2*omega2-omega1*omega1)/const*omega1*omega2*omega2
    endif
    if(itype==2)then
        consttemp=0.5*((1-dexp(-(sigma1+sigma2)*tw)*dcmplx(dcos((omega1-omega2)*tw),-dsin((omega1-omega2)*tw)))/dcmplx(sigma1+sigma2,omega1-omega2)-(1-dexp(-(sigma1+sigma2)*tw)*dcmplx(dcos((omega1+omega2)*tw),-dsin((omega1+omega2)*tw)))/dcmplx(sigma1+sigma2,omega1+omega2))
        c1=dreal(consttemp)*omega1*omega2
        c2=-dimag(consttemp)*omega1*omega2
    endif
    if(itype==1)then
        if(dabs(omega1-omega2)<1e-10)then
            c1=0.5*tw*omega1*omega2
            c2=0
        else
            c1=0.5*(dsin((omega1-omega2)*tw))/(omega1-omega2)*omega1*omega2
            c2=0.5*(dcos((omega1-omega2)*tw)-1)/(omega1-omega2)*omega1*omega2
        endif
    endif
    tp=dcmplx(c1,-c2)
    ampsum=ampsum+tp*amp*comp(3,3)
    deallocate(comp)
    return
    end
subroutine add(ll,ampsum,omega1,q1,domega,nps,uf)
    integer ll,i,j,k,nspec,iomegac1,lmax
    real(8) domega,pi,omega1,q1,sigma1,sigma0,taper,lrate
    real(8) omegan2,rhon,g
    complex(8) ampsum
    complex(8) uf(-nps:nps-1)
    lmax=900
    lrate=(1.0*ll)/(1.0*lmax)
    taper=1
    if(lrate>1)then
        taper=0
    elseif(lrate>0.9)then
        taper=erf(-(lrate-1.0)/(lrate-0.899))
    endif
    rhon=5515
    g=6.6723e-11
    pi=3.1415926535
    nspec=500
    iomegac1=int(omega1/domega)
    omegan2=rhon*pi*g
    sigma1=0.5*omega1/q1
    sigma0=domega
    if (iomegac1+nspec<nps-1)then
        do i=-nspec,nspec
            omega=domega*(i+iomegac1)
            uf(i+iomegac1)=uf(i+iomegac1)+(2.0*ll+1.0)*ampsum*taper/dcmplx((sigma1+sigma0),omega-omega1)/(4*pi*omegan2*omegan2)
        enddo
    endif
    end
subroutine GL_generator_MR(l,gl,sita,N1,N2)
!!according to Masters&Richards-Dinger(1998)
    integer l,N,m,i,N1,N2
    real(8) gl(-l:l,N1:N2)
    real(8), allocatable :: gld(:)
    real(8) s,c,sita,sumd,sumu,md,mu,dd,du,temp,temp1,temp2,ad,au,norm,mulu,muld,mul
    
    s=dsin(sita)
    c=dcos(sita)
    
    allocate(gld(-l:l))
    do N=N1,N2,1
        sumd=1
        sumu=1
        mc=int(N*c)
        !-l<lm<l
        !downward
        if(mod(l-N,2)==0)then
             gl(l,N)=1
             gld(l)=(l*c-N)/s
        else
            gl(l,N)=-1
            gld(l)=-(l*c-N)/s
        endif
        
        m=l
        norm=gl(l,N)
        do while (m>mc)
            gl(m-1,N)=-(gld(m)+(m*c-N)*norm/s)/sqrt((l+m)*(l-m+1.0))
            gld(m-1)=((m-1)*c-N)*gl(m-1,N)/s+sqrt((l+m)*(l-m+1.0))*norm
            sumd=sumd+gl(m-1,N)*gl(m-1,N)
            if (1<dabs(gl(m-1,N)))then
                temp=dabs(gl(m-1,N))
                !do i=m-1,l,1
                !    gl(i,N)=gl(i,N)/temp
                !enddo
                gld(m-1)=gld(m-1)/temp
                sumd=sumd/(temp*temp)
                if(gl(m-1,N)>0)then
                    norm=1
                else
                    norm=-1
                endif
            else
                norm=gl(m-1,N)
                temp=1
            endif
            
            m=m-1
        enddo
        muld=temp
        md=gl(mc,N)/temp
        dd=gld(mc)
        !upward
        gl(-l,N)=1
        gld(-l)=(l*c+N)/s
        temp=1
        m=-l
        norm=gl(-l,N)
        do while (m<mc)
            gl(m+1,N)=(gld(m)-(m*c-N)*norm/s)/sqrt((l-m)*(l+m+1.0))
            gld(m+1)=-((m+1)*c-N)*gl(m+1,N)/s-sqrt((l-m)*(l+m+1.0))*norm
            sumu=sumu+gl(m+1,N)*gl(m+1,N)
            if (1<dabs(gl(m+1,N)))then
                temp=dabs(gl(m+1,N))
                !do i=-l,m+1,1
                !    gl(i,N)=gl(i,N)/temp
                !enddo
                gld(m+1)=gld(m+1)/temp
                sumu=sumu/(temp*temp)
                if(gl(m+1,N)>0)then
                    norm=1
                else
                    norm=-1
                endif
            else
                norm=gl(m+1,N)
                temp=1
            endif
            m=m+1
        enddo
        mulu=temp
        mu=gl(mc,N)/temp
        du=gld(mc)
        
        temp1=sqrt(md*md*sumu+mu*mu*sumd-mu*mu*md*md)
        temp2=sqrt(dd*dd*sumu+du*du*sumd-mu*mu*dd*dd)
        if(temp1>=temp2)then
            ad=dabs(mu)/temp1
            au=dabs(md)/temp1
            gl(mc,N)=mu*dabs(md)/temp1
        else
            ad=dabs(du)/temp2
            au=dabs(dd)/temp2
            gl(mc,N)=dabs(du)*md/temp2
        endif
        
        mul=au/mulu
        do m=mc-1,-l,-1
            if(dabs(gl(m,N))>1)then
                mulu=dabs(gl(m,N))
            else
                mulu=1
            endif
            mul=mul/mulu
            gl(m,N)=gl(m,N)*mul
        enddo
        mul=ad/muld
        do m=mc+1,l,1
            if(dabs(gl(m,N))>1)then
                muld=dabs(gl(m,N))
            else
                muld=1
            endif
            mul=mul/muld
            gl(m,N)=gl(m,N)*mul
        enddo
        
    enddo
    deallocate(gld)
    return
    end
    
