        module array
	  complex*16,allocatable:: H_c(:,:),work(:),H(:,:),work1(:)
	  integer,allocatable::pi(:),pj(:),pci(:),pcj(:)
	  double precision, allocatable::m_s(:),th_s(:),ph_s(:),rwork(:),evl_c(:),dval(:)
	  double precision, allocatable::m_c(:),th_c(:),ph_c(:),p_m(:),rwork1(:),evl_s(:)
          double precision, allocatable::n_total(:),n_total_c(:),ion_dis(:),ion_dis_c(:)
	end module array

	module global
	integer ::temp,count1,d,dc,temp_max,MCSW,intrvl,seed
	double precision ::t1,t2,U1,filling,gama,gama_m
	double precision ::T,dp,sfc1_pp,sfc1_op,sfc1_po,ds,n_d_pls,n_d_mns
	double precision ::n_d_pls_c,n_d_mns_c
	end module global

	program tca_mf_mc_code
	use mtmod
	use array
	use global
	implicit none
	integer ::ref,accept,counter,i,j,ik,mc_count,l,iw,c,xc,yc,x,y,itr
	double precision::mu_c_p,mu_c_cur,eng_c_p,eng_c_cur,temp_m_s,temp_th_s,temp_ph_s
	double precision::p,eng_final,get_mu,r,mu_c_final,mc_av_eng,mc_av_mu,f,m,fs
	double precision::mu_sys,eng_sys,get_mu_s,wc,sit_av_mu,sit_av_eng,mc_eng_sys,strnth
	double precision::dp1,ds1,inp
	open (5,file='input.dat',status='unknown')
	do i=1,13
	read (5,*)inp
	if (i.eq.1)d=int(inp)
	if (i.eq.2)dc=int(inp)
	if (i.eq.3)t1=dble(inp)
	if (i.eq.4)t2=dble(inp)
	if (i.eq.5)temp_max=int(inp)
	if (i.eq.6)mcsw=int(inp)
	if (i.eq.7)intrvl=int(inp)
	if (i.eq.8)U1=dble(inp)
	if (i.eq.9)filling=dble(inp)
	if (i.eq.10)gama=dble(inp)
	if (i.eq.11)gama_m=dble(inp)
	if (i.eq.12)seed=int(inp)
if (i.eq.13)strnth=dble(inp)
	end do


	print*,"system size=",d
	print*,"cluster size=",dc
	print*,"nn hopping parameter=",t1
	print*,"nnn hopping parameter=",t2
	print*,"maximum no of temperature points=",temp_max
	print*,"total no of system sweeps MCSW=",MCSW
	print*,"interval in system sweeps steps to cal. observables=",intrvl
	print*,"the interaction value U1=",U1
	print*,"filling in the system=",filling
	print*,"the broadening of the lorenzian for cal. DOS=",gama
	print*,"the broadening of the lorentzian for cal. of distributin p(m) of m=",gama_m
	print*,"strnth of the ionic disorder=",strnth       



	allocate(n_total(d**2),n_total_c(dc**2),ion_dis(d**2),ion_dis_c(dc**2))
	allocate(m_s(d**2),th_s(d**2),ph_s(d**2))
	allocate(m_c(dc**2),th_c(dc**2),ph_c(dc**2))
	allocate(pi(d**2),Pj(d**2),pci(d**2),Pcj(d**2))
	allocate(H_c(2*dc**2,2*dc**2),p_m(100),dval(1000))
	allocate(evl_c(2*dc**2),work(2*(2*dc**2)-1),rwork(3*(2*dc**2)-2))
allocate(evl_s(2*d**2),work1(2*(2*d**2)-1),rwork1(3*(2*d**2)-2),H(2*d**2,2*d**2))
	dp=filling*dc**2      !! dp > average no of particles in the cluster
	ds=filling*d**2       !! ds > average no of particles in the cluster 
H_c=cmplx(0.0d0,0.0d0)
	call sgrnd(seed)     !!sgrand(seed)  > random number generator (Mersenne Twister)   ! mtfort.f
	call mysamp(0,1)     !! mysamp >  generating initial aux. field configuration 

	!do l=1,d**2
	!write(11,*)l,m_s(l),th_s(l),ph_s(l)     !printing the initial aux. field configaration
	!enddo
	!!m_s(l),th_s(l),ph_s(l) are system aux. fields and s for system
	!!***************************************************************************************
	!   Levels for the system
	!!***************************************************************************************
	do l=1,d**2                     !l -> levels of the lattice points
	do y=1,d
	do x=1,d
	if (d*(y-1)+x.eq.l) then 
	pi(l)=x                   !x co-ordinate of the system
	pj(l)=y                   !y co-ordinate of the system
	endif
	enddo
	enddo
	enddo
	!________________________________________________________________________________________

	!*****************************************************************************************
	!  Levels for the cluster
	!*****************************************************************************************
	do c=1,dc**2
	do yc=1,dc
	do xc=1,dc
	if (dc*(yc-1)+xc.eq.c) then
	pci(c)=xc                            ! c for cluster
	pcj(c)=yc
	endif
	enddo
	enddo
	enddo
	!_________________________________________________________________________________________  
	l=1
	do y=1,d
	do x=1,d
	ion_dis(l)=strnth*(-1)**(x+y)  !ion_dis(l) -> site dependent ionic disorder with strength "strnth"
!       write(31,*)l,ion_dis(l)
	l=l+1
	enddo
	enddo
	!--------------------------------------------------------------------------------------

	!      do l=1,d**2
	!       n_total(l)=filling       !fixing of initial density per site of the system
	!     enddo  
	!----------------------------------------------------------------------------------------   

	T=0.1750d0

	do temp=1,temp_max              !! Start of Temperature loop
	if(temp.le.3)T=T-0.025
	if((temp.gt.3).and.(temp.le.12))T=T-0.01
	if(temp.eq.13)T=0.005
	if(temp.eq.14)T=0.001
	!if(temp.le.13)T=T-0.025d0
	!if((temp.gt.13).and.(temp.le.22))T=T-0.01d0
	!if(temp.gt.22)T=0.005d0

	!print*,temp,T

	accept=0
	mc_av_eng=0.0d0
	mc_av_mu=0.0d0
	p_m=0.0d0
	sfc1_pp=0.0d0
	sfc1_op=0.0d0
	sfc1_po=0.0d0
	mc_count=0
	dval=0.0d0
	mc_eng_sys=0.0d0

	do count1=1,MCSW                     !! Start of Monte Carlo system sweep loop
	counter=0
	sit_av_eng=0.0d0
	sit_av_mu=0.0d0
	do ref=1,d**2                        !! Start of site loop 
	call mapping_levels(ref)             !! mapping_levels(ref) > maps system levels to cluster levels with refference point ref

	call clus_mat_gen                     !! clus_mat_gen > generating cluster matrix
	mu_c_p=get_mu(dble(dp))                !! mu_c_p > previous(before update) cluster chemical potential
	call energy(mu_c_p,eng_c_p)            !! eng_c_p > previous(before update) cluster energy

	temp_m_s=m_s(ref)
	temp_th_s=th_s(ref)
temp_ph_s=ph_s(ref)

	call mysamp(1,ref)                   !! mysamp(1,ref)  > called to update aux. f at site ref in the array of the system

	call mapping_levels(ref)             !!  maps again after the update at referece point ref from system to cluster levels
	call clus_mat_gen                     !! construct the cluster matrix again

	mu_c_cur=get_mu(dble(dp))        !! mu_c_cur > current cluster chemical potential after update

	call energy(mu_c_cur,eng_c_cur)      !!eng_c_cur > current cluster energy after update
	counter=counter+1
	!_______________________________________________________________________
	!Metropolis algorithm
	!_______________________________________________________________________
	if((eng_c_cur-eng_c_p).lt.0.0d0)then
	eng_final=eng_c_cur
	mu_c_final=mu_c_cur
	accept=accept+1
	!-----------------------------------------------------------------------
	! self-consistantly density is calcualated
	do itr=1,1                        !iteration loop for fixing density
call clus_avg_n(mu_c_final,dp)
	if(count1.eq.MCSW)then
write(93,*)T,n_d_pls_c/dble(dc**2),n_d_mns_c/dble(dc**2)
	endif
	call clus_mat_gen 
mu_c_final=get_mu(dble(dp)) 
	enddo        
	!-----------------------------------------------------------------------
	else
	p=exp(-((eng_c_cur-eng_c_p)/T))
call rannum(r)
	if (r.lt.p)then
	eng_final=eng_c_cur
	mu_c_final=mu_c_cur
	accept=accept+1

	!-----------------------------------------------------------------------
	! self-consistantly density is calcualated
	do itr=1,1                    !iteration loop for fixing density
call clus_avg_n(mu_c_final,dp)
	if(count1.eq.MCSW)then
write(93,*)T,n_d_pls_c/dble(dc**2),n_d_mns_c/dble(dc**2)
	endif
	call clus_mat_gen 
mu_c_final=get_mu(dble(dp)) 
	enddo   
	!-----------------------------------------------------------------------       
	else
	eng_final=eng_c_p                 !!eng_final > final cluster energy when the update is accepted
	mu_c_final=mu_c_p                  !!mu_c_final > final cluster chemical potential when the update is accepted
	m_s(ref)=temp_m_s
	th_s(ref)=temp_th_s
	ph_s(ref)=temp_ph_s
	endif

	endif

	!***************************************************************** !end of Metropolis algorithm
	f=0.0d0                             !!f > average no of particle in the cluster
	do i=1,2*dc**2
f=f+(1.0d0/(exp((evl_c(i)-mu_c_final)/T)+1.0d0))
	end do
	write(18,*)T,ref,f,mu_c_final


call mapping_levels_clus_sys(ref)
	!write(12,*)T,mu_c_final
	sit_av_eng=sit_av_eng+eng_final                         !! sit_av_eng > site average of cluster energy
	sit_av_mu=sit_av_mu+mu_c_final/(dble(d**2))                    !!  sit_av_mu > site average of cluster chemical potential

	enddo           !site levels



	print*,T,count1

	if(count1.eq.MCSW)then
	do l=1,d**2
write(51,*)T,l,n_total(l)
	enddo
	endif

	mc_av_mu=mc_av_mu+sit_av_mu                               !!mc_av_mu > sum of cluster chemical potn. over system sweeps
	!________________________________________________________________________
	!!calculation of outputs
	!________________________________________________________________________

	if((count1.ge.MCSW/2).and.(mod(count1,intrvl).eq.0))then !calculation of output
	mc_count=mc_count+1                                 !!mc_count > no of system sweep over which ovservables are calculated
	!print*,mc_count,sit_av_eng
	mc_av_eng=mc_av_eng+(sit_av_eng/(dble(d**2)))              !!mc_av_eng > average cluster energy
354 format (1x,i4,4f16.8)
	do l=1,d**2
	write(300+temp,354)l,m_s(l),th_s(l),ph_s(l),n_total(l)
flush(300+temp)
	enddo
	call cal_n
write(39,*)T,n_d_pls/dble(d**2),n_d_mns/dble(d**2)

	call distri_local_mnt              !! distribution of aux.f m
	call struc_factr                    !! structure factors

	call system_mat_gen


	mu_sys=get_mu_s(dble(ds))            !!mu_sys > system chemical potential
	if(count1.eq.MCSW)then
	do l=1,2*d**2
write(1600+temp,*)evl_s(l)
	enddo
	endif



	write(33,*)T,sit_av_mu,mu_sys
flush(33)

	call energy_s(sit_av_mu,eng_sys)         !!eng_sys > system energy
	mc_eng_sys=mc_eng_sys+eng_sys
	!write(19,*)T,eng_sys,mc_count
!flush(19)
	call dos
	endif

	enddo   !______________________________________________system sweep loop
	!*****************************************************************************************
	!do i=1,2*dc**2
	!write(23,*)evl_c(i)                    !! evl_c(i) > eigen values of the cluster Hamiltonian
	!enddo

	fs=0.0d0                                 !!fs > average no of particle in the system
	do i=1,2*d**2
fs=fs+(1.0d0/(exp((evl_s(i)-sit_av_mu)/T)+1.0d0))
	!! evl_s(i) > eigen values of the system Hamiltonian
	end do
	write(13,*)T,fs,sit_av_mu



	!      write(14,*)T,(mc_av_eng/dble(mc_count)),mc_count,accept

	!      write(15,*)T,mc_av_mu/dble(MCSW),f
	!      write(16,*)T,(mc_eng_sys/dble(mc_count))     !mc_eng_sys  > system average energy


	m=0.0d0
	do ik=1,100
	m=m+0.010d0
	write(600+temp,*)m,p_m(ik)/dble(mc_count)      !! p_m(ik) > distribution function of m aux.f
	enddo

	wc=-30.0d0
	do iw=1,1000
	wc=wc+0.06d0

write(800+temp,*)(wc-sit_av_mu),(dval(iw)/dble(mc_count))
	enddo

	write(30,*)T,(sfc1_pp/dble(mc_count)),(sfc1_op/dble(mc_count)),(sfc1_po/dble(mc_count))
flush(30)

	!!sfc1_pp > (pi,pi) structure factor ; sfc1_op > (0,pi) structure factor ;sfc1_po > (pi,o) structure factor 

	enddo !_________________________________________________temperature loop





	end !end of the main program

	!*******************************************************************
	!subroutine for calculation of random spin components
	!*******************************************************************
subroutine mysamp(flag,position)
	use array
	!	use input
	use global
	implicit none
	integer::i,position,flag,l
	double precision r,th1,ph1
	if(flag.eq.0)then
	do l=1,d**2
	call ang(th1,ph1)
call rannum(r)
	m_s(l)=1.50d0*r
	th_s(l)=th1
	ph_s(l)=ph1
	n_total(l)=filling
	enddo
	endif
	if(flag.eq.1)then
	call rannum(r)
call ang(th1,ph1)
	m_s(position)=1.50d0*r
	th_s(position)=th1
	ph_s(position)=ph1
	endif
	return
	end
	!*****************************************************************
	!random no generator
	!********************************************************************
subroutine rannum(r1)
	use mtmod
	implicit none
	double precision r1
r1=grnd()
	return
	end
	!********************************************************************
	subroutine ang(th,ph)      !for choosing angles
	use mtmod
	IMPLICIT NONE
	double precision x1,x2,norm,pi,ran2,ph,th,rm
	pi=acos(-1.0d0)
call rannum(rm)
	20     x1=1.0d0-2.0d0*rm
call rannum(rm)
	x2=1.0d0-2.0d0*rm
	norm = x1**2 + x2**2
	if(norm.ge.1.0d0) goto 20
	th=acos(1.0d0-2.0d0*norm)
ph=atan(x2/x1)
	if(x1.lt.0) ph=Pi-ph
	return
	end
	!***********************************************************************
	!subroutine for mapping of system to cluster levels 
	!*************************************************************************

subroutine mapping_levels(rf)
	!use input
	use array
	use global
	implicit none
	integer::ls,rf,q,c,p,i,j,pbcx,x,y,k,s,pbcy,yc,xc,l




	ls=rf
	q=0
	i=0
	x=pi(rf)               !x=pi(rf) > x levels in the system
	y=pj(rf)                !y=pj(rf) > y levels in the system
!write(44,*)rf,pi(rf),pj(rf)
	do c=1,dc**2
	xc=pci(c)
	yc=pcj(c)
!write(45,*)c,pci(c),pcj(c)
	q=q+1
	if(((d-(x-1)).lt.dc).and.((x-1+xc).gt.d))then
	pbcx=-d                                        !pbcx > periodic boundary condition along x
	else
	pbcx=0
	endif
	if(((d-(y-1)).lt.dc).and.((y-1+yc).gt.d))then
	pbcy=-d**2                                      !pbcy > periodic boundary condition along y
	else
	pbcy=0
	endif
	p=ls+(q-1)+pbcx+pbcy                       ! p > variable for mapping levels from system to cluster including PBC
	m_c(c)=m_s(p)                             !!m_c(c), th_c(c), ph_c(c) are aux.f in the cluster
	th_c(c)=th_s(p)
	ph_c(c)=ph_s(p)
	n_total_c(c)=n_total(p)
ion_dis_c(c)=ion_dis(p)
	if (mod(c,dc).eq.0)then
	q=0
	i=i+1
	ls=rf+d*i
	endif
	!print*,rf,c,pbcx,pbcy,p,dc
	enddo

	return
	end




	!***********************************************************************
	!generation of matrix for the cluster and diagonalization
	!***********************************************************************
	subroutine clus_mat_gen
	!	use input
	use array
	use global
	implicit none
	integer :: c,ii,id,ji,jd,k,a,b,i,j,info
	double precision r
	!call avg_n
H_c=cmplx(0.0d0,0.0d0)

	do c=1,dc**2             ! c > site index in the cluster 
	ii=1
	id=-1
	ji=1
	jd=-1
	i=pci(c)
	j=pcj(c)
!        write(46,*)c,pci(c),pcj(c)
	if (i.eq.1) id=-1+dc
	if (i.eq.dc) ii=1-dc
	if (j.eq.1) jd=-1+dc
	if (j.eq.dc) ji=1-dc
	do k=1,dc**2
	if(c.eq.k)then
	a=2*k-1
	b=2*k-1
	H_c(a,b)=(-abs(m_c(k))*cos(th_c(k))*(U1/2.0d0))+((n_total_c(k))*U1/2.0d0)+ ion_dis_c(k)
	H_c(a+1,b+1)=(abs(m_c(k))*cos(th_c(k))*(U1/2.0d0))+((n_total_c(k))*U1/2.0d0)+ ion_dis_c(k)
	H_c(a,b+1)=-abs(m_c(k))*sin(th_c(k))*cmplx(cos(ph_c(k)),-sin(ph_c(k)))*(U1/2.0d0)
H_c(a+1,b)=conjg(H_c(a,b+1))
	if((temp.eq.temp_max).and.(count1.eq.MCSW))then
write(34,*)real(H_c(a,b)),real(H_c(a+1,b+1))
	endif

	endif
	if (((pci(k).eq.(i+ii)) .and. (pcj(k).eq.j))&
			&.or. ((pci(k).eq.i) .and. (pcj(k).eq.(j+ji)))) then
	a=2*c-1
	b=2*k-1
	H_c(a,b)=t1
	H_c(a+1,b+1)=t1
	endif
	if (((pci(k).eq.(i+id)) .and. (pcj(k).eq.j))&
			&.or. ((pci(k).eq.i) .and. (pcj(k).eq.(j+jd)))) then
	a=2*c-1
	b=2*k-1
	H_c(a,b)=t1
	H_c(a+1,b+1)=t1
	endif
	if (((pci(k).eq.(i+ii)).and.(pcj(k).eq.(j+ji)))&
			&.or. ((pci(k).eq.(i+id)).and.(pcj(k).eq.(j+jd)))) then
	a=2*c-1
	b=2*k-1
	H_c(a,b)=t2
	H_c(a+1,b+1)=t2
	endif
	if (((pci(k).eq.(i+id)).and.(pcj(k).eq.(j+ji)))&
			&.or. ((pci(k).eq.(i+ii)).and.(pcj(k).eq.(j+jd)))) then
	a=2*c-1
	b=2*k-1
	H_c(a,b)=t2
	H_c(a+1,b+1)=t2
	endif

	enddo
	enddo
	!  do i=1,2*d**2
	!  do j=1,2*d**2
!   print*,i,j,H_c(i,j)
	!enddo
	!enddo
	!stop
	call zheev ('V','U',2*dc**2,H_c,2*dc**2,evl_c,work,2*(2*dc**2)-1,rwork,info)
	!do i=1,2*dc**2
	!print*,i,evl_c(i),info
	!enddo
	return
	end
	!********************************************************************************
!function for calculating of mu(chemical potential)
	!********************************************************************************
double precision function get_mu(fill)
	use array
	!       use input
	use global
	implicit none
	double precision f0, f, fL2, fR, mR, mL, rtmp,m_d
	integer i
	double precision fill
	mR = maxval(evl_c)       !right-side chemical potential
	fr=0.0d0
	do i=1,2*dc**2
fr=fr+(1.0d0/(exp((evl_c(i)-mR)/T)+1.0d0))
	end do
	mL = minval(evl_c)       !left-side chemical potential
	fL2=0.0d0
	do i=1,2*dc**2
fL2=fL2+(1.0d0/(exp((evl_c(i)-mL)/T)+1.0d0))
	end do
	m_d = 0.5d0*(mL+mR)    !middle chemical potential
	f=0.0d0
	do i=1,2*dc**2
f=f+(1.0d0/(exp((evl_c(i)-m_d)/T)+1.0d0))
	end do
	!print*,f,fill
	do while(abs(f-fill).ge.1e-8)
m_d = 0.5d0*(mL+mR)
	f=0.0d0
	do i=1,2*dc**2
f=f+(1.0d0/(exp((evl_c(i)-m_d)/T)+1.0d0))
	end do
	if(f.gt.fill)then
	!if middle filling is above target, make it the new right bound.
	mR = m_d
	fR = f
	elseif(f.lt.fill)then
	!if middle filling is below target, make it the new left bound.
	mL = m_d
	fR = f
	endif
	enddo
	!Return the middle value
	get_mu = m_d
	return
	end function get_mu
	!********************************************************************
	!subroutine for calculation of energy
	!********************************************************************
subroutine energy(EF,Ein)
	use array
	!	use input
	use global
	implicit none
	integer:: i
	double precision ::sum,sum1,Ein,EF   	


	sum=0.0d0
	!if(temp.le.9) then
	do i=1,2*dc**2
	sum=sum+evl_c(i)*0.5d0*(1.0d0+tanh((EF-evl_c(i))/(2.0d0*T)))
!/((exp(evl(i)-EF)/T)+1.0d0)
	enddo

	sum1=0.0d0
	do i=1,dc**2
sum1=sum1+(m_c(i)**2)*(U1/4.0d0)-((n_total_c(i))**2)*(U1/4.0d0)
	enddo
	Ein=sum+sum1
	return
	end
	!***********************************************************************
	!Calculation of distribution of local moment 
	!***********************************************************************
	subroutine distri_local_mnt
	!    use input
	use global
	use array
	implicit none
	integer::i,j
	double precision::p_m_sum,pp,m
	!    gama_m=0.1
pp=acos(-1.0d0)
	m=0.0d0
	do j=1,100
	m=m+0.010d0
	p_m_sum=0.0d0
	do i=1,d**2
	p_m_sum=p_m_sum+((gama_m/pp)/((m-m_s(i))**2+(gama_m)**2))   !!gama_m > broading of the lorentzian
	enddo
	p_m(j)=p_m(j)+p_m_sum

	enddo

	return
	end
	!!**********************************************************************
	! subroutine for magnetic structure factor
	!***********************************************************************
	subroutine struc_factr
	use array
	!	use input
	use global
	implicit none
	integer::i,j
	double precision ::sf_sum_pp,sf_sum_op,sf_sum_po
	double precision ::arg_pp,arg_op,arg_po,pp
	sf_sum_pp=0.0d0
	sf_sum_op=0.0d0
	sf_sum_po=0.0d0
pp=acos(-1.0d0)
	do i=1,d**2
	do j=1,d**2

arg_pp=(pp*(pi(i)-pi(j))+pp*(pj(i)-pj(j)))
	sf_sum_pp=sf_sum_pp+(m_s(i)*m_s(j))*(sin(th_s(i))*cos(ph_s(i))*sin(th_s(j))*cos(ph_s(j))&
			+sin(th_s(i))*sin(ph_s(i))*sin(th_s(j))*sin(ph_s(j))&
			+cos(th_s(i))*cos(th_s(j)))*cmplx(cos(arg_pp),sin(arg_pp))
arg_op=(0.0d0*(pi(i)-pi(j))+pp*(pj(i)-pj(j)))
	sf_sum_op=sf_sum_op+(m_s(i)*m_s(j))*(sin(th_s(i))*cos(ph_s(i))*sin(th_s(j))*cos(ph_s(j))&
			+sin(th_s(i))*sin(ph_s(i))*sin(th_s(j))*sin(ph_s(j))&
			+cos(th_s(i))*cos(th_s(j)))*cmplx(cos(arg_op),sin(arg_op))
arg_po=(pp*(pi(i)-pi(j))+0.0d0*(pj(i)-pj(j)))
	sf_sum_po=sf_sum_po+(m_s(i)*m_s(j))*(sin(th_s(i))*cos(ph_s(i))*sin(th_s(j))*cos(ph_s(j))&
			+sin(th_s(i))*sin(ph_s(i))*sin(th_s(j))*sin(ph_s(j))&
			+cos(th_s(i))*cos(th_s(j)))*cmplx(cos(arg_po),sin(arg_po))
	enddo
	enddo
	sfc1_pp=sfc1_pp+(sf_sum_pp/dble((d)**4))
	sfc1_op=sfc1_op+(sf_sum_op/dble((d)**4))
sfc1_po=sfc1_po+(sf_sum_po/dble((d)**4))
	!	write(89,*)T,sfc1_pp
	return
	end
	!********************************************************************************************
	!subroutine for matrix formation and diagonalization for the system
	!********************************************************************************************
	subroutine system_mat_gen
	!	use input
	use array
	use global
	implicit none
	integer :: l,ii,id,ji,jd,k,a,b,i,j,info
	double precision r
H=cmplx(0.0d0,0.0d0)

	do l=1,d**2
	ii=1
	id=-1
	ji=1
	jd=-1
	i=pi(l)
	j=pj(l)
!        write(47,*)l,pi(l),pj(l)
	if (i.eq.1) id=-1+d
	if (i.eq.d) ii=1-d
	if (j.eq.1) jd=-1+d
	if (j.eq.d) ji=1-d
	do k=1,d**2
	if(l.eq.k)then
	a=2*k-1
	b=2*k-1
	H(a,b)=(-abs(m_s(k))*cos(th_s(k))*(U1/2.0d0))+((n_total(k))*U1/2.0d0)+ ion_dis(k)
	H(a+1,b+1)=(abs(m_s(k))*cos(th_s(k))*(U1/2.0d0))+((n_total(k))*U1/2.0d0)+ ion_dis(k)
	H(a,b+1)=-abs(m_s(k))*sin(th_s(k))*cmplx(cos(ph_s(k)),-sin(ph_s(k)))*(U1/2.0d0)
H(a+1,b)=conjg(H(a,b+1))
	endif
	if (((pi(k).eq.(i+ii)) .and. (pj(k).eq.j))&
			&.or. ((pi(k).eq.i) .and. (pj(k).eq.(j+ji)))) then
	a=2*l-1
	b=2*k-1
	H(a,b)=t1
	H(a+1,b+1)=t1
	endif
	if (((pi(k).eq.(i+id)) .and. (pj(k).eq.j))&
			&.or. ((pi(k).eq.i) .and. (pj(k).eq.(j+jd)))) then
	a=2*l-1
	b=2*k-1
	H(a,b)=t1
	H(a+1,b+1)=t1
	endif
	if (((pi(k).eq.(i+ii)).and.(pj(k).eq.(j+ji)))&
			&.or. ((pi(k).eq.(i+id)).and.(pj(k).eq.(j+jd)))) then
	a=2*l-1
	b=2*k-1
	H(a,b)=t2
	H(a+1,b+1)=t2
	endif
	if (((pi(k).eq.(i+id)).and.(pj(k).eq.(j+ji)))&
			&.or. ((pi(k).eq.(i+ii)).and.(pj(k).eq.(j+jd)))) then
	a=2*l-1
	b=2*k-1
	H(a,b)=t2
	H(a+1,b+1)=t2
	endif
	enddo
	enddo
	call zheev ('N','U',2*d**2,H,2*d**2,evl_s,work1,2*(2*d**2)-1,rwork1,info)

	return
	end
	!****************************************************************************
	!calculation of chemical potential for the system
	!***************************************************************************
double precision function get_mu_s(fill)
	use array
	!    use input
	use global
	implicit none
	double precision f0, f, fL2, fR, mR, mL, rtmp,m_d
	integer i
	double precision fill
	mR = maxval(evl_s)       !right-side chemical potential
	fr=0.0d0
	do i=1,2*d**2
fr=fr+(1.0d0/(exp((evl_s(i)-mR)/T)+1.0d0))
	end do
	mL = minval(evl_s)       !left-side chemical potential
	fL2=0.0d0
	do i=1,2*d**2
fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/T)+1.0d0))
	end do
	m_d = 0.5d0*(mL+mR)    !middle chemical potential
	f=0.0d0
	do i=1,2*d**2
f=f+(1.0d0/(exp((evl_s(i)-m_d)/T)+1.0d0))
	end do
	!print*,f,fill
	do while(abs(f-fill).ge.1e-8)
m_d = 0.5d0*(mL+mR)
	f=0.0d0
	do i=1,2*d**2
f=f+(1.0d0/(exp((evl_s(i)-m_d)/T)+1.0d0))
	end do
	if(f.gt.fill)then
	!if middle filling is above target, make it the new right bound.
	mR = m_d
	fR = f
	elseif(f.lt.fill)then
	!if middle filling is below target, make it the new left bound.
	mL = m_d
	fR = f
	endif
	enddo
	!Return the middle value
	get_mu_s = m_d
	return
	end function get_mu_s
	!****************************************************************************
	!calculation of system energy
	!****************************************************************************

subroutine energy_s(EF,Ein)
	use array
	!	use input
	use global
	implicit none
	integer:: i
	double precision ::sum,sum1,Ein,EF   	


	sum=0.0d0
	!if(temp.le.9) then
	do i=1,2*d**2
	sum=sum+evl_s(i)*0.5d0*(1.0d0+tanh((EF-evl_s(i))/(2.0d0*T)))
!/((exp(evl(i)-EF)/T)+1.0d0)
	enddo

	sum1=0.0d0
	do i=1,d**2
sum1=sum1+(m_s(i)**2)*(U1/4.0d0)-((n_total(i))**2)*(U1/4.0d0)
	enddo
	Ein=sum+sum1
	return
	end
	!*************************************************************************
	!subroutine for density of states
	!*************************************************************************
	subroutine dos
	use array
	!	use input
	use global
	implicit none
	integer::jc,i,iw
	double precision::wc,d_sum,pp

pp=acos(-1.0d0)
	wc=-30.0d0
	do iw=1,1000
	wc=wc+0.06d0

	d_sum=0.0d0
	do i=1,2*d**2
	d_sum=d_sum+(gama/pp)/((wc-evl_s(i))**2+(gama)**2)  !gama > broading of the lorentzian
	enddo
dval(iw)=dval(iw)+(d_sum/dble(2*d**2))



	enddo

	return
	end

	!**************************************************************************
	!subroutine for calculation of total no electrons of the system
	!**************************************************************************
subroutine sys_avg_n(EF,d1_s)
	use array
	!	use input
	use global
	implicit none
double precision, allocatable::n_up(:),n_down(:)
	integer ::i,j
	double precision ::sum_nu,sum_nd,s_nup,s_ndown,EF,d1_s
allocate(n_up(2*d**2),n_down(2*d**2))
	n_up=0.0d0
	n_down=0.0d0
	s_nup=0.0d0
	do i=1,2*d**2,2
	sum_nu=0.0d0
	do j=1,2*d**2
	sum_nu=sum_nu+(H(i,j)*conjg(H(i,j))*(1.0d0/(exp((evl_s(j)-EF)/T)+1.0d0)))
s_nup=s_nup+(H(i,j)*conjg(H(i,j))*(1.0d0/(exp((evl_s(j)-EF)/T)+1.0d0)))
	enddo                               !sum_nu=total no of up electrons
	n_up(i)=sum_nu
	enddo

	s_ndown=0.0d0
	do i=2,2*d**2,2
	sum_nd=0.0d0
	do j=1,2*d**2
	sum_nd=sum_nd+(H(i,j)*conjg(H(i,j))*(1.0d0/(exp((evl_s(j)-EF)/T)+1.0d0)))
s_ndown=s_ndown+(H(i,j)*conjg(H(i,j))*(1.0d0/(exp((evl_s(j)-EF)/T)+1.0d0)))
	enddo
	n_down(i)=sum_nd           !sum_nd=total no of down electrons
	enddo
	!	d1_s= s_ndown+s_nup
	write(38,*)d1_s,s_ndown,s_nup
	do i=1,d**2
n_total(i)=n_up(2*i-1)+n_down(2*i)
	enddo
	return
	end

	!**************************************************************************
	!subroutine for calculation of total no electrons of the cluster
	!**************************************************************************
subroutine clus_avg_n(EF,d1_c)
	use array
	!	use input
	use global
	implicit none
double precision, allocatable::n_up(:),n_down(:)
	integer ::i,j,p,x,y,l,ll
	double precision ::sum_nu,sum_nd,s_nup,s_ndown,EF,d1_c
allocate(n_up(2*d**2),n_down(2*d**2))
	n_up=0.0d0
	n_down=0.0d0

	s_nup=0.0d0
	ll=0
	do i=1,2*dc**2,2
	ll=ll+1
	sum_nu=0.0d0
	do j=1,2*dc**2
	sum_nu=sum_nu+(H_c(i,j)*conjg(H_c(i,j))*(1.0d0/(exp((evl_c(j)-EF)/T)+1.0d0)))
s_nup=s_nup+(H_c(i,j)*conjg(H_c(i,j))*(1.0d0/(exp((evl_c(j)-EF)/T)+1.0d0)))
	enddo                               !sum_nu=total no of up electrons
	n_up(ll)=sum_nu
	!	  n_up(i)=sum_nu
	enddo

	s_ndown=0.0d0
	ll=0
	do i=2,2*dc**2,2
	ll=ll+1
	sum_nd=0.0d0
	do j=1,2*dc**2
	sum_nd=sum_nd+(H_c(i,j)*conjg(H_c(i,j))*(1.0d0/(exp((evl_c(j)-EF)/T)+1.0d0)))
s_ndown=s_ndown+(H_c(i,j)*conjg(H_c(i,j))*(1.0d0/(exp((evl_c(j)-EF)/T)+1.0d0)))
	enddo
	n_down(ll)=sum_nd           !sum_nd=total no of down electrons
	!         n_down(i)=sum_nd  
	enddo
	!	d1_c= s_ndown+s_nup
	write(37,*)s_ndown,s_nup,s_ndown+s_nup,s_ndown-s_nup
	do i=1,dc**2
!	 n_total_c(i)=n_up(2*i-1)+n_down(2*i)

n_total_c(i)=n_up(i)+n_down(i)
	enddo

	n_d_pls_c=0.0d0
	n_d_mns_c=0.0d0
	l=1
	do y=1,dc
	do x=1,dc
	p=(-1)**(x+y)
	if(p.eq.1)n_d_pls_c=n_d_pls_c+n_up(l)+n_down(l)  
if(p.eq.-1)n_d_mns_c=n_d_mns_c+n_up(l)+n_down(l)
	l=l+1
	enddo
	enddo



	return
	end
	!***********************************************************************
	!subroutine for mapping back of cluster to system levels 
	!*************************************************************************

subroutine mapping_levels_clus_sys(rf)
	!use input
	use array
	use global
	implicit none
	integer::ls,rf,q,c,p,i,j,pbcx,x,y,k,s,pbcy,yc,xc,l




	ls=rf
	q=0
	i=0
	x=pi(rf)               !x=pi(rf) > x levels in the system
	y=pj(rf)                !y=pj(rf) > y levels in the system
!write(44,*)rf,pi(rf),pj(rf)
	do c=1,dc**2
	xc=pci(c)
	yc=pcj(c)
!write(45,*)c,pci(c),pcj(c)
	q=q+1
	if(((d-(x-1)).lt.dc).and.((x-1+xc).gt.d))then
	pbcx=-d                                        !pbcx > periodic boundary condition along x
	else
	pbcx=0
	endif
	if(((d-(y-1)).lt.dc).and.((y-1+yc).gt.d))then
	pbcy=-d**2                                      !pbcy > periodic boundary condition along y
	else
	pbcy=0
	endif
	p=ls+(q-1)+pbcx+pbcy                       ! p > variable for mapping levels from system to cluster including PBC
	n_total(p)=n_total_c(c)
!print*,c,p,m_c(c),th_c(c),ph_c(c)
	if (mod(c,dc).eq.0)then
	q=0
	i=i+1
	ls=rf+d*i
	endif
!write(32,*)c,p,n_total(p)

	enddo

	return
	end
	!***********************************************************************
	!subroutine for calculation of sub-lattice density
	subroutine cal_n
	use array
	use global
	implicit none
	integer i

	n_d_pls=0.0d0
	n_d_mns=0.0d0
	do i=1,(d**2)/2
	n_d_pls=n_d_pls+n_total(2*i-1)
n_d_mns=n_d_mns+n_total(2*i)
	enddo

	return
	end

	!***********************************************************************
	!subroutine for calculating density in the cluster

	!	subroutine clus_density
	!	use array
	!	use global
	!	implicit none
	!	integer
	!	double precision

	!	do i=1,dc**2
	!	 do a=1,2*dc**2
!	  nup=nup+(H_c(2*i-1,a)*conjg(H_c(2*i-1,a))*(1.0d0/(exp((evl_c(j)-EF)/T)+1.0d0)))
	!	 

	!	 enddo

	!	enddo



	!return
	!end
