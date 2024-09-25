CC_use_mod_instead_c --- ------------------------------------------------------------------
CC_use_mod_instead_c --- common blocks related to the dynamic/thermodynamic core
CC_use_mod_instead_c --- ------------------------------------------------------------------
CC_use_mod_instead_c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  viceold,vicenew
      common /ice1/ viceold, vicenew 
c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
CC_use_mod_instead_     .  u,v,           ! velocity components
CC_use_mod_instead_     .  dp,            ! layer thickness
CC_use_mod_instead_     .  dpold,         ! layer thickness at old time level
CC_use_mod_instead_     .  dpu,dpv,       ! layer thickness at u- and v-points
CC_use_mod_instead_     .  temp,          ! temperature
CC_use_mod_instead_     .  saln,          ! salinity
CC_use_mod_instead_     .  sigma          ! potential density
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
CC_use_mod_instead_     .  p,             ! interface pressure
CC_use_mod_instead_     .  pu,pv,         ! interface pressure at u- and v-points
CC_use_mod_instead_     .  phi            ! interface geopotential
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
CC_use_mod_instead_     .  sigmar,        ! reference potential density
CC_use_mod_instead_     .  temmin,        ! minimum temperature allowed in an isopycnic layer
CC_use_mod_instead_     .  dpuold,dpvold, ! layer thickness at u- and v-points at old time level
CC_use_mod_instead_     .  told,          ! temperature at old time level
CC_use_mod_instead_     .  sold,          ! salinity at old time level
CC_use_mod_instead_     .  diaflx         ! time integral of diapycnal flux
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
CC_use_mod_instead_     .  corioq,        ! coriolis parameter at q-point
CC_use_mod_instead_     .  coriop,        ! coriolis parameter at p-point
CC_use_mod_instead_     .  betafp,        ! latitudinal variation of the coriolis param. at p-point
CC_use_mod_instead_     .  potvor         ! potential vorticity
CC_use_mod_instead_c
CC_use_mod_instead_      common /micom1/ u,v,dp,dpold,dpu,dpv,temp,saln,sigma,p,pu,pv,phi,
CC_use_mod_instead_     .                sigmar,temmin,dpuold,dpvold,told,sold,
CC_use_mod_instead_     .                diaflx,corioq,coriop,betafp,potvor
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
CC_use_mod_instead_     .  uflx,vflx,     ! horizontal mass fluxes
CC_use_mod_instead_     .  utflx,vtflx,   ! horizontal heat fluxes
CC_use_mod_instead_     .  usflx,vsflx    ! horizontal salt fluxes
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
CC_use_mod_instead_     .  umfltd,vmfltd, ! horizontal mass fluxes due to thickness diffusion
CC_use_mod_instead_     .  utfltd,vtfltd, ! horizontal heat fluxes due to thickness diffusion
CC_use_mod_instead_     .  utflld,vtflld, ! horizontal heat fluxes due to lateral diffusion
CC_use_mod_instead_     .  usfltd,vsfltd, ! horizontal salt fluxes due to thickness diffusion
CC_use_mod_instead_     .  usflld,vsflld  ! horizontal salt fluxes due to lateral diffusion
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3) ::
CC_use_mod_instead_     .  ubflxs,vbflxs  ! barotropic mass flux sums
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
CC_use_mod_instead_     .  pb,            ! bottom pressure
CC_use_mod_instead_     .  ubflx,vbflx,   ! barotropic mass fluxes
CC_use_mod_instead_     .  pb_mn,         ! bottom pressure
CC_use_mod_instead_     .  ubflx_mn,vbflx_mn, ! barotropic mass fluxes
CC_use_mod_instead_     .  pbu,pbv,       ! bottom pressure at velocity points
CC_use_mod_instead_     .  ub,vb,         ! barotropic velocity components
CC_use_mod_instead_     .  ubflxs_p,vbflxs_p, ! predicted barotropic mass flux sums
CC_use_mod_instead_     .  pvtrop         ! potential vorticity of barotropic flow
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
CC_use_mod_instead_     .  pb_p,              ! predicted bottom pressure
CC_use_mod_instead_     .  pbu_p,pbv_p,       ! predicted bottom pressure at velocity points
CC_use_mod_instead_     .  pvtrop_o,      ! potential vorticity of barotropic flow at old time lev.
CC_use_mod_instead_     .  ubcors_p,vbcors_p, ! predicted sums of barotropic coriolis terms
CC_use_mod_instead_     .  defor1,defor2, ! deformation components
CC_use_mod_instead_     .  utotm,vtotm,   ! total (barotropic+baroclinic) ...
CC_use_mod_instead_     .  utotn,vtotn,   ! ... velocities at 2 time levels
CC_use_mod_instead_     .  uflux,vflux,   ! horizontal mass fluxes
CC_use_mod_instead_     .  uflux2,vflux2, ! more mass fluxes
CC_use_mod_instead_     .  uflux3,vflux3  ! more mass fluxes
CC_use_mod_instead_c
CC_use_mod_instead_      common /micom2/ uflx,vflx,utflx,vtflx,usflx,vsflx,umfltd,vmfltd,
CC_use_mod_instead_     .                utfltd,vtfltd,utflld,vtflld,usfltd,vsfltd,usflld,
CC_use_mod_instead_     .                vsflld,ubflxs,vbflxs,pb,ubflx,vbflx,pb_mn,
CC_use_mod_instead_     .                ubflx_mn,vbflx_mn,pbu,pbv,ub,vb,ubflxs_p,vbflxs_p,
CC_use_mod_instead_     .                pvtrop,pb_p,pbu_p,pbv_p,pvtrop_o,ubcors_p,
CC_use_mod_instead_     .                vbcors_p,defor1,defor2,utotm,vtotm,utotn,vtotn,
CC_use_mod_instead_     .                uflux,vflux,uflux2,vflux2,uflux3,vflux3
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
CC_use_mod_instead_     .  pgfx,pgfy      ! horizontal pressure gradient force
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
CC_use_mod_instead_     .  pgfxo,pgfyo    ! horizontal pressure gradient force at old time level
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
CC_use_mod_instead_     .  pgfxm,pgfym,   ! PGF terms not dependent on barotropic bottom pressure
CC_use_mod_instead_     .  xixp,xixm,     ! PGF terms dependent on barotropic bottom pressure
CC_use_mod_instead_     .  xiyp,xiym      ! PGF terms dependent on barotropic bottom pressure
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
CC_use_mod_instead_     .  pgfxm_o,pgfym_o,   ! PGF terms not dep. on barotr. bot. pres. at old ti.
CC_use_mod_instead_     .  xixp_o,xixm_o,     ! PGF terms dep. on barotr. bot. pres. at old tim. l.
CC_use_mod_instead_     .  xiyp_o,xiym_o,     ! PGF terms dep. on barotr. bot. pres. at old tim. l.
CC_use_mod_instead_     .  util1,util2,   ! arrays for temporary storage
CC_use_mod_instead_     .  util3,util4,   ! arrays for temporary storage
CC_use_mod_instead_     .  scqx,scqy,     ! mesh size at q-points in x,y direction
CC_use_mod_instead_     .  scpx,scpy,     ! mesh size at p-points in x,y direction
CC_use_mod_instead_     .  scux,scuy,     ! mesh size at u-points in x,y direction
CC_use_mod_instead_     .  scvx,scvy,     ! mesh size at v-points in x,y direction
CC_use_mod_instead_     .  scq2,scp2,     ! grid box size at q- and p-points
CC_use_mod_instead_     .  scu2,scv2,     ! grid box size at u- and v-points
CC_use_mod_instead_     .  scq2i,scp2i,   ! inverses of scq2,scp2
CC_use_mod_instead_     .  scuxi,scvyi,   ! inverses of scux,scvy
CC_use_mod_instead_     .  scuyi,scvxi,   ! inverses of scuy,scvx
CC_use_mod_instead_     .  umax,vmax,     ! maximum allowable velocities
CC_use_mod_instead_     .  depths         ! water depth
CC_use_mod_instead_c
CC_use_mod_instead_      common /micom3/ pgfx,pgfy,pgfxo,pgfyo,pgfxm,pgfym,xixp,xixm,
CC_use_mod_instead_     .                xiyp,xiym,pgfxm_o,pgfym_o,xixp_o,xixm_o,
CC_use_mod_instead_     .                xiyp_o,xiym_o,util1,util2,util3,util4,scqx,scqy,
CC_use_mod_instead_     .                scpx,scpy,scux,scuy,scvx,scvy,scq2,scp2,scu2,scv2,
CC_use_mod_instead_     .                scq2i,scp2i,scuxi,scvyi,scuyi,scvxi,umax,vmax,
CC_use_mod_instead_     .                depths
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
CC_use_mod_instead_     .  difint,        ! layer interface diffusivity
CC_use_mod_instead_     .  difiso,        ! isopycnal diffusivity
CC_use_mod_instead_     .  difdia         ! diapycnal diffusivity
CC_use_mod_instead_c
CC_use_mod_instead_      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
CC_use_mod_instead_     .  uja,ujb,       ! velocities at lateral ...
CC_use_mod_instead_     .  via,vib,       ! ... neighbor points
CC_use_mod_instead_     .  difmxp,        ! maximum lateral diffusivity at p-points
CC_use_mod_instead_     .  difmxq,        ! maximum lateral diffusivity at q-points
CC_use_mod_instead_     .  difwgt,        ! eddy diffusivity weight
CC_use_mod_instead_     .  sealv,         ! sea surface height
CC_use_mod_instead_     .  surflx,        ! surface thermal energy flux
CC_use_mod_instead_     .  surrlx,        ! surface relaxation thermal energy flux
CC_use_mod_instead_     .  sswflx,        ! surface solar energy flux
CC_use_mod_instead_     .  salflx,        ! surface salinity flux
CC_use_mod_instead_     .  brnflx,        ! surface brine flux
CC_use_mod_instead_     .  salrlx,        ! surface relaxation salinity flux
CC_use_mod_instead_     .  taux,tauy,     ! surface stress components
CC_use_mod_instead_     .  ustar,         ! surface friction velocity
CC_use_mod_instead_     .  ustarb,        ! bottom friction velocity
CC_use_mod_instead_     .  buoyfl,        ! surface buoyancy flux
CC_use_mod_instead_     .  twedon,        ! tidal wave energy diffipation over buoyancy frequency
CC_use_mod_instead_     .  pbrnda         ! brine plume pressure depth
CC_use_mod_instead_c
CC_use_mod_instead_      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
CC_use_mod_instead_     .  kfpla          ! index of first physical layer
CC_use_mod_instead_c
CC_use_mod_instead_      common /micom4/ difint,difiso,difdia,uja,ujb,via,vib,difmxp,
CC_use_mod_instead_     .                difmxq,difwgt,sealv,surflx,surrlx,sswflx,salflx,
CC_use_mod_instead_     .                brnflx,salrlx,taux,tauy,ustar,ustarb,buoyfl,
CC_use_mod_instead_     .                twedon,pbrnda,kfpla
CC_use_mod_instead_c
CC_use_mod_instead_      real time,delt1,dlt,area,avgbot
CC_use_mod_instead_      integer nstep,nstep1,nstep2,lstep
CC_use_mod_instead_c
CC_use_mod_instead_      common /varbls/ time,delt1,dlt,area,avgbot,
CC_use_mod_instead_     .                nstep,nstep1,nstep2,lstep
CC_use_mod_instead_c
CC_use_mod_instead_c --- 'baclin' = baroclinic time step
CC_use_mod_instead_c --- 'batrop' = barotropic time step
CC_use_mod_instead_c --- 'mdv2hi' = Laplacian diffusion velocity (cm/s) for momentum dissipation
CC_use_mod_instead_c --- 'mdv2lo' = same as mdv2hi but used when Rossby radius is resolved
CC_use_mod_instead_c --- 'mdv4hi' = Biharmonic diffusion velocity (cm/s) for momentum dissipation
CC_use_mod_instead_c --- 'mdv4lo' = same as mdv4hi but used when Rossby radius is resolved
CC_use_mod_instead_c --- 'mdc2hi' = Laplacian diffusivity (cm**2/s) for momentum dissipation
CC_use_mod_instead_c --- 'mdc2lo' = same as mdc2hi but used when Rossby radius is resolved
CC_use_mod_instead_c --- 'vsc2hi' = parameter used in deformation-dependent Laplacian viscosity
CC_use_mod_instead_c --- 'vsc2lo' = same as vsc2hi but used when Rossby radius is resolved
CC_use_mod_instead_c --- 'vsc4hi' = parameter used in deformation-dependent Biharmonic viscosity
CC_use_mod_instead_c --- 'vsc4lo' = same as vsc4hi but used when Rossby radius is resolved
CC_use_mod_instead_c --- slip = +1  for free-slip boundary cond., slip = -1  for non-slip cond.
CC_use_mod_instead_c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
CC_use_mod_instead_c --- 'cb'     = coefficient of quadratic bottom friction
CC_use_mod_instead_c --- 'cwbdts' = coastal wave breaking damping resiprocal time scale (1/s)
CC_use_mod_instead_c --- 'cwbdls' = coastal wave breaking damping length scale (m)
CC_use_mod_instead_c --- 'wuv1/2' = weights for time smoothing of u,v field
CC_use_mod_instead_c --- 'wts1/2' = weights for time smoothing of t,s field
CC_use_mod_instead_c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
CC_use_mod_instead_c --- 'wpgf'   = weight for time averaging of pressure gradient force
CC_use_mod_instead_c --- 'mltmin' = minimum mixed-layer thickness (m)
CC_use_mod_instead_c --- 'thktop' = thickness of top layer (m)
CC_use_mod_instead_c --- 'thkbot' = thickness of bottom boundary layer (pressure units)
CC_use_mod_instead_c --- 'acurcy' = permissible roundoff error in column integral calc.
CC_use_mod_instead_c --- 'egc'    = the parameter c in the Eden and Greatbatch (2008)
CC_use_mod_instead_c ---            parameterization
CC_use_mod_instead_c --- 'eggam'  = the parameter gamma in the Eden and Greatbatch (2008)
CC_use_mod_instead_c ---            parameterization [].
CC_use_mod_instead_c --- 'egmndf' = minimum diffusivity in the Eden and Greatbatch (2008)
CC_use_mod_instead_c ---            parameterization [cm**2/s].
CC_use_mod_instead_c --- 'egmxdf' = maximum diffusivity in the Eden and Greatbatch (2008)
CC_use_mod_instead_c ---            parameterization [cm**2/s]
CC_use_mod_instead_c --- 'egidfq' = factor relating the isopycnal diffusivity to the layer
CC_use_mod_instead_c ---            interface diffusivity in the Eden and Greatbatch (2008)
CC_use_mod_instead_c ---            parameterization. egidfq=difint/difiso
CC_use_mod_instead_c --- 'csdiag' = if set to .true., then output check sums
CC_use_mod_instead_c --- 'cnsvdi' = if set to .true., then output conservation diagnostics
CC_use_mod_instead_c
CC_use_mod_instead_      real baclin,batrop,mdv2hi,mdv2lo,mdv4hi,mdv4lo,mdc2hi,mdc2lo,
CC_use_mod_instead_     .     vsc2hi,vsc2lo,vsc4hi,vsc4lo,slip,cbar,cb,cwbdts,cwbdls,
CC_use_mod_instead_     .     wuv1,wuv2,wts1,wts2,wbaro,wpgf,mltmin,thktop,thkbot,acurcy,
CC_use_mod_instead_     .     egc,eggam,egmndf,egmxdf,egidfq
CC_use_mod_instead_      logical csdiag,cnsvdi
CC_use_mod_instead_c
CC_use_mod_instead_      common /parms1/ baclin,batrop,mdv2hi,mdv2lo,mdv4hi,mdv4lo,
CC_use_mod_instead_     .                mdc2hi,mdc2lo,vsc2hi,vsc2lo,vsc4hi,vsc4lo,slip,
CC_use_mod_instead_     .                cbar,cb,cwbdts,cwbdls,wuv1,wuv2,wts1,wts2,wbaro,
CC_use_mod_instead_     .                wpgf,mltmin,thktop,thkbot,acurcy,egc,eggam,
CC_use_mod_instead_     .                egmndf,egmxdf,egidfq,csdiag,cnsvdi
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'g'      = gravity acceleration
c --- 'rearth' = radius of the earth
c --- 'spcifh' = specific heat of sea water (j/g/deg)
c --- 'rhoa_r' = reference air density (g/cm**3)
c --- 'cd_r'   = reference transfer coefficient of momentum
c --- 'ch_r'   = reference transfer coefficient of sensible heat
c --- 'ce_r'   = reference transfer coefficient of tent heat
c --- 'wg2_r'  = reference gustiness squared (cm**2/s**2)
c --- 'alpha0' = reference value of specific volume (cm**3/g)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c --- 'raddep' = maximum depth of light penetration (m)
c --- 'redfac' = red fraction of light aborbed in mixed layer (jerlov 1)
c --- 'betabl' = blue light extinction coefficient (m) (jerlov 1)
c
      real tenm,onem,tencm,onecm,onemm,g,rearth,spcifh,rhoa_r,cd_r,ch_r,
     .     ce_r,wg2_r,alpha0,epsil,raddep,redfac,betabl,huge,radian,pi
c
      common /consts/ tenm,onem,tencm,onecm,onemm,g,rearth,spcifh,
     .                rhoa_r,cd_r,ch_r,ce_r,wg2_r,alpha0,epsil,raddep,
     .                redfac,betabl,huge,radian,pi
c
c --- grid point where detailed diagnostics are desired:
c
      integer itest,jtest,ptest
c
      common /testpt/ itest,jtest,ptest
c
      character*80 path,path1,path2,runid
      integer path_len,path1_len,path2_len,runid_len
c
      common /iovars/ path,path1,path2,runid,
     .                path_len,path1_len,path2_len,runid_len
c
c
c> Revision history:
c>
c-----------------------------------------------------------------------------
c
