import numpy as np
import random
import scipy.integrate as integ

#------ Some constants ------
clgt=2.99792458e8
unitM=9.11e-31
unitQ=1.602e-19
eps0=8.854e-12
Nbuf=1
F_T=4.0/3.0
T_T=2.0/3.0
#-------------------------------------



#----- load(): generate plasma particles  ------
# _n0: density in [m^-3], _sp: superParticle, _dx: mesh size [m], 
# _ncll: number of cells, 
# _prof: profile function - argument x is [m] 
# _T: temperature [eV], _m: mass [kg]. 
# Density at x is n0*prof(x), where x goes from 0 thru _ncll
# velocity distribution is Gaussian with temperature _T: ~exp[-v**2/vth**2]
# where 0.5*m*vth**2=kB*T, up to +-3*vth (3-sigma).
# The return is arrays of dx-normalized x's of particles and 
# initial position xi=x (it is used for current calculation in EM sim),
# velocities in x,y,z directions, corresponding gamma (relativistic), and 
# magnitude of velocity, and finally number of created particles.    
# the velocity is normalized to speed of light
def load(_n0,_sp,_dx,_ncll,_prof,_T,_m):
	_np = (int)(_n0*_dx/_sp*_ncll); _npr=0 
	_xx=[]; _xxi=[]; _uu1=[]; _uu2=[]; _uu3=[];  _gm=[];  _um=[];
	for _i in range(_np):
		_x = _ncll*random.random(); _j = (int)(_x)
		_den=_prof(_x)
		if _den > random.random():
			_xx.append(_x);  _xxi.append(_j); _uu=thermal(_T,_m)
			_uu1.append(_uu[0])
			_uu2.append(_uu[1])
			_uu3.append(_uu[2])
			_uusq = _uu[0]**2+_uu[1]**2+_uu[2]**2
			_um.append(np.sqrt(_uusq))
			_gm.append(1.0/np.sqrt(1.0-_uusq))  
			_npr +=1
	_xx=np.array(_xx); _xxi=np.array(_xxi); 
	_uu1=np.array(_uu1); _uu2=np.array(_uu2); _uu3=np.array(_uu3);
	_gm = np.array(_gm)
	_um = np.array(_um)
	return _xx,_xxi,_uu1,_uu2,_uu3,_gm,_um,_npr


#---------------- Creat particles for EM 1D simulation ---------------
def loadEM1(_m,_n0,_prof,_T,_npc,_dx,_ncll,_bf):
	_np = _npc*_ncll;  
	_nMm = int(_np*_bf)
	_x  =np.zeros(_nMm); _xo =np.zeros(_nMm); _xi =np.zeros(_nMm);
	_ux =np.zeros(_nMm); _uy =np.zeros(_nMm); _uz =np.zeros(_nMm); 
	_u  =np.zeros(_nMm); _gm =np.zeros(_nMm);  
	_exp = np.zeros(_nMm); _eyp = np.zeros(_nMm); _ezp = np.zeros(_nMm);
	_bxp = np.zeros(_nMm); _byp = np.zeros(_nMm); _bzp = np.zeros(_nMm);

	_xi = [0]*_nMm; _xi = np.array(_xi) 

	_mm = _m*unitM

	_id=0
	for _i in range(_np):
		_xx = _ncll*random.random(); _j = (int)(_xx)
		_den=_prof(_xx*_dx)

		if _den > random.random():
			_x[_id] = _xx-_j;  _xo[_id] = _x[_id]; _xi[_id] = _j; 
			_uu=thermal(_T,_mm)
			_ux[_id] =_uu[0]; _uy[_id] =_uu[1]; _uz[_id] =_uu[2]
			_usq = _uu[0]**2+_uu[1]**2+_uu[2]**2
			_u[_id] = np.sqrt(_usq); _gm[_id] = np.sqrt(1.0+_usq)
			_id +=1

	return _x,_xi,_xo,_ux,_uy,_uz,_gm,_u,_exp,_eyp,_ezp,_bxp,_byp,_bzp,_id


#---------------------- Creat arrays of fields -----------------------
def createFld1(_ncll):
	_nMm = _ncll + 1+ 2*Nbuf
	_ex  =np.zeros(_nMm); _ey  =np.zeros(_nMm); _ez  =np.zeros(_nMm); 
	_bx  =np.zeros(_nMm); _by  =np.zeros(_nMm); _bz  =np.zeros(_nMm); 
	_bxa =np.zeros(_nMm); _bya =np.zeros(_nMm); _bza =np.zeros(_nMm); 
	_jx  =np.zeros(_nMm); _jy  =np.zeros(_nMm); _jz  =np.zeros(_nMm); 
	return _ex,_ey,_ez,_bx,_by,_bz,_bxa,_bya,_bza,_jx,_jy,_jz


#----- generate plasma particles over a specified range  ------
# The same as load(), but only over a specified range [_rs,_re]
# and add them to the original particle list
# _rs, _re: lower and upper limit of the range (mesh unit)
# _p: position [m] corresponding to _rs 
# _npo: num. particles of original particle list
# _x,...,_gm: position and velocities of original particle list
def addPtcl(_n0,_sp,_dx,_p,_rs,_re,_prof,_T,_m,_npo,_xx,_xxo,_xxi,_uu,_ux,_uy,_uz,_gm):
	_l = _re-_rs;	
	_np = (int)(_n0*_dx/_sp*_l); _npr=0 
	for _i in range(_np):
		_x = _l*random.random(); _j = (int)(_x); _j += _rs
		_den=_prof(_x*_dx+_p)
		if _den > random.random():
			_id = _i+_npo
			_xx[_id]= _x;  _xxi[_id]=_j; _uv=thermal(_T,_m)
			_ux[_id] = _uv[0]
			_uy[_id] = _uv[1]
			_uz[_id] = _uv[2]
			_usq = _uv[0]**2+_uv[1]**2+_uv[2]**2
			_gm[_id] = 1.0/np.sqrt(1.0-_usq)
			_uu[_id] = np.sqrt(_usq)
			_npr +=1
	return _npr+_npo

#------ thermal(): returns Gaussian velocity ------
def thermal(_T,_m):
	if _T==0: return np.array([0.0, 0.0, 0.0])
	_vth=np.sqrt(2.0*unitQ*_T/_m)/clgt	
	_st=0
	while _st !=1:
		_eta=3.0*random.random()
		_pr=4.0*np.exp(-_eta**2)*(_eta)**2/np.sqrt(np.pi)
		if random.random() <= _pr: _st=1
	_eta *=_vth
	_v1=2.0*(random.random()-0.5)
	_v2=2.0*(random.random()-0.5)
	_v3=2.0*(random.random()-0.5)
	_sm=np.sqrt(_v1**2 + _v2**2 + _v3**2)
	_v1 *=_eta/_sm; _v2 *=_eta/_sm; _v3 *=_eta/_sm
	return np.array([_v1, _v2, _v3])


#--- linear interpolate field exerting on particles in EM simulations ---
# _np: number of simulation particles, 
# _xx, _xxi: dx-normalized position (Note: _x+_xxi is the position).
# _ex,... _bz: fields on meshes
# _exp ... _bzp: fields interpolated on particles. Return values
def fieldEM1lin(_np,_xx,_xxi,_ex,_ey,_ez,_bx,_by,_bz,_exp,_eyp,_ezp,_bxp,_byp,_bzp):
	for _id in range(_np):
		_x=_xx[_id];  _j=_xxi[_id] 
		_j+= Nbuf;   # lift index by amount of buffer  

		if _x<=0.5: _xa=_x+0.5; _is=0
		else:       _xa=_x-0.5; _is=1

		_p0=1.0-_xa; _p1=_xa; _q0=1.0-_x;  _q1=_x

		_jj=_j+_is

		_exp[_id] = _p0*_ex[_jj-1]+_p1*_ex[_jj] 
		_byp[_id] = _p0*_by[_jj-1]+_p1*_by[_jj] 
		_bzp[_id] = _p0*_bz[_jj-1]+_p1*_bz[_jj]
		_eyp[_id] = _q0*_ey[_j]+_q1*_ey[_j+1]
		_ezp[_id] = _q0*_ez[_j]+_q1*_ez[_j+1]
		_bxp[_id] = _q0*_bx[_j]+_q1*_bx[_j+1]
	return


#-------- Add external field to particles ---------------
def addExtFld(_np,_dx,_x,_xi,_exp,_eyp,_ezp,_bxp,_byp,_bzp,_ex0,_ey0,_ez0,_bx0,_by0,_bz0):
	for _id in range(_np):
		_xx = _dx*(_x[_id]+_xi[_id])
		_exp[_id] += _ex0(_xx);
		_eyp[_id] += _ey0(_xx);
		_ezp[_id] += _ez0(_xx);
		_bxp[_id] += _bx0(_xx);
		_byp[_id] += _by0(_xx);
		_bzp[_id] += _bz0(_xx);
	return


#----- push particles by Lorentz equation in EM simulations ------- 
# Input
# _np, _ncll: num. particles and cells
# _q,_m: charge [C] and mass[kg] of a single (not super) particle
# _w,_dx,_dt: EM angular freq [rad/s], mesh size [m], time step [s]
# _x,_xo,_xi: relative position, previous step value of _x, mesh index. 
#    Note: _x+_xi is the current position in unit of dx
# _gmm,_uu,_uxx,_uyy,_uzz: relativistic gamma factor, speed and velocity[c]
# _ex,... _bzz: electromagnetic fields interpolated on particles. 
# _optMW: option for moving window. 0 for off, 1 for on at the right edge
#  When it is on, merely the particles xi>ncll are preserved instead of 
#  being removed. 
# Output
# _gmm,_uu,_uxx,_uyy,_uzz: updated by Lorentz equation 
# _x,_xo,_xi: updated by dx/dt=u/gam 
def move(_np,_ncll,_q,_m,_w,_dx,_dt,_x,_xo,_xi,_gmm,_uu,_uxx,_uyy,_uzz,_ex,_ey,_ez,_bxx,_byy,_bzz,_optMW):
	qmhdt = 0.5*(_q/unitQ)*(unitM/_m)*(_w*_dt)
	_al = _dt*clgt/_dx
	_id=0;  
	_xo = _x
	while _id < _np:
		_ax = qmhdt*_ex[_id]; _ay =qmhdt*_ey[_id]; _az =qmhdt*_ez[_id]
		_ux = _uxx[_id]; _uy = _uyy[_id]; _uz = _uzz[_id]
		_bx = _bxx[_id]; _by = _byy[_id]; _bz = _bzz[_id];

		_ux += _ax;    _uy += _ay;    _uz += _az
		_gm =np.sqrt(1.0+_ux**2 +_uy**2 +_uz**2)

		_qdt_gm = qmhdt/_gm 
		_tx = _qdt_gm*_bx; _ty = _qdt_gm*_by; _tz = _qdt_gm*_bz
		_t2p1 = 2.0/(1.0 + _tx**2 +_ty**2 +_tz**2)

		_upx = _ux + _uy*_tz-_uz*_ty
		_upy = _uy + _uz*_tx-_ux*_tz
		_upz = _uz + _ux*_ty-_uy*_tx

		_ux += _t2p1*(_upy*_tz-_upz*_ty)
		_uy += _t2p1*(_upz*_tx-_upx*_tz)
		_uz += _t2p1*(_upx*_ty-_upy*_tx)

		_ux += _ax; _uy += _ay; _uz += _az
		_usq = _ux**2+ _uy**2 + _uz**2
		_gm = np.sqrt(1.0+ _usq)

		_gmm[_id] = _gm; _uxx[_id]=_ux; _uyy[_id]=_uy; _uzz[_id]=_uz
		_uu[_id] = np.sqrt(_usq)
		_x[_id] += (_ux/_gm)*_al 

		if _x[_id]>1:  _x[_id] -=1;  _xo[_id] -=1;  _xi[_id] +=1
		elif _x[_id]<0:  _x[_id] +=1;  _xo[_id] +=1;  _xi[_id] -=1

		if _xi[_id]<0:
			if _id <_np-1:
				_x[_id] = _x[_np-1];  _xi[_id] = _xi[_np-1]
				_uxx[_id]=_uxx[_np-1]; _uyy[_id]=_uyy[_np-1]
				_uzz[_id]=_uzz[_np-1]; _gmm[_id]=_gmm[_np-1]
				_xo[_id]=_xo[_np-1]; _np -=1
				_uu[_id] = _uu[_np-1]
			else:   _np -=1;  return _np
		else: _id+=1

		if _xi[_id]>=_ncll and _optMW==0:
			if _id <_np-1:
				_x[_id] = _x[_np-1];  _xi[_id] = _xi[_np-1]
				_uxx[_id]=_uxx[_np-1]; _uyy[_id]=_uyy[_np-1]
				_uzz[_id]=_uzz[_np-1]; _gmm[_id]=_gmm[_np-1]
				_xo[_id]=_xo[_np-1]; _np -=1
				_uu[_id] = _uu[_np-1]
			else:   _np -=1;  return _np
		else: _id+=1
	return _np


#----- moveES(): push particles electrostatic environment  ------- 
# Difference from move() is just that moveES() does not track xo 
# (old position), which is used for J in electromagnetic sim. 
# _np: num. particles, _ncll: num. cells, _q,_m: charge and mass of particles
# _wp: plasma freq. used for time normalization, _dx: mesh size [m] - before
# normalization, _dt: time step [s] - before normalization, 
# _x, _gmm, _uu, _uxx, _uyy,_uzz: arrays of position (dx-normalized), gamma
# factor, velocity magnitude and components (c-normalized),
# _ex,_ey,_ez,_bxx,_byy,_bzz: arrays of fields, 
# _lBcOpt, _rBcOpt: boundary condtion at left and right. 0 means reflection,
# 1 means absorbing 
def moveES(_np,_ncll,_q,_m,_wp,_dx,_dt,_x,_gmm,_uu,_uxx,_uyy,_uzz,_ex,_ey,_ez,_bxx,_byy,_bzz,_lBcOpt,_rBcOpt):
	qmhdt = 0.5*(_q/unitQ)*(unitM/_m)*(_wp*_dt)
	_al = _dt*clgt/_dx
	_id=0;  
	while _id < _np:
		_ax = qmhdt*_ex[_id]; _ay =qmhdt*_ey[_id]; _az =qmhdt*_ez[_id]
		_ux = _uxx[_id]; _uy = _uyy[_id]; _uz = _uzz[_id]
		_bx = _bxx[_id]; _by = _byy[_id]; _bz = _bzz[_id];

		_ux += _ax;    _uy += _ay;    _uz += _az
		_gm =np.sqrt(1.0+_ux**2 +_uy**2 +_uz**2)

		_qdt_gm = qmhdt/_gm 
		_tx = _qdt_gm*_bx; _ty = _qdt_gm*_by; _tz = _qdt_gm*_bz
		_t2p1 = 2.0/(1.0 + _tx**2 +_ty**2 +_tz**2)

		_upx = _ux + _uy*_tz-_uz*_ty
		_upy = _uy + _uz*_tx-_ux*_tz
		_upz = _uz + _ux*_ty-_uy*_tx

		_ux += _t2p1*(_upy*_tz-_upz*_ty)
		_uy += _t2p1*(_upz*_tx-_upx*_tz)
		_uz += _t2p1*(_upx*_ty-_upy*_tx)

		_ux += _ax; _uy += _ay; _uz += _az
		_usq = _ux**2 + _uy**2 + _uz**2
		_gm = np.sqrt(1.0+ _usq)

		_gmm[_id] = _gm; _uxx[_id]=_ux; _uyy[_id]=_uy; _uzz[_id]=_uz
		_uu[_id] = np.sqrt(_usq)
		_x[_id] += (_ux/_gm)*_al 
		
		if _x[_id]<0:
			if _lBcOpt ==0:   # absorb at the left
				if _id <_np-1:
					_x[_id] = _x[_np-1]
					_uxx[_id]=_uxx[_np-1]
					_uyy[_id]=_uyy[_np-1]
					_uzz[_id]=_uzz[_np-1]
					_gmm[_id]=_gmm[_np-1]; _np -=1
					_uu[_id] = _uu[_np-1]
				else:   _np -=1;  return _np
			elif _lBcOpt==1:  # reflection   
				_x[_id] *= -1.0;  _uxx[_id] *= -1.0; _id +=1
			elif _lBcOpt==2: # appears at the other side (periodic)
				_x[_id] += _ncll;  _id +=1
			else: _id +=1

		elif _x[_id]>=_ncll:
			if _rBcOpt==0:  # absorbe at the right
				if _id <_np-1:
					_x[_id] = _x[_np-1]
					_uxx[_id]=_uxx[_np-1]
					_uyy[_id]=_uyy[_np-1]
					_uzz[_id]=_uzz[_np-1]
					_gmm[_id]=_gmm[_np-1]; _np -=1
					_uu[_id] = _uu[_np-1]
				else:   _np -=1;  return _np
			elif _rBcOpt==1:  # reflection 
				_x[_id]=2.0*_ncll-_x[_id]; _uxx[_id]*= -1.0; _id+=1
			elif _rBcOpt==2: # appears at the other side (periodic)
				_x[_id] -= _ncll;  _id +=1
			else: _id +=1
		else: _id+=1
	return _np


#----------------- Move particles by analytic field ------------------
# Written by MSHur, 20210307
# Push particles by analytic fields, which are given by functions.
# Useful for single particle simulations under prescribed field. 
# _np: num. simulation ptcls. Multiple particles can be traced simultaneously.
# _q,_m: q/e, m/me, i.e. normalized charge and mass of a single particle.
# _dt: w-nrm time step
# _x,_y,_z: k-nrm particle position (or its array when multiple ptcls)
# _gmm,_uu,_uxx,_uyy,_uzz: relativistic gamma and c-nrm speed and vel.
# _ef, _bf: functions for electric and magnetic fields, respectively. 
#           their arguments should be (x,y,z) (k-nrm) and return
#           (Ex,Ey,Ez) and (Bx,By,Bz) respectively.
def moveAF(_np,_q,_m,_dt,_x,_y,_z,_gmm,_uu,_uxx,_uyy,_uzz,_ef,_bf):
	qmhdt = 0.5*(_q/_m)*_dt
	_id=0;  
	while _id < _np:
		_ex,_ey,_ez = _ef(_x[_id],_y[_id],_z[_id])
		_bx,_by,_bz = _bf(_x[_id],_y[_id],_z[_id])
		#_ex,_ey,_ez = _ef(_x,_y,_z)
		#_bx,_by,_bz = _bf(_x,_y,_z)

		_ux = _uxx[_id]; _uy = _uyy[_id]; _uz = _uzz[_id]

		_ax = qmhdt*_ex; _ay =qmhdt*_ey; _az =qmhdt*_ez

		_ux += _ax;    _uy += _ay;    _uz += _az
		_gm =np.sqrt(1.0+_ux**2 +_uy**2 +_uz**2)

		_qdt_gm = qmhdt/_gm 
		_tx = _qdt_gm*_bx; _ty = _qdt_gm*_by; _tz = _qdt_gm*_bz
		_t2p1 = 2.0/(1.0 + _tx**2 +_ty**2 +_tz**2)

		_upx = _ux + _uy*_tz-_uz*_ty
		_upy = _uy + _uz*_tx-_ux*_tz
		_upz = _uz + _ux*_ty-_uy*_tx

		_ux += _t2p1*(_upy*_tz-_upz*_ty)
		_uy += _t2p1*(_upz*_tx-_upx*_tz)
		_uz += _t2p1*(_upx*_ty-_upy*_tx)

		_ux += _ax; _uy += _ay; _uz += _az
		_usq = _ux**2+ _uy**2 + _uz**2
		_gm = np.sqrt(1.0+ _usq)

		_gmm[_id] = _gm; _uxx[_id]=_ux; _uyy[_id]=_uy; _uzz[_id]=_uz
		_uu[_id] = np.sqrt(_usq)

		_x[_id] +=(_ux/_gm)*_dt 
		_y[_id] +=(_uy/_gm)*_dt 
		_z[_id] +=(_uz/_gm)*_dt 

		_id +=1

	return _np


#----------------- Magnetic field by a circular loop -----------------
# Written by MSHur, 20210307
# Calculate the magnetic fields generated by current in a circular loop.
# Not exactly satisfies Biot-Savart law, but use approximate expression
# (with A=1) in the appendix of note PP2-2.
# _r: k-nrm radius from axis, which is always the z-axis.
# _z: k-nrm z, relative to the center of the loop _z0. 
# _z0: z (k-nrm) of the cntr of the loop. Note that r-position is always 0. 
# _a: radius of the loop (k-nrm).
# Br, Bz should be always used together to assure div.B=curl.B=0.
def Br(_r,_z,_z0,_a):
	_br=1.5*(_z-_z0)*_r/((_z-_z0)**2+_a**2)**2.5
	return _br

def Bz(_r,_z,_z0,_a):
	_zz = _z-_z0
	_bz=1/(_zz**2+_a**2)**1.5
	_fc2 = 1.0 - 5.0*_zz**2/(_zz**2+_a**2)
	_fc =  1.0 + 0.75*_r**2*_fc2/(_zz**2+_a**2)
	_bztot = _bz*_fc 
	return _bztot


#---------------- Magnetic field by an infinite wire  ----------------
# Written by MSHur, 20210307
# Calculate the magnetic field generated from current in an infinite wire.
# Simply use the expression that B_\phi ~ 1/r.
# _x,_y,_z: positions to evaluate the field (k-nrm)
# _xc,_yc,_zc: a point which the wire passes through. k-nrm
# _ax: direction of the wire passing thru (_xc,_yc,_zc). 0,1,2 means
#      x, y, and z direction, respectively. No oblique case allowed.  
# Returns Bx,By,Bz (Cartesian).
def BwireXYZ(_x,_y,_z,_xc,_yc,_zc, _ax):

	if _ax ==0:   # wire is along x
		_r = np.sqrt((_y -_yc)**2 + (_z-_zc)**2); _b = 1.0/_r
		_bx = 0.0;  _by = -_b*(_z/_r);  _bz = _b*(_y/_r)

	elif _ax ==1:   # wire is along y
		_r = np.sqrt((_z -_zc)**2 + (_x-_xc)**2); _b = 1.0/_r
		_bx = _b*(_z/_r);  _by = 0.0;  _bz = -_b*(_x/_r) 

	elif _ax ==2:   # wire is along z
		_r = np.sqrt((_x -_xc)**2 + (_y-_yc)**2); _b = 1.0/_r
		_bx = -_b*(_y/_r);  _by = _b*(_x/_r);  _bz = 0.0

	else:  print("Axis should be one of 0,1,2 in BwireXYZ")

	return _bx, _by, _bz


#------------ Electric field by an infinite line chargee  ------------
# Written by MSHur, 20210307
# Calculate the electric field generated from a line charge.
# Simply use the expression that E_r ~ 1/r.
# _x,_y,_z: positions to evaluate the field (k-nrm)
# _xc,_yc,_zc: a point which the wire passes through. k-nrm
# _ax: direction of the wire passing thru (_xc,_yc,_zc). 0,1,2 means
#      x, y, and z direction, respectively. No oblique case allowed.  
# Returns Ex,Ey,Ez (Cartesian).
def ElineXYZ(_x,_y,_z,_xc,_yc,_zc, _ax):

	if _ax ==0:   # line is along x
		_r = np.sqrt((_y -_yc)**2 + (_z-_zc)**2); _e = 1.0/_r
		_ex = 0.0;  _ey = _e*(_y/_r);  _ez = _e*(_z/_r)

	elif _ax ==1:   # line is along y
		_r = np.sqrt((_z -_zc)**2 + (_x-_xc)**2); _e = 1.0/_r
		_ex = _e*(_x/_r);  _ey = 0.0;  _ez = _e*(_z/_r) 

	elif _ax ==2:   # line is along z
		_r = np.sqrt((_x -_xc)**2 + (_y-_yc)**2); _e = 1.0/_r
		_ex = _e*(_x/_r);  _ey = _e*(_y/_r);  _ez = 0.0

	else:  print("Axis should be one of 0,1,2 in ElineXYZ")

	return _ex, _ey, _ez


#-------------- Field line constructor from x,y,z fields -------------
# Written by MSHur, 20210307
# Construct a field line starting from a point, for given analytic fields.  
# _x0,_y0,_z0: a point which the field line begins from.
# _ds: step size to solve the ODE, which gets the field line. 
# _smax: max. step to solve the ODE. The end position of the field line is
#        determined by this.
# _fxyz: analytic function that represents the field, i.e. F\vec(x,y,z).
#       Its returning values are Fx,Fy,Fz (vector fields) 
def fLineXYZ(_x0,_y0,_z0,_ds, _smax,_fxyz):
	def _f(_rr,_tt):
		_xx,_yy,_zz=_rr
		_xxp,_yyp,_zzp = _fxyz(_xx,_yy,_zz)
		return [_xxp,_yyp,_zzp]

	_r0=[_x0, _y0, _z0]
	_s = np.arange(0, _smax*_ds, _ds)
	sol = integ.odeint(_f,_r0,_s)	

	return sol


#----- density()  ------		
# Calculate number density from given array of position 
# _den: density calculation is stored here. - num per length [m^-1], 
# not [m^-3], because it is 1D. _ncll: num. cells, 
# _np: num. particles, _sp: superParticle
# _dx: mesh size [m], _xx: array of particle positions (dx-normalized)
def density(_den,_ncll,_np,_sp,_dx,_xx):
	_den *=0.0
	for _i in range(_np): 
		_j = (int)(_xx[_i])
		if _j<_ncll and _j>=0:
			_dw = _xx[_i]-_j 
			_den[_j] += (1.0-_dw)
			_den[_j+1] += _dw
	_den *= (_sp/_dx)
	return


#----- Poisson()  ------		
# solve Poisson equation by directly inverting tridiagonal matrix
# _nmsh: num. meshes = _ncll+1,  _rh: charge density, _phi: potential
# It just solves _phi[j+1]-2.0*_phi[j]+_phi[j-1]=-_rh[j]
# How to scale _rh and _phi is up to the user
def Poisson(_nmsh,_rh,_phi):
	_A=np.full(_nmsh,-2.0)
	for _i in range(_nmsh-3):  _A[_i+2] -= 1.0/_A[_i+1]
	_v0 = _phi[0]; _vN= _phi[_nmsh-1]
	_rh[1] += _v0;  _rh[_nmsh-2] += _vN
	for _i in range(_nmsh-3):  _rh[_i+2] -= _rh[_i+1]/_A[_i+1]
	_i=_nmsh-2; _phi[_i] = -_rh[_i]/_A[_i]
	while _i>1:  
		_phi[_i-1] = (-_rh[_i-1]- _phi[_i])/_A[_i-1];  _i -=1
	return



#---------  getPP(): get second derivative for cubicspline ---------
def getPP(_n,_y,_ypp):
	B=np.full(_n,4.0)
	for _i in range(_n-3): B[_i+2] -= 1.0/B[_i+1]  

	_ypp[0] = -_y[4]+3.0*_y[3]-3.0*_y[2]+_y[1]
	_ypp[_n-1] = _y[_n-2]-3.0*_y[_n-3]+3.0*_y[_n-4]-_y[_n-5]
	for _i in range(_n-2):
		_ypp[_i+1] = 6.0*(_y[_i+2]-2.0*_y[_i+1]+_y[_i])
	_ypp[1] -= _ypp[0]
	_ypp[_n-2] -= _ypp[_n-1]
	for _i in range(_n-3):  _ypp[_i+2] -= _ypp[_i+1]/B[_i+1]
	_i=_n-2; _ypp[_i] /=B[_i]
	while _i >1:	
		_ypp[_i-1] = (_ypp[_i-1]-_ypp[_i])/B[_i-1]; _i-=1
	#_ypp[0] = 3.0*_ypp[1]-3.0*_ypp[2]+_ypp[3]
	#_ypp[_n-1] = 3.0*_ypp[_n-2]-3.0*_ypp[_n-3]+_ypp[_n-4]
	return


#--------- staticE(): Get -dphi/dx=Ex, by differentiating CS-phi  ---------
# Based on cubic spline (see NRC or Hur's note of CP), get Ex on each particle
# from differentiation of phi, i.e. Ex=-d phi/dx. It just differentiate
# given array of phi (_y here) so scale of Ex depends on that of phi, which 
# is defined by user. 
# _np: num. particles, _x: array of particle position, _n: num. mesh = _ncll+1
# _y: potential (phi), _ypp: array to store 2nd derivative of _y. It is
# calculated inside the routine by calling getPP(). _ncll: num. cells
# _ex: array of E-field exerting on each particle
def staticE(_np,_x,_n,_y,_ypp,_ncll,_ex):
	getPP(_n,_y,_ypp)
	for _i in range(_np):	
		_j=(int)(_x[_i])
		_dw = _x[_i] - _j
		if _j<_ncll and _j>=0: 
			_ex[_i] = _y[_j+1]+(3.0*_dw**2-1.0)*_ypp[_j+1]/6.0-_y[_j]-(3.0*(1.0-_dw)**2-1.0)*_ypp[_j]/6.0
			_ex[_i] *= -1.0
		else: _ex[_i]=0.0
	return



#----- collide(): rearrange the velocities by elastic collision -----
def collide(_nu,_dt,_N,_u,_ux,_uy,_uz):
	_prb = _nu*_dt
	for _i in range(_N):
		if random.random() <_prb: 
			_q1 = np.pi*random.random()
			_q2 = 2.0*np.pi*random.random()
			_ux[_i] = _u[_i]*np.cos(_q1) 
			_ut= _u[_i]*np.sin(_q1)
			_uy[_i] = _ut*np.cos(_q2)
			_uz[_i] = _ut*np.sin(_q2)
	return


#------------------ 1D wave solve on Yee mesh -----------------------
# _w (angular freq of EM wave),_dx, _dt are in SI unit, i.e. [rad/s],
# [m], [s]. All other variables including dt should be normalized 
# by EM normalization convention (see ppt)
# Note that Nbuf is assumed to be >=1. Otherwise the solver does not work.
def maxwell(_w,_dx,_dt,_ex,_ey,_ez,_bx,_by,_bz,_jx,_jy,_jz):
	_al = clgt*_dt/_dx;   _dt = _w*_dt
	_a = Nbuf;  _b = -Nbuf 
	_by[_a:_b-1] +=  _al*(_ez[_a+1:_b]-_ez[_a:_b-1])
	_bz[_a:_b-1] += -_al*(_ey[_a+1:_b]-_ey[_a:_b-1]) 
	_ex[_a:_b-1] += -_jx[_a:_b-1]*_dt
	_ey[_a+1:_b-1] +=-_al*(_bz[_a+1:_b-1]-_bz[_a:_b-2]) -_jy[_a+1:_b-1]*_dt 
	_ez[_a+1:_b-1] += _al*(_by[_a+1:_b-1]-_by[_a:_b-2]) -_jz[_a+1:_b-1]*_dt

	#_by[:-1] +=  _al*(_ez[1:]-_ez[:-1])
	#_bz[:-1] += -_al*(_ey[1:]-_ey[:-1]) 
	#_ey[1:-1] +=-_al*(_bz[1:-1]-_bz[:-2]) 
	#_ez[1:-1] += _al*(_by[1:-1]-_by[:-2])
	return


#------------- calculation of current density (3rd order) --------------
def J3(_np,_w,_dx,_dt,_q,_xx,_xxi,_xxo,_ux,_uy,_uz,_gam,_jx,_jy,_jz):
	_jx *=0.0; _jy *=0.0; _jz *=0.0  # initialize J arrays by zero
	_nc = _w**2*eps0*unitM/unitQ**2
	_nrmF = 1.0/(unitQ*_nc*clgt*_dt)
	_nrmFp = 1.0/(unitQ*_nc*_dx)

	for _i in range(_np): 
		_xn = _xx[_i];  _xo = _xxo[_i];   _gm = _gam[_i]  
		_pL2=(_xn-_xo)**2;  _ns=1;  _im = _imm=0;  _j = _xxi[_i]+Nbuf

		while _ns <=2:
			_ip=0.0;  _imm=0
			if _xo<=1 and _xo >=0: _ns=100
			elif _xo<0: _imm=-1; _ip=0.999999; _xo=0.0
			elif _xo>1: _imm=1; _ip=-0.999999; _xo=1.0
			else: print('wrong xo'); quit()
			_wb=0.5*(_xn+_xo); _dw=_xn-_xo; _wbc=1.0-_wb
			_wb1=_wb+1.0; _wbc1=_wbc+1.0;
			_WF = _nrmF*_q*_dw
			_jx[_j+_im-1] += 0.5*_wbc**2*_WF 
			_jx[_j+_im] += (0.5+_wb*_wbc)*_WF 
			_jx[_j+_im+1] += 0.5*_wb**2*_WF

			_fr=_dw*_dw
			if _pL2 > _fr+1e-12: _fr=np.sqrt(_fr/_pL2)
			else: _fr=1.0
			_fc1 =F_T-_wb1*(2.0-_wb1*(1.0-_wb1/6.0))
			_fc  =T_T-_wb**2*(1.0-0.5*_wb)
			_fcc =T_T-_wbc**2*(1.0-0.5*_wbc)
			_fcc1=F_T-_wbc1*(2.0-_wbc1*(1.0-_wbc1/6.0))

			_WF=_fr*_nrmFp*_q*(_uy[_i]/_gm)
			_jy[_j+_im-1] += _fc1*_WF 
			_jy[_j+_im]   += _fc*_WF
			_jy[_j+_im+1] += _fcc*_WF
			_jy[_j+_im+2] += _fcc1*_WF

			_WF=_fr*_nrmFp*_q*(_uz[_i]/_gm)
			_jz[_j+_im-1] += _fc1*_WF 
			_jz[_j+_im]   += _fc*_WF
			_jz[_j+_im+1] += _fcc*_WF
			_jz[_j+_im+2] += _fcc1*_WF

			_im +=_imm; _xn = _xo+_ip; _xo = _xxo[_i] - _im; 
			_ns+=1
	return


#-------------- calculation of J (1st order) --------------		
# Input
# _np: num. particles
# _w,_dx,_dt: angular freq. EM wave [rad/s], mesh size [m], time step [s]
# _q: charge of a superParticle
# _xx,_xxi,_xxo: same as in move()
# _ux,_uy,_uz,_gam: velocity [c] and relativistic gamma
# Output
# _jx,_jy,_jz: current density 
def J1(_np,_w,_dx,_dt,_q,_xx,_xxi,_xxo,_ux,_uy,_uz,_gam,_jx,_jy,_jz):
	_jx *=0.0; _jy *=0.0; _jz *=0.0  # initialize J arrays by zero
	_nc = _w**2*eps0*unitM/unitQ**2
	_nrmF = 1.0/(unitQ*_nc*clgt*_dt)
	_nrmFp = 1.0/(unitQ*_nc*_dx)

	for _i in range(_np): 
		_xn = _xx[_i];  _xo = _xxo[_i];   _gm = _gam[_i]  
		_pL2=(_xn-_xo)**2;  _ns=1;  _im = _imm=0;  _j = _xxi[_i]+Nbuf

		while _ns <=2:
			_ip=0.0;  _imm=0
			if _xo<=1 and _xo >=0: _ns=100
			elif _xo<0: _imm=-1; _ip=0.999999; _xo=0.0
			elif _xo>1: _imm=1; _ip=-0.999999; _xo=1.0
			else: print('wrong xo'); quit()

			_wb=0.5*(_xn+_xo); _dw=_xn-_xo; _WF = _nrmF*_q*_dw
			_jx[_j+_im] += _WF 

			_fr=_dw*_dw
			if _pL2 > _fr+1e-12: _fr=np.sqrt(_fr/_pL2)
			else: _fr=1.0

			_WF=_fr*_nrmFp*_q*(_uy[_i]/_gm)
			_jy[_j+_im]   += (1.0-_wb)*_WF
			_jy[_j+_im+1] += _wb*_WF

			_WF=_fr*_nrmFp*_q*(_uz[_i]/_gm)
			_jz[_j+_im]   += (1.0-_wb)*_WF
			_jz[_j+_im+1] += _wb*_WF

			_im +=_imm; _xn = _xo+_ip; _xo = _xxo[_i] - _im; 
			_ns+=1
	return


#-------------------- J 1D linear interpolation ----------------------
# _np: number of simulation particles to trace
# _n0: n0/nc, where n0 is the peak plasma density
# _q: q/e
# _Nc: number of simulation particles per cell at the peak density
# _dt,_dx: wdt, kdx, i.e. normalized time step and mesh size
# _xx,_xxi,_xxo: dx-normalized position variables
# _ux,_uy,_uz,_gam: c-normalized relativistic velocities and gamma factor
# _jx,_jy,_jz: ecnc-normalized current density. Returning values.
def J1lin(_np,_n0,_q,_Nc,_dt,_dx,_xx,_xxi,_xxo,_ux,_uy,_uz,_gam,_jx,_jy,_jz):
	_jx *=0.0; _jy *=0.0; _jz *=0.0  # initialize J arrays by zero
	_gpFc = _q*_n0/_Nc
	_gFc  = _gpFc*(_dx/_dt)

	for _i in range(_np): 
		_xn = _xx[_i];  _xo = _xxo[_i];   _gm = _gam[_i]  
		_pL2=(_xn-_xo)**2;  _ns=1;  _im = _imm=0;  _j = _xxi[_i]+Nbuf

		while _ns <=2:
			_ip=0.0;  _imm=0
			if _xo<=1 and _xo >=0: _ns=100
			elif _xo<0: _imm=-1; _ip=0.999999; _xo=0.0
			elif _xo>1: _imm=1; _ip=-0.999999; _xo=1.0
			else: print('wrong xo'); quit()

			_wb=0.5*(_xn+_xo); _dw=_xn-_xo; _WF = _gFc*_dw
			_jx[_j+_im] += _WF 

			_fr=_dw*_dw
			if _pL2 > _fr+1e-12: _fr=np.sqrt(_fr/_pL2)
			else: _fr=1.0

			_WF=_fr*_gpFc*(_uy[_i]/_gm)
			_jy[_j+_im]   += (1.0-_wb)*_WF
			_jy[_j+_im+1] += _wb*_WF

			_WF=_fr*_gpFc*(_uz[_i]/_gm)
			_jz[_j+_im]   += (1.0-_wb)*_WF
			_jz[_j+_im+1] += _wb*_WF

			_im +=_imm; _xn = _xo+_ip; _xo = _xxo[_i] - _im; 
			_ns+=1
	return


#-------------------- move 1D for EM simulation ----------------------
# Boris mover, with tracking old position and mesh indices of particles
# _np: number of simulation particles to push
# _q, _m: q/e and m/me, i.e. normalized charge and mass of a single particle
# _dt,_dx: wdt and kdx, i.e. normalized time step and mesh size
# _x,_xo,_xi: dx-normalized particle position arrays, Updated values 
# _gmm: gamma factor. Return value
# _uu,_uxx,_uyy,_uzz: c-normalized relativistic velocities. Updated values
# _ex,...,_bzz: normalized fields on each particle 
def moveEM1(_np,_q,_m,_dt,_dx,_x,_xo,_xi,_gmm,_uu,_uxx,_uyy,_uzz,_ex,_ey,_ez,_bxx,_byy,_bzz):
	qmhdt = 0.5*_dt*_q/_m
	_al = _dt/_dx
	_id=0;  
	_xo[:_np] = _x[:_np]  #save current position to _xo, i.e. old position 
	while _id < _np:
		_ax = qmhdt*_ex[_id]; _ay =qmhdt*_ey[_id]; _az =qmhdt*_ez[_id]
		_ux = _uxx[_id]; _uy = _uyy[_id]; _uz = _uzz[_id]
		_bx = _bxx[_id]; _by = _byy[_id]; _bz = _bzz[_id];

		_ux += _ax;    _uy += _ay;    _uz += _az
		_gm =np.sqrt(1.0+_ux**2 +_uy**2 +_uz**2)

		_qdt_gm = qmhdt/_gm 
		_tx = _qdt_gm*_bx; _ty = _qdt_gm*_by; _tz = _qdt_gm*_bz
		_t2p1 = 2.0/(1.0 + _tx**2 +_ty**2 +_tz**2)

		_upx = _ux + _uy*_tz-_uz*_ty
		_upy = _uy + _uz*_tx-_ux*_tz
		_upz = _uz + _ux*_ty-_uy*_tx

		_ux += _t2p1*(_upy*_tz-_upz*_ty)
		_uy += _t2p1*(_upz*_tx-_upx*_tz)
		_uz += _t2p1*(_upx*_ty-_upy*_tx)

		_ux += _ax; _uy += _ay; _uz += _az
		_usq = _ux**2+ _uy**2 + _uz**2
		_gm = np.sqrt(1.0+ _usq)

		_gmm[_id] = _gm; _uxx[_id]=_ux; _uyy[_id]=_uy; _uzz[_id]=_uz
		_uu[_id] = np.sqrt(_usq)
		_x[_id] += (_ux/_gm)*_al 

		if _x[_id]>1:  _x[_id] -=1;  _xo[_id] -=1;  _xi[_id] +=1
		elif _x[_id]<0:  _x[_id] +=1;  _xo[_id] +=1;  _xi[_id] -=1

		_id+=1
	return


#---------------- handle particles in OB for EM case -----------------
# _np: number of simulation particles
# _ncll: number of cells
# _x,_xo,_xi: dx-normalized positions
# _gm,_u,_ux,_uy,_uz: rel-factor, c-normalized rel speed and velocities
# _lbc, _rbc: left and right boundary options
#    0: simple remove OB particles
#    1: reflected
#    2: packman-like
# Returns updated value of _np
def treatOBpEM1(_np,_ncll,_x,_xo,_xi,_gm,_u,_ux,_uy,_uz,_lbc,_rbc):
	_id=0;  

	while _id < _np:
		if _xi[_id]<0:
			if _lbc==0:
				_li=_np-1;  _np -=1
				if _id < _li:
					_x[_id]  =_x[_li]; _xi[_id] =_xi[_li]
					_xo[_id] =_xo[_li]; _u[_id]  =_u[_li]
					_ux[_id] =_ux[_li]; _uy[_id] =_uy[_li]
					_uz[_id] =_uz[_li]; _gm[_id] =_gm[_li]
				else:  return _np 
			elif _lbc==1:
				_xi[_id] = 0;  _x[_id] = 1.0 - _x[_id]
				_xo[_id] = 1.0 - _xo[_id]; _ux[_id] *= -1.0
				_id +=1
			elif _lbc==2:  _xi[_id] += _ncll;  _id +=1
			else:  print("Wrong _lbc in treatOBpEM1")

		elif _xi[_id]>=_ncll:
			if _rbc==0:
				_li=_np-1;  _np -=1
				if _id < _li:
					_x[_id]  =_x[_li]; _xi[_id] =_xi[_li]
					_xo[_id] =_xo[_li]; _u[_id]  =_u[_li]
					_ux[_id] =_ux[_li]; _uy[_id] =_uy[_li]
					_uz[_id] =_uz[_li]; _gm[_id] =_gm[_li]
				else:  return _np  
			elif _rbc==1:
				_xi[_id] = _ncll-1;  _x[_id] = 1.0 - _x[_id]
				_xo[_id] = 1.0 - _xo[_id]; _ux[_id] *= -1.0
				_id +=1
			elif _rbc==2:  _xi[_id] -= _ncll;  _id +=1
			else:  print("Wrong _rbc in treatOBpEM1")

		else: _id+=1
	return _np


#------------- Shift simulation window in EM simulations -------------
# Written by MSHur, 20210314
# Shift field and particle position backward by one mesh. That yields
# the same effect as moving simulation window one step forward.
# Note that for particles, xi has only to be shifted.
#
# _np: num. particles
# _xi: array of particle's mesh index 
# _ex,..,_bz: electromagnetic field on meshes
def shiftWindow(_np,_xi,_ex,_ey,_ez,_bx,_by,_bz):
	_a = Nbuf;  _b = -Nbuf 

	#_ex[_a:_b-1] = _ex[_a+1:_b]
	#_ey[_a:_b-1] = _ey[_a+1:_b]
	#_ez[_a:_b-1] = _ez[_a+1:_b]
	#_bx[_a:_b-1] = _bx[_a+1:_b]
	#_by[_a:_b-1] = _by[_a+1:_b]
	#_bz[_a:_b-1] = _bz[_a+1:_b]

	#_ex[_b-1]=_ey[_b-1]=_ez[_b-1]=0.0
	#_bx[_b-1]=_by[_b-1]=_bz[_b-1]=0.0

	_ex[:-1] = _ex[1:]
	_ey[:-1] = _ey[1:]
	_ez[:-1] = _ez[1:]
	_bx[:-1] = _bx[1:]
	_by[:-1] = _by[1:]
	_bz[:-1] = _bz[1:]


	for _id in range(_np): _xi[_id] -=1

	return


#------------------ Load particles in only one-cell ------------------ 
# Written by MSHur, 20210314
# Mostly used combined with shiftWindow, to load new particles from 
# the window edge. 
# _m, _n0: m/me and n0/nc, respectively 
# _prof: density profile function. Single argument is k-nrm position
# _T: temperature in [eV]
# _npc: num. ptcls per cell.
# _dx: k*dx
# _icll: index of the cell to load ptcls into.
# _orgn: displacement of the origin (x=0) in unit of cell. It is required
#        when shiftwindow is running, to find actual position for _prof.
# _np: current num. ptcls in ptcl arrays (_x,_xo,...,_gm).
# _x,...,_gm: arrays of ptcls, to which the new ptcls are added. Assumes 
#             they have enough extra spaces to add new ones. 
def loadSlice(_m,_n0,_prof,_T,_npc,_dx,_icll,_orgn,_np,_x,_xo,_xi,_u,_ux,_uy,_uz,_gm):
	_mm = _m*unitM;  _id=0
	for _i in range(_npc):
		_xx = random.random(); 
		_den=_prof((_xx+_icll+_orgn)*_dx)

		if _den > random.random():
			_ii = _id + _np 
			_x[_ii] = _xx;  _xo[_ii] = _xx; _xi[_ii] = _icll; 
			_uu=thermal(_T,_mm)
			_ux[_ii] =_uu[0]; _uy[_ii] =_uu[1]; _uz[_ii] =_uu[2]
			_usq = _uu[0]**2+_uu[1]**2+_uu[2]**2
			_u[_ii] = np.sqrt(_usq); _gm[_ii] = np.sqrt(1.0+_usq)
			_id +=1

	return _np+_id


#--------------- Field line constructor from r-z fields --------------
def fLineRZ(_r0,_z0,_ds, _smax,_frz):
	def _f(_yy,_tt):
		_rr,_zz=_yy
		_rrp,_zzp = _frz(_rr,_zz)
		return [_rrp,_zzp]

	_y0=[_r0, _z0]
	_s = np.arange(0, _smax*_ds, _ds)
	sol = integ.odeint(_f,_y0,_s)	

	return sol




#------------------- Faraday and Ampere-Maxwell laws -----------------
def faraday1D(_dt,_dx,_ey,_ez,_by,_bz,_bya,_bza):
	_al = _dt/_dx;   
	_a = Nbuf;  _b = -Nbuf 

	_bya[_a:_b-1] = _by[_a:_b-1]
	_bza[_a:_b-1] = _bz[_a:_b-1]
	_by[_a:_b-1] +=  _al*(_ez[_a+1:_b]-_ez[_a:_b-1])
	_bz[_a:_b-1] += -_al*(_ey[_a+1:_b]-_ey[_a:_b-1]) 
	_bya[_a:_b-1] = 0.5*(_bya[_a:_b-1] + _by[_a:_b-1])
	_bza[_a:_b-1] = 0.5*(_bza[_a:_b-1] + _bz[_a:_b-1])

	return

def ampereMaxwell1D(_dt,_dx,_ex,_ey,_ez,_bx,_by,_bz,_jx,_jy,_jz):
	_al = _dt/_dx;  
	_a = Nbuf;  _b = -Nbuf 
	_ex[_a:_b-1] += -_jx[_a:_b-1]*_dt
	_ey[_a+1:_b-1] +=-_al*(_bz[_a+1:_b-1]-_bz[_a:_b-2]) -_jy[_a+1:_b-1]*_dt
	_ez[_a+1:_b-1] += _al*(_by[_a+1:_b-1]-_by[_a:_b-2]) -_jz[_a+1:_b-1]*_dt

	return


#------------ Calculate density from X, Xi data in EM sim. -----------
def densityEM(_np,_ncll,_npc,_x,_xi,_n):
	_n *= 0.0 
	for _id in range(_np):
		_j = _xi[_id]
		_n[_j]   += 1.0 - _x[_id]  
		_n[_j+1] += _x[_id]  
	_n /= _npc

	return


#-------------- Artificial beam to generate wakefield ----------------
def artBeam(_t,_E,_ncll,_dx,_E0,_v,_sgm,_c0):
	_c = _v*_t+_c0
	for _i in range(_ncll+1):
		_x = _dx*_i
		_xx = _x-_c
		_E[_i+Nbuf] = _E0*(-2.0*_xx/_sgm)*np.exp(-(_xx/_sgm)**2) 
		#_E[_i+Nbuf] = _E0*np.exp(-(_xx/_sgm)**2) 
	return
#---------------------------------------------------------------------
