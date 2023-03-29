import numpy as np
#import matplotlib.pyplot as plt
_pi = np.pi;  _clgt=2.99792458e8

def emitEMwave(t,E,dX,mthd,pol,phi,theta,frq,psi0,chrp,fl,r0,dur,Apk): 
	_dm = E.ndim
	if _dm > 1: (_szh, _szv) = E.shape
	else: (_szh,_szv) = (len(E),1)
	_dh = _dv = dX[0]
	if len(dX)>1:  _dv = dX[1]
	_r0h = _r0v = r0[0]
	if len(r0)>1:  _r0v = r0[1]
	_ch = _dh* (_szh-1)/2;  _cv = _dv* (_szv-1)/2
	_sf = np.sin(phi);  _sq = np.sin(theta)
	_cf = np.cos(phi);  _cq = np.cos(theta)
#	_lmda = _clgt/frq;   _w = 2.0*_pi*frq;   _k = _w/_clgt
	_lmda = 1.0/frq;   _w = 2.0*_pi*frq;   _k = _w
	_zh = _pi*_r0h**2/_lmda;  _zv = _pi*_r0v**2/_lmda
	_A0 = Apk*np.exp(-((t-2.0*dur)/dur)**2)

	for _i in range(_szh):  
		for _j in range(_szv):
			_h = _i*_dh-_ch;  _v = _j*_dv-_cv
			_np = _h*_sf*_cq + _v*_sq
			_hp = _h*_cf;   _vp = -_h*_sf*_sq + _v*_cq 
			_zf = _np-fl
			_ah = -_zf/_zh;  _av = -_zf/_zv
			_bh = 0.5*np.arctan(_ah);  _bv = 0.5*np.arctan(_av)
			_rh = _r0h*np.sqrt(1+_ah**2);  _rv = _r0v*np.sqrt(1+_av**2)
			_Dh = (_hp/_rh)**2;  _Dv = (_vp/_rv)**2
			_psi = -_k*_np + _w*t + _Dh*_ah + _Dv*_av - _bh -_bv +psi0
			_S = np.sqrt((_r0h/_rh)*(_r0v/_rv))*np.exp(-_Dh-_Dv)
			_Ep = _A0*_S*np.cos(_psi)
			if pol == 'p':
				_Enp = (_lmda/_pi)*_Ep*(_hp/_rh**2)*(_ah*np.cos(_psi)+np.sin(_psi))
			else: 
				_Enp = (_lmda/_pi)*_Ep*(_vp/_rv**2)*(_av*np.cos(_psi)+np.sin(_psi))
			_fld = _Enp*_cq*_sf+_Ep*_cf
			if mthd=='sum':  E[_i] += _fld
			else:  E[_i] = _fld 
	return
