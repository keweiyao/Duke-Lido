#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>

// for vecotr data
struct VectorT{
    double t, x, y, z;
    friend inline std::ostream& operator<<(std::ostream& os, const VectorT& a) {
      os << a.t << ' ' << a.x << ' ' << a.y << ' ' << a.z;
      return os;
    }
    friend inline VectorT operator+(VectorT a, const VectorT& b) {
        a.t += b.t;
        a.x += b.x;
		a.y += b.y;
		a.z += b.z;
		return a;
    }
	VectorT operator* (double u) const{
		return VectorT{t*u, x*u, y*u, z*u};
	}
};

// for tensor data
struct TensorT{
    double tt, tx, ty, tz;
    double xt, xx, xy, xz;
    double yt, yx, yy, yz;
    double zt, zx, zy, zz;
    friend inline std::ostream& operator<<(std::ostream& os, const TensorT& a) {
      os << a.tt << ' ' << a.tx << ' ' << a.ty << ' ' << a.tz << std::endl 
         << a.xt << ' ' << a.xx << ' ' << a.xy << ' ' << a.xz << std::endl 
		 << a.yt << ' ' << a.yx << ' ' << a.yy << ' ' << a.yz << std::endl 
		 << a.zt << ' ' << a.zx << ' ' << a.zy << ' ' << a.zz;
      return os;
    }
    friend inline TensorT operator+(TensorT a, const TensorT& b) {
        a.tt += b.tt; a.xt += b.xt; a.yt += b.yt; a.zt += b.zt;
        a.tx += b.tx; a.xx += b.xx; a.yx += b.yx; a.zx += b.zx;
		a.ty += b.ty; a.xy += b.xy; a.yy += b.yy; a.zy += b.zy;
		a.tz += b.tz; a.xz += b.xz; a.yz += b.yz; a.zz += b.zz;
		return a;
    }
	TensorT operator* (double u) const{
		return TensorT{ tt*u, tx*u, ty*u, tz*u,
						xt*u, xx*u, xy*u, xz*u,
						yt*u, yx*u, yy*u, yz*u,
						zt*u, zx*u, zy*u, zz*u};
	}
};
#endif
