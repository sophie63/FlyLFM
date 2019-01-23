/*
	Numerical Integration over 2D sphere (disk) by product of 
	Gauss-type Quadratures. High-precision abscissas and weights are used.

	Project homepage: http://www.holoborodko.com/pavel/?page_id=1879
	Contact e-mail:   pavel@holoborodko.com

	Copyright (c)2007-2010 Pavel Holoborodko
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:
	
	1. Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.
	
	2. Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.
	
	3. Redistributions of any form whatsoever must retain the following
	acknowledgment:
	"
         This product includes software developed by Pavel Holoborodko
         Web: http://www.holoborodko.com/pavel/
         e-mail: pavel@holoborodko.com
	
	"

	4. This software cannot be, by any means, used for any commercial 
	purpose without the prior permission of the copyright holder.
	
	Any of the above conditions can be waived if you get permission from 
	the copyright holder. 

	THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
	FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
	OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
	HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
	LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
	OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
	SUCH DAMAGE.
*/

#ifndef __GAUSS_2D_SPHERE_H__
#define __GAUSS_2D_SPHERE_H__

#ifdef __cplusplus
extern "C"
{
#endif

	/* 
	 2D Numerical computation of int(f(x,y), (x-Xc)^2+(y-Yc)^2<=R^2) by Gauss-type high precision quadratures. 
	 2D integration over the sphere (disk) derived as a product of two Gaussian rules for: 
	 1. int(|r|*f(r), r=-1..1)
	 2. int((1-t^2)^(-1/2)*f(t),t=-1..1) - Gauss-Chebyshev quadrature

	 See for details: 
	 A.H. Stroud Approximate Calculation of Multiple Integrals. pages 32-35.

	 Cubature of n-th order has n^2 nodes and exact on polynomials of total degree at most 2*n-1.
	 All nodes and weights are precise up to 25 digits (some even more)

     Parameters:
		[in]n       - quadrature order, n = 1..50, 64, 128, 192, 256, 320, 384, 448, 512, 576, 640, 704, 768, 832, 896, 960, 1024
		[in]f       - integrand
		[in]data    - pointer on user-defined data which will 
					  be passed to f every time it called (as third parameter).
		[in]R       - radius of the sphere
		[in](Xc,Yc) - center of the sphere. Unit disk is R=1, Xc=Yc=0. 

	 return:
			-computed integral value or -1.0 if n order quadrature is not supported
	*/
	double gauss_product_2D_sphere(int n, double (*f)(double,double,void*), void* data, double R, double Xc, double Yc);

#ifdef __cplusplus
}
#endif

#endif /* __GAUSS_2D_SPHERE_H__ */


