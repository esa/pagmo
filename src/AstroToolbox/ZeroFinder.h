/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef ZEROFINDER_H
#define ZEROFINDER_H

#define NUMERATOR(dab,dcb,fa,fb,fc) fb*(dab*fc*(fc-fb)-fa*dcb*(fa-fb))
#define DENOMINATOR(fa,fb,fc) (fc-fb)*(fa-fb)*(fa-fc)


namespace ZeroFinder
{
	
	//// Class for one dimensional functions
	//  with some parameters
	/*  The ()-operator with one double argument
	 *  has to be overloaded for a derived class
	 *  The return value is the ordinate computed for
	 *  the abscissa-argument.
	 */
	class Function1D
	{
	public:
		Function1D(const double &a, const double &b):p1(a),p2(b) {}
		// parameters
		const double p1,p2;
	};
	
	class Function1D_7param
	{
	public:
		Function1D_7param(const double &a, const double &b, const double &c,
				  const double &d, const double &e, const double &f,
				  const double &g):p1(a),p2(b),p3(c),p4(d),p5(e),p6(f),p7(g) {}
		// parameters
		const double p1,p2,p3,p4,p5,p6,p7;
	};
	
	
	class FZero
	{
	private:
		double a, c; // lower and upper bounds
		
	public:
		FZero(const double &a_, const double &c_):a(a_),c(c_) {} // constructor
		// fzero procedure
		template <class Functor>
				double FindZero(const Functor &f);
	};
	
	//-------------------------------------------------------------------------------------//
	// This part is an adaption of the 'Amsterdam method', which is an inversre quadratic  //
	// interpolation - bisection method													   //
	// See http://mymathlib.webtrellis.net/roots/amsterdam.html							   //
	//
	//-------------------------------------------------------------------------------------//
	template <class Functor>
			double FZero::FindZero(const Functor &f)
	{
		static const int max_iterations = 500;
		static const double tolerance = 1e-15;
		
		double fa = f(a), b = 0.5 * ( a + c ), fc = f(c), fb = fa * fc;
		double delta, dab, dcb;
		int i;
		
		// If the initial estimates do not bracket a root, set the err flag and //
		// return.  If an initial estimate is a root, then return the root.     //
		
		if ( fb >= 0.0 ) {
			if ( fb > 0.0 )  { return 0.0; }
			else {return ( fa == 0.0 ) ? a : c;}
		}
		
		// Insure that the initial estimate a < c. //
		
		if ( a > c ) {
			delta = a; a = c; c = delta; delta = fa; fa = fc; fc = delta;
		}
		
		// If the function at the left endpoint is positive, and the function //
		// at the right endpoint is negative.  Iterate reducing the length    //
		// of the interval by either bisection or quadratic inverse           //
		// interpolation.  Note that the function at the left endpoint        //
		// remains nonnegative and the function at the right endpoint remains //
		// nonpositive.                                                       //
		
		if ( fa > 0.0 )
			for ( i = 0; i < max_iterations; i++) {
			
			// Are the two endpoints within the user specified tolerance ?
			
			if ( ( c - a ) < tolerance ) return 0.5 * ( a + c);
			
			// No, Continue iteration.
			
			fb = f(b);
			
			// Check that we are converging or that we have converged near //
			// the left endpoint of the inverval.  If it appears that the  //
			// interval is not decreasing fast enough, use bisection.      //
			if ( ( c - a ) < tolerance ) return 0.5 * ( a + c);
			if ( ( b - a ) < tolerance ) {
				if ( fb > 0 ) {
					a = b; fa = fb; b = 0.5 * ( a + c ); continue;
				}
				else {
					return b;
				}
			}
			
			// Check that we are converging or that we have converged near //
			// the right endpoint of the inverval.  If it appears that the //
			// interval is not decreasing fast enough, use bisection.      //
			
			if ( ( c - b ) < tolerance ) {
				if ( fb < 0 ) {
					c = b; fc = fb; b = 0.5 * ( a + c ); continue;
				}
				else {
					return b;
				}
			}
			
			// If quadratic inverse interpolation is feasible, try it. //
			
			if (  ( fa > fb ) && ( fb > fc ) ) {
				delta = DENOMINATOR(fa,fb,fc);
				if ( delta != 0.0 ) {
					dab = a - b;
					dcb = c - b;
					delta = NUMERATOR(dab,dcb,fa,fb,fc) / delta;
					
					// Will the new estimate of the root be within the   //
					// interval?  If yes, use it and update interval.    //
					// If no, use the bisection method.                  //
					
					if ( delta > dab && delta < dcb ) {
						if ( fb > 0.0 ) { a = b; fa = fb; }
						else if ( fb < 0.0 )  { c = b; fc = fb; }
						else return b;
						b += delta;
						continue;
					}
				}
			}
			
			// If not, use the bisection method. //
			
			fb > 0 ? ( a = b, fa = fb ) : ( c = b, fc = fb );
			b = 0.5 * ( a + c );
		}
		else
			
			// If the function at the left endpoint is negative, and the function //
			// at the right endpoint is positive.  Iterate reducing the length    //
			// of the interval by either bisection or quadratic inverse           //
			// interpolation.  Note that the function at the left endpoint        //
			// remains nonpositive and the function at the right endpoint remains //
			// nonnegative.                                                       //
			
			for ( i = 0; i < max_iterations; i++) {
			if ( ( c - a ) < tolerance ) return 0.5 * ( a + c);
			fb = f(b);
			
			if ( ( b - a ) < tolerance ) {
				if ( fb < 0 ) {
					a = b; fa = fb; b = 0.5 * ( a + c ); continue;
				} else {
					return b;
				}
			}
			
			if ( ( c - b ) < tolerance ) {
				if ( fb > 0 ) {
					c = b; fc = fb; b = 0.5 * ( a + c ); continue;
				} else {
					return b;
				}
			}
			
			if (  ( fa < fb ) && ( fb < fc ) ) {
				delta = DENOMINATOR(fa,fb,fc);
				if ( delta != 0.0 ) {
					dab = a - b;
					dcb = c - b;
					delta = NUMERATOR(dab,dcb,fa,fb,fc) / delta;
					if ( delta > dab && delta < dcb ) {
						if ( fb < 0.0 ) { a = b; fa = fb; }
						else if ( fb > 0.0 )  { c = b; fc = fb; }
						else return b;
						b += delta;
						continue;
					}
				}
			}
			fb < 0 ? ( a = b, fa = fb ) : ( c = b, fc = fb );
			b = 0.5 * ( a + c );
		}
		return  b;
	}
}

#undef NUMERATOR
#undef DENOMINATOR

#endif
