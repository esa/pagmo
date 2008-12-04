// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#ifndef ZEROFINDER_H
#define ZEROFINDER_H

/// Namespace
namespace ZeroFinder
{

	/// Class for one dimensional functions
	//  with some parameters
	/**  The ()-operator with one double argument
	 *  has to be overloaded for a derived class
	 *  The return value is the ordinate computed for
	 *  the abscissa-argument.
	 */
	class Function1D
	{
	public:
		//virtual double Compute(double x)=0;
		virtual double operator()(double x)=0;
		// parameters
		double p1,p2;
		void SetParameters(double a, double b);
		virtual ~Function1D() {}
	};

	class Function1D_7param
	{
	public:
		//virtual double Compute(double x)=0;
		virtual double operator()(double x)=0;
		// parameters
		double p1,p2,p3,p4,p5,p6,p7;
		void SetParameters(double a, double b, double c, double d, double e, double f, double g);
		virtual ~Function1D_7param() {}
	};


	class FZero
	{
		private:
			double a, c; // lower and upper bound

		public:
			FZero(double a, double b); // constructor
			// fzero procedure
			double FindZero(Function1D& f);
			double FindZero7(Function1D_7param& f);
			void SetBounds(double a, double b);
	};


}
#endif
