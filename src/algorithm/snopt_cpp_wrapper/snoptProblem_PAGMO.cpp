#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "../../problem/base.h"
#include "../../algorithm/snopt.h"
#include "snopt_PAGMO.h"
#include "snfilewrapper_PAGMO.h"
#include "snoptProblem_PAGMO.h"

using namespace std;

snoptProblem_PAGMO::snoptProblem_PAGMO(const pagmo::problem::base& problema, pagmo::algorithm::snopt::preallocated_memory *di_comodo):
  iSpecs(0), iSumm(0), iPrint(0), initCalled(0),cu((char *)&problema),lencu(500),ru((doublereal *)di_comodo), lenru(500)
{
  init2zero();

  // Nothing incremented yet
  fortranStyleObj = 0;
  fortranStyleAG  = 0;

  //  iSpecs =  4;
  //  if (screen) iSumm  = 6;
  //  iPrint = 15;

  // Create temporary memory for the call to sninit_.
  // Lengths must all be >= 500.
  lencw  = 500;
  leniw  = 500;
  lenrw  = 500;
  this->alloc( 500, 500, 500 );
  // sninit_ "undefines" the optional parameters

  this->init();
}

snoptProblem_PAGMO::~snoptProblem_PAGMO()
{
  //Close print and spec files if necessary.
  if (iPrint != 0 ) {
    snclose_( &iPrint );
  }
  if (iSpecs != 0 ) {
    snclose_( &iSpecs );
  }

  //Delete work arrays.
  delete [] rw;
  delete [] iw;
  delete [] cw;
}

void snoptProblem_PAGMO::init2zero()
{
  // Data that must be set by user.

  n = 0; neF = 0;

  ObjRow =  0;
  ObjAdd =  0;

  iAfun = 0; jAvar = 0;  A = 0;
  iGfun = 0; jGvar = 0;

  x = 0; xlow = 0; xupp = 0; xmul = 0;
  F = 0; Flow = 0; Fupp = 0; Fmul = 0;

  xstate  = 0; Fstate  = 0;
  nxnames = 0; nFnames = 0;

  usrfun = 0;

  neA  = -1;  // Indicate that neA is yet to be assigned
  neG  = -1;
  lenA = -1;
  lenG = -1;
}

void snoptProblem_PAGMO::userDataSet()
{
  if ( n    == 0)  errMsgExit( "n"  );
  if ( neF  == 0)  errMsgExit( "neF");

  if ( x    == 0 ) errMsgExit( "x"    );
  if ( xlow == 0 ) errMsgExit( "xlow" );
  if ( xupp == 0 ) errMsgExit( "xupp" );
  if ( xmul == 0 ) errMsgExit( "xmul" );

  if ( F    == 0 ) errMsgExit( "F"    );
  if ( Flow == 0 ) errMsgExit( "Flow" );
  if ( Fupp == 0 ) errMsgExit( "Fupp" );
  if ( Fmul == 0 ) errMsgExit( "Fmul" );

  if ( xnames  == 0 ) errMsgExit( "xnames"  );
  if ( Fnames  == 0 ) errMsgExit( "Fnames"  );
  if ( nxnames == 0 ) errMsgExit( "nxnames" );
  if ( nFnames == 0 ) errMsgExit( "nFnames" );

  if ( usrfun ==  0 ) errMsgExit( "usrfun" );
  if ( lenA   == -1 ) errMsgExit( "lenA" );
  if ( lenG   == -1 ) errMsgExit( "lenG" );

  if ( (neA > 0) && (iAfun == 0) ) errMsgExit( "iAfun" );
  if ( (neA > 0) && (jAvar == 0) ) errMsgExit( "jAvar" );
  if ( neA > 0 && A == 0 ) errMsgExit( "A"     );

  if ( neG > 0 && iGfun == 0 ) errMsgExit( "iGfun" );
  if ( neG > 0 && jGvar == 0 ) errMsgExit( "jGvar" );
}

void snoptProblem_PAGMO::errMsgExit( const char *var )
{
  cerr << "****************************************************\n";
  cerr << "Error: " << var << " must be set prior to call to " << endl
       << "snoptProblem_PAGMO::solve() or snoptProblem_PAGMO::computeJac()!\n";
  cerr << "****************************************************\n";
  exit(1);
}

void snoptProblem_PAGMO::setMemory()
{
  int memoryGuess;
  memoryGuess = this->snmema(mincw, miniw, minrw);
  if ( mincw > lencw || miniw > leniw || minrw > lenrw ) {
    // Reallocate memory while retaining the values set in sninit_
    this->realloc( mincw, miniw, minrw );
    // Save the lengths of the new work arrays.
    this->setIntParameter("Total real workspace   ", lenrw );
    this->setIntParameter("Total integer workspace", leniw );

    // Did we have to guess values of neA and neG for snmema()
    if ( memoryGuess == 1 ) {
      this->computeJac();
      memoryGuess = this->snmema(mincw, miniw, minrw);
      assert( memoryGuess == 0 );
      this->realloc( mincw, miniw, minrw );
      this->setIntParameter("Total real workspace   ", lenrw );
      this->setIntParameter("Total integer workspace", leniw );
    }
  }
}

void snoptProblem_PAGMO::alloc( integer alencw, integer aleniw, integer alenrw )
{
  // Reset work array lengths.
  lencw = alencw;
  leniw = aleniw;
  lenrw = alenrw;

  // Allocate new memory for work arrays.
  cw = new char[8*lencw];
  iw = new integer[leniw];
  rw = new doublereal[lenrw];
}

void snoptProblem_PAGMO::realloc(  integer alencw, integer aleniw, integer alenrw )
{
  // Call to this->alloc will overwrite these values => must save.
  integer tlencw = lencw;
  integer tleniw = leniw;
  integer tlenrw = lenrw;

  // Call to this->alloc will create new values for cw, iw, rw => must save.
  char       *tcw = cw;
  integer    *tiw = iw;
  doublereal *trw = rw;

  // Allocate new memory
  this->alloc   ( alencw, aleniw, alenrw );
  // Copy in old values, previously set
  this->memcpyIn( tcw, tiw, trw, tlencw, tleniw, tlenrw );

  // Delete temporary work arrays
  delete [] tcw;
  delete [] tiw;
  delete [] trw;
}

void snoptProblem_PAGMO::memcpyIn( char *tcw, integer *tiw, doublereal *trw,
			     integer tlencw, integer tleniw, integer tlenrw )
{
  integer mlencw = lencw < tlencw ? lencw : tlencw;
  integer mleniw = leniw < tleniw ? leniw : tleniw;
  integer mlenrw = lenrw < tlenrw ? lenrw : tlenrw;

  memcpy( cw, tcw, 8*mlencw*sizeof( char ) );
  memcpy( iw, tiw, mleniw*sizeof(integer));
  memcpy( rw, trw, mlenrw*sizeof(doublereal));
}

void snoptProblem_PAGMO::memcpyOut( char *tcw, integer *tiw, doublereal *trw,
			     integer tlencw, integer tleniw, integer tlenrw )
{
  integer mlencw = lencw < tlencw ? lencw : tlencw;
  integer mleniw = leniw < tleniw ? leniw : tleniw;
  integer mlenrw = lenrw < tlenrw ? lenrw : tlenrw;

  memcpy( tcw, cw, 8*mlencw*sizeof( char ) );
  memcpy( tiw, iw, mleniw*sizeof(integer));
  memcpy( trw, rw, mlenrw*sizeof(doublereal));
}

void snoptProblem_PAGMO::increment()
{
  if( !fortranStyleObj ) {
    //Increment row indicator.
    ObjRow++;
    fortranStyleObj = 1;
  }

  if( !fortranStyleAG ) {
    //Increment A indices.
    for( int k = 0; k < neA; k++ ) {
      iAfun[k]++; jAvar[k]++;
    }
    //Increment G indices.
    for( int k = 0; k < neG; k++ ) {
      iGfun[k]++; jGvar[k]++;
    }
    fortranStyleAG = 1;
  }
}

void snoptProblem_PAGMO::decrement()
{
  if( fortranStyleObj ) {
    //Decrement row indicator.
    ObjRow--;
    fortranStyleObj = 0;
  }

  if (fortranStyleAG) {
    //Decrement A indices.
    for( int k = 0; k < neA; k++ ) {
      iAfun[k]--; jAvar[k]--;
    }
    //Decrement G indices.
    for( int k = 0; k < neG; k++ ) {
      iGfun[k]--; jGvar[k]--;
    }
    fortranStyleAG = 0;
  }
}

void snoptProblem_PAGMO::computeJac()
{
  //Ensures all user data has been initialized.
  userDataSet();
  this->snmema( mincw, miniw, minrw );
  if ( mincw > lencw || miniw > leniw || minrw > lenrw ) {
    // Reallocate memory while retaining the values set in sninit_
    this->realloc( mincw, miniw, minrw );
    // Save the lengths of the new work arrays.
    this->setIntParameter("Total real workspace   ", lenrw );
    this->setIntParameter("Total integer workspace", leniw );
  }
  snjac_( &inform, &neF, &n, usrfun,
	  iAfun, jAvar, &lenA, &neA, A,
	  iGfun, jGvar, &lenG, &neG,
	  x, xlow, xupp, &mincw, &miniw, &minrw,
          cu, &lencu, iw, &leniw, ru, &lenru, //User Workspaces!!! Contain the Pointer Structure that will be reconstructed inside the static snopt_function_
	  cw, &lencw, iw, &leniw, rw, &lenrw,
	  8*500, 8*500 );
  //snjac_ will generate fortran style arrays.
  fortranStyleAG = 1;
}

int snoptProblem_PAGMO::snmema
   ( integer &amincw, integer &aminiw, integer &aminrw)
{
  int memoryGuess = 0;
  integer nxname = 1; integer nfname = 1;
  if ( neA < 0 ) {
    neA = n*neF;
    memoryGuess = 1;
  }
  if ( neG < 0 ) {
    neG = n*neF;
    memoryGuess = 1;
  }
  snmema_( &inform, &neF, &n, &nxname, &nfname, &neA, &neG,
	   &amincw, &aminiw, &aminrw, cw, &lencw, iw, &leniw,
	   rw, &lenrw, 8*500 );
  return memoryGuess;
}

void snoptProblem_PAGMO::init()
{
  initCalled = 1;
  sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
}

void snoptProblem_PAGMO::setParameter( char *stropt )
{
  assert( initCalled == 1 );

  integer iPrt   = 0; // suppresses printing
  integer iSum   = 0;
  integer stropt_len = strlen(stropt);
  snset_( stropt, &iPrt, &iSum, &inform, cw, &lencw, iw, &leniw,
	  rw, &lenrw, stropt_len, 8*500 );
}

void snoptProblem_PAGMO::getParameter( char *stroptin, char *stroptout )
{
  assert( initCalled == 1 );
  integer stroptin_len  = strlen(stroptin);
  integer stroptout_len = strlen(stroptout);
  sngetc_( stroptin, stroptout, &inform, cw, &lencw, iw, &leniw,
	   rw, &lenrw, stroptin_len, stroptout_len, 8*500 );
}

void snoptProblem_PAGMO::setIntParameter( const char *stropt, integer opt )
{
  assert( initCalled == 1 );

  integer iPrt   = 0; // suppresses printing
  integer iSum   = 0;
  integer stropt_len = strlen(stropt);
  snseti_( (char*)(void*)stropt, &opt, &iPrt, &iSum, &inform,
	   cw, &lencw, iw, &leniw, rw, &lenrw, stropt_len, 8*500 );
}

void snoptProblem_PAGMO::getIntParameter( const char *stropt, integer &opt )
{
  assert( initCalled == 1 );
  integer stropt_len = strlen(stropt);
  sngeti_( (char*)(void*)stropt, &opt, &inform, cw, &lencw, iw, &leniw,
	   rw, &lenrw, stropt_len, 8*500 );
}

void snoptProblem_PAGMO::setRealParameter( const char *stropt, doublereal opt )
{
  assert( initCalled == 1 );

  integer iPrt   = 0; // suppresses printing
  integer iSum   = 0;
  integer stropt_len = strlen(stropt);
  snsetr_( (char*)(void*)stropt, &opt, &iPrt, &iSum, &inform,
	   cw, &lencw, iw, &leniw, rw, &lenrw, stropt_len, 8*500 );
}

void snoptProblem_PAGMO::getRealParameter( const char *stropt, doublereal &opt )
{
  assert( initCalled == 1 );
  integer stropt_len = strlen(stropt);
  sngetr_( (char*)(void*)stropt, &opt, &inform, cw, &lencw, iw, &leniw,
	   rw, &lenrw, stropt_len, 8*500 );
}

void snoptProblem_PAGMO::solve( integer starttype )
{
  assert( initCalled == 1 );
  //Ensures all user data initialized.
  userDataSet();
  //Unlike snjac_ we also need neA and neG to be set.
  if ( neA == -1 || neG == -1 ) {
    cerr << "Warning: neA and neG must be set before calling"
         << "snoptProblem_PAGMO::solve()\n";
    exit(1);
  }
  integer npname = strlen(Prob);
  integer nS, nInf;
  doublereal sInf;
  this->increment(); //Convert array entries to Fortran style
  this->setMemory();
  snopta_( &starttype, &neF, &n, &nxnames,
	   &nFnames,
	   &ObjAdd, &ObjRow, Prob,
	   usrfun, iAfun, jAvar, &lenA, &neA, A,
	   iGfun, jGvar, &lenG, &neG,
	   xlow, xupp, xnames, Flow,
	   Fupp, Fnames, x, xstate,
	   xmul, F, Fstate, Fmul,
	   &inform, &mincw, &miniw, &minrw, &nS, &nInf, &sInf,
           cu, &lencu, iw, &leniw, ru, &lenru,      //User Workspaces!!! Contain the Pointer Structure that will be reconstructed inside the static snopt_function_
           cw, &lencw, iw, &leniw, rw, &lenrw,
	   npname, 8*(nxnames), 8*(nFnames), 8*500, 8*500);
  this->decrement();  //Convert array entries to C style
}

void snoptProblem_PAGMO::setPrintFile( const char aprintname[] )
{
    assert( initCalled = 1 );
    if (iPrint != 0 ) {
      snclose_( &iPrint );
    }
    iPrint = 15;
    strcpy( printname, aprintname );  prnt_len = strlen(printname);
    snopenappend_( &iPrint, printname,   &inform, prnt_len );
    this->setIntParameter("Print file", iPrint);
}

void snoptProblem_PAGMO::setSpecFile( char aspecname[] )
{
  assert( initCalled == 1 );
  if (iSpecs != 0 ) {
    snclose_( &iSpecs );
  }
  iSpecs = 4;
  strcpy( specname, aspecname );    spec_len = strlen(specname);
  snfilewrapper_( specname, &iSpecs, &inform, cw, &lencw,
		  iw, &leniw, rw, &lenrw, spec_len, 8*lencw);
  if( inform != 101 ){
    printf("Warning: unable to find specs file %s \n", specname);
  }
}

void snoptProblem_PAGMO::setProblemSize( integer an, integer aneF )
{
  //  checkSet = checkSet+2;
  n   = an;
  neF = aneF;
}

void snoptProblem_PAGMO::setObjective( integer aObjRow, doublereal aObjAdd )
{
  //checkSet = checkSet+2;
  ObjRow = aObjRow;
  ObjAdd = aObjAdd;
}

void snoptProblem_PAGMO::setA( integer alenA, integer *aiAfun,
			 integer *ajAvar, doublereal *aA )
{
  //checkSet = checkSet+4;
  lenA  = alenA;
  iAfun = aiAfun;
  jAvar = ajAvar;
  A     = aA;
}

void snoptProblem_PAGMO::setNeA( integer aneA )
{
  //checkSet = checkSet+1;
  neA = aneA;
}
void snoptProblem_PAGMO::setG( integer alenG, integer *aiGfun,
			 integer *ajGvar)
{
  //checkSet = checkSet+3;
  lenG  = alenG;
  iGfun = aiGfun;
  jGvar = ajGvar;
}

void snoptProblem_PAGMO::setNeG( integer aneG )
{
  //checkSet = checkSet+1;
  neG = aneG;
}

void snoptProblem_PAGMO::setX( doublereal *ax, doublereal *axlow, doublereal *axupp,
			 doublereal *axmul, integer *axstate )
{
  //checkSet = checkSet+5;
  x      = ax;
  xlow   = axlow;
  xupp   = axupp;
  xmul   = axmul;
  xstate = axstate;
}

void snoptProblem_PAGMO::setF( doublereal *aF, doublereal *aFlow, doublereal *aFupp,
			 doublereal *aFmul, integer *aFstate )
{
  //checkSet = checkSet+5;
  F      = aF;
  Flow   = aFlow;
  Fupp   = aFupp;
  Fmul   = aFmul;
  Fstate = aFstate;
}

void snoptProblem_PAGMO::setXNames( char *axnames, integer anxnames )
{
  //checkSet = checkSet+2;
  xnames  = axnames;
  nxnames = anxnames;
}

void snoptProblem_PAGMO::setFNames( char *aFnames, integer anFnames )
{
  //checkSet = checkSet+2;
  Fnames  = aFnames;
  nFnames = anFnames;
}

void snoptProblem_PAGMO::setProbName( const char *aProb )
{
  //checkSet = checkSet+1;
  sprintf(Prob, "%8s", aProb );
}

void snoptProblem_PAGMO::setUserFun( My_fp ausrfun )
{
  //checkSet = checkSet+1;
  usrfun = ausrfun;
}
