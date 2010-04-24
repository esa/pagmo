#ifndef SNOPTPROBLEM_PAGMO
#define SNOPTPROBLEM_PAGMO

#include "snopt_PAGMO.h"
#include "../../problem/base.h"
#include "../../algorithm/snopt.h"


// class snoptProblem_PAGMO performs problem set-up, initialization,
// and problem-specific variable storage.
class snoptProblem_PAGMO {
protected:
  //************************************************************
  // The variables below must be set by the user with initialize
  char       Prob[200];

  integer    n, neF;
  integer    ObjRow;
  doublereal ObjAdd;

  integer    lenA, neA;
  integer    *iAfun, *jAvar;
  doublereal *A;
  integer    lenG, neG;
  integer    *iGfun, *jGvar;

  doublereal *x, *xlow, *xupp, *xmul;
  doublereal *F, *Flow, *Fupp, *Fmul;

  integer *xstate, *Fstate;

  char *xnames, *Fnames;
  integer nxnames, nFnames;

  My_fp usrfun;
  //***********************************************************
  //***********************************************************
  void userDataSet();
  void errMsgExit( const char *var );

  integer    inform, fortranStyleObj, fortranStyleAG;
  integer    minrw, miniw, mincw;
  integer    lenrw, leniw, lencw;
  doublereal *rw;
  integer    *iw;
  char       *cw;

  integer    iSpecs, iSumm, iPrint;
  integer    initCalled;
  char       specname[200], printname[200];
  integer    spec_len, prnt_len;

  //added for PAGMO compatibility
  char       *cu;
  integer    lencu;
  double       *ru;
  integer     lenru;

  void init2zero();
  void setMemory();
  void alloc    ( integer lencw, integer leniw, integer lenrw );
  void realloc  ( integer lencw, integer leniw, integer lenrw );
  void memcpyIn ( char *tcw, integer *tiw, doublereal *trw,
		  integer tlencw, integer tleniw, integer tlenrw);
  void memcpyOut( char *tcw, integer *tiw, doublereal *trw,
		  integer tlencw, integer tleniw, integer tlenrw);
public:
  snoptProblem_PAGMO(const pagmo::problem::base& problem, pagmo::algorithm::snopt::preallocated_memory *di_comodo);
  ~snoptProblem_PAGMO();
  void increment();
  void decrement();
  void computeJac();

  int snmema( integer &mincw, integer &miniw, integer &minrw);

  void init();
  void setParameter    ( char *stroptin );
  void getParameter    ( char *stroptin, char *stroptout );
  void setIntParameter ( const char *stropt,   integer     opt );
  void getIntParameter ( const char *stropt,   integer    &opt );
  void setRealParameter( const char *stropt,   doublereal  opt );
  void getRealParameter( const char *stropt,   doublereal &opt );
  void solve           ( integer starttype );
  void setPrintFile    ( const char printname[] );
  void setSpecFile     ( char specname[] );

  // Functions that set up the problem data:
  void setProblemSize( integer n, integer neF );
  void setObjective  ( integer ObjRow, doublereal ObjAdd );
  void setA          ( integer lenA, integer *iAfun, integer *jAvar,
		       doublereal *A );
  void setNeA        ( integer neA );
  void setG          ( integer lenG, integer *iGfun, integer *jGvar );
  void setNeG        ( integer neG );
  integer getNeG        (){return neG;}
  void setX          ( doublereal *x, doublereal *xlow, doublereal *xupp,
		       doublereal *xmul, integer *xstate );
  void setF          ( doublereal *F, doublereal *Flow, doublereal *Fupp,
		       doublereal *Fmul, integer *Fstate );
  void setXNames     ( char *xnames, integer nxnames );
  void setFNames     ( char *Fnames, integer nFnames );
  void setProbName   ( const char *Prob );
  void setUserFun    ( My_fp usrfun );
  void activate_screen_output() {iSumm = 6;}
  void deactivate_screen_output() {iSumm = 0;}
};

#endif
