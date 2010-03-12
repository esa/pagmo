/* Mike Gertz - 2-Aug-98 */

#ifndef SNFILEWRAPPER
#define SNFILEWRAPPER
#include "snopt_PAGMO.h"
#ifdef __cplusplus
extern "C" {
#endif

  void snopenappend_ (integer *iunit, char *name,
		      integer *inform, ftnlen name_len);

  void snfilewrapper_(char *name__, integer *ispec, integer *inform__,
		      char *cw, integer *lencw, integer *iw,
		      integer *leniw, doublereal *rw, integer *lenrw,
		      ftnlen name_len, ftnlen cw_len);

  void snclose_      (integer *iunit);

  void snopenread_   (integer *iunit, char *name, integer *inform,
		      ftnlen name_len);

#ifdef __cplusplus
}
#endif

#endif
