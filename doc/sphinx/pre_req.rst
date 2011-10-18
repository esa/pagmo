Preparing your system to compile PyGMO
======================================


We assume you have root access to your unix machine (or you are administering your Windows machine)
and a working python installation. As PyGMO has some dependencies, you will need to install
a few packages before being able to install PyGMO on your machine. All are supported by most package managers,
so that a simple apt-get or emerge should suffice to prepare your system.


* Install `git <http://git-scm.com>`_
* Install `CMake <http://www.cmake.org>`_ with its ccmake utility
* Install `boost libraries <http://www.boost.org>`_ both headers and compiled libraries

Other packages are actually optional and they will enhance PyGMO functionalities if present:

* SNOPT (proprietary)
* `IPOPT <https://projects.coin-or.org/Ipopt>`_
* `SciPy <http://www.scipy.org/>`_
* `NLOPT <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_ (compiled with the c++ flag activated)
* `GSL <http://www.gnu.org/s/gsl/>`_ (version 1.15 required)
* `PyKEP <http://keptoolbox.sourceforge.net/>`_ (version 1.15 required)

These packages need to be compiled in such a way as to allow PyGMO 1) to find them 2) tho use them.
See the dedicated pages on how to make sure this happens.


