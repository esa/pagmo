.. _howtoinstall:

Install PyGMO
======================

Compiling and Installing under Unix
-----------------------------------

Assuming you have just downloaded the source code following the instructions given, see :ref:`howtodownload`, you will have 
created a directory pagmo in your current directory, move there::

  cd pagmo

You will now need to create a build directory where to build the source code, so::

  mkdir build

You can now move there::

  cd build

and have ccmake help you select the options that are most suitable for you::

  ccmake ../

After pressing c once, a typical ccmake screen will look like 

.. image:: images/ccmake1.png

note that the PyGMO option is selcted as well as other options requiring external libraries ....

At this point (after pressing c again) in case all required system libraries are found
you should be seeing something like this on the screen:

.. image:: images/ccmake2.png


You can now press 'g' to generate a make file and exit ccmake utility. You are back to the prompt where you can now type::

  make

and::

  sudo make install

Watch carefully the message in the terminal where the installation path is given to check 
that the correct python dist-packages or site-packages directory has been located

Here is a typical example of the output obtained (gentoo system):

.. image:: images/install1.png

Compiling and Installing under Windows
--------------------------------------

Same as under Unix, just make sure that

* You have compiled the boost libraries correctly (i.e invoking bjam with the option toolset=gcc link=shared). 
* Place the whole boost directory where the CMake script can find it (e.g. in C:/boost). This may also require renaming the folder from boost_x_xx_xx to boost)
* Check, when running CMake, that all libraries are found correctly
* When running a make install, Windows will probably put your PyGMO directory under Program Files/pagmo,
  move it to the correct place (e.g. C:/Python27/Lib/site-packages/)
* Put all dll in PyGMO/core
* Hope for the best (honestly, if it works for you just download the binaries.... it is easier)