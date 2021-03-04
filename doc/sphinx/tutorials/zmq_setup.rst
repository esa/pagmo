.. _zeromq_setup:

=========
The Setup
=========

In order for ZeroMQ islands to find their peers and announce when they connect to and disconnect from the network, we must first set up the broker server. Our broker software of choice is Redis. This server acts as a bootstrap and communication channel. 

There are packages available in most distributions of Linux, but in this tutorial we will use the `debian:latest` Docker image. It's not necessary to install it in a container, but it makes it easier to reproduce the setup steps. 

First we will update the package list and install the Redis server:

.. code-block:: bash

    root@dbc834f0d0aa:/# apt-get -qq update
    root@dbc834f0d0aa:/# apt-get -qq install redis-server -y
    (output omitted...)
    root@dbc834f0d0aa:/#

Now we'll install the compilation tools and necessary libraries (this will take a while)

.. code-block:: bash
    
    root@dbc834f0d0aa:/# apt-get -qq install build-essential cmake git libboost-dev libboost-system-dev libboost-serialization-dev libboost-thread-dev libboost-python-dev libhiredis-dev libev-dev libzmq5-dev -y
    (output omitted...)

We also need to compile and install `redox` because Debian does not provide a binary package.

.. code-block:: bash

    root@dbc834f0d0aa:~# git clone https://github.com/hmartiro/redox /tmp/redox
    (output omitted...)
    root@dbc834f0d0aa:~# cd /tmp/redox/
    root@dbc834f0d0aa:/tmp/redox# mkdir build && cd build
    root@dbc834f0d0aa:/tmp/redox/build# cmake .. -DLIB_SUFFIX="/"
    (output omitted...)
    root@dbc834f0d0aa:/tmp/redox/build# make && make install && ldconfig

Next, we'll clone the PaGMO source tree and compile it with PyGMO and ZeroMQ enabled. The `make` command takes around 15-20 minutes to complete (could take longer), so this is the perfect time to grab a caffeinated beverage.

.. code-block:: bash

    root@dbc834f0d0aa:/# git clone https://github.com/jdiez17/pagmo /tmp/pagmo
    root@dbc834f0d0aa:/# cd /tmp/pagmo
    root@dbc834f0d0aa:/tmp/pagmo# mkdir build && cd build 
    root@dbc834f0d0aa:/tmp/pagmo/build# cmake .. -DENABLE_ZMQ=ON -DBUILD_PYGMO=ON
    root@dbc834f0d0aa:/tmp/pagmo/build# make && make install && ldconfig

If everything went well, you should be able to import PyGMO and zmq_island from Python!

.. code-block:: bash

    root@dbc834f0d0aa:/tmp/pagmo/build# python
    Python 2.7.9 (default, Mar  1 2015, 12:57:24) 
    [GCC 4.9.2] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import PyGMO
    >>> 'zmq_island' in dir(PyGMO)
    True
    >>> 

Finally, we'll start a Redis server and daemonize it so it runs in the background.

.. code-block:: bash

    root@dbc834f0d0aa:/# redis-server --daemonize yes

At this point, you are ready to run ZMQ islands!
