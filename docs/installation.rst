Instalation
***********

Dependencies
============

MIRZA
-----

Download and install `MIRZA <http://www.clipz.unibas.ch/index.php?r=tools/mirza/Submission/index>`_.
It is statically compiled so you do not need to compile it again

CONTRAfold
----------

Download and install `CONTRAfold <http://contra.stanford.edu/contrafold/download.html>`_.
It might be that you experience an error when compiling CONTRAfold. Something like this:

.. code-block:: bash

    In file included from LBFGS.hpp:52:0,
                     from InnerOptimizationWrapper.hpp:12,
                     from OptimizationWrapper.hpp:12,
                     from Contrafold.cpp:16:
    LBFGS.ipp: In instantiation of ‘Real LBFGS<Real>::Minimize(std::vector<T>&) [with Real = double]’:
    OptimizationWrapper.ipp:260:9:   required from ‘void OptimizationWrapper<RealT>::LearnHyperparameters(std::vector<int>, std::vector<T>&) [with RealT = double]’
    Contrafold.cpp:451:9:   required from ‘void RunTrainingMode(const Options&, const std::vector<FileDescription>&) [with RealT = double]’
    Contrafold.cpp:68:54:   required from here
    LBFGS.ipp:112:105: error: ‘DoLineSearch’ was not declared in this scope, and no declarations were found by argument-dependent lookup at the point of instantiation [-fpermissive]
    LBFGS.ipp:112:105: note: declarations in dependent base ‘LineSearch<double>’ are not found by unqualified lookup
    LBFGS.ipp:112:105: note: use ‘this->DoLineSearch’ instead
    make: *** [Contrafold.o] Error 1



To fix it:

* add -fpermissive flag to CSXXFLAGS in Makefile:

.. code-block:: c

    CXXFLAGS = -O3 -DNDEBUG -W -pipe -Wundef -Winline --param large-function-growth=100000 -Wall -fpermissive

    instead of
    CXXFLAGS = -O3 -DNDEBUG -W -pipe -Wundef -Winline --param large-function-growth=100000 -Wall

* add in Utilities.hpp:

.. code-block:: c

  #include <limits.h>

Jobber
------

Download and setup Jobber python library for workflow managment.

.. code-block:: bash

    pip install Jobber

After installation start the Jobber daemon:

.. code-block:: bash

    $ nohup jobber_server > jobber.log 2>&1 &

.. note::

    If you installed Jobber as user you might not have an access to the jobber_server. By
    default the binary location is $HOME/.local/bin and you have to export it in bash:

    .. code-block:: python

        $ export PATH="$HOME/.local/bin:$PATH"


    or add this statement to .bashrc file.

    jobber_server produces ~/.jobber/jobber.pid file that indicates whether the Jobber is already
    running. If the file exists one cannot start new instance of the jobber_server. This file is
    not clean when jobber_server is killed - only when it was stopped with stop command. Thus,
    after some crash one have to remove this file in order to start jobber_server again.


This will automatically create a ~/.jobber and ~/jobber/log directories and
it will put there config.py and executers.py files. Look at them and adjust
according to your needs.

Python
------

Install python modules:
 * Jobber (see upper paragraph)
 * drmaa (if you are going to submit it to the cluster)
 * statsmodels
 * pandas
 * BioPython
 * dendropy
 * numpy
 * scipy

