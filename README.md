Project Template
================

A project template for building shared libraries and binaries with cmake. It also installs headers and the proper cmake files so that you can use `find_package(this_project)` in subsequent cmake builds.

Getting started
---------------

Here we install to your home directory (instead of /usr/local):

```bash
git clone https://github.com/whit2333/project_template.git
mkdir project_build
cd project_build
cmake ../project_template/. -DCMAKE_INSTALL_PREFIX=$HOME
make && make install
export PATH=$HOME/bin:$PATH                         # if it is not included already
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH   # if it is not included already
```

Make sure the last two lines are in your bashrc if not a standard location.


Send questions to Whitney Armstrong ( whit@jlab.org )

