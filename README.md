Project Template
================

A project template for building shared libraries and binaries with cmake. It 
also installs headers and the proper cmake files so that you can use 
`find_package(this_project)` in subsequent cmake builds.

Getting started
---------------

Here we install to your home directory (instead of /usr/local):

```bash
git clone https://github.com/whit2333/project_template.git
mkdir project_build
cd project_build
cmake ../project_template/. -DCMAKE_INSTALL_PREFIX=$HOME
make && make install
```

```
export PATH=$HOME/bin:$PATH                         # if it is not included already
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH   # if it is not included already
```

Make sure the last two lines are in your bashrc if not a standard location.

Comments
--------

This template is not meant to be the best example and it is probably not using 
the most current features of cmake. So please don't use it as a sole source for 
learning cmake because cmake has wonderful documentation.

It is meant to provide a quick starting point. 


Please feel free to send questions or comments to Whitney Armstrong 
(whit@jlab.org). Also pull requests are greatly appreciated!

Todo
----

1. Add renaming script
2. Clean up cmake files
3. Document features/convention


