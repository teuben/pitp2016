# pitp2016
Useful notes for PITP 2016 software (Linux and Mac)

# Package Management Reminders

For most homeworks you will need to install programs and python modules.
Commands such as "make", "gcc", "gfortran", "autoconf", "cmake", "ipython"
should all work. A frequent source of nuisance is that you think you have
a library, but when you try and link it with your program, it fails. Many
package managers split the runtime from the compile time.


## Linux Debian-style (apt)

  sudo apt install libgsl-dev       # install package(s)
  dpkg -S /usr/bin/ls               # which packate does /bin/ls belong to?
  dpkg -L coreutils                 # lists files in a package
  dpkg --list                       # lists all packages you have


## Linux Redhat-style (yum)

  sudo yum install NAME1            # install packages

  rpm -qf /bin/ls                   # which packate does /bin/ls belong to?
  rpm -ql coreutils                 # lists files in a package
  rpm -qa                           # lists all packages you have

## Mac HomeBrew-style

A virgin mac doesn't have a package manager.
Before anything, you will first need to install Xcode.


## Mac MacPorts-style

A virgin mac doesn't have a package manager. 
Before anything, you will first need to install Xcode.

## Mac Fink-style

A virgin mac doesn't have a package manager. 
Before anything, you will first need to install Xcode.




