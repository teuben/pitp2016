dnl ######################################################################
dnl
dnl File:	distro.m4
dnl
dnl Purpose:	Determine unix distro name
dnl
dnl Version:	$Id: distro.m4,v 1.4 2009/06/29 17:22:49 pteuben Exp $
dnl
dnl
dnl ######################################################################


AC_MSG_CHECKING(for Unix name) 

DISTRO="Unknown"

case "$host" in

  *-darwin*)
    dnl Figure out that cat name (10.2 jaguar; 10.3 panther; 10.4 tiger 10.5 leopard )
    DISTRO_NAME=`sw_vers -productVersion`
    DISTRO="Mac"
    ;;

  *-cygwin*)
    DISTRO_NAME="Unknown"
    DISTRO="Cygwin"
    ;;

  *-linux*)
    linuxver=`uname -r | sed 's/-.*$//' | sed 's/\.[0-9]*$//'`
    linuxminver=`echo $linuxver | sed 's/^.*\.//'`
    linuxmajver=`echo $linuxver | sed 's/\..*$//'`
    for f in /etc/redhat-release /etc/debian_version /etc/issue; do
    	if test -f $f ; then
    	   DISTRO_NAME=`cat $f`
           break
  	fi
    done
    dnl  Ubuntu: 7.04 (Feisty Fawn) 7.10 (Gutsy Gibbon) 
    dnl          8.04 (Hardy Heron) 8.10 (Intrepid Ibex) 
    dnl          9.04 (Jaunty Jackalope) 
    dnl  Debian: ?   squeeze
    dnl  Debian: 5.0 lenny/sid
    dnl  Debian: 4.0 etch
    dnl  Debian: 3.1 sarge
    dnl  Debian: 3.0 woody
    dnl  Debian: 2.2 potato
    dnl  Debian: 2.1 slink
    dnl  Debian: 2.0 hamm

    dnl figure out which exists: /etc/redhat-release /etc/debian_version
    DISTRO="Linux"
    ;;

esac

AC_MSG_RESULT($DISTRO)

AC_MSG_CHECKING(for Distro name) 
AC_MSG_RESULT($DISTRO_NAME)

AC_SUBST(DISTRO)
AC_SUBST(DISTRO_NAME)



dnl echo done
