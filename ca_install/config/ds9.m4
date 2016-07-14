dnl ######################################################################
dnl
dnl File:	ds9.mp4
dnl
dnl Purpose:	Determine where the ds9 and xpa tools are.  
dnl
dnl Version:	$Id: ds9.m4,v 1.1 2009/05/28 13:37:25 pteuben Exp $
dnl
dnl
dnl ######################################################################


AC_PATH_PROGS(DS9,   [ds9] ,     no)
AC_PATH_PROGS(XPA,   [xpaset] ,  no)

dnl ##



