dnl ######################################################################
dnl
dnl File:	cvs.m4
dnl
dnl Purpose:	Determine where cvs and svn are installed
dnl
dnl Version:	$Id: cvs.m4,v 1.1 2009/06/02 20:35:59 pteuben Exp $
dnl
dnl
dnl ######################################################################


AC_PATH_PROGS(CVS,   [cvs] ,  no)
AC_PATH_PROGS(SVN,   [svn] ,  no)

dnl ##



