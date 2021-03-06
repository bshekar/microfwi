/* This file is automatically generated. DO NOT EDIT! */

#ifndef _smth2d_h
#define _smth2d_h


/* Prototype for function used internally */
/*****************************************************************************
*****************************************************************************/
void tripd(float *d, float *e, float *b, int n);
/*<****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b
Input:
d[]     the diagonal of A 
e[]     the superdiagonal of A
b[]     the rhs of Ax=b

Output:
b[]     b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
****************************************************************************>*/


/***************************************************************************
***************************************************************************/
void smooth2d(float **v, int n1, int n2, float r1, float r2, float rw);
/*<**************************************************************************
	int n1;		 number of points in x1 (fast) dimension
	int n2;		 number of points in x1 (fast) dimension 
        float **v0;       array of input velocities 
        float r1;        smoothing parameter for x1 direction
        float r2;        smoothing parameter for x2 direction
        float rw;        smoothing parameter for window
" Notes:								",
" Larger r1 and r2 result in a smoother data. Recommended ranges of r1 	", 
" and r2 are from 1 to 20.						",
"									",
" Smoothing can be implemented in a selected window. The range of 1st   ",
" dimension for window is from win[0] to win[1]; the range of 2nd   	",
" dimension is from win[2] to win[3]. 					",
"									",
" Smoothing the window function (i.e. blurring the edges of the window)	",
" may be done by setting a nonzero value for rw, otherwise the edges	",
" of the window will be sharp.						",
" 									",
**************************************************************************>*/

#endif
