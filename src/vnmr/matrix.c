/*
 * Copyright (C) 2015  University of Oregon
 *
 * You may distribute under the terms of either the GNU General Public
 * License or the Apache License, as specified in the LICENSE file.
 *
 * For more information, see the LICENSE file.
 */

/******************************************************************
*  matrix.c  -  symmetric and non symmetric matrix transposition  *
*               either in-place or out-of-place                   *
******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnmrsys.h"			/* To define UNIX (or VMS)	*/

#pragma GCC optimize ("O3")

#define COMPLETE	0
#define ERROR		1
#define REAL		1		/* real word length		*/
#define COMPLEX		2		/* complex word length		*/
#define HYPERCOMPLEX	4		/* hypercomplex word length	*/
#define BUF_INIT	0x1		/* initialize transpose buffer	*/
#define BUF_ALLOCATE	0x2		/* allocate transpose buffer	*/
#define BUF_RELEASE	0x4		/* release transpose buffer	*/

struct _fcomplex
{
   float	re;
   float	im;
};

struct _hypercomplex
{
   float	rere;
   float	reim;
   float	imre;
   float	imim;
};

struct xposebuf
{
   int          buffersize;
   float        *bufferpntr;
};

typedef struct _fcomplex	fcomplex;
typedef struct _hypercomplex	hypercomplex;

static struct xposebuf  xpose;

static int xposebufalloc(int, int);


/*---------------------------------------
|					|
|	     initxposebuf()/0		|
|					|
+--------------------------------------*/
void initxposebuf()
{
   xposebufalloc(0, BUF_INIT);
}


/*---------------------------------------
|					|
|	     closexposebuf()/0		|
|					|
+--------------------------------------*/
void closexposebuf()
{
   xposebufalloc(0, BUF_RELEASE);
}


/*---------------------------------------
|					|
|	       invert4()/5		|
|					|
+--------------------------------------*/
static void invert4(float *dmatrix, int offset, int size,
                    int inum, int jnum, int max2)
{
   register int		m,
			l,
			k,
			i,
			j;
   register float	s,
			*i1,
			*i2,
			*i3,
			*i4,
			*submatrix;


   submatrix = dmatrix + offset;
   m = max2;

   for (j = 0; j < (size - 1); j++)
   {
      for (i = (j + 1); i < size; i++)
      {
         i3 = submatrix + (i*inum) + (max2*j);
	 i4 = submatrix + j + (max2*i*jnum);
	 k = inum;
	 while (k--)
	 {
            i1 = i3++;
	    i2 = i4++;
	    l = jnum;
	    while (l--)
	    {
               s = *i1;
               *i1 = *i2;
               *i2 = s;
	       i1 += m;
               i2 += m;
            }
	 }
      }
   }
}


/*---------------------------------------
|					|
|	       invert8()/5		|
|					|
+--------------------------------------*/
static void invert8(float *dmatrix, int offset, int size,
                    int inum, int jnum, int max2)
{
   register int		m,
			l,
			k,
			i,
			j;
   register float	s;
   register fcomplex	*i1,
			*i2,
			*i3,
			*i4,
			*submatrix;


   submatrix = ((fcomplex *) (dmatrix)) + offset;
   m = max2;

   for (j = 0; j < (size - 1); j++)
   {
      for (i = (j + 1); i < size; i++)
      {
         i3 = submatrix + (i*inum) + (max2*j);
	 i4 = submatrix + j + (max2*i*jnum);
	 k = inum;
	 while (k--)
	 {
            i1 = i3++;
	    i2 = i4++;
	    l = jnum;
	    while (l--)
	    {
               s = i1->re;
               i1->re = i2->re;
               i2->re = s;
	       s = i1->im;
               i1->im = i2->im;
               i2->im = s;
	       i1 += m;
               i2 += m;
            }
	 }
      }
   }
}


/*---------------------------------------
|					|
|	      invert16()/5		|
|					|
+--------------------------------------*/
static void invert16(float *dmatrix, int offset, int size,
                     int inum, int jnum, int max2)
{
   register int			m,
				l,
				k,
				i,
				j;
   register float		s;
   register hypercomplex	*i1,
				*i2,
				*i3,
				*i4,
				*submatrix;


   submatrix = ((hypercomplex *) (dmatrix)) + offset;
   m = max2;

   for (j = 0; j < (size - 1); j++)
   {
      for (i = (j + 1); i < size; i++)
      {
         i3 = submatrix + (i*inum) + (max2*j);
	 i4 = submatrix + j + (max2*i*jnum);
	 k = inum;
	 while (k--)
	 {
            i1 = i3++;
	    i2 = i4++;
	    l = jnum;
	    while (l--)
	    {
               s = i1->rere;
               i1->rere = i2->rere;
               i2->rere = s;
	       s = i1->reim;
               i1->reim = i2->reim;
               i2->reim = s;
               s = i1->imre;
               i1->imre = i2->imre;
               i2->imre = s;
               s = i1->imim;
               i1->imim = i2->imim;
               i2->imim = s;
	       i1 += m;
               i2 += m;
            }
	 }
      }
   }
}


/*---------------------------------------
|					|
|	       itrans()/4		|
|					|
+--------------------------------------*/
static void itrans(float *matrix, int max2, int max1, int datatype)
{
   int	ii,
	kk,
	matsize,
	blocks,
	length1,
	length2;
   void (*invertfunc)();	/* pointer to a function */


   if (datatype == HYPERCOMPLEX)
   {
      invertfunc = invert16;
   }
   else if (datatype == COMPLEX)
   {
      invertfunc = invert8;
   }
   else
   {
      invertfunc = invert4;
   }


   if (max1 >= max2)
   {
      blocks = max1/max2; 
      matsize=2;
      for (length1 = max1/2; length1 >= max2; length1 /= 2)
      {
         for (length2 = 1; length2 <= max2/2; length2 *= 2) 
         {
            for (ii = 0; ii < max2/(matsize*length2); ii++)
            {
               for (kk = 0; kk < max1/(matsize*length1); kk++)
	       {
                  (*invertfunc) (matrix,
			       matsize * (ii*length2 + max2*kk*length1),
			       matsize, length2, length1, max2);
               }
            }
         }
      }

      matsize = max2;
      for (ii = 0; ii < blocks; ii++)
         (*invertfunc) (matrix, max2*matsize*ii, matsize, 1, 1, max2);
   }
   else
   {
      blocks = max2/max1; 
      matsize = max1;
      for (ii = 0; ii < blocks; ii++)
         (*invertfunc) (matrix, matsize*ii, matsize, 1, 1, max2);

      matsize=2;
      for (length2 = max1; length2 <= max2/2; length2 *= 2)
      {
         for (length1 = max1/2; length1 >= 1; length1 /= 2)
         {
            for (ii = 0; ii < max1/(matsize*length1); ii++)
            {
               for (kk = 0; kk < max2/(matsize*length2); kk++)
	       {
                  (*invertfunc) (matrix,
			       matsize * (kk*length2 + max2*ii*length1),
			       matsize, length2, length1, max2);
	       }
            }
         }
      }
   }
}


/*-----------------------------------------------
|                                               |
|              xposebufalloc()/2                |
|                                               |
|   This function allocates a block of memory   |
|   for out-of-place transposition of a non-    |
|   symmetric matrix.                           |
|                                               |
+----------------------------------------------*/
int xposebufalloc(int nbytes, int bufstatus)
{
   if (bufstatus & BUF_INIT)
   {
      xpose.buffersize = 0;
      xpose.bufferpntr = NULL;
   }

   if (bufstatus & BUF_RELEASE)
   {
      xpose.buffersize = 0;
      if (xpose.bufferpntr != NULL)
      {
         free((char *)xpose.bufferpntr);
         xpose.bufferpntr = NULL;
      }
   }
 
   if (bufstatus & BUF_ALLOCATE)
   {
      if (xpose.buffersize != nbytes)
      {
         if (xpose.bufferpntr != NULL)
         {
            free((char *)xpose.bufferpntr);
            xpose.bufferpntr = NULL;
            xpose.buffersize = 0;
         }

         xpose.buffersize = nbytes;
         xpose.bufferpntr = (float *) (malloc( (unsigned) (nbytes) ));
         if (xpose.bufferpntr == NULL)
         {
            xpose.buffersize = 0;
            return(ERROR);
         }
      }
   }

   return(COMPLETE);
}   
 
 
/*-----------------------------------------------
|                                               |
|                  symxpose()/3                 |
|                                               |
|   This function performs an in-place trans-   |
|   position of a real, complex, or hyper-      |
|   complex half-transformed 2D data set.       |
|                                               |
+----------------------------------------------*/

/* SIMD-friendly blocks for symxpose's inner loop */
#define RBS 8
#define CBS 4
 
static inline __attribute__((always_inline)) void t8x8(const float s[RBS][RBS], float d[RBS][RBS])
{ for (int i=0;i<RBS;i++) for (int j=0;j<RBS;j++) d[j][i]=s[i][j]; }
 
static inline __attribute__((always_inline)) void t4x4(const double s[CBS][CBS], double d[CBS][CBS])
{ for (int i=0;i<CBS;i++) for (int j=0;j<CBS;j++) d[j][i]=s[i][j]; }
 
static void symxpose_real(float *mat, int n)
{
   int i, j;
   for (i = 0; i + RBS <= n; i += RBS) {
      float in[RBS][RBS], out[RBS][RBS];
      for (int k=0;k<RBS;k++) memcpy(in[k], mat+(size_t)(i+k)*n+i, RBS*sizeof(float));
      t8x8(in, out);
      for (int k=0;k<RBS;k++) memcpy(mat+(size_t)(i+k)*n+i, out[k], RBS*sizeof(float));
      for (j = i + RBS; j + RBS <= n; j += RBS) {
         float a_in[RBS][RBS], a_out[RBS][RBS], b_in[RBS][RBS], b_out[RBS][RBS];
         for (int k=0;k<RBS;k++) {
            memcpy(a_in[k], mat+(size_t)(i+k)*n+j, RBS*sizeof(float));
            memcpy(b_in[k], mat+(size_t)(j+k)*n+i, RBS*sizeof(float));
         }
         t8x8(a_in, a_out); t8x8(b_in, b_out);
         for (int k=0;k<RBS;k++) {
            memcpy(mat+(size_t)(i+k)*n+j, b_out[k], RBS*sizeof(float));
            memcpy(mat+(size_t)(j+k)*n+i, a_out[k], RBS*sizeof(float));
         }
      }
   }
   for (int r = 0; r < n; r++) {
      int cstart = (r >= i) ? 0 : i;
      for (int c = (r+1 > cstart ? r+1 : cstart); c < n; c++) {
         float t = mat[r*n+c]; mat[r*n+c] = mat[c*n+r]; mat[c*n+r] = t;
      }
   }
}
 
static void symxpose_complex(float *matf, int n)
{
   double *mat = (double *)matf;
   int i, j;
   for (i = 0; i + CBS <= n; i += CBS) {
      double in[CBS][CBS], out[CBS][CBS];
      for (int k=0;k<CBS;k++) memcpy(in[k], mat+(size_t)(i+k)*n+i, CBS*sizeof(double));
      t4x4(in, out);
      for (int k=0;k<CBS;k++) memcpy(mat+(size_t)(i+k)*n+i, out[k], CBS*sizeof(double));
      for (j = i + CBS; j + CBS <= n; j += CBS) {
         double a_in[CBS][CBS], a_out[CBS][CBS], b_in[CBS][CBS], b_out[CBS][CBS];
         for (int k=0;k<CBS;k++) {
            memcpy(a_in[k], mat+(size_t)(i+k)*n+j, CBS*sizeof(double));
            memcpy(b_in[k], mat+(size_t)(j+k)*n+i, CBS*sizeof(double));
         }
         t4x4(a_in, a_out); t4x4(b_in, b_out);
         for (int k=0;k<CBS;k++) {
            memcpy(mat+(size_t)(i+k)*n+j, b_out[k], CBS*sizeof(double));
            memcpy(mat+(size_t)(j+k)*n+i, a_out[k], CBS*sizeof(double));
         }
      }
   }
   for (int r = 0; r < n; r++) {
      int cstart = (r >= i) ? 0 : i;
      for (int c = (r+1 > cstart ? r+1 : cstart); c < n; c++) {
         double t = mat[r*n+c]; mat[r*n+c] = mat[c*n+r]; mat[c*n+r] = t;
      }
   }
}
 
static void symxpose_hyper(float *mat, int n)
{
   const int TILE = 32;
   for (int bi = 0; bi < n; bi += TILE) {
      int ih = bi+TILE<n ? TILE : n-bi;
      for (int bj = bi; bj < n; bj += TILE) {
         int jw = bj+TILE<n ? TILE : n-bj;
         for (int i = 0; i < ih; i++) {
            int j0 = (bj==bi) ? i+1 : 0;
            for (int j = j0; j < jw; j++) {
               int r = bi+i, c = bj+j;
               if (r >= c && bj == bi) continue;
               float tmp[4];
               memcpy(tmp, mat+((size_t)r*n+c)*4, 16);
               memcpy(mat+((size_t)r*n+c)*4, mat+((size_t)c*n+r)*4, 16);
               memcpy(mat+((size_t)c*n+r)*4, tmp, 16);
            }
         }
      }
   }
}

/* npoints      number of "datatype" points			*/
/* datatype     1 = real    2 = complex    4 = hypercomplex	*/
static void symxpose(float *matrix, int npoints, int datatype)
{
  if (datatype == HYPERCOMPLEX)      symxpose_hyper(matrix, npoints);
  else if (datatype == COMPLEX)      symxpose_complex(matrix, npoints);
  else if (datatype == REAL)         symxpose_real(matrix, npoints);
}

static int is_pow2(int n)
{
   return (n > 0) && ((n & (n - 1)) == 0);
}
 
/*-----------------------------------------------
|                                               |
|               nonsymxpose()/4                 |
|                                               |
|   This function performs an out-of-place      |
|   transpose of a non-symmetric matrix.  The   |
|   transposed matrix is stored back in the     |
|   original matrix.  All temporary storage is  |
|   handled within this routine.                |
|                                               |
+----------------------------------------------*/

/* SIMD-friendly blocks for nonsymxpose's copy-into-buffer */
static inline __attribute__((always_inline)) void tile_real_np(const float *src, int src_ld, float *dst, int dst_ld)
{
    float in[RBS][RBS], out[RBS][RBS];
    for (int i=0;i<RBS;i++) memcpy(in[i], src+(size_t)i*src_ld, RBS*sizeof(float));
    t8x8(in, out);
    for (int i=0;i<RBS;i++) memcpy(dst+(size_t)i*dst_ld, out[i], RBS*sizeof(float));
}
 
static inline __attribute__((always_inline)) void tile_complex_np(const double *src, int src_ld, double *dst, int dst_ld)
{
    double in[CBS][CBS], out[CBS][CBS];
    for (int i=0;i<CBS;i++) memcpy(in[i], src+(size_t)i*src_ld, CBS*sizeof(double));
    t4x4(in, out);
    for (int i=0;i<CBS;i++) memcpy(dst+(size_t)i*dst_ld, out[i], CBS*sizeof(double));
}
 
static void scalar_block_real(const float *src, float *dst, int bi, int bj, int bh, int bw, int src_cols, int dst_cols)
{
    for (int i=0;i<bh;i++) for (int j=0;j<bw;j++)
        dst[(bj+j)*dst_cols+(bi+i)] = src[(bi+i)*src_cols+(bj+j)];
}
 
static void scalar_block_complex(const double *src, double *dst, int bi, int bj, int bh, int bw, int src_cols, int dst_cols)
{
    for (int i=0;i<bh;i++) for (int j=0;j<bw;j++)
        dst[(bj+j)*dst_cols+(bi+i)] = src[(bi+i)*src_cols+(bj+j)];
}
 
static void xpose_into_buffer_real(const float *matrix, float *buffer, int ncols, int nrows)
{
    int i, j;
    for (i = 0; i + RBS <= nrows; i += RBS) {
        for (j = 0; j + RBS <= ncols; j += RBS)
            tile_real_np(matrix+(size_t)i*ncols+j, ncols, buffer+(size_t)j*nrows+i, nrows);
        if (j < ncols) scalar_block_real(matrix, buffer, i, j, RBS, ncols-j, ncols, nrows);
    }
    if (i < nrows) scalar_block_real(matrix, buffer, i, 0, nrows-i, ncols, ncols, nrows);
}
 
static void xpose_into_buffer_complex(const float *matrixf, float *bufferf, int ncols, int nrows)
{
    const double *matrix = (const double *)matrixf;
    double *buffer = (double *)bufferf;
    int i, j;
    for (i = 0; i + CBS <= nrows; i += CBS) {
        for (j = 0; j + CBS <= ncols; j += CBS)
            tile_complex_np(matrix+(size_t)i*ncols+j, ncols, buffer+(size_t)j*nrows+i, nrows);
        if (j < ncols) scalar_block_complex(matrix, buffer, i, j, CBS, ncols-j, ncols, nrows);
    }
    if (i < nrows) scalar_block_complex(matrix, buffer, i, 0, nrows-i, ncols, ncols, nrows);
}
 
static void xpose_into_buffer_hyper(const float *matrix, float *buffer, int ncols, int nrows)
{
    const int TILE = 32;
    for (int bi = 0; bi < nrows; bi += TILE) {
        int ih = bi+TILE<nrows ? TILE : nrows-bi;
        for (int bj = 0; bj < ncols; bj += TILE) {
            int jw = bj+TILE<ncols ? TILE : ncols-bj;
            for (int i = 0; i < ih; i++)
                for (int j = 0; j < jw; j++)
                    memcpy(buffer+((size_t)(bj+j)*nrows+(bi+i))*4,
                           matrix+((size_t)(bi+i)*ncols+(bj+j))*4, 4*sizeof(float));
        }
    }
}

static int nonsymxpose(float *matrix, int ncols, int nrows, int datatype)
{
   int                  nbytes;
   register int         i,
                        j,
                        skip1,
                        skip2;
 

   skip1 = datatype;
   skip2 = skip1 * (nrows - 1);
   nbytes = sizeof(float) * nrows * ncols * skip1;
 
   if (xposebufalloc(nbytes, BUF_ALLOCATE))
   {
/* itrans()'s block-recursive in-place algorithm
 * only produces correct output when BOTH dimensions 
 * are individually powers of two  */
	  if (!is_pow2(nrows) || !is_pow2(ncols))
         return(ERROR);
      itrans(matrix, ncols, nrows, datatype);
      return(COMPLETE);
   }

   if (skip1 == REAL)
      xpose_into_buffer_real(matrix, xpose.bufferpntr, ncols, nrows);
   else if (skip1 == COMPLEX)
      xpose_into_buffer_complex(matrix, xpose.bufferpntr, ncols, nrows);
   else if (skip1 == HYPERCOMPLEX)
      xpose_into_buffer_hyper(matrix, xpose.bufferpntr, ncols, nrows);
 
   memcpy(matrix, xpose.bufferpntr, nbytes);

   return(COMPLETE);
}


/*-----------------------------------------------
|                                               |
|                 transpose()/4                 |
|                                               |
|   This function calls the appropriate xpose   |
|   routine.                                    |
|                                               |
+----------------------------------------------*/
int transpose(void *matrix, int ncols, int nrows, int datatype)
{
   if (nrows == ncols)
   {
      symxpose((float *)matrix, nrows, datatype);
   }
   else
   {
      if (nonsymxpose((float *)matrix, ncols, nrows, datatype))
         return(ERROR);
   }
 
   return(COMPLETE);
}
