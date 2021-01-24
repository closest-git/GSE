#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "..\pch\grus_def.h"
#include "..\pch\grus_version.h"
#include "..\pch\grus_version_lib.h"

#define GRUS_GSE_DLL __declspec(dllexport)
#include "GSE_Krylov.h"

#define NEV_select 20

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static IMPLEMENT_ONE_ZERO( one_,fuyi_,zero_ );	
//DECLARE_GRUS_TRANS( )
//DECLARE_ONE_ZERO( one_,fuyi_,zero_ )

/*
	7/10/2007
		v0.1	cys
*/
double CS_eigen_resi( int dim,int *ptr,int *ind,double *val,DoubleComplex lenda,DoubleComplex* v )	{
	DoubleComplex *z=malloc( sizeof(DoubleComplex)*dim ),s,t;
	int i,j;
	double off=0.0,real,imag;

	memset( z,0x0,sizeof(DoubleComplex)*dim );
	for( i = 0; i < dim; i++ )	{
		ASSIGN( s,v[i] );
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			SET( t,val[j],0.0 );
			MULT_ADD( z[ind[j]],t,s );	
		}
	}

	for( i = 0; i < dim; i++ )	{
		MULT_SUB( z[i],lenda,v[i] );
		real=REAL_COMPONENT(z[i]),		imag=IMAG_COMPONENT(z[i]);
		off += (real*real+imag*imag);
	}
	off = sqrt( off );

	return off;
}

/*
	call dgeev to computes the eigenvalues and right eigenvectors 

	7/9/2007
		v0.1	cys
	7/16/2007		增加直接读入满阵的功能
		v0.2	cys
*/
void GSE_QR_D( int nCol,int *ptr,int *ind,double *val,int flag,double sft[2] )	{
	int i,j,k,n,ldv,lwork=-1,info,no_min=-1,*select,*NO;
	double *A,*wr,*wi,*work,*vl=0x0,*vr=0x0,off,off_min=1.e100,res;
	char V[1]={'V'},N[1]={'N'},sPath[80];	//only eigenvalue,no schur vector
	DoubleComplex lenda[NEV_select],*ev;

	n=nCol;	
	if( 1 )	{
		A=malloc( sizeof(double)*n*(n+2) );
		memset( A,0x0,sizeof(double)*n*n );
		for( i = 0; i < n; i++ )	{
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				A[i*n+ind[j]]=val[j]; 
			}
		}
	}else	{
//		DGE_load( &n,&A,"H:\\GSP\\GSE\\test\\H.dat",0x0 );
	}
	wr=malloc( sizeof(double)*2),			wi=wr+n;
	ldv=n;
	ev=malloc( sizeof(DoubleComplex)*NEV_select*n );

	lwork=-1;
	work=wr;
	dgeev(N, N,&n,A,&n, wr,wi, vl,&ldv, vr, &ldv, work,&lwork,&info);
	ASSERT( info==0 );
	lwork=work[0];
	work=malloc( sizeof(double)*(lwork+n*n) );		vr=work+lwork;
	dgeev(N, V,&n,A,&n, wr,wi, vl,&ldv, vr, &ldv, work,&lwork,&info);

	select = malloc( sizeof(int)*(n+NEV_select) );	NO=select+n;
	for( i = 0; i < n; i++ )	select[i]=0;
	for( i = 0; i < NEV_select; i++ )	{
		no_min=-1;		off_min=1.e100;		
		for( j = 0; j < n; j++ )	{
			if( select[j]!=0 )
				continue;
			off=sqrt( (wr[j]-sft[0])*(wr[j]-sft[0]) + (wi[j]-sft[1])*(wi[j]-sft[1]) );
			if( off < off_min )	
			{	off_min=off;	no_min=j;		}			
		}
		select[no_min]=i+1;
		if( wi[no_min] != 0.0 )	{
			select[no_min+1]=-i-1;
		}
	}

	for( i = 0; i < n; i++ )	{
		if( select[i]<=0 )
			continue;
		NO[select[i]-1]=i;
	}
	printf( "--NO--------------------LENDA-----------|Ax-lenda*x|----|lenda-shift|--\r\n" );
	for( i = 0; i < NEV_select; i++ )	{
		j=NO[i];
		off=sqrt( (wr[j]-sft[0])*(wr[j]-sft[0]) + (wi[j]-sft[1])*(wi[j]-sft[1]) );
		SET( lenda[i],wr[j],wi[j] );
		if( wi[j]==0.0 )	{
			for( k = 0; k < n; k++ )	
				SET( ev[i*n+k],vr[j*n+k],0.0 );
		}else	{
			for( k = 0; k < n; k++ )	
				SET( ev[i*n+k],vr[j*n+k],vr[(j+1)*n+k] );
		}
		res = CS_eigen_resi( nCol,ptr,ind,val,lenda[i],ev+i*n );
		printf( "%3d: (%.12g,%.12g)\t%g\t%g\r\n",i+1,wr[j],wi[j],res,off );
	}
	sprintf( sPath,"F:\\IEEE\\ORIGIN\\dgeev_%d_out.info",n );
	GSE_outputEV( NEV_select,n,lenda,ev,sPath );

	free( A );			
	free( wr );
	free( select );

	free( ev );
}

/*
	call zgeev to computes the eigenvalues and right eigenvectors 

	7/16/2007		
		v0.1	cys
*/
void GSE_QR_Z( int nCol,int *ptr,int *ind,DoubleComplex *val,int flag,DoubleComplex sft )	{
	int checkShift=0;
	int i,j,k,n,ldv,lwork=-1,info,no_min=-1,*select,*NO,type=1,ilo=1,ihi,mm,m,*ifaill,*ifailr;
	double off,off_min=1.e100,res,*rwork,a_max,a_min;
	DoubleComplex *A,*B,*w,*work,*vl=0x0,*vr=0x0,lenda[NEV_select],*ev,s,*z=0x0;
	char E[1]={'E'},V[1]={'V'},N[1]={'N'},R[1]={'R'},sPath[80];	//only eigenvalue,no schur vector

	n=nCol;	
	if( type==1 )	{
		A=malloc( sizeof(DoubleComplex)*n*n );
		memset( A,0x0,sizeof(DoubleComplex)*n*n );
		for( i = 0; i < n; i++ )	{
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				ASSIGN( A[i*n+ind[j]],val[j] ); 
			}
		}
	}else	{
		ZGE_load( &n,&A,"H:\\GSP\\GSE\\test\\H.dat",0x0 );
	}
	a_max=0.0,		a_min=1.0e100;
	for( i = 0; i < n; i++ )	{
	for( j = 0; j < n; j++ )	{
		ABS( off,A[i*n+j] );
		if( off==0.0 )		continue;
		a_max=MAX( a_max,off );		a_min=MIN( a_min,off );
	}
	}
	w=malloc( sizeof(DoubleComplex)*n );
	ldv=n;
	ev=malloc( sizeof(DoubleComplex)*NEV_select*n );
/*
	lwork=-1;
	work=w;
	rwork = malloc( sizeof(double)*2*n );
	zgeev(N, N,&n,A,&n, w, vl,&ldv, vr, &ldv, work,&lwork,rwork,&info);
	ASSERT( info==0 );
	lwork=(int)(REAL_COMPONENT(work[0]));
	work=malloc( sizeof(DoubleComplex)*(lwork+n*n) );		vr=work+lwork;
	B=malloc( sizeof(DoubleComplex)*n*n );
	memcpy( B,A,sizeof(DoubleComplex)*n*n );
	zgeev(N, V,&n,B,&n, w, vl,&ldv, vr, &ldv, work,&lwork,rwork,&info);
	free(B);
	ASSERT( info==0 );
	free( rwork );
*/	
	lwork=2*n;			ihi=n;			mm = n;
	work=malloc( sizeof(DoubleComplex)*(lwork+n*n) );		vr=work+lwork;
	B=malloc( sizeof(DoubleComplex)*n*n );
	memcpy( B,A,sizeof(DoubleComplex)*n*n );
	zhseqr( E,N,&n,&ilo,&ihi,B,&n,w,z,&n, work,&lwork,&info );
	ASSERT( info==0 );
	rwork=work+n*n;		
	vr=B;
	select=(int*)malloc( sizeof(int)*(3*n) );
	for( i = 0; i < n; i++ )	select[i]=1;
	ifaill=select+n;		ifailr=ifaill+n;
	select[0]=0;		mm=19;
	zhsein ( R,N,N,select,&n,A,&n,w,vl,&n,vr,&n,&mm,&m,work,rwork,ifaill,ifailr,&info );
	ASSERT( info==0 && m==mm );

	select = malloc( sizeof(int)*(n+NEV_select) );	NO=select+n;
	for( i = 0; i < n; i++ )	select[i]=0;
	if( checkShift )	{
		for( i = 0; i < NEV_select; i++ )	{
			no_min=-1;		off_min=1.e100;		
			for( j = 0; j < n; j++ )	{
				if( select[j]!=0 )
					continue;
				OFF( off,w[j],sft );
				if( off < off_min )	
				{	off_min=off;	no_min=j;		}			
			}
			select[no_min]=i+1;
		}
		for( i = 0; i < n; i++ )	{
			if( select[i]<=0 )
				continue;
			NO[select[i]-1]=i;
		}
	}else	{
		for( i = 0; i < NEV_select; i++ )	
		{	select[i]=i;		NO[i]=i;		}
	}
	printf( "--NO--------------------LENDA-----------|Ax-lenda*x|----|lenda-shift|--\r\n" );
	for( i = 0; i < NEV_select; i++ )	{
		j=NO[i];
		OFF( off,w[j],sft );
		ASSIGN( lenda[i],w[j] );
		for( k = 0; k < n; k++ )	
			ASSIGN( ev[i*n+k],vr[j*n+k] );	
		if( type==0 )	{
			res = CS_eigen_resi( nCol,ptr,ind,val,lenda[i],ev+i*n );
		}else	{		//|lenda*x-Ax|
//			BLAS_GEMV( BLAS_NOTRANS,n,n,fuyi_,A,n,ev+i*n,inc_1,zero_,work,inc_1 );	
			BLAS_GEMV( BLAS_N,n,n,fuyi_,A,n,ev+i*n,inc_1,zero_,work,inc_1 );	
			BLAS_AXPY( n,lenda[i],ev+i*n,inc_1,work,inc_1 );
			res = BLAS_NRM2( n,work,inc_1 );	
		}
		printf( "%3d: (%.12lf,%.12lf)\t%g\t%g\r\n",i+1,
			REAL_COMPONENT(w[j]),IMAG_COMPONENT(w[j]),res,off );
	}
	free( A );			free( w );
	free( select );

	sprintf( sPath,"F:\\IEEE\\ORIGIN\\dgeev_%d_out.info",n );
	GSE_outputEV( NEV_select,n,lenda,ev,sPath );
	free( ev );
}

/*
	call zgeev to computes the eigenvalues and right eigenvectors 

	7/9/2007
		v0.1	cys
*/
void GSE_QR_main_1( int nCol,int *ptr,int *ind,double *val,int flag,double sft[2] )	{
	int i,j,n=nCol,ldv=n,info,ilo=1,ihi=n,no_min=-1,
		ldvl=n,mm,m,*ifaill=0x0,ifailr[NEV_select],*select,lwork=2*n;
	GRUS_FLOATING *A=malloc( sizeof(GRUS_FLOATING)*n*(n+2) ),*w=A+n*n,*work=w+n,s,
		*z=0x0,shift,*vl=0x0,*vr,*rwork;
	char E[1]={'E'},N[1]={'N'},R[1]={'R'},S[1]={'S'};	//only eigenvalue,no schur vector
	double off,off_min=1.e100,res;

	memset( A,0x0,sizeof(GRUS_FLOATING)*n*n );
	for( i = 0; i < n; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			SET( A[i*n+ind[j]],val[j],0.0 ); 
		}
	}
	SET( shift,sft[0],sft[1] );
	for( i = 0; i < n; i++ )	{
		DECREMENT( A[i*n+i],shift );
	}

	w=malloc( sizeof(GRUS_FLOATING)*n*(NEV_select+1) );
	vr=w+n;
	work=malloc( sizeof(GRUS_FLOATING)*n*(n+1) );	rwork=work+n*n;			
	zgeev(N, N,&n,A,&n, w, vl,&ldvl, vl,&ldvl, work,&lwork, rwork,&info);
	ASSERT( info==0 );
	for( i = 0; i < NEV_select; i++ )	{
		no_min=-1;		off_min=1.e100;		
		for( j = i; j < n; j++ )	{
			ABS( off,w[j] );
			if( off < off_min )	
			{	off_min=off;	no_min=j;		}			
		}
		ASSIGN( s,w[i] );		ASSIGN( w[i],w[no_min] );		ASSIGN( w[no_min],s );
	}
//	zgees(N, S, select,&n, a,&n, &m, w, vs, ldvs, work,&lwork,rwork, bwork, info)
/*	zhseqr(E,N,&n,&ilo,&ihi,A,&n,w,z,&n, work,&lwork,&info);

	select = malloc( sizeof(int)*n );
	for( i = 0; i < n; i++ )	select[i] = i<NEV_select ? 1 : 0;
	mm=NEV_select;
	zhsein ( R,N,N,select,&n,A,&n,w,vl,&ldvl,
		vr,&n,&mm,&m,work,rwork,ifaill,ifailr,&info );
	ASSERT( info==0 && m<=mm );
	free( select );*/
	free( work );
	free( A );

	printf( "--NO--------------------LENDA-----------|Ax-lenda*x|----|lenda-shift|--\r\n" );
	for( i = 0; i < NEV_select; i++ )	{
		ASSEMBLE( w[i],shift );
		res = CS_eigen_resi( nCol,ptr,ind,val,w[i],vr+i*n );
		SUB( s,w[i],shift );
		ABS( off,s );
		printf( "%3d: (%.12lf,%.12lf)\t%g\t%g\r\n",i+1,REAL_COMPONENT(w[i]),IMAG_COMPONENT(w[i]),res,off );
	}
	free( w );

}

