/*
	Residual Krylov method

	Start:	
		9/6/2007		cys

*/
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "..\pch\grus_def.h"
#include "..\pch\grus_version.h"
#include "..\pch\grus_version_lib.h"

#define GRUS_GSE_DLL __declspec(dllexport)
#include "GSE_Krylov.h"
#include "GSE_solver.h"
#include "GSE_util.h"

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static IMPLEMENT_ONE_ZERO( one_,fuyi_,zero_ );	

/*
	z=A*(v1,v2,...vj)*H(:,j)		ÆäÖÐz=V(:,j+1)=A*vj
	1¡¢orthogonalization
	2¡¢update H

	v0.1	cys
		9/6/2007
*/
void GSE_RatioK_expand_( GSE_RATIOK* hRatK,int p,int m  )	{
	int i,j,dim=hRatK->dim,ldh=hRatK->ldh;
	double *D=hRatK->D,norm;
	GRUS_FLOATING s,*V=hRatK->V,*t=hRatK->work;

	hRatK->V_m = m;
	for( j =p; j < m; j++ )	{
		BLAS_GEMV( BLAS_N,dim,m,one_,V,dim,t,inc_1,zero_,V+(j+1)*dim,inc_1 );	//V[:,j+1]=V[:,1:j]*t
		if( hRatK->balance_mode != GSE_BALANCE_NO )
		{	for( i = 0; i< dim; i++ )	SCALE_RECIP( V[j*dim+i],D[i] );	}
		GRUS_driver_switch( V+j*dim,V+(j+1)*dim,GRUS_DRIVER_OPX,0 );
		if( hRatK->balance_mode != GSE_BALANCE_NO )	{
			for( i = 0; i < dim; i++ )	SCALE_DIV( V[j*dim+i],D[i] );
			for( i = 0; i < dim; i++ )	SCALE_DIV( V[(j+1)*dim+i],D[i] );
		}
		if( hRatK->AV != GRUS_NULL )	
			BLAS_COPY( dim,V+(j+1)*dim,inc_1,hRatK->AV+j*dim,inc_1 );
		GSE_Krylov_update_( hRatK,j,0x0 );		
		ABS( norm,hRatK->H[j*ldh+j+1] );	
		ASSERT(  norm > hRatK->tole );
	}
}

/*
	Rational Krylov method.	
	[ref]:	Rational Krylov, a practical algorithm for large sparse nonsymmetric matrix pencils
	
	v0.1	cys
		10/31/2007
*/
unsigned __stdcall RatioK_ThreadProc(  void* pArguments  )	{
	GSE_RATIOK *hRatK=(GSE_KRYLOV *)pArguments;
	int i,j,k,dim=hRatK->dim,ret=GSE_OK,m=hRatK->nGV,ldh=hRatK->ldh,K=hRatK->nEV,t;
	int n,ldq=hRatK->ldq,info,lwork,*select=hRatK->i_temp,*permut,isConverge,base;
	double norm,res,F_norm,*D=hRatK->D;
	GRUS_FLOATING s,*w=hRatK->w_H,*vl=0x0,*vr=hRatK->ev_H,*H,*Q,*work=hRatK->work,beta;
	clock_t start=clock( );
	char N[1]={'N'},F[1]={'F'},V[1]={'V'};
		
	if( hRatK->init_V0 == 0 )
		GSE_random_1( dim,hRatK->V );
	if( hRatK->balance_mode == GSE_BALANCE_ONESIDE )
//		GSE_Krylov_Balance( hRatK,3,0x0 );
		GSE_Balance_1( dim,3,hRatK->D,hRatK->d_temp,hRatK->work,0x0 );
	else if( hRatK->balance_mode == GSE_BALANCE_TWOSIDE )
		GSE_Krylov_Balance_2( hRatK,3,0x0 );

	while( 1 )	{
		ASSERT( hRatK->nCurEV<=hRatK->V_p );
		GSE_RatioK_expand_( hRatK,hRatK->V_p,m );
		n=m;			t=hRatK->nCurEV;
		F_norm = LAPACK_LANGE( F,&n,&n,hRatK->H,&ldh,hRatK->d_temp );
		H=hRatK->H;			Q=hRatK->Q;
//		GSE_MG_evs( n,t,H,ldh,w,Q,ldq,vr,ldh,work,hRatK->d_temp,&info );
		permut = select+n;
		GSE_Select( n-t,w+t,K-t,select,permut+t,zero_,0.00,GSE_ORDER_MAX );
		for( i = 0; i < n; i++ )		select[i]=i<t ? 1 : 0;
		for( i = t; i < K; i++ )	{
			base = permut[i]+t;
			MULT( s,beta,hRatK->ev_H[base*ldh+n-1] );		
			ABS( res,s );
			isConverge=1;
			switch( hRatK->stop_mode )	{
			case GSE_STOP_KRYLOV:
				ABS( norm,hRatK->w_H[base] );
				if( res >= MAX( mch_eps*F_norm,hRatK->tole*norm ) )		
					isConverge=0;
				break;//		
			default:
				if( res < hRatK->tole )	{
				}
				else
				{	isConverge=0;			}
				break;
			}
			if( isConverge==1 )				{
				select[base] = 1;
				hRatK->nCurEV++;
			}	
		}	
		info = GSE_MS_sort( n-t,K-t,0,H+t*ldh+t,ldh,Q+t*ldq+t,ldq,w+t,select+t,0x0 );
//		info = GSE_MS_sort( n,K,t,H,ldh,Q,ldq,hRatK->w_H,select,0x0 );
		ASSERT( info == 0 );
		KSchur_reduce_2( hRatK,beta,t );
		if( hRatK->nCurEV == hRatK->nEV )
			break;
//		ASSERT( KScuhr_Verify_Resid( hRatK,K )==0  );
		if( hRatK->nRestart>=hRatK->MAX_RESTART )	
			break;
	}

//	ZGE_output( n,vr,ldh,"F:\\GSP\\GSE\\test\\H_1.dat",GSE_OUTPUT_SHORT );
	for( i = 0; i < hRatK->nCurEV; i++ )	{
		OFF( norm,hRatK->w_H[i],hRatK->H[i*ldh+i] );
		ASSERT( norm == 0.0 );
		ARE_trans_lenda( hRatK,&(hRatK->w_H[i]) );
	}	
	
END:
	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	GRUS_driver_switch( GRUS_NULL,GRUS_NULL,GRUS_DRIVER_FINISH,0 );
	 _endthreadex( 0 );
	return ret;
}