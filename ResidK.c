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
//DECLARE_GRUS_TRANS( )
//DECLARE_ONE_ZERO( one_,fuyi_,zero_ )

/*
	校验U'AU=H的偏差

	v0.1	cys
		8/1/2007
*/
int ResidK_Verify_Resid( GSE_KRYLOV *hResK )	{
	int i,j,k,dim=hResK->dim,ldh=hResK->ldh,m=hResK->V_m,ret=0;
	GRUS_FLOATING *work=hResK->work,*H=hResK->H,*V=hResK->V,*AV=hResK->AV;
	double off,off_max=0.0;
	clock_t start=clock( );

//	BLAS_GEMM( BLAS_CONJTRANS,BLAS_NOTRANS,m,m,dim,one_,V,dim,AV,dim,zero_,work,ldh );	
	BLAS_GEMM_0( BLAS_C,BLAS_N,m,m,dim,one_,V,dim,AV,dim,zero_,work,ldh );	
	for( i = 0; i < m; i++ )	{
	for( j = 0; j < m; j++ )	{
		k = i*ldh+j;
		OFF( off,work[k],H[k] );
		off_max = MAX( off_max,off );
		if( off > 1.0e-5 )	{		//unreasonable large error
			ret=-(i+1);		goto END;
		}
	}
	}

END:
	GRUS_SE_INFO[GRUS_SE_MAX_OFF_KRYLOV]=MAX( GRUS_SE_INFO[GRUS_SE_MAX_OFF_KRYLOV],off_max );
	GRUS_SE_INFO[GRUS_SE_VERIFY_TIME] += clock( )-start;
	return ret;
}

/*
	v0.1	cys
		9/6/2007
*/
void ResidK_Expand_( GSE_KRYLOV *hResK,GRUS_FLOATING *resi ){ 
	int i,m=hResK->V_m,dim=hResK->dim,ldh=hResK->ldh,k;
	double *D=hResK->D,norm;
	GRUS_FLOATING *V=hResK->V,*AV=hResK->AV,*H=hResK->H;

//	ASSERT( ResidK_Verify_Resid(hResK)==0 );
	norm = BLAS_NRM2( dim,resi,inc_1 );
//	BLAS_COPY( dim,resi,inc_1,V+m*dim,1 );
	BLAS_COPY( dim,resi,inc_1,V+m*dim,inc_1 );
	GSE_Krylov_update_( hResK,m-1,1 );		//其中对H的update是多余的
	if( hResK->balance_mode != GSE_BALANCE_NO )
	{	for( i = 0; i< dim; i++ )	SCALE_RECIP( V[m*dim+i],D[i] );	}
	GRUS_driver_switch( V+m*dim,AV+m*dim,GRUS_DRIVER_OPX,0 );
	if( hResK->balance_mode != GSE_BALANCE_NO )	{
		for( i = 0; i < dim; i++ )	SCALE_DIV( V[m*dim+i],D[i] );
		for( i = 0; i < dim; i++ )	SCALE_DIV( AV[m*dim+i],D[i] );
	}	
//Expand H
//	BLAS_GEMV( BLAS_CONJTRANS,dim,m,one_,V,dim,AV+m*dim,inc_1,zero_,H+ldh*m,inc_1 );	
	BLAS_GEMV( BLAS_C,dim,m,one_,V,dim,AV+m*dim,inc_1,zero_,H+ldh*m,inc_1 );	
	k=1;
//	BLAS_GEMM( BLAS_CONJTRANS,BLAS_NOTRANS,k,m,dim,one_,V+m*dim,dim,AV,dim,zero_,H+m,ldh );		
	BLAS_GEMM_0( BLAS_C,BLAS_N,k,m,dim,one_,V+m*dim,dim,AV,dim,zero_,H+m,ldh );		
//	BLAS_DOT_c(dim,V+m*dim,inc_1,AV+m*dim,inc_1,H[ldh*m+m] );
	cblas_zdotc_sub(dim,V+m*dim,inc_1,AV+m*dim,inc_1,H+ldh*m+m );
	ABS( norm,hResK->H[ldh*m+m] );	
//	GRUS_DOT(dim,V+m*dim,inc_1,AV+m*dim,inc_1,H+ldh*m+m,1 );
//	ABS( norm,hResK->H[ldh*m+m] );	
//	ASSERT(  norm > hResK->tole );
	hResK->V_m++;
	ASSERT( ResidK_Verify_Resid(hResK)==0 );
}

/*
	v0.1	cys
		9/7/2007
*/
int  ResidK_restart( GSE_KRYLOV *hKSchur,int t )	{
	int i,ret=GSE_OK,dim=hKSchur->dim,n,m=hKSchur->V_m,K=hKSchur->nEV;
	int ldq=hKSchur->ldq,ldh=hKSchur->ldh,info,*permut,lwork=1,sdim=0,*bwork=0x0;
	GRUS_FLOATING *Q=hKSchur->Q,*work=hKSchur->work,*w,*Z,*H,s;
	char N[1]={'N'},V[1]={'V'};
	double *rwork=hKSchur->d_temp,a;

	Z = GRUS_alloc( sizeof(GRUS_FLOATING)*dim*m );

	n=m;
	ASSERT(	ZU_verify( n,n,Q,ldq,0x0 )==0x0 );
	H = hKSchur->H;		Q=hKSchur->Q;
//	BLAS_GEMM( BLAS_NOTRANS,BLAS_NOTRANS,dim,K,n,one_,hKSchur->V,dim,Q,ldq,zero_,Z,dim );
	BLAS_GEMM_0( BLAS_N,BLAS_N,dim,K,n,one_,hKSchur->V,dim,Q,ldq,zero_,Z,dim );
//	BLAS_COPY( dim*K,Z,1,hKSchur->V,1 );
	i=dim*K;		BLAS_COPY( i,Z,inc_1,hKSchur->V,inc_1 );
	if( hKSchur->AV!=GRUS_NULL )	{
//		BLAS_GEMM( BLAS_NOTRANS,BLAS_NOTRANS,dim,K,n,one_,hKSchur->AV,dim,Q,ldq,zero_,Z,dim );
		BLAS_GEMM_0( BLAS_N,BLAS_N,dim,K,n,one_,hKSchur->AV,dim,Q,ldq,zero_,Z,dim );
//		BLAS_COPY( dim*K,Z,1,hKSchur->AV,1 );
		i=dim*K;		BLAS_COPY( i,Z,inc_1,hKSchur->AV,inc_1 );
	}
	for( i = 0; i < K; i++ )	{		//b'=beta*Q[m,:]
		if( i < hKSchur->nCurEV )
		{	CLEAR( H[i*ldh+K] );	}
	}
	if( t > 0 )	{
		n -= t;
//		BLAS_GEMM( BLAS_NOTRANS,BLAS_NOTRANS,t,n,n,one_,hKSchur->H+t*ldh,ldh,Q+t*ldh+t,ldq,zero_,Z,t );
		BLAS_GEMM_0( BLAS_N,BLAS_N,t,n,n,one_,hKSchur->H+t*ldh,ldh,Q+t*ldh+t,ldq,zero_,Z,t );
		for( i = 0; i < n; i++ )	{		
			BLAS_COPY( t,Z+i*t,inc_1,hKSchur->H+(t+i)*ldh,inc_1 );
		}
		n += t;
	}
	hKSchur->V_p=K;
//	BLAS_COPY( dim,hKSchur->V+dim*m,1,hKSchur->V+dim*K,1 );
	BLAS_COPY( dim,hKSchur->V+dim*m,inc_1,hKSchur->V+dim*K,inc_1 );
	ASSERT(	ZU_verify( dim,K,hKSchur->V,dim,0x0 )==0x0 );

	GRUS_free( Z );

	i = GRUS_SE_INFO[GRUS_SE_MV_COUNT];
	ARE_trace( "\r\trestart=%d\tConverged=%d\tOP*x=%d",
		hKSchur->nRestart+1,hKSchur->nCurEV,i );

	hKSchur->nRestart++;

	return ret;
}

/*
	[REF]:	Improving Eigenvectors in Arnoldi’s Method
	ISSUE-NEEDED:
		how to get beta?

	v0.1	cys
		10/30/2007
*/
void ResidK_Refined_RV( GSE_KRYLOV *hSolver,int ev,GRUS_FLOATING *resi,GRUS_FLOATING *rrv,double res_2 )	{
	int i,ret=GSE_OK,dim=hSolver->dim,ldh=hSolver->ldh,no=-1,m=hSolver->V_m;
	GRUS_FLOATING *work=hSolver->work,*x=work,*ev_H=hSolver->ev_H,*V=hSolver->V,lenda,alpha,beta,eta,s;
	double norm,ksi,gama,a,eta_2;

	GRUS_driver_switch( V+dim*m,x,GRUS_DRIVER_OPX,0 );
	ASSIGN( alpha,hSolver->w_H[ev] );		SCALE_RECIP( alpha,-1 );	
	BLAS_AXPY( dim,alpha,V+dim*m,inc_1,x,inc_1 );		//x=(A-lenda*I)*V[:,m+1]
	ksi=BLAS_NRM2( dim,x,inc_1 );		
	ksi = ksi*ksi;
	GRUS_DOT( dim,V+dim*m,inc_1,x,inc_1,&eta,1 );		//eta=|V[:,m+1]'*(A-lenda*I)*V[:,m+1]|
	ABS( eta_2,eta );			eta_2=eta_2*eta_2;
	gama=sqrt( (res_2-ksi)*(res_2-ksi)+4.0*res_2*eta_2 );
	ASSERT( ksi-gama < res_2 );
	a = (ksi-res_2+gama)/sqrt( (ksi-res_2+gama)*(ksi-res_2+gama)+4*res_2*eta_2 );
	SET( alpha,a,0.0 );
	a = 2.0*a/(res_2-ksi-gama);
	MULT( s,hSolver->H[(m-1)*ldh+m],ev_H[m-1] );
	MULT_CONJ( beta,s,eta );
	SCALE_RECIP( beta,a );
//rrv = alpha*V*y+beta*V[:,m+1]
	for( i = 0; i < dim; i++ )	CLEAR( rrv[i] );
	BLAS_GEMV( BLAS_N,dim,m,one_,V,dim,ev_H+ldh*ev,inc_1,zero_,x,inc_1 );	//x=V*y	
	BLAS_AXPY( dim,alpha,x,inc_1,rrv,inc_1 );		
	BLAS_AXPY( dim,beta,V+dim*m,inc_1,rrv ,inc_1 );	

//verify
	ABS( a,alpha );					ABS( norm,beta );
//	ASSERT( fabs(a*a+norm*norm-1.0)<DBL_EPSILON );
	DIV( s,beta,alpha );			ABS( a,s );
	a -= 2*sqrt(res_2*eta_2)/(ksi-res_2+gama);
//	ASSERT( fabs(a)<DBL_EPSILON );

	return ;
}

/*
	Residual Krylov method.	
	[ref]:	Residual Arnoldi Methods Theory, Package, and Experiments
	
	v0.1	cys
		9/6/2007
*/
unsigned __stdcall ResidK_ThreadProc(  void* pArguments  )	{
	GSE_KRYLOV *hResK=(GSE_KRYLOV *)pArguments;
	int i,t,dim=hResK->dim,ret=GSE_OK,m=hResK->nGV,ldh=hResK->ldh,K=hResK->nEV,res_no,nConverg=0;
	int n,ldq=hResK->ldq,info,*select=hResK->i_temp,*permut,isConverge,base,ldh_2=ldh*ldh;
	double norm,res,F_norm,*D=hResK->D;
	GRUS_FLOATING s,*w=hResK->w_H,*vl=0x0,*vr=hResK->ev_H,*H,*Q,*work=hResK->work,beta,*resi,*H_temp;
	clock_t start=clock( );
	char N[1]={'N'},F[1]={'F'},V[1]={'V'};
		
	if( hResK->init_V0 == 0 )
		GSE_random_1( dim,hResK->V );
	if( hResK->balance_mode == GSE_BALANCE_ONESIDE )
//		GSE_Krylov_Balance( hResK,3,0x0 );
		GSE_Balance_1( dim,3,hResK->D,hResK->d_temp,hResK->work,0x0 );
	else if( hResK->balance_mode == GSE_BALANCE_TWOSIDE )
		GSE_Krylov_Balance_2( hResK,3,0x0 );

	resi = GRUS_alloc( sizeof(GRUS_FLOATING)*(dim+ldh*ldh) );
	H_temp = resi+dim;
	GSE_Krylov_expand_( hResK,0,m );
	H=hResK->H;		Q=hResK->Q;

	while( 1 )	{
//		ASSERT( ResidK_Verify_Resid(hResK)==0 );
		n=hResK->V_m;			t=hResK->nCurEV;
//		ASSIGN(beta,hResK->H[(m-1)*ldh+m] );
		F_norm = LAPACK_LANGE( F,&n,&n,hResK->H,&ldh,hResK->d_temp );
//		BLAS_COPY( ldh*ldh,hResK->H,inc_1,H_temp,inc_1 );			
		BLAS_COPY( ldh_2,hResK->H,inc_1,H_temp,inc_1 );			
		GSE_MG_evs( n,t,H_temp,ldh,w,Q,ldq,vr,ldh,work,hResK->d_temp,&info );
		permut = select+n;
		GSE_Select( n-t,w+t,K-t,select,permut+t,zero_,0.00,GSE_ORDER_MAX );
		for( i = 0; i < n; i++ )		select[i]=i<t ? 1 : 0;
		res_no=-1;				nConverg=0;
		for( i = t; i < K; i++ )	{
			base = permut[i]+t;
//			MULT( s,beta,hResK->ev_H[base*ldh+n-1] );		
//			ABS( res,s );
			res = GSE_Krylov_Resid( hResK,base,base,res_no==-1 ?resi:GRUS_NULL,GSE_KRYLOV_ONAV );		
			isConverge=1;
			switch( hResK->stop_mode )	{
			case GSE_STOP_KRYLOV:
				ABS( norm,hResK->w_H[base] );
				if( res >= MAX( mch_eps*F_norm,hResK->tole*norm ) )		
					isConverge=0;
				break;//		
			default:
				break;
			}
			if( isConverge==1 )				{
				select[base] = 1;
				hResK->nCurEV++;		nConverg++;
			}else if( res_no==-1 )	{
				res_no=i;
				ResidK_Refined_RV( hResK,base,resi,resi,res*res );
			}
		}	
//		if( res_no != -1 && hResK->nCurEV!=hResK->nEV )
//			GSE_Krylov_Resid( hResK,res_no,res_no,resi,GSE_KRYLOV_ONAV );		
		if( hResK->V_m == hResK->nGV || nConverg>0 )	{
//			BLAS_COPY( ldh*ldh,H_temp,inc_1,hResK->H,inc_1 );			
			BLAS_COPY( ldh_2,H_temp,inc_1,hResK->H,inc_1 );			
			info = GSE_MS_sort( n-t,K-t,0,H+t*ldh+t,ldh,Q+t*ldq+t,ldq,w+t,select+t,0x0 );
			ASSERT( info == 0 );
			ResidK_restart( hResK,t );
			hResK->V_m = hResK->V_p;
		}
		if( hResK->nCurEV == hResK->nEV )
			break;
		ResidK_Expand_( hResK,resi );
		if( hResK->nRestart>=hResK->MAX_RESTART )	
			break;
	}

//	ZGE_output( n,vr,ldh,"F:\\GSP\\GSE\\test\\H_1.dat",GSE_OUTPUT_SHORT );
	for( i = 0; i < hResK->nCurEV; i++ )	{
		OFF( norm,hResK->w_H[i],hResK->H[i*ldh+i] );
		ASSERT( norm == 0.0 );
		ARE_trans_lenda( hResK,&(hResK->w_H[i]) );
	}	

	GRUS_free( resi );
	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	GRUS_driver_switch( GRUS_NULL,GRUS_NULL,GRUS_DRIVER_FINISH,0 );
	 _endthreadex( 0 );
	return ret;
}