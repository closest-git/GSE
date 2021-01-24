/*
	Krylov-Schur method

*/
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
	校验AV[:,1:m]=V[:,1:m]*S+V[:,m+1]*b'的偏差

	v0.1	cys
		8/1/2007
*/
int KScuhr_Verify_Resid( GSE_KRYLOV *hKSchur,int K )	{
	int ret=0;
	int i,t=hKSchur->nCurEV,dim=hKSchur->dim,ldh=hKSchur->ldh,no=-1/*,m=hKSchur->V_m*/;
	GRUS_FLOATING *work=hKSchur->work,*vs=work+dim,*H=hKSchur->H,*V=hKSchur->V;
	double norm,max_dot=0.0,*off=GRUS_SE_INFO+GRUS_SE_MAX_OFF_KRYLOV;
	clock_t start=clock( );

	for( i = 0; i < K; i++ )	{
//		BLAS_GEMV( BLAS_NOTRANS,dim,K,one_,V,dim,H+ldh*i,inc_1,zero_,vs,inc_1 );		
		BLAS_GEMV( BLAS_N,dim,K,one_,V,dim,H+ldh*i,inc_1,zero_,vs,inc_1 );		
		GRUS_driver_switch( V+i*dim,work,GRUS_DRIVER_OPX,0 );
		BLAS_SCAL( dim,fuyi_,work,inc_1 );
//		BLAS_AXPY( dim,one_,vs,1,work,1 );
		BLAS_AXPY( dim,one_,vs,inc_1,work,inc_1 );
//		norm=BLAS_NRM2(dim,work,1 );
		norm=BLAS_NRM2(dim,work,inc_1 );
//		BLAS_AXPY( dim,H[ldh*i+K],hKSchur->V+K*dim,1,work,1 );
		BLAS_AXPY( dim,H[ldh*i+K],hKSchur->V+K*dim,inc_1,work,inc_1 );
//		norm=BLAS_NRM2(dim,work,1 );
		norm=BLAS_NRM2(dim,work,inc_1 );
		*off=MAX( *off,norm );
		if( norm > 1.0e-5 )	{		//unreasonable large error
			ret=-(i+1);		goto END;
		}
	}
END:
	GRUS_SE_INFO[GRUS_SE_VERIFY_TIME] += clock( )-start;
	return ret;
}

/*
	Transform B to Schur form. Then reorder

	v0.1	cys
		7/30/2007
*/
int KSchur_reduce( GSE_KRYLOV *hKSchur )	{
	int i,j,ret=GSE_OK,dim=hKSchur->dim,n,m=hKSchur->V_m,t=hKSchur->nCurEV,K,K_extra=0;
	int ldq=hKSchur->ldq,ldh=hKSchur->ldh,info,*select=0x0,*permut,lwork=1,sdim=0,*bwork=0x0,ifst,ilst;
	GRUS_FLOATING *Q=hKSchur->Q,*work=hKSchur->work,*w=hKSchur->w_H+t,beta,*Z,*H,s;
	char N[1]={'N'},V[1]={'V'};
	double *rwork=hKSchur->d_temp,sep,norm;

//	ASSERT(	KScuhr_Verify_Resid( hKSchur,m )==0x0 );
	ASSIGN( beta,hKSchur->H[(m-1)*ldh+m] );
	select = malloc( sizeof(int)*m*2 ),permut=select+m;
	K=hKSchur->nEV-t;
	n=m-t;		//ifst=1;		ilst=n;
	lwork=m*m;
	H = hKSchur->H+t*ldh+t;
	zgees( V,N, select,&n,H,&ldh,&sdim,w,Q,&ldq,work,&lwork,rwork,bwork,&info);
	ASSERT( info==0 && sdim==0 );
/*	for( i = 0; i < K; i++ )	{
		double a,b;
		ifst=i+1;		ilst=i+1;
		ABS( a,H[i*ldh+i] );
		for( j = i+1; j < n; j++ )	{
			ABS( b,H[j*ldh+j] );
			if( b>a )
			{	a=b;	ifst=j+1;	}
		}		
		if( ifst!=ilst )	{
			ztrexc ( V,&n,H,&ldh,Q,&ldq,&ifst,&ilst,&info );
			ASSERT( info==0 );
//			ASSIGN( s,w[ifst-1] );	ASSIGN( w[ifst-1],w[ilst-1] );	ASSIGN( w[ilst-1],s );
		}
		ASSIGN( w[i],H[i*ldh+i] );
	}
	for( i = K; i < n; i++ )	ASSIGN( w[i],H[i*ldh+i] );	*/
	K_extra=GSE_Select( n,w,K,select,permut,zero_,0.00,GSE_ORDER_MAX );
	ASSERT( K_extra>=0 && K_extra<=n-K );
	if( K_extra==n-K )
		return GSE_FAIL_CLUSTER;
	else if( K_extra != 0 )
		{	K += K_extra;	}
	ztrsen( N,V, select,&n,H,&ldh,Q,&ldq,w,&sdim,&s,&sep, work, &lwork,&info);
	ASSERT( info==0 && sdim==K );
	ASSERT( ZV_select_verify(n,w,K,0x0)== 0 );
	free( select );

	ASSERT(	ZU_verify( n,n,Q,ldq,0x0 )==0x0 );

		Z = GRUS_alloc( sizeof(GRUS_FLOATING)*dim*m );
	if( 0 )	{
//		BLAS_GEMM( BLAS_NOTRANS,BLAS_NOTRANS,dim,m,m,one_,hKSchur->V,dim,Q,ldq,zero_,Z,dim );	
		BLAS_GEMM_0( BLAS_N,BLAS_N,dim,m,m,one_,hKSchur->V,dim,Q,ldq,zero_,Z,dim );	
//		BLAS_COPY( dim*m,Z,1,hKSchur->V,1 );
		i = dim*m;		BLAS_COPY( i,Z,inc_1,hKSchur->V,inc_1 );
		for( i = 0; i < m; i++ )			//b'=beta*Q[m,:]
			MULT( hKSchur->H[i*ldh+m],beta,Q[i*ldq+m-1] );
	}else	{
//		BLAS_GEMM( BLAS_NOTRANS,BLAS_NOTRANS,dim,K,n,one_,hKSchur->V+t*dim,dim,Q,ldq,zero_,Z,dim );
		BLAS_GEMM_0( BLAS_N,BLAS_N,dim,K,n,one_,hKSchur->V+t*dim,dim,Q,ldq,zero_,Z,dim );
//		BLAS_COPY( dim*K,Z,1,hKSchur->V+t*dim,1 );
		i = dim*K;		BLAS_COPY( i,Z,inc_1,hKSchur->V+t*dim,inc_1 );
		for( i = 0; i < K; i++ )			//b'=beta*Q[m,:]
			MULT( H[i*ldh+K],beta,Q[i*ldq+n-1] );
		if( t > 0 )	{
			ASSERT( n+t==m );
//			BLAS_GEMM( BLAS_NOTRANS,BLAS_NOTRANS,t,n,n,one_,hKSchur->H+t*ldh,ldh,Q,ldq,zero_,Z,t );
			BLAS_GEMM_0( BLAS_N,BLAS_N,t,n,n,one_,hKSchur->H+t*ldh,ldh,Q,ldq,zero_,Z,t );
			for( i = 0; i < n; i++ )	{		
				BLAS_COPY( t,Z+i*t,inc_1,hKSchur->H+(t+i)*ldh,inc_1 );
			}
		}
	}
		GRUS_free( Z );

	hKSchur->V_p=hKSchur->nEV+K_extra;
	return ret;
}

/*
	Transform B to Schur form. Then reorder

	ISSUE-NEEDED:	
		1、b[i]<tole的处理

	v0.1	cys
		8/9/2007
*/
int KSchur_reduce_2( GSE_KRYLOV *hKSchur,GRUS_FLOATING beta,int t )	{
	int i,ret=GSE_OK,dim=hKSchur->dim,n,m=hKSchur->V_m,K=hKSchur->nEV;
	int ldq=hKSchur->ldq,ldh=hKSchur->ldh,info,*permut,lwork=1,sdim=0,*bwork=0x0;
	GRUS_FLOATING *Q=hKSchur->Q,*work=hKSchur->work,*w,*Z,*H,s;
	char N[1]={'N'},V[1]={'V'};
	double *rwork=hKSchur->d_temp,a;

	Z = GRUS_alloc( sizeof(GRUS_FLOATING)*dim*m );

	if( hKSchur->nCurEV > 0 )	{		//K动态调整 参见ARPACK之naup2.f	"adjust the size of NEV"
		K = K+MIN( hKSchur->nCurEV,(m-hKSchur->nEV)/2 );
	}

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
		i=dim*K;	BLAS_COPY( i,Z,inc_1,hKSchur->AV,inc_1 );
	}
	for( i = 0; i < K; i++ )	{		//b'=beta*Q[m,:]
		MULT( H[i*ldh+K],beta,Q[i*ldq+n-1] );
		ABS( a,H[i*ldh+K] );
		if( i < hKSchur->nCurEV )
		{	CLEAR( H[i*ldh+K] );	}
//		else if( i < hKSchur->nEV )
//			ASSERT( fabs(a)>hKSchur->tole );
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
	ASSERT(	ZU_verify( dim,K+1,hKSchur->V,dim,0x0 )==0x0 );

	GRUS_free( Z );

	i = GRUS_SE_INFO[GRUS_SE_MV_COUNT];
	ARE_trace( "\r\trestart=%d\tConverged=%d\tOP*x=%d",
		hKSchur->nRestart+1,hKSchur->nCurEV,i );

	hKSchur->nRestart++;

	return ret;
}

/*
	v0.1	cys
		8/2/2007
*/
int KSchur_deflate( GSE_KRYLOV *hKSchur )	{
	int i,ret=GSE_OK,t=hKSchur->nCurEV,dim=hKSchur->dim,K;
	int *select,n,ldh=hKSchur->ldh,ldvl=1,ldq=hKSchur->ldq,mm,nSelect,info,need_order=0,ifst,ilst;
	int isConverge=0,isBreak=0;
	GRUS_FLOATING *z,lenda,*H,*vl=0x0,*vr,s,*w,dot;
	double F_norm,norm,res,epsi=hKSchur->tole*hKSchur->a_norm,ratio=1.0e100,res_B=0.0;
	char R[1]={'R'},S[1]={'S'},F[1]={'F'};

	K=hKSchur->nEV;
	n = hKSchur->V_m;			nSelect=K;			mm=n;
	w = hKSchur->w_H;
	vr= hKSchur->ev_H;
	H = hKSchur->H;
	select = malloc( sizeof(int)*n ); 
	for( i = 0; i < n; i++ )	select[i] = i < K ? 1 : 0;
	ztrevc( R,S,select,&n,H,&ldh,vl,&ldvl,vr,&ldh,&mm,&nSelect,hKSchur->work,hKSchur->d_temp,&info);
	ASSERT( info==0 && nSelect==K );
	free( select );
//	ZGE_output( n,vr,ldh,"F:\\GSP\\GSE\\test\\H_0.dat",GSE_OUTPUT_SHORT );

	for( i = hKSchur->V_p; i < hKSchur->V_m;i++ )	{
		DIV( s,hKSchur->w_H[t],hKSchur->w_H[i] );
		ABS( norm,s );
		ratio=MIN( ratio,norm );
	}
//	BLAS_COPY( ldh*n,hKSchur->H,inc_1,hKSchur->work,inc_1 );
	F_norm = LAPACK_LANGE( F,&n,&n,hKSchur->H,&ldh,hKSchur->d_temp );
	for( i = t; i < K; i++ )	{
//		ABS( norm,H[i*ldh+hKSchur->V_p] );
		isConverge=1;
//		norm = BLAS_NRM2( n,vr+i*ldh,1 );							//1、normalize ev_H
		norm = BLAS_NRM2( n,vr+i*ldh,inc_1 );							//1、normalize ev_H
		SET( s,1.0/norm,0.0 );
//		BLAS_SCAL( n,s,vr+i*ldh,1 );
		BLAS_SCAL( n,s,vr+i*ldh,inc_1 );
//		BLAS_DOT_u( n,hKSchur->Q+n-1,ldq,vr+i*ldh,inc_1,dot );		//2、get krylov-resid
		GRUS_DOT( n,hKSchur->Q+n-1,ldq,vr+i*ldh,inc_1,dot,0 );		//2、get krylov-resid
		MULT( s,H[i*ldh+hKSchur->V_p],dot );		ABS( norm,s );
		res_B = norm;	//sqrt(  res_B*res_B+norm*norm );		

		z=hKSchur->rv_A+i*dim;										//3、get ritz vector
//		BLAS_GEMV( BLAS_NOTRANS,dim,n,one_,hKSchur->V,dim,vr+i*ldh,inc_1,zero_,z,inc_1 );		
		BLAS_GEMV( BLAS_N,dim,n,one_,hKSchur->V,dim,vr+i*ldh,inc_1,zero_,z,inc_1 );		
		switch( hKSchur->stop_mode )	{
		case GSE_STOP_KRYLOV:
			ABS( norm,hKSchur->w_H[i] );
			if( res_B >= MAX( mch_eps*F_norm,hKSchur->tole*norm ) )		
				isConverge=0;
			break;//		
		default:
			ASSIGN( lenda,w[i] );
			ARE_trans_lenda( hKSchur,&lenda );	
			res=Eigen_Defect( hKSchur,lenda,z,hKSchur->work );
			if( res < epsi )	{
				if( need_order==1 )	
				{	ifst=i+1;		ilst=i+1;	}				
			}
			else
			{	isConverge=0;	need_order=1;		}
			break;
		}
		if( isConverge==1 ){
			hKSchur->nCurEV++;
			CLEAR( hKSchur->H[i*ldh+hKSchur->V_p] );
		}			
		else	{
			break;
		}
	}
	GRUS_SE_INFO[GRUS_SE_RESTART_RATIO] = ratio;
	res_B=0.0;
	for( i = t; i < K; i++ )	{
		ABS( norm,H[i*ldh+hKSchur->V_p] );
		res_B = sqrt(  res_B*res_B+norm*norm );
	}
	i = GRUS_SE_INFO[GRUS_SE_MV_COUNT];
	ARE_trace( "\r\trestart=%d\tConverged=%d\tOP*x=%d\tratio=%.2g,res_B=%.3g",
		hKSchur->nRestart+1,hKSchur->nCurEV,i,ratio,res_B );

	if( hKSchur->nCurEV==hKSchur->nEV && ZV_select_verify( hKSchur->V_m,w,hKSchur->V_p,0x0 )!=0 )
	{	GRUS_SE_INFO[GRUS_SE_STATUS]=GSE_MISS_ENGENVALUE;		return GSE_MISS_ENGENVALUE;	}

	return ret;
}

/*

	v0.1	cys
		7/30/2007
*/
int KSchur_restart( GSE_KRYLOV *hKSchur )	{
	int i,j,ret=GSE_OK,dim=hKSchur->dim,m=hKSchur->V_m,K;
	double norm;

	hKSchur->nRestart++;

	K=hKSchur->V_p;
//	BLAS_COPY( dim,hKSchur->V+dim*m,1,hKSchur->V+dim*K,1 );
	BLAS_COPY( dim,hKSchur->V+dim*m,inc_1,hKSchur->V+dim*K,inc_1 );
	ASSERT(	ZU_verify( dim,K+1,hKSchur->V,dim,0x0 )==0x0 );

//	ASSERT(	KScuhr_Verify_Resid( hKSchur,K )==0x0 );
/*
	for( i = 0; i < K; i++ )	{		//b'=beta*Q[m,:]
		MULT( hKSchur->H[i*ldh+K],beta,Q[i*ldq+m-1] );
	}*/
	return ret;

}

/*
	对w排序,然后变换H,V	

	v0.1	cys
		8/8/2007
*/
void KSchur_Spectrum_Order( GSE_KRYLOV *hKSchur,int flag  )	{
	int i,K,*ipiv,ifst,ilst,info,m=hKSchur->nGV,ldh=hKSchur->ldh,ldq=ldh;
	char V[1]={'V'};
	GRUS_FLOATING s,*Q=hKSchur->work,*w=hKSchur->w_H;
	double a;

	K = hKSchur->nCurEV;
	ipiv = GRUS_alloc( sizeof(int)*K );
//	ZV_order_x( K,w,zero_,ipiv,flag|GSE_ORDER_SWAP );

	for( i = 0; i < K; i++ )	{
		ASSERT( ipiv[i] >= i );
		if( ipiv[i] == i )		continue;
		ifst=ipiv[i]+1;		ilst=i+1;
		ztrexc ( V,&m,hKSchur->H,&ldh,Q,&(ldq),&ifst,&ilst,&info );
		ASSERT( info==0 );
		OFF(a,w[i],hKSchur->H[i*ldh+i] );
		ASSERT( a== 0.0 );
	}
	GRUS_free( ipiv );
}

/*
	注意：
		code决定流程
	v0.1	cys
		7/30/2007
*/
int KSchur_core( void* hSolver,GRUS_FLOATING **u_0,GRUS_FLOATING **u_1,int *code,int flag )	{
	GSE_KRYLOV *hKSchur=hSolver;
	int i,j,k,dim=hKSchur->dim,ret=GSE_OK,m=hKSchur->nGV,ldh=hKSchur->ldh;
	double norm,res;
	GRUS_FLOATING s;
	clock_t start=clock( );
		
	switch( *code )	{
	case GRUS_DRIVER_START:
		if( hKSchur->init_V0 == 0 )
			GSE_random_1( dim,hKSchur->V );
		break;
	case GRUS_DRIVER_OPX:
		goto OPX_NEXT;
		break;
	default:
		break;
	}

	while( 1 )	{
		hKSchur->V_m = m;
		k = MAX( hKSchur->nCurEV,hKSchur->V_p );
		for( j = k; j < m; j++ )	{
			*u_1 = hKSchur->V+(j+1)*dim;		
			*u_0 = hKSchur->V+j*dim;
			*code=GRUS_DRIVER_OPX;			sprintf( hKSchur->enviro_str,"%d",j );
			goto END;
//			ARE_trans_V( hKSchur,u_0,u_1 );
OPX_NEXT:
			sscanf( hKSchur->enviro_str,"%d",&j );
			GSE_Krylov_update_( hKSchur,j,0x0 );		
			ABS( norm,hKSchur->H[j*ldh+j+1] );		
//			ASSERT(  norm > mch_eps );
			ASSERT(  norm > hKSchur->tole );
		}
		if( KSchur_reduce( hKSchur )!= GSE_OK )
			break;
		KSchur_deflate( hKSchur );
		if( hKSchur->nCurEV == hKSchur->nEV )
			break;
		KSchur_restart( hKSchur );
		if( hKSchur->nRestart>=hKSchur->MAX_RESTART )	
			break;
	}

	for( i = 0; i < hKSchur->nCurEV; i++ )	{
		OFF( norm,hKSchur->w_H[i],hKSchur->H[i*ldh+i] );
		ASSERT( norm == 0.0 );
		ARE_trans_lenda( hKSchur,&(hKSchur->w_H[i]) );
		ASSIGN( (*u_0)[i],hKSchur->w_H[i] );
	}	
	*code=GRUS_DRIVER_FINISH;
END:
	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	return ret;
}


/*
	注意：
		code决定流程
	v0.1	cys
		7/30/2007

int KSchur_core_1( void* hSolver,GRUS_FLOATING **u_0,GRUS_FLOATING **u_1,int *code,int flag )	{
	GSE_KRYLOV *hKSchur=hSolver;
	int i,j,k,dim=hKSchur->dim,ret=GSE_OK,m=hKSchur->nGV,ldh=hKSchur->ldh,K=hKSchur->nEV,t;
	int n,ldq=hKSchur->ldq,info,lwork,*select=hKSchur->i_temp,*permut,isConverge,base;
	double norm,res,F_norm;
	GRUS_FLOATING s,*w=hKSchur->w_H,*vl=0x0,*vr=hKSchur->ev_H,*H,*Q,*work=hKSchur->work,beta;
	clock_t start=clock( );
	char N[1]={'N'},F[1]={'F'},V[1]={'V'};
		
	switch( *code )	{
	case GRUS_DRIVER_START:
		if( hKSchur->init_V0 == 0 )
			GSE_random_1( dim,hKSchur->V );
		break;
	case GRUS_DRIVER_OPX:
		goto OPX_NEXT;
		break;
	default:
		break;
	}

	while( 1 )	{
		hKSchur->V_m = m;
		k = MAX( hKSchur->nCurEV,hKSchur->V_p );
		for( j = k; j < m; j++ )	{
			*u_1 = hKSchur->V+(j+1)*dim;		
			*u_0 = hKSchur->V+j*dim;
			*code=GRUS_DRIVER_OPX;			sprintf( hKSchur->enviro_str,"%d",j );
			goto END;
OPX_NEXT:
			sscanf( hKSchur->enviro_str,"%d",&j );
			GSE_Krylov_update_( hKSchur,j,0x0 );		
			ABS( norm,hKSchur->H[j*ldh+j+1] );		
			ASSERT(  norm > hKSchur->tole );
		}		
		ASSIGN(beta,hKSchur->H[(m-1)*ldh+m] );

		n=m;			t=hKSchur->nCurEV;
		F_norm = LAPACK_LANGE( F,&n,&n,hKSchur->H,&ldh,hKSchur->d_temp );
		H=hKSchur->H;			Q=hKSchur->Q;
		GSE_MG_evs( n,t,H,ldh,w,Q,ldq,vr,ldh,work,hKSchur->d_temp,&info );
		permut = select+n;
		GSE_Select( n-t,w+t,K-t,select,permut+t,zero_,0.00,GSE_ORDER_MAX );
		for( i = 0; i < n; i++ )		select[i]=i<t ? 1 : 0;
		for( i = t; i < K; i++ )	{
			base = permut[i]+t;
			MULT( s,H[(n-1)*ldh+n],hKSchur->ev_H[base*ldh+n-1] );		
			ABS( res,s );
			isConverge=1;
			switch( hKSchur->stop_mode )	{
			case GSE_STOP_KRYLOV:
				ABS( norm,hKSchur->w_H[base] );
				if( res >= MAX( mch_eps*F_norm,hKSchur->tole*norm ) )		
					isConverge=0;
				break;//		
			default:
				if( res < hKSchur->tole )	{
				}
				else
				{	isConverge=0;			}
				break;
			}
			if( isConverge==1 )				{
				select[base] = 1;
				hKSchur->nCurEV++;
			}	
		}	
		info = GSE_MS_sort( n-t,K-t,0,H+t*ldh+t,ldh,Q+t*ldq+t,ldq,w+t,select+t,0x0 );
//		info = GSE_MS_sort( n,K,t,H,ldh,Q,ldq,hKSchur->w_H,select,0x0 );
		ASSERT( info == 0 );
		if( hKSchur->nCurEV == hKSchur->nEV )
			break;
		KSchur_reduce_2( hKSchur,beta,t );
		if( hKSchur->nRestart>=hKSchur->MAX_RESTART )	
			break;
	}

	ZGE_output( n,vr,ldh,"F:\\GSP\\GSE\\test\\H_1.dat",GSE_OUTPUT_SHORT );
	for( i = 0; i < hKSchur->nCurEV; i++ )	{
		OFF( norm,hKSchur->w_H[i],hKSchur->H[i*ldh+i] );
		ASSERT( norm == 0.0 );
		ARE_trans_lenda( hKSchur,&(hKSchur->w_H[i]) );
		ASSIGN( (*u_0)[i],hKSchur->w_H[i] );
	}	
	*code=GRUS_DRIVER_FINISH;
END:
	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	return ret;
}
*/

/*
	v0.1	cys
		8/13/2007
*/
unsigned __stdcall KSchur_ThreadProc(  void* pArguments  )	{
	GSE_KRYLOV *hKSchur=(GSE_KRYLOV *)pArguments;
	int i,j,k,dim=hKSchur->dim,ret=GSE_OK,m=hKSchur->nGV,ldh=hKSchur->ldh,K=hKSchur->nEV,t;
	int n,ldq=hKSchur->ldq,info,lwork,*select=hKSchur->i_temp,*permut,isConverge,base;
	double norm,res,F_norm,*D=hKSchur->D;
	GRUS_FLOATING s,*w=hKSchur->w_H,*vl=0x0,*vr=hKSchur->ev_H,*H,*Q,*work=hKSchur->work,beta;
	clock_t start=clock( );
	char N[1]={'N'},F[1]={'F'},V[1]={'V'};
		
	if( hKSchur->init_V0 == 0 )
		GSE_random_1( dim,hKSchur->V );
	if( hKSchur->balance_mode == GSE_BALANCE_ONESIDE )
//		GSE_Krylov_Balance( hKSchur,3,0x0 );
		GSE_Balance_1( dim,3,hKSchur->D,hKSchur->d_temp,hKSchur->work,0x0 );
	else if( hKSchur->balance_mode == GSE_BALANCE_TWOSIDE )
		GSE_Krylov_Balance_2( hKSchur,3,0x0 );

	while( 1 )	{
		ASSERT( hKSchur->nCurEV<=hKSchur->V_p );
		GSE_Krylov_expand_( hKSchur,hKSchur->V_p,m );
/*		for( j = hKSchur->V_p; j < m; j++ )	{
			if( hKSchur->balance_mode != GSE_BALANCE_NO )
			{	for( i = 0; i< dim; i++ )	SCALE_RECIP( hKSchur->V[j*dim+i],D[i] );	}
			GRUS_driver_switch( hKSchur->V+j*dim,hKSchur->V+(j+1)*dim,GRUS_DRIVER_OPX,0 );
			if( hKSchur->balance_mode != GSE_BALANCE_NO )	{
				for( i = 0; i < dim; i++ )	SCALE_DIV( hKSchur->V[j*dim+i],D[i] );
				for( i = 0; i < dim; i++ )	SCALE_DIV( hKSchur->V[(j+1)*dim+i],D[i] );
			}
			GSE_Krylov_update_( hKSchur,j,0x0 );		
			ABS( norm,hKSchur->H[j*ldh+j+1] );		
			ASSERT(  norm > hKSchur->tole );
		}	*/	
		ASSIGN(beta,hKSchur->H[(m-1)*ldh+m] );
//		ASSERT( KScuhr_Verify_Resid( hKSchur,m )==0  );

		n=m;			t=hKSchur->nCurEV;
		F_norm = LAPACK_LANGE( F,&n,&n,hKSchur->H,&ldh,hKSchur->d_temp );
		H=hKSchur->H;			Q=hKSchur->Q;
		GSE_MG_evs( n,t,H,ldh,w,Q,ldq,vr,ldh,work,hKSchur->d_temp,&info );
		permut = select+n;
		GSE_Select( n-t,w+t,K-t,select,permut+t,zero_,0.00,GSE_ORDER_MAX );
		for( i = 0; i < n; i++ )		select[i]=i<t ? 1 : 0;
		for( i = t; i < K; i++ )	{
			base = permut[i]+t;
			MULT( s,beta,hKSchur->ev_H[base*ldh+n-1] );		
			ABS( res,s );
			isConverge=1;
			switch( hKSchur->stop_mode )	{
			case GSE_STOP_KRYLOV:
				ABS( norm,hKSchur->w_H[base] );
				if( res >= MAX( mch_eps*F_norm,hKSchur->tole*norm ) )		
					isConverge=0;
				break;//		
			default:
				if( res < hKSchur->tole )	{
				}
				else
				{	isConverge=0;			}
				break;
			}
			if( isConverge==1 )				{
				select[base] = 1;
				hKSchur->nCurEV++;
			}	
		}	
		info = GSE_MS_sort( n-t,K-t,0,H+t*ldh+t,ldh,Q+t*ldq+t,ldq,w+t,select+t,0x0 );
//		info = GSE_MS_sort( n,K,t,H,ldh,Q,ldq,hKSchur->w_H,select,0x0 );
		ASSERT( info == 0 );
		KSchur_reduce_2( hKSchur,beta,t );
		if( hKSchur->nCurEV == hKSchur->nEV )
			break;
//		ASSERT( KScuhr_Verify_Resid( hKSchur,K )==0  );
		if( hKSchur->nRestart>=hKSchur->MAX_RESTART )	
			break;
		if( BIT_TEST( hKSchur->sft_mode,GSE_SHIFT_DYNAMIC ) )	{
			GSE_Krylov_DynShift( hKSchur );
		}

	}

//	ZGE_output( n,vr,ldh,"F:\\GSP\\GSE\\test\\H_1.dat",GSE_OUTPUT_SHORT );
	for( i = 0; i < hKSchur->nCurEV; i++ )	{
		OFF( norm,hKSchur->w_H[i],hKSchur->H[i*ldh+i] );
		ASSERT( norm == 0.0 );
		ARE_trans_lenda( hKSchur,&(hKSchur->w_H[i]) );
	}	
	
END:
	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	GRUS_driver_switch( GRUS_NULL,GRUS_NULL,GRUS_DRIVER_FINISH,0 );
	 _endthreadex( 0 );
	return ret;
}

/*
	v0.1	cys
		8/13/2007
*/
unsigned __stdcall KSchur_ThreadProc_1(  void* pArguments  )	{
	GSE_KRYLOV *hKSchur=(GSE_KRYLOV *)pArguments;
	int i,j,k,dim=hKSchur->dim,ret=GSE_OK,m=hKSchur->nGV,ldh=hKSchur->ldh,K=hKSchur->nEV,t;
	int n,ldq=hKSchur->ldq,info,lwork,*select=hKSchur->i_temp,*permut,isConverge,base,nConverg;
	double norm,res,F_norm,*D=hKSchur->D;
	GRUS_FLOATING s,*w=hKSchur->w_H,*vl=0x0,*vr=hKSchur->ev_H,*H,*Q,*work=hKSchur->work,beta,*H_temp;
	clock_t start=clock( );
	char N[1]={'N'},F[1]={'F'},V[1]={'V'};
		
	if( hKSchur->init_V0 == 0 )
		GSE_random_1( dim,hKSchur->V );
	if( hKSchur->balance_mode == GSE_BALANCE_ONESIDE )
//		GSE_Krylov_Balance( hKSchur,3,0x0 );
		GSE_Balance_1( dim,3,hKSchur->D,hKSchur->d_temp,hKSchur->work,0x0 );
	else if( hKSchur->balance_mode == GSE_BALANCE_TWOSIDE )
		GSE_Krylov_Balance_2( hKSchur,3,0x0 );

	H_temp = GRUS_alloc( sizeof(GRUS_FLOATING)*(dim+ldh*ldh) );
	while( 1 )	{
		ASSERT( hKSchur->nCurEV<=hKSchur->V_p );
		for( j = hKSchur->V_p; j < m; j++ )	{
			if( hKSchur->balance_mode != GSE_BALANCE_NO )
			{	for( i = 0; i< dim; i++ )	SCALE_RECIP( hKSchur->V[j*dim+i],D[i] );	}
			GRUS_driver_switch( hKSchur->V+j*dim,hKSchur->V+(j+1)*dim,GRUS_DRIVER_OPX,0 );
			if( hKSchur->balance_mode != GSE_BALANCE_NO )	{
				for( i = 0; i < dim; i++ )	SCALE_DIV( hKSchur->V[j*dim+i],D[i] );
				for( i = 0; i < dim; i++ )	SCALE_DIV( hKSchur->V[(j+1)*dim+i],D[i] );
			}
			GSE_Krylov_update_( hKSchur,j,0x0 );		
			ABS( norm,hKSchur->H[j*ldh+j+1] );		
			ASSERT(  norm > hKSchur->tole );
			if( j < K )		continue;

			ASSIGN(beta,hKSchur->H[j*ldh+j+1] );
			n=j+1;			t=hKSchur->nCurEV;
			F_norm = LAPACK_LANGE( F,&n,&n,hKSchur->H,&ldh,hKSchur->d_temp );
			H=hKSchur->H;			Q=hKSchur->Q;
//			BLAS_COPY( ldh*ldh,hKSchur->H,inc_1,H_temp,inc_1 );			
			i=ldh*ldh;		BLAS_COPY( i,hKSchur->H,inc_1,H_temp,inc_1 );			
			GSE_MG_evs( n,t,H_temp,ldh,w,Q,ldq,vr,ldh,work,hKSchur->d_temp,&info );
			permut = select+n;	
			GSE_Select( n-t,w+t,K-t,select,permut+t,zero_,0.00,GSE_ORDER_MAX );
			nConverg = 0;
			for( i = 0; i < n; i++ )		select[i]=i<t ? 1 : 0;
			for( i = t; i < K; i++ )	{
				base = permut[i]+t;
				MULT( s,beta,hKSchur->ev_H[base*ldh+n-1] );		
				ABS( res,s );
				isConverge=1;
				switch( hKSchur->stop_mode )	{
				case GSE_STOP_KRYLOV:
					ABS( norm,hKSchur->w_H[base] );
					if( res >= MAX( mch_eps*F_norm,hKSchur->tole*norm ) )		
						isConverge=0;
					break;//		
				default:
					if( res < hKSchur->tole )	{
					}
					else
					{	isConverge=0;			}
					break;
				}
				if( isConverge==1 )				{
					select[base] = 1;
					hKSchur->nCurEV++;			nConverg++;
				}	
			}	
			if( nConverg > 0 || j==m-1 )	{
//				BLAS_COPY( ldh*ldh,H_temp,inc_1,hKSchur->H,inc_1 );
				i=ldh*ldh;		BLAS_COPY( i,H_temp,inc_1,hKSchur->H,inc_1 );			
				info = GSE_MS_sort( n-t,K-t,0,H+t*ldh+t,ldh,Q+t*ldq+t,ldq,w+t,select+t,0x0 );
				ASSERT( info == 0 );
				KSchur_reduce_2( hKSchur,beta,t );
			}
			if( hKSchur->nCurEV == hKSchur->nEV )
				break;		
		}		
//		ASSERT( KScuhr_Verify_Resid( hKSchur,m )==0  );


//		ASSERT( KScuhr_Verify_Resid( hKSchur,K )==0  );
		if( hKSchur->nRestart>=hKSchur->MAX_RESTART )	
			break;
	}

//	ZGE_output( n,vr,ldh,"F:\\GSP\\GSE\\test\\H_1.dat",GSE_OUTPUT_SHORT );
	for( i = 0; i < hKSchur->nCurEV; i++ )	{
		OFF( norm,hKSchur->w_H[i],hKSchur->H[i*ldh+i] );
		ASSERT( norm == 0.0 );
		ARE_trans_lenda( hKSchur,&(hKSchur->w_H[i]) );
	}	
	
END:
	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	GRUS_driver_switch( GRUS_NULL,GRUS_NULL,GRUS_DRIVER_FINISH,0 );
	 _endthreadex( 0 );
	return ret;
}