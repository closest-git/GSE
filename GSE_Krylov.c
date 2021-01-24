/*
	COMMON FUNCTION FOR GSE_KRYLOV_SOLVER

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
#include "GSE_win.h"

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static IMPLEMENT_ONE_ZERO( one_,fuyi_,zero_ );	
//DECLARE_GRUS_TRANS( )
//DECLARE_ONE_ZERO( one_,fuyi_,zero_ )
extern double GRUS_SE_INFO[GRUS_SE_INFO_ITEM];

double dlamch( char* );

/*
	v0.1	cys
		5/2/2007
	v0.2	cys
		8/4/2007		取消(ptr,ind,val)参数

	ISSUE-NEEDED:
		a_norm的确定

	注意
		1、对于m阶arnoldi分解，V,H的rank为m+1
*/
void *GSE_Krylov_Init( int dim,double *param,double *shift )	{
	int dimQH,i,ret,s_mode;
	ARE_R_SOLVER *hAR_R=GRUS_NULL;
	char cmach[1]={'e'};

	GRUS_init_heapinfo( );
	for( i = 0; i < GRUS_SE_INFO_ITEM; i++ )	
		GRUS_SE_INFO[i]=0.0;
	GRUS_SE_INFO[GRUS_SE_MAX_OFF_ORTHONORMAL]=-1.0;
	GRUS_SE_INFO[GRUS_SE_MAX_INVERT_RES]=-1.0;
	GRUS_SE_INFO[GRUS_SE_MAX_OFF_H]=-1.0;
	GRUS_SE_INFO[GRUS_SE_MAX_OFF_KRYLOV]=-1.0;

	hAR_R=GRUS_alloc(sizeof(ARE_R_SOLVER) );

	hAR_R->dim=dim;		hAR_R->nnz=-1;
	hAR_R->A_ptr=NULL;	hAR_R->A_ind=NULL;		hAR_R->A=NULL;		
	hAR_R->ptr=NULL;	hAR_R->ind=NULL;		hAR_R->val=NULL;
	hAR_R->hLU=NULL;

	hAR_R->method=param[GSE_METHOD];
	hAR_R->A_type=param[GSE_DATA_TYPE];
	hAR_R->val_type=GSE_DATA_COMPLEX;		//val暂定为COMPLEX
	hAR_R->nEV=param[GSE_PARAM_NEV];
	hAR_R->nGV=param[GSE_PARAM_NGV];	
	if( hAR_R->nGV > dim )
		hAR_R->nGV = dim;
	if( hAR_R->nEV >= hAR_R->nGV )
		hAR_R->nEV = hAR_R->nGV-1;
	if( hAR_R->nEV<=0 )
		return NULL;

	hAR_R->ivt_mode=param[GSE_INVERT_MODE];	
	hAR_R->sft_mode=param[GSE_SHIFT_MODE];	
	hAR_R->stop_mode=param[GSE_STOP_MODE];	
	hAR_R->restart_mode=param[GSE_RESTART_MODE];	
	hAR_R->MAX_RESTART=param[GSE_PARAM_MAX_RESTART];
	SET( hAR_R->delta,shift[0],shift[1] );
	s_mode = hAR_R->sft_mode&GSE_SHIFT_TYPE;
	switch( s_mode )	{
		case GSE_SHIFT_MA:
			SET( hAR_R->shift,(shift[0]+shift[1]),0.0 );
			break;
		case GSE_SHIFT_NORMAL:
			SET( hAR_R->shift,shift[0],shift[1] );
			break;
		default:
			SET( hAR_R->shift,0.0,0.0 );
			break;
	}
/*---------------internal parameter----------------------*/
	hAR_R->checkOffBase = 0;
	hAR_R->balance_mode = GSE_BALANCE_ONESIDE;	//TWOSIDE;ONESIDE;NO		
	hAR_R->adapt_Vp = 0;					//adaptive set vp for restart
/*---------------internal parameter----------------------*/

	dimQH=hAR_R->ldh=hAR_R->ldq=hAR_R->nGV+1;
	hAR_R->work=GRUS_alloc( sizeof(GRUS_FLOATING)*( (dim+1)*2+dimQH*(dimQH+2) ) );	//see LAPACK_HSEIN
	hAR_R->d_temp=GRUS_alloc( sizeof(double)*dim*7 );	//see LAPACK_HSEIN,LAPACK_HEEVX
	hAR_R->i_temp=GRUS_alloc( sizeof(double)*MAX(dim,dimQH*2) );
	hAR_R->V=GRUS_alloc( sizeof(GRUS_FLOATING)*dim*dimQH );	
	hAR_R->AV=GRUS_NULL;	
	hAR_R->nRV=hAR_R->nEV+1;	
	hAR_R->rv_A=GRUS_alloc( sizeof(GRUS_FLOATING)*dim*(hAR_R->nRV) );		//ritz vector
	hAR_R->H=GRUS_alloc( sizeof(GRUS_FLOATING)*dimQH*dimQH );
	hAR_R->Q=GRUS_alloc( sizeof(GRUS_FLOATING)*dimQH*dimQH );
	hAR_R->ev_H=GRUS_alloc( sizeof(GRUS_FLOATING)*dimQH*dimQH );		//ritz vector
	hAR_R->w_H=GRUS_alloc( sizeof(GRUS_FLOATING)*dimQH );	//real part of eigenvector of H
	if( hAR_R->checkOffBase )
		hAR_R->Base=GRUS_alloc( sizeof(GRUS_FLOATING)*dim );	
	if( hAR_R->balance_mode != GSE_BALANCE_NO )
		hAR_R->D=GRUS_alloc( sizeof(double)*dim );

	memset( hAR_R->H,0x0,sizeof(GRUS_FLOATING)*dimQH*dimQH );

	hAR_R->nRestart = 0;
	hAR_R->nCurEV = 0;
	hAR_R->V_p = 0;
	hAR_R->init_V0 = 0;
	
	hAR_R->a_norm = 0.0;
	mch_eps = dlamch( cmach );
	if( hAR_R->stop_mode==GSE_STOP_DIRECT )
		hAR_R->tole = 100*MAX(param[GSE_CONVERGE_TOL],mch_eps);		//The convergence tolerance
	else
		hAR_R->tole = param[GSE_CONVERGE_TOL];
	hAR_R->ortho_thresh = param[GSE_REORTHO_THRESH];	//sin(param[GSE_REORTHO_THRESH]/180.0*PI);		

//#ifdef	GSE_VIRTUAL_DRIVER
	switch( hAR_R->method )	{
		case GSE_METHOD_KSCHUR:
			ret = GRUS_driver_init( &KSchur_ThreadProc,(void*)hAR_R );
			break;
		case GSE_METHOD_RESK:
			ret = GRUS_driver_init( &ResidK_ThreadProc,(void*)hAR_R );
			hAR_R->AV=GRUS_alloc( sizeof(GRUS_FLOATING)*dim*dimQH );	
			hAR_R->tole = 10*MAX(param[GSE_CONVERGE_TOL],mch_eps);		//适当放大
			break;
		default:
			break;
	}
	if(  ret != GSE_OK )	{
		GSE_Krylov_Clear( hAR_R );
		hAR_R=GRUS_NULL;
	}
//#endif*/
	return (void*)hAR_R;
}

/*
	v0.1	cys
		5/2/2007

	注意
		1、对于m阶arnoldi分解，V,H的rank为m+1
*/
void *GSE_Krylov_Init_ccs( int dim,int *ptr,int *ind,void *val,double *param,double *shift )	{
	ARE_R_SOLVER *hAR_R;
	hAR_R = GSE_Krylov_Init( dim,param,shift );
	if( hAR_R != GRUS_NULL )
	{	hAR_R->A_ptr=ptr;	hAR_R->A_ind=ind;		hAR_R->A=val;	}	
	hAR_R->a_norm = CCS_DLA_F( dim,hAR_R->A_ptr,hAR_R->A_ind,hAR_R->A );
	return (void*)hAR_R;
}


/*
	v0.1	cys
		5/4/2007
*/
void GSE_Krylov_Clear( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;;

	GRUS_driver_clear( );

	GRUS_free( hAR_R->ev_H );
	GRUS_free( hAR_R->work );
	GRUS_free( hAR_R->d_temp );
	GRUS_free( hAR_R->i_temp );
	GRUS_free( hAR_R->V ); 
	GRUS_free( hAR_R->rv_A );
	GRUS_free( hAR_R->H );
	GRUS_free( hAR_R->Q );
	GRUS_free( hAR_R->w_H );
	if( hAR_R->AV!=GRUS_NULL )
		GRUS_free( hAR_R->AV );	
	if( hAR_R->checkOffBase )
		GRUS_free( hAR_R->Base );	

	if( hAR_R->val != hAR_R->A )	{
		GRUS_free( hAR_R->val );
		GRUS_free( hAR_R->ptr );
		GRUS_free( hAR_R->ind );
	}
	if( hAR_R->balance_mode != GSE_BALANCE_NO )
		GRUS_free ( hAR_R->D );

	if( hAR_R->hLU != NULL )	{
		switch( hAR_R->ivt_mode )	{
		case GSE_INVERT_GSS:
			mf_numeric_clear( hAR_R->hLU );
			GRUS_free( hAR_R->hLU );
			break;
		case GSE_INVERT_MF_1:
//			mf_clear_1( hAR_R->hLU );
//			GRUS_free( hAR_R->hLU );
			hAR_R->hLU=NULL;
			break;
		default:
			break;
		}
	}
	GRUS_free( hAR_R );

	GRUS_exit_heapinfo( );

	return ;
}

/*
	采用normwise backward error
	支持返回特征值,特征向量,schur向量等多种功能
	缺省情况下按特征值排序

	v0.1	cys
		5/11/2007
	v0.2	cys
		8/7/2007
*/
int GSE_Krylov_Post( GSE_KRYLOV* hSolver,GRUS_FLOATING *w,GRUS_FLOATING *V,int *nEV,int *ret )	{
	int i,j,dim=hSolver->dim,K,n,ldh=hSolver->ldh,ldq=hSolver->ldq,flag=*ret;
	int *select,*permut,nSelect,ldvl=1,mm,o_flag,cur,s_mode=hSolver->sft_mode&GSE_SHIFT_TYPE;
	double *info=GRUS_SE_INFO,norm,*D=hSolver->D;
	GRUS_FLOATING s,*vl=0x0,*vr=hSolver->ev_H,*z,*work=hSolver->work;
	char F_norm[1]={'F'},R[1]={'R'},S[1]={'S'};		//Frobenius norm

	*nEV = hSolver->nCurEV;

	n=hSolver->V_m;
	K=hSolver->nCurEV;
	select = malloc( sizeof(int)*n*2 );		permut=select+n;

	if( hSolver->a_norm != 0.0 )	{
		ARE_trace( "\r\n\t||A|%g\r\n",hSolver->a_norm );
	}
	if( hSolver->V_p > hSolver->nEV ) 
		ARE_trace( "\tCLUSERT index=%d\r\n",hSolver->V_p-hSolver->nEV );

	if( s_mode != GSE_SHIFT_NO )
	{	ASSIGN(s,hSolver->shift);	o_flag=GSE_ORDER_MIN;	}	
	else	
	{	ASSIGN(s,zero_);			o_flag=GSE_ORDER_MAX;	}
	GSE_Select( K,hSolver->w_H,K,select,permut,s,0.0,o_flag );

	if( BIT_TEST(flag,GSE_EIGEN_VALUE) )	{
		ASSERT( w!=GRUS_NULL );
		for( i = 0; i < K; i ++ )	{
			cur = permut[i];
			ASSIGN( w[i],hSolver->w_H[cur] );
		}
		BIT_RESET(flag,GSE_EIGEN_VALUE);
	}

	if( BIT_TEST(flag,GSE_RIGHT_EIGENVECTOR) )	{
		for( i = 0; i < n; i++ )	select[i] = i < K ? 1 : 0;
		mm=n;
		ztrevc( R,S,select,&n,hSolver->H,&ldh,vl,&ldvl,vr,&ldh,&mm,&nSelect,work,hSolver->d_temp,&info);
		if( info!=0 || nSelect!=K )
			goto END;
		for( i = 0; i < K; i++ )	{
			cur = permut[i];
			norm = BLAS_NRM2( n,vr+cur*ldh,inc_1 );	//normalize ev_H
			SET( s,1.0/norm,0.0 );
//			BLAS_SCAL( n,s,vr+cur*ldh,1 );
			BLAS_SCAL_0( n,s,vr+cur*ldh,inc_1 );
			z=V+i*dim;		
//			BLAS_GEMV( BLAS_NOTRANS,dim,n,one_,hSolver->V,dim,vr+cur*ldh,inc_1,zero_,z,inc_1 );	
			BLAS_GEMV( BLAS_N,dim,n,one_,hSolver->V,dim,vr+cur*ldh,inc_1,zero_,z,inc_1 );	
			if( hSolver->balance_mode != GSE_BALANCE_NO )
				for( j = 0; j < dim; j++ )	SCALE_RECIP( z[j],D[j] );
		}
		BIT_RESET(flag,GSE_RIGHT_EIGENVECTOR);
	}
	if( BIT_TEST(flag,GSE_SCHUR_VECTOR) )	{

	}
END:
	free( select );

	return *ret=flag;
}

/*
	输出一些统计信息
	 
	v0.1	cys
		7/8/2007
*/
void GSE_Krylov_Dump( GSE_KRYLOV* hSolver,int flag )	{
	double *info=GRUS_SE_INFO,s=1.0/CLOCKS_PER_SEC;
	
	ARE_trace( "\tMAX_OFF: " );
	ARE_trace( "KRYLOV=%g(%g),U=%g,SOLVER=%g\r\n",info[GRUS_SE_MAX_OFF_KRYLOV],hSolver->ortho_thresh,
		info[GRUS_SE_MAX_OFF_ORTHONORMAL],info[GRUS_SE_MAX_INVERT_RES] );

	ARE_trace( "\tTime : Total =%g\r\n",info[GRUS_SE_TOTAL_TIME]*s );
	ARE_trace( "\tMV=%g(%g)\tInvert=%g\tUpdate=%g(Reortho=%4.3g%%)\r\n",info[GRUS_SE_MV_TIME]*s,
		info[GRUS_SE_MV_COUNT],info[GRUS_SE_INVERT_TIME]*s,info[GRUS_SE_UPDATE_TIME]*s,
		info[GRUS_SE_REORTHO_TIME]*100.0/info[GRUS_SE_UPDATE_TIME] );
	ARE_trace( "\tRestart Time=%g(%g)\tVerify Time=%g\r\n",info[GRUS_SE_RESTART_TIME]*s,
		info[GRUS_SE_RESTART_COUNT],info[GRUS_SE_VERIFY_TIME]*s );
}

/*
	z=A*(v1,v2,...vj)*H(:,j)		其中z=V(:,j+1)=A*vj
	1、orthogonalization
	2、update H

	v0.1	cys
		9/6/2007
*/
void GSE_Krylov_update_( GSE_KRYLOV* hSolver,int j,int flag  )	{
	GRUS_FLOATING *h,*z,*c,*u,s;
	int dim=hSolver->dim,ldh=hSolver->ldh,n=j+1,i,t;
	double norm,norm_0,norm_1;
	clock_t start=clock( );

	z=hSolver->V+(j+1)*dim;
	norm_0 = BLAS_NRM2( dim,z,inc_1 );
	if( j+1==hSolver->nCurEV || flag==1 )	//如何证明	必相当于deflation过程
		h=hSolver->work+ldh;
	else
		h=hSolver->H+j*ldh;
//	orthogonal:		z = z-V*(V,z)
//	BLAS_GEMV( BLAS_CONJTRANS,dim,n,one_,hSolver->V,dim,z,inc_1,zero_,h,inc_1 );	//h=V'z
	BLAS_GEMV( BLAS_C,dim,n,one_,hSolver->V,dim,z,inc_1,zero_,h,inc_1 );	//h=V'z
//	BLAS_GEMV( BLAS_NOTRANS,dim,n,fuyi_,hSolver->V,dim,h,inc_1,one_,z,inc_1 );	
	BLAS_GEMV( BLAS_N,dim,n,fuyi_,hSolver->V,dim,h,inc_1,one_,z,inc_1 );	
//DGKS correction if needed
	norm = BLAS_NRM2( dim,z,inc_1 );	//normalize z
	if( norm < hSolver->ortho_thresh*norm_0 )	{
		clock_t s_0=clock( );
		norm_0 = norm;
		c = hSolver->work;
//		BLAS_GEMV( BLAS_CONJTRANS,dim,n,one_,hSolver->V,dim,z,inc_1,zero_,c,inc_1 );	//h=-V'z
		BLAS_GEMV( BLAS_C,dim,n,one_,hSolver->V,dim,z,inc_1,zero_,c,inc_1 );	//h=-V'z
//		BLAS_GEMV( BLAS_NOTRANS,dim,n,fuyi_,hSolver->V,dim,c,inc_1,one_,z,inc_1 );	
		BLAS_GEMV( BLAS_N,dim,n,fuyi_,hSolver->V,dim,c,inc_1,one_,z,inc_1 );	
		i=j+1;			BLAS_AXPY( i,one_,c,inc_1,h,inc_1 );	
		norm = BLAS_NRM2( dim,z,inc_1 );	//normalize z
		if( norm < 0.707*norm_0 )	{	//z lies in the span of V_{j} numerically
			printf("\r\n!!!z lies in the span of V_{j} numerically!!!\r\n");
//			for( i = 0; i < dim; i++ )	CLEAR( z[i] );
		}
		GRUS_SE_INFO[GRUS_SE_REORTHO_TIME] += clock( )-s_0;	
	}
	SET(h[j+1],norm,0.0);
	SET(s,1.0/norm,0.0);
//	BLAS_SCAL( dim,s,z,1 );
	BLAS_SCAL_0( dim,s,z,inc_1 );

	GRUS_SE_INFO[GRUS_SE_UPDATE_TIME] += clock( )-start;
}

/*
	z=A*(v1,v2,...vj)*H(:,j)		其中z=V(:,j+1)=A*vj
	1、orthogonalization
	2、update H

	v0.1	cys
		9/6/2007
*/
void GSE_Krylov_expand_( GSE_KRYLOV* hSolver,int p,int m  )	{
	int i,j,dim=hSolver->dim,ldh=hSolver->ldh;
	double *D=hSolver->D,norm;
	GRUS_FLOATING s,*V=hSolver->V;

	hSolver->V_m = m;
	for( j =p; j < m; j++ )	{
		if( hSolver->balance_mode != GSE_BALANCE_NO )
		{	for( i = 0; i< dim; i++ )	SCALE_RECIP( V[j*dim+i],D[i] );	}
		GRUS_driver_switch( V+j*dim,V+(j+1)*dim,GRUS_DRIVER_OPX,0 );
		if( hSolver->balance_mode != GSE_BALANCE_NO )	{
			for( i = 0; i < dim; i++ )	SCALE_DIV( V[j*dim+i],D[i] );
			for( i = 0; i < dim; i++ )	SCALE_DIV( V[(j+1)*dim+i],D[i] );
		}
		if( hSolver->AV != GRUS_NULL )	
			BLAS_COPY( dim,V+(j+1)*dim,inc_1,hSolver->AV+j*dim,inc_1 );
		GSE_Krylov_update_( hSolver,j,0x0 );		
		ABS( norm,hSolver->H[j*ldh+j+1] );	
		ASSERT(  norm > hSolver->tole );
	}
}

/*
	res = Ax-lenda*x;		x=Vy
	必要时返回residual vector

	v0.1	cys
		9/6/2007
*/
double GSE_Krylov_Resid( GSE_KRYLOV *hSolver,int ev_0,int ev_1,GRUS_FLOATING *resi,int flag )	{
	int i,ret=GSE_OK,dim=hSolver->dim,ldh=hSolver->ldh,no=-1,ev,m=hSolver->V_m;
	GRUS_FLOATING *work=hSolver->work,*x=work,*ev_H=hSolver->ev_H,*V=hSolver->V,lenda;
	double norm;

	if( resi==GRUS_NULL )	{
		resi=x+dim;
	}
	ev = ev_0;
//	BLAS_GEMV( BLAS_NOTRANS,dim,m,one_,V,dim,ev_H+ldh*ev,inc_1,zero_,x,inc_1 );		
	BLAS_GEMV( BLAS_N,dim,m,one_,V,dim,ev_H+ldh*ev,inc_1,zero_,x,inc_1 );		
	if( 1 )		{		//for shift-invert residual krylov
		GRUS_driver_switch( x,resi,GRUS_DRIVER_OPX,0 );
		ASSIGN( lenda,hSolver->w_H[ev] );
//		ARE_trans_lenda( hSolver,&(lenda) );
		SCALE_RECIP( lenda,-1 );
//		DIV(lenda,fuyi_,hSolver->w_H[ev] );
		BLAS_AXPY( dim,lenda,x,inc_1,resi,inc_1 );
	}	else	{
		MULT( lenda,fuyi_,hSolver->w_H[ev] );
		if( BIT_TEST(flag,GSE_KRYLOV_ONAV) && hSolver->AV != GRUS_NULL )	
//		{	BLAS_GEMV( BLAS_NOTRANS,dim,m,one_,hSolver->AV,dim,ev_H+ldh*ev,inc_1,zero_,resi,inc_1 );	}
		{	BLAS_GEMV( BLAS_N,dim,m,one_,hSolver->AV,dim,ev_H+ldh*ev,inc_1,zero_,resi,inc_1 );	}
		else
			GRUS_driver_switch( x,resi,GRUS_DRIVER_OPX,0 );
		BLAS_AXPY( dim,lenda,x,inc_1,resi,inc_1 );
	}
	norm=BLAS_NRM2(dim,resi,inc_1 );

	return norm;
}

/*
	返回shift为中心的rank阶子空间的代表向量

	v0.1	cys
		8/6/2007
*/
int GSE_Krylov_getrv( GSE_KRYLOV *hSolver,double *shift,GRUS_FLOATING *z,int rank )	{
	int i,base,ldh=hSolver->ldh,info,m=hSolver->V_m,mm=m,ldvl=1,nSelect;
	int ret=GSE_OK,dim=hSolver->dim,*select,*permut,extra;
	GRUS_FLOATING *w=hSolver->work,*z_0,*z_1,sft,*ev=hSolver->ev_H,s,*vl=0x0;
	double norm;
	char R[1]={'R'},S[1]={'S'};

	SET( sft,shift[0],shift[1] );
/*	BLAS_COPY( m,hSolver->w_H,inc_1,w,inc_1 );
	for( i = 0; i < m; i++ )	{
		DECREMENT( w[i],sft );
	}*/
	select = GRUS_alloc( sizeof(int)*m*2 );		permut=select+m;
	extra = GSE_Select( m,w,rank,select,permut,sft,0.00,GSE_ORDER_MIN );
	ztrevc( R,S,select,&m,hSolver->H,&ldh,vl,&ldvl,ev,&ldh,&mm,&nSelect,hSolver->work,hSolver->d_temp,&info);
	ASSERT( info==0 && nSelect==rank+extra );
	z_0=w+m;			z_1=z_0+dim;
	memset( z_1,0x0,sizeof(GRUS_FLOATING)*dim );
	for( i = 0; i < rank; i++ )	{
		base=permut[i];		ASSERT( select[base]==1 );
//		BLAS_GEMV( BLAS_NOTRANS,dim,m,one_,hSolver->V,dim,ev+i*ldh,inc_1,zero_,z_0,inc_1 );	
		BLAS_GEMV( BLAS_N,dim,m,one_,hSolver->V,dim,ev+i*ldh,inc_1,zero_,z_0,inc_1 );	
		BLAS_AXPY( dim,one_,z_0,inc_1,z_1,inc_1 );		
	}
	norm = BLAS_NRM2( dim,z_1,inc_1 );
	SET(s,1.0/norm,0.0);
//	BLAS_SCAL( dim,s,z_1,1 );
	BLAS_SCAL_0( dim,s,z_1,inc_1 );
	GRUS_free( select );

	BLAS_COPY( dim,z_1,inc_1,z,inc_1 );

	return ret;
}

/*
	flag==0,则统计信息清零

	v0.1	cys
		8/5/2007
*/
void GSE_Krylov_Shift( void *hSol,double *shift,int flag )	{
	GSE_KRYLOV *hSolver=hSol;
	int dimQH=hSolver->ldh,i,isReset=flag==0,s_mode=hSolver->sft_mode&GSE_SHIFT_TYPE;

	if( isReset==0 )	{
		for( i = 0; i < GRUS_SE_INFO_ITEM; i++ )	
			GRUS_SE_INFO[i]=0.0;
		GRUS_SE_INFO[GRUS_SE_MAX_OFF_ORTHONORMAL]=-1.0;
		GRUS_SE_INFO[GRUS_SE_MAX_INVERT_RES]=-1.0;
		GRUS_SE_INFO[GRUS_SE_MAX_OFF_H]=-1.0;
		GRUS_SE_INFO[GRUS_SE_MAX_OFF_KRYLOV]=-1.0;
	}

	hSolver->init_V0 = GSE_Krylov_getrv( hSolver,shift,hSolver->V,hSolver->nEV )==GSE_OK;

	SET( hSolver->delta,shift[0],shift[1] );
	switch( s_mode )	{
		case GSE_SHIFT_MA:
			SET( hSolver->shift,(shift[0]+shift[1]),0.0 );
			break;
		case GSE_SHIFT_NORMAL:
			SET( hSolver->shift,shift[0],shift[1] );
			break;
		default:
			SET( hSolver->shift,0.0,0.0 );
			break;
	}
/*---------------internal parameter----------------------*/
	hSolver->checkOffBase = 0;
/*---------------internal parameter----------------------*/
	memset( hSolver->H,0x0,sizeof(GRUS_FLOATING)*dimQH*dimQH );
	if( isReset )	{
		hSolver->nRestart = 0;
	}
	hSolver->nCurEV = 0;
	hSolver->V_p = 0;

	return ;
}

void dlarnv( int *,int *, int *, double * );
/*
	[REF]:	Balancing Sparse Matrices for Computing Eigenvalues
	v0.1	cys
		8/14/2007
*/
int GSE_Krylov_Balance( GSE_KRYLOV *hSolver,int loop,int flag )	{
	int i,j,ret=GSE_OK,dim=hSolver->dim,idist=2,iseed[4]={1,3,5,7};
	double *x=hSolver->d_temp,a,*D=hSolver->D,d_max=0.0,d_min=1.0e100,norm;
	GRUS_FLOATING *p=hSolver->work,*z=p+dim,s;

	for( j = 0; j < dim; j++ ) D[j]=1.0;
	for( i = 0; i < loop; i++ )	{
		dlarnv( &idist,iseed,&dim,x );
		for( j = 0; j < dim; j++ )	{
			if( x[j]<=0.0 )
			{	SET( z[j],-1.0,0.0 );	}
			else
			{	SET( z[j],1.0,0.0 );	}
		}
		if( i > 0 )	{
			for( j = 0; j < dim; j++ )	{
				SCALE_RECIP( z[j],D[j] );
			}
		}
		GRUS_driver_switch( z,p,GRUS_DRIVER_OPX,0 );
		if( i > 0 )	{
			for( j = 0; j < dim; j++ )	{
				SCALE_DIV( p[j],D[j] );
			}
		}		
		for( j = 0; j < dim; j++ )	{
			ABS( a,p[j] );
			if( a==0.0 )
				a = 1.0;	
		/*	D[j] *= a;*/
			if( i==0 )
			{	D[j]=a;	}
			else
			{	D[j]/=a;	}
		}
	}
	norm=0.0;
	for( j = 0; j < dim; j++ )		{
		norm += D[j]*D[j];
	}
	norm = sqrt( norm );
	for( j = 0; j < dim; j++ )		{
		D[j] /= norm;
		d_max = MAX( d_max,D[j] );
		d_min = MIN( d_min,D[j] );
	}

	return ret;
}

/*
	[REF]:	LAPACK\dgebal.f
	D=sqrt( |q|/|p| )

	v0.1	cys
		8/17/2007
*/
int PQ_Balance_D( int dim,GRUS_FLOATING *p,GRUS_FLOATING *q,double *D,double sfmin1,double thresh ){
	double radix=2.0,sfmax1,sfmin2,sfmax2,r,c,s,f,g,factor=0.95,a;
	int i,isConv=1;

	sfmax1 = 1.0/sfmin1;
	sfmin2 = sfmin1*radix;					sfmax2 = 1.0/sfmin2;	
	for( i = 0; i < dim; i++ )	{
		ABS( r,p[i] );			ABS( c,q[i] );		
		if( r<thresh || c==0.0)			//取a=1.0
			continue;	
		a = sqrt( c/r );
		s = r+c;		g = r/radix;		f=1.0;
		while( !(c>=g || MAX(f,c)>=sfmax2 || MIN(r,g)<=sfmin2) )	{
			f*=radix;		c*=radix;
			r/=radix;		g/=radix;
		}
		g = c / radix;
		while( !(g<r || r>=sfmax2 || MIN( f,MIN(c,g))<=sfmin2) )	{
			f/=radix;		c/=radix;			g/=radix;
			r*=radix;		
		}

		if( c+r >= factor*s )
			continue;
		if( f < 1.0 && D[i] < 1.0 && f*D[i] < sfmin1 )
			continue;
		if( f > 1.0 && D[i] > 1.0 && D[i] > sfmax1/f )
			continue;
		a = D[i]*f;
//		if( a > 1.0e6 || a < 1.0e-6 )
//			continue;
		D[i] = a;
	//	g=1.0/f;
	//	D[i] *= g;

		isConv = 0;
	}

	return isConv;
}

/*
	[REF]:	Balancing Sparse Matrices for Computing Eigenvalues

	v0.1	cys
		8/17/2007
*/
int GSE_Krylov_Balance_2( GSE_KRYLOV *hSolver,int loop,int flag )	{
	int i,j,ret=GSE_OK,dim=hSolver->dim,idist=2,iseed[4]={1,3,5,7},checkErr=0,isConv;
	double *x=hSolver->d_temp,a,b,*D=hSolver->D,d_max=0.0,d_min=1.0e100,norm,thresh,off=0.0,sfmin1;
	GRUS_FLOATING *q,*p,*z=hSolver->work,*ax=z+dim,s;
	char S[1]={'S'},P[1]={'P'};

	sfmin1 = dlamch( S ) / dlamch( P );		

	p=GRUS_alloc( sizeof(GRUS_FLOATING)*dim*2 );		q=p+dim;
	for( j = 0; j < dim; j++ )				D[j]=1.0;
	thresh = 0.0;
	for( i = 0; i < loop; i++ )	{
		for( j = 0; j < 3; j++ )	iseed[i] *=(i+1);
		dlarnv( &idist,iseed,&dim,x );
		for( j = 0; j < dim; j++ )	{
			if( x[j]<=0.0 )
			{	SET( z[j],-1.0,0.0 );	}
			else
			{	SET( z[j],1.0,0.0 );	}
		}
		memcpy( ax,z,sizeof(GRUS_FLOATING)*dim );
		if( i > 0 )	{
			for( j = 0; j < dim; j++ )	{
				SCALE_RECIP( z[j],D[j] );
			}
		}
		GRUS_driver_switch( z,p,GRUS_DRIVER_OPX,0 );
		if( checkErr==1 )	{
			GRUS_driver_switch( p,ax,GRUS_DRIVER_AX,0 );
			BLAS_AXPY( dim,fuyi_,z,inc_1,ax,inc_1 );
			norm = BLAS_NRM2( dim,ax,inc_1 );
			off=MAX( off,norm );
		}		
		if( i > 0 )	{
			for( j = 0; j < dim; j++ )	{
				SCALE_DIV( p[j],D[j] );
			}
		}	

//to get q
		memcpy( z,ax,sizeof(GRUS_FLOATING)*dim );
		if( i > 0 )	{
			for( j = 0; j < dim; j++ )	{
				SCALE_DIV( z[j],D[j] );
			}
		}		
		GRUS_driver_switch( z,q,GRUS_DRIVER_OPTX,0 );
		if( checkErr==1 )	{
			GRUS_driver_switch( q,ax,GRUS_DRIVER_ATX,0 );
			BLAS_AXPY( dim,fuyi_,z,inc_1,ax,inc_1 );
			norm = BLAS_NRM2( dim,ax,inc_1 );
			off=MAX( off,norm );
		}		
		if( i > 0 )	{
			for( j = 0; j < dim; j++ )	{
				SCALE_RECIP( q[j],D[j] );
			}
		}	
		
		isConv = PQ_Balance_D( dim,p,q,D,sfmin1,thresh );
		if( isConv )
			break;
	}
	ASSERT( off < 1.0e-6 );
	GRUS_SE_INFO[GRUS_SE_MAX_INVERT_RES]=MAX( GRUS_SE_INFO[GRUS_SE_MAX_INVERT_RES],off );
	norm=0.0;
	for( j = 0; j < dim; j++ )		{
		norm += D[j]*D[j];
	}
	norm = sqrt( norm );
	for( j = 0; j < dim; j++ )		{
//		D[j] /= norm;
		d_max = MAX( d_max,D[j] );
		d_min = MIN( d_min,D[j] );
	}
	GRUS_free( p );

	return ret;
}

/*
	v0.1	cys
		11/1/2007
*/
void GSE_Krylov_DynShift( GSE_KRYLOV* hSolver )	{
	int i,m=hSolver->V_m;
	GRUS_FLOATING *work=hSolver->work,*w=hSolver->w_H,s,s_0;

/*	CLEAR( s );
	for( i = 0; i < m; i++ )	{
		ASSIGN( s_0,w[i] );
		ARE_trans_lenda( hSolver,&s_0 );
		ASSEMBLE( s,s_0 );
	}
	SCALE_DIV( s,m );*/

	ASSIGN( s,w[0] );
	ARE_trans_lenda( hSolver,&s );
	ASSIGN( work[0],s );
	GRUS_driver_switch( work,GRUS_NULL,GRUS_DRIVER_SHIFT_INVERT,0 );

	BIT_RESET( hSolver->sft_mode,GSE_SHIFT_DYNAMIC );
}