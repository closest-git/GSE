#pragma comment( user, "Compiled on " __DATE__ " at " __TIME__ ) 

//import GSS lib
#ifdef _USE_GSS_DLL_6_
		#pragma comment( lib, "GSS_DLL_6" )
		extern "C" void GSS_demo_6( int dim,int *ptr,int *ind,double *val );
#else
	#ifdef _DEBUG
		#pragma comment( lib, "GSS_LIB_ZD" )
	#else
		#pragma comment( lib, "GSS_LIB_ZD.lib" )
	#endif
#endif

//import BLAS lib
#ifdef _USE_GOTO_BLAS_
	#pragma comment( lib, "goto_p4" )
#else
	#ifdef _USE_ATLAS_BLAS_
	//	#pragma comment( lib, "atlas_P4" )
		#pragma comment( lib, "cblas" )
		#pragma comment( lib, "atlas" )
		#pragma comment( lib, "lapack" )
	#else
		#pragma comment( lib, "mkl_c" )
		#pragma comment( lib, "libguide" )
		void MKL_FreeBuffers(void);
	#endif
#endif

/*import METIS lib
*/
#pragma comment( lib, "mkl_solver" )

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
//IMPLEMENT_GRUS_TRANS( )
//IMPLEMENT_ONE_ZERO( one_,fuyi_,zero_ )

double GRUS_SE_INFO[GRUS_SE_INFO_ITEM];

double mch_eps;

/*
	generate random vector with norm-1

	v0.1	cys
		5/25/2007
*/
void	GSE_random_1( int dim,GRUS_FLOATING *V )	{
	int iseed[4]={1,3,5,7},idist=4,i;
	double norm;
	GRUS_FLOATING s;

	LAPACK_LARNV( &idist,iseed,&dim,V );
//	for( i = 0; i < dim; i++ )	{	SET(V[i],1,0.0);	}
	if( 0 ) {		//failed test
		int i;
		double *x=GRUS_alloc( sizeof(double)*dim );
		dlarnv( &idist,iseed,&dim,x );
		for( i = 0; i < dim; i++ )	{
			if( x[i]<=0.0 )
			{	SET( V[i],-1.0,0.0 );	}
			else
			{	SET( V[i],1.0,0.0 );	}
		}
		GRUS_driver_switch( V,V,GRUS_DRIVER_OPX,0 );
		for( i = 0; i < dim; i++ )
			RECIPROCAL( V[i] );
		GRUS_free( x );
	}
	norm = BLAS_NRM2( dim,V,inc_1 );
	SET( s,1.0/norm,0.0 );
//	BLAS_SCAL( dim,s,V,inc_1 );
	BLAS_SCAL_0( dim,s,V,inc_1 );
}

/*
	z=(A-shift*I)*v		A is the matrix before transformation

	注意：A与val的区别

	v0.1	cys
		5/12/2007
*/
void ARE_AV_0( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v,GRUS_FLOATING *z )	{
	GRUS_FLOATING s,t,*val;
	int i,j,dim=hAR_R->dim,*ptr=hAR_R->A_ptr,*ind=hAR_R->A_ind;
	double *d_val,*x=hAR_R->d_temp;

	for( i = 0; i < dim; i++ )	CLEAR(z[i]);
	
	if( hAR_R->A_type==GSE_DATA_DOUBLE )	{	//double version
		d_val = hAR_R->A;
		for( i = 0; i < dim; i++ )	{
			ASSIGN( s,v[i] );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				SET( t,d_val[j],0.0 );
				MULT_ADD( z[ind[j]],t,s );	
			}
		}
	}else	{									//complex version
		ASSERT( hAR_R->A_type==GSE_DATA_COMPLEX );
		val = hAR_R->A;
		for( i = 0; i < dim; i++ )	{
			ASSIGN( s,v[i] );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				MULT_ADD( z[ind[j]],val[j],s );	
			}
		}
	}	
}

/*
	z=(A-shift*I)*a*v+z

	v0.1	cys
		5/25/2007
*/
void ARE_AVPY(  ARE_R_SOLVER *hAR_R,GRUS_FLOATING a,GRUS_FLOATING *v,GRUS_FLOATING *z )	{	
	GRUS_FLOATING s,t,*val;
	int i,j,dim=hAR_R->dim,*ptr=hAR_R->A_ptr,*ind=hAR_R->A_ind;
	double *d_val;

//	MULT( s,a,hAR_R->shift );
//	SCALE_RECIP( s,-1.0 );
//	BLAS_AXPY( dim,s,v,1,z,1 );	//
	
	if( hAR_R->A_type==GSE_DATA_DOUBLE )	{	//double version
		d_val = hAR_R->A;
		for( i = 0; i < dim; i++ )	{
//			ASSIGN( s,v[i] );
			MULT( s,a,v[i] );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				SET( t,d_val[j],0.0 );
				MULT_ADD( z[ind[j]],t,s );	
			}
			MULT_SUB( z[i],hAR_R->shift,s );	
		}
	}else	{									//complex version
		ASSERT( hAR_R->A_type==GSE_DATA_COMPLEX );
		val = hAR_R->A;
		for( i = 0; i < dim; i++ )	{
			ASSIGN( s,v[i] );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				MULT_ADD( z[ind[j]],val[j],s );	
			}
		}
	}	
}

/*
	res=|lenda*e-Ae|	
	注意：r=A*e-lenda*e

	v0.1	cys
		7/18/2007
*/
double Eigen_Defect( ARE_R_SOLVER *hAR_R,GRUS_FLOATING lenda,GRUS_FLOATING *ev,GRUS_FLOATING *r )	{
	int dim=hAR_R->dim;
	double res;
	GRUS_FLOATING s;

	ASSIGN( s,lenda );		SCALE_RECIP( s,-1.0 );	
	ARE_AV_0( hAR_R,ev,r );						
//	BLAS_SCAL( dim,fuyi_,work,1 );
	BLAS_AXPY( dim,s,ev,inc_1,r,inc_1 );
	res = BLAS_NRM2( dim,r,inc_1 );

	return res;
}

/*
	v0.1	cys
		7/29/2007
*/
void ARE_mode2str( ARE_R_SOLVER *hAR_R,char *restart )	{
	switch( hAR_R->restart_mode )	{
		case GSE_RESTART_SORENSEN:
			strcpy( restart,"SORENSEN" );
			break;
		default:
			break;
	}
}

/*
	v0.1	cys
		5/16/2007
*/
void SQRT( GRUS_FLOATING *x )	{
	double real=REAL_COMPONENT(*x),imag=IMAG_COMPONENT(*x),r,xita;
	r = sqrt( real*real+imag*imag );
	r = sqrt( r );
	xita = atan2( imag,real )/2.0;

	real = r*cos(xita );		imag = r*sin(xita );
	REAL_COMPONENT(*x)=real,	IMAG_COMPONENT(*x)=imag;
}

/*
	solve z=inv(A)*v with GSS 

	v0.1	cys
		5/12/2007
*/
void ARE_gss_solve( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v,GRUS_FLOATING *z )	{
	int i,dim=hAR_R->dim;
	double *x=hAR_R->d_temp;

	for( i = 0; i < dim; i++ )	x[i]=REAL_COMPONENT(v[i]);
	ASSERT(0);
//	GSS_solve_udi( hAR_R->hLU,dim,hAR_R->ptr,hAR_R->ind,hAR_R->val,x );
	for( i = 0; i < dim; i++ )	REAL_COMPONENT(z[i])=x[i];
	for( i = 0; i < dim; i++ )	x[i]=IMAG_COMPONENT(v[i]);
	ASSERT(0);
//	GSS_solve_udi( hAR_R->hLU,dim,hAR_R->ptr,hAR_R->ind,hAR_R->val,x );
	for( i = 0; i < dim; i++ )	IMAG_COMPONENT(z[i])=x[i];
}

/*
	z=spectrum transform(v)

	v0.1	cys
		5/12/2007
*/
void ARE_trans_V( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v,GRUS_FLOATING *z )	{
	GRUS_FLOATING s,t,*val,*z_1,c;
	int i,j,dim=hAR_R->dim,*ptr=hAR_R->ptr,*ind=hAR_R->ind,info=0;
	double *d_val,*x=hAR_R->d_temp,d,res;
	GRUS_FLOATING *b=hAR_R->work+dim;
	clock_t start=clock( );

	for( i = 0; i < dim; i++ )	CLEAR(z[i]);

	if( hAR_R->ivt_mode > 0 )	{	//inv(A)*v
		switch( hAR_R->ivt_mode )	{
		case GSE_INVERT_GSS:
			if( hAR_R->sft_mode == GSE_SHIFT_MA )	{
				z_1 = hAR_R->work;		ASSERT( z_1 != z );
				ARE_gss_solve( hAR_R,v,z_1 );
				ARE_gss_solve( hAR_R,z_1,z );
				SUB( c,hAR_R->shift,hAR_R->delta );
				ABS( d,c );							SET( c,-d*d,0.0 );
				BLAS_SCAL( dim,c,z,inc_1 );
				d = REAL_COMPONENT(hAR_R->delta)-REAL_COMPONENT(hAR_R->shift);
				SET( c,2*d,0.0 );
				BLAS_AXPY( dim,c,z_1,inc_1,z,inc_1 );
			}else
				ARE_gss_solve( hAR_R,v,z );
			break;
		case GSE_INVERT_MF_1:
//			mf_solve_1( BLAS_NOTRANS,dim,hAR_R->hLU,v,z,info );
//			mf_solve_1( BLAS_N,dim,hAR_R->hLU,v,z,info );
			break;
		default:		//z=Av
			break;
		}	
	}else	{
		/*if( hAR_R->data_type==GSE_DATA_DOUBLE )	{	//double version
			d_val = hAR_R->val;
			for( i = 0; i < dim; i++ )	{
				ASSIGN( s,v[i] );
				for( j = ptr[i]; j < ptr[i+1]; j++ )	{
					SET( t,d_val[j],0.0 );
					MULT_ADD( z[ind[j]],t,s );	
				}
			}
		}else*/	{									//complex version
			val = hAR_R->val;
			for( i = 0; i < dim; i++ )	{
				ASSIGN( s,v[i] );
				for( j = ptr[i]; j < ptr[i+1]; j++ )	{
					MULT_ADD( z[ind[j]],val[j],s );	
				}
			}
		}
	}
	GRUS_SE_INFO[GRUS_SE_MV_TIME] += clock( )-start;
	GRUS_SE_INFO[GRUS_SE_MV_COUNT] ++;
}

/*
	v0.1	cys
		5/16/2007
*/
void ARE_trans_lenda( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *lenda )	{
	GRUS_FLOATING a,b,c,x,y,z;
	double d,s_a,r_l=REAL_COMPONENT(*lenda),i_l=IMAG_COMPONENT(*lenda);
	GRUS_FLOATING lenda_1,lenda_2;
	int s_mode = hAR_R->sft_mode&GSE_SHIFT_TYPE;
	if( s_mode == GSE_SHIFT_MA )	{
		s_a=REAL_COMPONENT( hAR_R->shift );
		ASSIGN( a,*lenda );
		REAL_COMPONENT( b ) = -2*REAL_COMPONENT(hAR_R->delta)	+(1.0-r_l)*2*s_a;
		IMAG_COMPONENT( b ) =									 (   -i_l)*2*s_a;
//		ABS( d,hAR_R->delta );
//		REAL_COMPONENT( c ) = d*d	-(1.0-r_l)*s_a*s_a;
//		IMAG_COMPONENT( c ) =		-(   -i_l)*s_a*s_a;
		MULT_CONJ( c,hAR_R->delta,hAR_R->delta );
		REAL_COMPONENT( c ) += -(1.0-r_l)*s_a*s_a;
		IMAG_COMPONENT( c ) += -(   -i_l)*s_a*s_a;
		MULT( x,b,b );			
		MULT( y,a,c );			SCALE_RECIP( y,4.0 );
		DECREMENT( x,y );		//b*b-4*a*c 
		SQRT( &x );
		ASSIGN( y,a );			SCALE_RECIP( y,2.0 );
		SUB( z,x,b );			DIV( lenda_1,z,y );	
		ASSIGN( z,x);			SCALE_RECIP( z,-1.0 );		DECREMENT( z,b );				
		DIV( lenda_2,z,y );
		ASSIGN( *lenda,lenda_2 );			
	}else if( s_mode == GSE_SHIFT_NORMAL ){
		if( hAR_R->ivt_mode > 0 )	{	
			RECIPROCAL( *lenda );
		}
		ASSEMBLE( *lenda,hAR_R->shift );

	}else if( hAR_R->ivt_mode > 0 )	{	
		RECIPROCAL( *lenda );
	}
}

/*
	v0.1	cys
		5/10/2007

int ARE_invert( void* hSolver,int flag )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int nRet=GSE_OK,dim=hAR_R->dim,info=0x0,*ptr=hAR_R->ptr,*ind=hAR_R->ind;
	clock_t start=clock( );

	hAR_R->hLU = NULL;
	ASSERT( hAR_R->ivt_mode > 0 );
	switch( hAR_R->ivt_mode )	{
	case GSE_INVERT_GSS:
		ASSERT( hAR_R->val_type==GSE_DATA_DOUBLE );
		memset( GRUS_MF_CONTROL,0x0,sizeof(double)*GRUS_MF_CONTROL_ITEM );
		memset( GRUS_MF_INFO,0x0,sizeof(double)*GRUS_MF_INFO_ITEM );
		info = mf_G_InitDefault_udi( dim,dim,hAR_R->ptr,hAR_R->ind,hAR_R->val,GRUS_MF_UNI,GRUS_MF_AUTO );
		if( info != GRUS_OK )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}		
		hAR_R->hLU = mf_symbol_udi( dim,dim,hAR_R->ptr,hAR_R->ind,hAR_R->val );	
		if( hAR_R->hLU == NULL )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}		
		mf_numeric_udi( hAR_R->val,hAR_R->hLU );
		if( GRUS_MF_INFO[GRUS_MF_STATUS] != GRUS_OK )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}
		info = 0x0;
		break;
	case GSE_INVERT_MF_1:
		ASSERT( hAR_R->val_type==GSE_DATA_COMPLEX );
		hAR_R->hLU = mf_symbol_1( dim,ptr[dim],ptr,ind );			
		if( hAR_R->hLU == NULL )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}
		mf_numeric_1( dim,hAR_R->val,hAR_R->hLU,info );		
		break;
	default:		//z=Av
		break;
	}	

END:
	if( info != 0x0 )	{
		printf( "\r\n\tERROR at INVERT. ERROR CODE:%g\r\n",info );
		nRet=GSE_INVERT_FAIL;
	}
	GRUS_SE_INFO[GRUS_SE_INVERT_TIME] += clock( )-start;
	return nRet;
}
*/

/*
	A -= shift*I

	假设	val_type==DoubleComplex
		
	v0.1	cys
		5/16/2007
	v0.2	cys
		7/10/2007
*/
int ARE_shift( void* hSolver,int flag )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,j,nRet=GSE_OK,dim=hAR_R->dim,nnz=hAR_R->nnz,*ptr,*ind,nX=0,nz,isFind;
	double real=REAL_COMPONENT(hAR_R->shift),imag=IMAG_COMPONENT(hAR_R->shift),*d_val;
	GRUS_FLOATING *val;

	ASSERT( hAR_R->sft_mode > 0 || hAR_R->val_type!=hAR_R->A_type /*&& !IS_ZERO(hAR_R->shift)*/ );

	ASSERT( hAR_R->val_type==GSE_DATA_COMPLEX );
	ptr=hAR_R->A_ptr;		ind=hAR_R->A_ind;
	for( i = 0; i < dim; i++ )	{
		isFind = 0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			if( ind[j]==i )
				isFind=1;
		}
		if( isFind == 0 )
			nX++;
	}

	hAR_R->val=GRUS_alloc( sizeof(GRUS_FLOATING)*(nnz+nX) );		
	hAR_R->ptr=GRUS_alloc( sizeof(int)*(dim+1) );		
	hAR_R->ind=GRUS_alloc( sizeof(int)*(nnz+nX) );		
	val = (GRUS_FLOATING *)(hAR_R->val);
	nz=0;
	hAR_R->ptr[0]=0;
	for( i = 0; i < dim; i++ )	{
		isFind=0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			if( hAR_R->A_type==GSE_DATA_DOUBLE )	{	//double=>complex
				d_val = hAR_R->A;
				if( ind[j]!=i )		{
					if( ind[j]>i && isFind==0 )	{
						SET( val[nz],-real,-imag );
						hAR_R->ind[nz]=i;		isFind=1;
						nz++;
					}
					hAR_R->ind[nz]=ind[j];
					SET( val[nz],d_val[j],0.0 );
				}else	{
					isFind=1;
					hAR_R->ind[nz]=ind[j];
					SET( val[nz],d_val[j]-real,-imag );
				}		
			}else	{									//complex version
				ASSERT( 0 );
			}
			nz++;
		}
		if( isFind == 0 )	{
			SET( val[nz],-real,-imag );
			hAR_R->ind[nz]=i;
			nz++;
		}
		hAR_R->ptr[i+1]=nz;
	}
	ASSERT( nnz+nX==nz );
END:
	return nRet;
}

/*
	v0.1	cys
		5/2/2007
*/
void ai_expand( void* hSolver,GRUS_FLOATING*z,int j )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int dim=hAR_R->dim,i;
	GRUS_FLOATING *v,c;
	double res,*max_res=GRUS_SE_INFO+GRUS_SE_MAX_INVERT_RES;

	v = hAR_R->V+j*dim;
	ARE_trans_V( hAR_R,v,z );
/*	if( hAR_R->ivt_mode > 0 ) {		//	| v-(A-delta*I)*z |
		ARE_AV_0( hAR_R,z,hAR_R->work );		
		if( !IS_ZERO(hAR_R->shift) )	{
			SET( c,-REAL_COMPONENT(hAR_R->shift),-IMAG_COMPONENT(hAR_R->shift) );
			BLAS_AXPY( dim,c,z,1,hAR_R->work,1 );	//-delta*z
		}
		BLAS_AXPY( dim,fuyi_,v,1,hAR_R->work,1 );	//Az-v
		res = BLAS_NRM2( dim,hAR_R->work,1 );		
		*max_res = MAX( *max_res,res );
	}	*/
}

/*
	z=A*(v1,v2,...vj)*H(:,j)		其中z=V(:,j+1),在ai-loop中等于A*vj

	v0.1	cys
		5/2/2007
*/
void ai_update_( ARE_R_SOLVER *hAR_R,int j  )	{
	GRUS_FLOATING *h,*z,*c,*u,s;
	int dim=hAR_R->dim,ldh=hAR_R->ldh,n=j+1,i,t;
	double norm,norm_0,norm_1;
	clock_t start=clock( );

	z=hAR_R->V+(j+1)*dim;
	norm_0 = BLAS_NRM2( dim,z,inc_1 );
	if( j+1==hAR_R->nCurEV )	//如何证明	必相当于deflation过程
		h=hAR_R->work+ldh;
	else
		h=hAR_R->H+j*ldh;
//	orthogonal:		z = z-V*(V,z)
//	BLAS_GEMV( BLAS_CONJTRANS,dim,n,one_,hAR_R->V,dim,z,inc_1,zero_,h,inc_1 );	//h=V'z
	BLAS_GEMV( BLAS_C,dim,n,one_,hAR_R->V,dim,z,inc_1,zero_,h,inc_1 );	//h=V'z
//	BLAS_GEMV( BLAS_NOTRANS,dim,n,fuyi_,hAR_R->V,dim,h,inc_1,one_,z,inc_1 );	
	BLAS_GEMV( BLAS_N,dim,n,fuyi_,hAR_R->V,dim,h,inc_1,one_,z,inc_1 );	
//DGKS correction if needed
	norm = BLAS_NRM2( dim,z,inc_1 );	//normalize z
	if( norm < hAR_R->ortho_thresh*norm_0 )	{
		clock_t s_0=clock( );
		norm_0 = norm;
		c = hAR_R->work;
//		BLAS_GEMV( BLAS_CONJTRANS,dim,n,one_,hAR_R->V,dim,z,inc_1,zero_,c,inc_1 );	//h=-V'z
		BLAS_GEMV( BLAS_C,dim,n,one_,hAR_R->V,dim,z,inc_1,zero_,c,inc_1 );	//h=-V'z
//		BLAS_GEMV( BLAS_NOTRANS,dim,n,fuyi_,hAR_R->V,dim,c,inc_1,one_,z,inc_1 );	
		BLAS_GEMV( BLAS_N,dim,n,fuyi_,hAR_R->V,dim,c,inc_1,one_,z,inc_1 );	
//		BLAS_AXPY( j+1,one_,c,inc_1,h,inc_1 );	
		i=j+1;		BLAS_AXPY( i,one_,c,inc_1,h,inc_1 );	
		norm = BLAS_NRM2( dim,z,inc_1 );	//normalize z
		if( norm < 0.707*norm_0 )	{	//z lies in the span of V_{j} numerically
			printf("\r\n!!!z lies in the span of V_{j} numerically!!!\r\n");
//			for( i = 0; i < dim; i++ )	CLEAR( z[i] );
		}
		GRUS_SE_INFO[GRUS_SE_REORTHO_TIME] += clock( )-s_0;	
	}
	SET(h[j+1],norm,0.0);
	SET(s,1.0/norm,0.0);
	BLAS_SCAL( dim,s,z,inc_1 );

	GRUS_SE_INFO[GRUS_SE_UPDATE_TIME] += clock( )-start;
}

/*
	v0.1	cys
		5/29/2007

void GSE_spectrum_order( ARE_R_SOLVER *hAR_R,int flag )	{
	int i,j,n=hAR_R->V_m;
	double a,b;
	GRUS_FLOATING t;

	for( i = hAR_R->nCurEV; i < n; i++ )	{		//order the spectrum
		ABS(a,hAR_R->w_H[i]);
		for( j = i+1; j < n; j++ )	{
			ABS(b,hAR_R->w_H[j]);
			if( a < b )	{
				ASSIGN(t,hAR_R->w_H[i]);	ASSIGN(hAR_R->w_H[i],hAR_R->w_H[j]);	ASSIGN(hAR_R->w_H[j],t);
				a = b;
			}
		}
	}
}
*/
/*
	compute right ritz pari of A;	and normalize ev_H

	v0.1	cys
		5/5/2007
*/
int GSE_H_eigen( void* hSolver,int H_type ){
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,j,k,n=hAR_R->V_m-hAR_R->nCurEV,dim=hAR_R->dim,ldvl=n,ldvr,*select,mm,m,*i_temp,*ifaill,*ifailr;
	int ldh=hAR_R->ldh,ldz=1,lwork=4*n,info=-1,ilo=1,ihi=n;
	char E[1]={'E'},V[1]={'V'},I[1]={'I'},N[1]={'N'},U[1]={'U'},R[1]={'R'},eigsrc[1]={'N'},initv[1]={'N'};	//only eigenvalue,no schur vector
	GRUS_FLOATING *z=0x0,*vl=0x0,*rwork,s,*H,*w_H,*ev_H=hAR_R->ev_H;
	double norm,*v_u=0x0,*v_l=0x0,abstol=0.0,norm_0;

	mm=hAR_R->nRV;
	i_temp=(int*)malloc( sizeof(int)*(mm*2+5*n) );
	select=i_temp;			
	ifaill=select+n;		ifailr=ifaill+mm;
	ldvr=hAR_R->ldh;
//	for( i = 0; i < n; i++ )	select[i]=1;		
	for( i = 0; i < hAR_R->V_m; i++ )	select[i]=i<hAR_R->nCurEV?0:1;			//7/6/2007
	H=hAR_R->work+lwork;
//	H= malloc( sizeof(GRUS_FLOATING)*ldh*ldh );
//	BLAS_COPY( ldh*ldh,hAR_R->H,inc_1,H,inc_1 );
	i=ldh*ldh;			BLAS_COPY( i,hAR_R->H,inc_1,H,inc_1 );
	H += hAR_R->nCurEV*ldh+hAR_R->nCurEV;		
	switch( H_type )	{
	case GSE_MAT_HESSENBERG:
#ifdef _G_COMPLEX_
		w_H = hAR_R->w_H+hAR_R->nCurEV;
		zhseqr(E,N,&n,&ilo,&ihi,H,&ldh,w_H,z,&ldz, hAR_R->work,&lwork,&info);
		ASSERT( info==0 );
//		GSE_spectrum_order( hAR_R,0x0 );
		ZV_order_x( hAR_R->V_m,hAR_R->w_H,zero_,GSE_ORDER_MAX );
		rwork=hAR_R->work+n*n;		
		mm = n;
		n=hAR_R->V_m;		
		ev_H=hAR_R->ev_H+hAR_R->nCurEV*ldh;
		zhsein ( R,eigsrc,initv,select,&n,hAR_R->H,&ldh,hAR_R->w_H,vl,
				&ldvl,ev_H,&ldvr,&mm,&m,hAR_R->work,rwork,ifaill,ifailr,&info );
		ASSERT( info==0 && m==mm );
#else
		rv=hAR_R->d_temp,		iv=rv+n;
		dhseqr(job,compz,&n,&ilo,&ihi,H,&ldh,rv,iv,z,&ldz, hAR_R->work,&lwork,&info);
		for( i = 0; i < n; i++ )	
		{	REAL_COMPONETN(hAR_R->w_H)=rv[i];	COMPLEX_COMPONETN(hAR_R->w_H)=iv[i];	}
//		GSE_spectrum_order(hAR_R);
		select[1]=hAR_R->wi_H[0]==0.0 ? 0 : 1;
		rv=hAR_R->d_temp,		iv=rv+n;
		for( i = 0; i < n; i++ )	
		{	rv[i]=REAL_COMPONETN(hAR_R->w_H);	iv[i]=COMPLEX_COMPONETN(hAR_R->w_H);	}
		H = hAR_R->H+ldh*hAR_R->nCurEV+hAR_R->nCurEV;
		dhsein ( job,eigsrc,initv,select,&n,H,&ldh,rv,iv,vl,
				&ldvl,ev_H,&ldvr,&mm,&m,hAR_R->work,ifaill,ifailr,&info );
#endif
		break;
	case GSE_MAT_HERMITIAN:
#ifdef ZINT
		ldz=n;		
		ifaill=i_temp+5*n;
		v_u=v_l=hAR_R->d_temp;
		zheevx(V, I, U,&n,H,&ldh, v_l,v_u, &ihi,&ihi,&abstol,&m,hAR_R->w_H,
			ev_H,&ldz,hAR_R->work,&lwork,hAR_R->d_temp,i_temp,ifaill,&info);
		IMAG_COMPONENT(hAR_R->w_H[0]) = 0.0;
		ASSERT( m==1 && info==0 );
//		zheev(V, U,&n,H,&ldh,hAR_R->w_H,hAR_R->work, &lwork, hAR_R->d_temp,&info);
//		m=1;
#else
		ASSERT( 0 );
#endif		
		break;
	case GSE_MAT_GENERAL:
		n=hAR_R->V_m;		
		H=hAR_R->work+lwork;			w_H=hAR_R->w_H;			ev_H=hAR_R->ev_H;
		zgeev(N, V,&n,H,&ldh, w_H, vl,&ldvl, ev_H, &ldh, hAR_R->work,&lwork,hAR_R->d_temp,&info);
		ASSERT( info==0 );
		for( i = 0; i < n; i++ )	{
			ABS( norm_0,w_H[i] );	ihi=-1;
			for( j = i; j < n; j++ )	{
				ABS( norm,w_H[j] )
				if( norm>norm_0 )	
				{	norm_0=norm;	ihi=j;		}
			}
			if( ihi!=-1 )	{
				j = ihi;
				for( k = 0; k<n; k++ )
				{	ASSIGN(s,ev_H[i*ldh+k]);ASSIGN(ev_H[i*ldh+k],ev_H[j*ldh+k]);ASSIGN(ev_H[j*ldh+k],s);}
				ASSIGN(s,w_H[i]);	ASSIGN(w_H[i],w_H[j]);	ASSIGN(w_H[j],s)
			}
		}
		break;
	default:
		break;
	}
	free( i_temp );

	ev_H=hAR_R->ev_H;
	for( i = hAR_R->nCurEV; i < n; i++ )	{
		norm = BLAS_NRM2( n,ev_H+i*ldh,inc_1 );	//normalize ev_H
		SET( s,1.0/norm,0.0 );
		BLAS_SCAL( n,s,ev_H+i*ldh,inc_1 );
	}


	return hAR_R->info;
}

/*
	subspace extraction 

	注意：		
		V orthonormal，|rv|=|ev_H|,必然等于1

	v0.1	cys
		5/4/2007
*/
int GSE_extraction( void* hSolver,int H_type ){
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,n=hAR_R->V_m,dim=hAR_R->dim,t=hAR_R->nCurEV;
	GRUS_FLOATING *rv=hAR_R->rv_A+t*dim,*w=hAR_R->work,*ev_H=hAR_R->ev_H+t*hAR_R->ldh,s,*z;
	double res,norm;

#ifdef _DEBUG
//	ZGE_output( n,hAR_R->H,hAR_R->ldh,"H:\\GSP\\GSE\\test\\H.dat",1 );
#endif
	ASSERT( n <= hAR_R->nGV );
	GSE_H_eigen( hAR_R,H_type );
//	BLAS_COPY( dim,rv,1,w,1);
//	BLAS_GEMV( BLAS_NOTRANS,dim,n,one_,hAR_R->V,dim,ev_H,inc_1,zero_,rv,inc_1 );		//
	BLAS_GEMV( BLAS_N,dim,n,one_,hAR_R->V,dim,ev_H,inc_1,zero_,rv,inc_1 );		//
	norm = BLAS_NRM2( dim,rv,inc_1 );	
	if( fabs(norm-1.0) > 1.0e-12 )	{
		ARE_trace( "\tERROR!!! code=%d,at EV-%d",GSE_LOST_ORTHOGONAL,t+1 );
		GRUS_SE_INFO[GRUS_SE_STATUS]=GSE_LOST_ORTHOGONAL;
		return GSE_LOST_ORTHOGONAL;	//
	}

	return GSE_OK;
}

/*
	m-step arnoldi factorization

	注意：
		需寻找ai的收敛标准,即H[j,j+1]<epsi的值

	[ref] An Arnoldi code for computing selected eigenvalues of sparse real unsymmetric matrices (1995)
	v0.1	cys
		5/2/2007
*/
int ARE_ai_core( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int j,k,m=hAR_R->nGV,dim=hAR_R->dim,ldh=hAR_R->ldh,ret=ARE_RET_NORMAL;
	GRUS_FLOATING *z,s;
	double norm,epsi=0.0;//	//mch_eps/*hAR_R->tole*hAR_R->a_norm;*/;

	hAR_R->V_m = m;
	k = MAX( hAR_R->nCurEV,hAR_R->V_p );
	for( j = k; j < m; j++ )	{
		z = hAR_R->V+(j+1)*dim;		
		ai_expand( hAR_R,z,j );			//expand krylov space
		ai_update_( hAR_R,j );		
		ABS( norm,hAR_R->H[j*ldh+j+1] );
/*		norm = BLAS_NRM2( dim,z,1 );	//normalize z
		SET( hAR_R->H[j*ldh+j+1],norm,0.0);
		SET(s,1.0/norm,0.0);
		BLAS_SCAL( dim,s,z,1 );*/
		ASSERT(  norm > mch_eps );
		if( norm < epsi )		{	//find invariant space		7/4/2007
			hAR_R->V_m = j+1;
			ret=ARE_AI_CONVERGE;	goto END;		
		}
//		if( AR_isConverge( )==1 )
//			break;
	}
END:
	return ret;
}

/*
	Arnoldi校验,校验Ax=xita*x+beta*V[:,m+1]的偏差,与restart密切相关

	v0.1	cys
		7/14/2007
*/
int ARE_Verify_Resid( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,t=hAR_R->nCurEV,t_1=hAR_R->nEV,dim=hAR_R->dim,ldh=hAR_R->ldh,no=-1,m=hAR_R->V_m;
	GRUS_FLOATING beta,s,*work=hAR_R->work,*ev_H=hAR_R->ev_H,*ev=work+dim;
	double norm,max_dot=0.0,*off=GRUS_SE_INFO+GRUS_SE_MAX_OFF_KRYLOV,err_3,err_h=0.0;
	clock_t start=clock( );

/*	for( i = t; i < m; i++ )	{
		BLAS_GEMV( BLAS_NOTRANS,m,m,one_,hAR_R->H,ldh,hAR_R->ev_H+i*ldh,inc_1,zero_,work,inc_1 );	
		MULT( s,fuyi_,hAR_R->w_H[i] );
		BLAS_AXPY( m,s,hAR_R->ev_H+i*ldh,1,work,1 );
		norm = BLAS_NRM2( m,work,1 );	
		err_h = MAX( err_h,norm );
	}*/

	for( i = t; i < t_1; i++ )	{
//		BLAS_GEMV( BLAS_NOTRANS,dim,m,one_,hAR_R->V,dim,ev_H+ldh*i,inc_1,zero_,ev,inc_1 );		
		BLAS_GEMV( BLAS_N,dim,m,one_,hAR_R->V,dim,ev_H+ldh*i,inc_1,zero_,ev,inc_1 );		
		BLAS_COPY( dim,ev,inc_1,work,inc_1 );
		ARE_trans_V( hAR_R,ev,work );			//A*x
//		ai_expand( hAR_R,work,i );			
		ASSIGN( s,hAR_R->w_H[i] );			SCALE_RECIP( s,-1 );
		BLAS_AXPY( dim,s,ev,inc_1,work,inc_1 );
		norm=BLAS_NRM2(dim,work,inc_1 );
		MULT( beta,hAR_R->H[(m-1)*ldh+m],ev_H[ldh*i+m-1] );		//
		SCALE_RECIP( beta,-1 );
		BLAS_AXPY( dim,beta,hAR_R->V+m*dim,inc_1,work,inc_1 );
		norm=BLAS_NRM2(dim,work,inc_1 );
		*off=MAX( *off,norm );
		if( *off > 1.0e-5 )	{		//unreasonable large error
//			BLAS_GEMV( BLAS_NOTRANS,m,m,one_,hAR_R->H,ldh,hAR_R->ev_H+i*ldh,inc_1,zero_,work,inc_1 );	
			BLAS_GEMV( BLAS_N,m,m,one_,hAR_R->H,ldh,hAR_R->ev_H+i*ldh,inc_1,zero_,work,inc_1 );	
			MULT( s,fuyi_,hAR_R->w_H[i] );
			BLAS_AXPY( m,s,hAR_R->ev_H+i*ldh,inc_1,work,inc_1 );
			err_3 = BLAS_NRM2( m,work,inc_1 );	//normalize z
		}
	}

	GRUS_SE_INFO[GRUS_SE_VERIFY_TIME] += clock( )-start;

	return hAR_R->info;
}

/*
	v0.1	cys
		6/10/2007
	注意：
		1、z0,z1,z2借用V空间，因此V_m>=nEV+2
		2、m_0=nEV,不同于nCurEV

	ref:	[The Arnoldi-Tchebycheff method for solving large complex non hermitian generalized eigenproblems.]

void GSE_CHEBYCHEV_0( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v )	{
	int m_0=hAR_R->nEV,m_1=hAR_R->V_m,i,j,isLarge,dim=hAR_R->dim,K;
	double x_1,x_0,y_1,y_0,x_d,y_d,a=0.0,b,c,t;
	GRUS_FLOATING *w_H=hAR_R->w_H,delta_0,delta_1,delta_i,*z0,*z1,*z2,s,s_0,s_1,s_2;
	GRUS_FLOATING *work=hAR_R->d_temp,*lenda=w_H+m_0-1;		

	ASSERT( m_1>m_0 );
	x_1=x_0=REAL_COMPONENT(w_H[m_0]);		
	y_1=y_0=IMAG_COMPONENT(w_H[m_0]);
	for( i = m_0+1;  i < m_1; i++ )	{	//[m_0,m_1) is the unwanted part of the spectrum
		x_1=MAX( x_1,REAL_COMPONENT(w_H[i]) );		x_0=MIN( x_0,REAL_COMPONENT(w_H[i]) );
		y_1=MAX( y_1,IMAG_COMPONENT(w_H[i]) );		y_0=MIN( y_0,IMAG_COMPONENT(w_H[i]) );
	}
	x_d=(x_1+x_0)/2.0;			y_d=(y_1+y_0)/2.0;
	a = (x_1-x_d)*1.001;		//适当放大一点点
	isLarge=1;
	while( isLarge==1 )	{
		b = 0.0;
		for( i =  m_0;  i < m_1; i++ )	{
			x_0=REAL_COMPONENT(w_H[i]);		y_0=IMAG_COMPONENT(w_H[i]);
			ASSERT( fabs(a-(x_0-x_d)) > 1.0e-16 );
			t = a*fabs(y_0-y_d)/sqrt( a*a-(x_0-x_d)*(x_0-x_d) );
			b= MAX( b,t );
		}
		isLarge = 0;
		for( i =  m_0;  i < m_1; i++ )	{
			x_0=REAL_COMPONENT(w_H[i]);		y_0=IMAG_COMPONENT(w_H[i]);
			t = (x_0-x_d)*(x_0-x_d)/(a*a)+(y_0-y_d)*(y_0-y_d)/(b*b);
			if( t > 1.0 )	{
				a += (REAL_COMPONENT(*lenda)-x_1)/2.0;
				isLarge = 1;		break;
			}
		}
	}
	if( a>b )
		c = sqrt( a*a-b*b );
	else
		c = sqrt( b*b-a*a );

	z0=v,				z1=z0+dim;
	SET( delta_0,c,0.0 );	
	SET( s,REAL_COMPONENT(*lenda)-x_d,IMAG_COMPONENT(*lenda)-y_d );
	DIV_1( delta_0,s );		//delta=e/(lenda_1-c)
//	printf("delta_0=(%g,%g)\r\n",REAL_COMPONENT(delta_0),IMAG_COMPONENT(delta_0));
	ARE_trans_V( hAR_R,z0,z1 );
	SET( s,-x_d/c,-y_d/c );			MULT( s_0,delta_0,s );
	SET( s_1,REAL_COMPONENT(delta_0)/c,IMAG_COMPONENT(delta_0)/c );			
	BLAS_SCAL( dim,s_1,z1,1 );
	for( j = 0; j < dim; j++ )	{		//z1=delta/e*z1-delta/e*c*z0
		MULT_ADD( z1[j],s_0,z0[j] );
	}	

	K = hAR_R->V_m;
	SET( delta_i,REAL_COMPONENT(delta_0)/2.0,IMAG_COMPONENT(delta_0)/2.0 );			RECIPROCAL( delta_i );
	for( i = 0; i < K; i ++ )	{
		SUB( delta_1,delta_0,delta_i );			RECIPROCAL( delta_1 );
		z0=v+(i%3)*dim,		z1=v+(i+1)%3*dim,		z2=v+(i+2)%3*dim;
		ARE_trans_V( hAR_R,z1,z2 );
		ASSIGN( s_2,delta_1 );			SCALE_RECIP( s_2,2.0/c );	//s_2
		SET(s,-x_d,-y_d);				MULT( s_1,s_2,s );			//s_1=2.0*delta_1/e*(-c)
		MULT( s_0,delta_1,delta_0 );	SCALE_RECIP( s_0,-1.0 );	//s_0=-delta_0*delta_1
		BLAS_SCAL( dim,s_2,z2,1 );
//	printf("s_1=(%g,%g),\ts_2=(%g,%g)\r\n",REAL_COMPONENT(s_1),IMAG_COMPONENT(s_1),REAL_COMPONENT(s_2),IMAG_COMPONENT(s_2));
		for( j = 0; j < dim; j++ )	{
			MULT_ADD( z2[j],s_1,z1[j] );
			MULT_ADD( z2[j],s_0,z0[j] );
		}
		ASSIGN( delta_0,delta_1 );
	}

	if( z2 != v )
		BLAS_COPY( dim,z2,1,v,1 );

}
*/
/*
	v0.1	cys
		6/10/2007
	注意：
		1、z0,z1,z2借用V空间，因此V_m>=nEV+2
		2、m_0=nEV,不同于nCurEV
*/
void GSE_CHEBYCHEV_0( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v )	{
	int m_0=hAR_R->nEV,m_1=hAR_R->V_m,i,j,isLarge,dim=hAR_R->dim,K;
	GRUS_FLOATING *w_H=hAR_R->w_H,delta_0,delta_1,delta_i,*z0,*z1,*z2,s,s_0,s_1,s_2;
	GRUS_FLOATING *lenda=w_H+m_0-1,e,c,*cheby,*c2,*c1,*c0;		
	double p,p_max=0.0,tao;

	if( m_1<=m_0 )		//没有unwanted set		7/10/2007
		return;

	tao=Optimal_Ellipse_0( m_0,m_1,w_H,&e,&c );	//e is foci and c is the center

	z0=v,				z1=z0+dim;
	ASSIGN( delta_0,e );				//SET( delta_0,c,0.0 );	
	SUB(s,*lenda,c );					//SET( s,REAL_COMPONENT(*lenda)-x_d,IMAG_COMPONENT(*lenda)-y_d );
	DIV_1( delta_0,s );		//delta=e/(lenda_1-c)
//	printf("delta_0=(%g,%g)\r\n",REAL_COMPONENT(delta_0),IMAG_COMPONENT(delta_0));
	ARE_trans_V( hAR_R,z0,z1 );
	DIV(s,c,e);		SCALE_RECIP( s,-1.0 );		//SET( s,-x_d/c,-y_d/c );			
	MULT( s_0,delta_0,s );
	DIV( s_1,delta_0,e);						//SET( s_1,REAL_COMPONENT(delta_0)/c,IMAG_COMPONENT(delta_0)/c );			
	BLAS_SCAL( dim,s_1,z1,inc_1 );
	for( j = 0; j < dim; j++ )	{		//z1=delta/e*z1-delta/e*c*z0
		MULT_ADD( z1[j],s_0,z0[j] );
	}	

	cheby=malloc( sizeof(GRUS_FLOATING)*m_1*3 );
	for( j = m_0-1; j < m_1; j++ )	{
		SET( cheby[3*j],1.0,0.0 );		ASSIGN( cheby[3*j+1],w_H[j] ); 
	}
	K = 5;//hAR_R->V_m;
	SET( delta_i,REAL_COMPONENT(delta_0)/2.0,IMAG_COMPONENT(delta_0)/2.0 );			RECIPROCAL( delta_i );
	for( i = 0; i < K; i ++ )	{
		SUB( delta_1,delta_0,delta_i );			RECIPROCAL( delta_1 );
		z0=v+(i%3)*dim,		z1=v+(i+1)%3*dim,		z2=v+(i+2)%3*dim;
		ARE_trans_V( hAR_R,z1,z2 );
		ASSIGN( s_2,delta_1 );		SCALE_RECIP(s_2,2.0);		DIV_1(s_2,e);		//SCALE_RECIP( s_2,2.0/c );	//s_2
		SET(s,-REAL_COMPONENT(c),-IMAG_COMPONENT(c));			MULT( s_1,s_2,s );			//s_1=2.0*delta_1/e*(-c)
		MULT( s_0,delta_1,delta_0 );	SCALE_RECIP( s_0,-1.0 );	//s_0=-delta_0*delta_1
		BLAS_SCAL( dim,s_2,z2,inc_1 );
//	printf("s_1=(%g,%g),\ts_2=(%g,%g)\r\n",REAL_COMPONENT(s_1),IMAG_COMPONENT(s_1),REAL_COMPONENT(s_2),IMAG_COMPONENT(s_2));
		for( j = 0; j < dim; j++ )	{
			MULT_ADD( z2[j],s_1,z1[j] );
			MULT_ADD( z2[j],s_0,z0[j] );
		}
		ASSIGN( delta_0,delta_1 );
		for( j = m_0-1; j < m_1; j++ )	{		//C(k+1)=2*lenda*C(k)-C(k-1)
			c0 = cheby+3*j+(i)%3;	c1 = cheby+3*j+(i+1)%3;	c2 = cheby+3*j+(i+2)%3;
			MULT( *c2,w_H[j],*c1 );			SCALE_RECIP(*c2,2.0);	
			DECREMENT( *c2,*c0 );			
		}	
	}

	c1 = cheby+3*(m_0-1)+(i-1)%3;
	for( j = m_0; j < m_1; j++ )	{
		c0 = cheby+3*j+(i-1)%3;	
		DIV(s,*c0,*c1 );
		ABS( p,s );
		p_max=MAX( p_max,p );
	}
	printf( "\tCHEBYCHEV:\tK=%d,tao=%g,p-MAX=%g\r\n",K,tao,p_max );
	free( cheby );

	if( z2 != v )
		BLAS_COPY( dim,z2,inc_1,v,inc_1 );
}

/*
	f=v'*A*v

	v0.1	cys
		7/28/2007
*/
void ARE_field( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v,GRUS_FLOATING *f )	{
	int dim=hAR_R->dim;
	GRUS_FLOATING *z=hAR_R->work;

	ARE_trans_V( hAR_R,v,z );
	BLAS_DOT_c( dim,v,inc_1,z,inc_1,*f );
}

/*
	t=epsi*M'*u-M'*r

	v0.1	cys
		7/18/2007
*/
void ARE_Restart_Jacobi( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *v )	{
	int t=hAR_R->nCurEV,info=0,dim=hAR_R->dim;
	GRUS_FLOATING *r=hAR_R->work,*x_v=r+dim,*x_r=hAR_R->d_temp,es_u,es_r,s,epsi;
	double res_0,res_1;

	ASSIGN( s,hAR_R->w_H[t] );
	ARE_trans_lenda( hAR_R,&s );

	res_0 = Eigen_Defect( hAR_R,s,v,r );

//	mf_solve_1( BLAS_NOTRANS,dim,hAR_R->hLU,v,x_v,info );
//	mf_solve_1( BLAS_N,dim,hAR_R->hLU,v,x_v,info );
//	mf_solve_1( BLAS_NOTRANS,dim,hAR_R->hLU,r,x_r,info );
//	mf_solve_1( BLAS_N,dim,hAR_R->hLU,r,x_r,info );
	BLAS_DOT_c( dim,v,inc_1,x_v,inc_1,es_u );
	BLAS_DOT_c( dim,v,inc_1,x_r,inc_1,es_r );
	DIV( epsi,es_r,es_u );
	BLAS_SCAL( dim,fuyi_,x_r,inc_1 );
	BLAS_AXPY( dim,epsi,x_v,inc_1,x_r,inc_1 );

	BLAS_AXPY( dim,one_,x_r,inc_1,v,inc_1 );
	res_1 = Eigen_Defect( hAR_R,s,v,r );
	ASSERT( res_1 < res_0 );
}
/*
	REF:
		IMPLICITLY RESTARTED ARNOLDI LANCZOS METHODS FOR LARGE SCALE EIGENVALUE CALCULATIONS

	ISSUE-NEEDED:
		1、k>=hAR_R->nEV的动态选择	

	v0.1	cys
		7/19/2007
*/
void ARE_Restart_Sorensen( ARE_R_SOLVER *hAR_R )	{
	int i,j,k,n,dim=hAR_R->dim,t=hAR_R->nCurEV,nEV=hAR_R->nEV,m=hAR_R->V_m,nP=m-nEV,ldh=hAR_R->ldh,ldq=m;
	GRUS_FLOATING *w=hAR_R->w_H,*shift,*V=hAR_R->V+t*dim,*H=hAR_R->H,*resid,*Q,*workL,*workD,*v,s;
	double norm,off=0.0,ratio,r_max;

	shift=w+nEV;		
	resid = V+dim*m;
	BLAS_SCAL( dim,H[(m-1)*ldh+m],resid,inc_1 );
	workD=hAR_R->work;		workL=workD+2*dim;		Q=workL+m;		//[m*(m+1)+2*dim]	
	k=hAR_R->nEV-t;
/*	r_max = 0.0;
	for( i = hAR_R->nEV; i < MIN( hAR_R->nEV*2,m-1 ); i++ )	{
		DIV( s,hAR_R->w_H[i-1],hAR_R->w_H[i] );
		ABS( ratio,s );
		if( ratio > r_max )
		{	r_max=ratio;		k=i-t;		nP=m-i;	}
	}
	GRUS_SE_INFO[GRUS_SE_RESTART_RATIO] = r_max;
	
	for( i = hAR_R->nEV+1; i < m-1; i++ )	{
		DIV( s,hAR_R->w_H[i],hAR_R->w_H[hAR_R->nEV] );
		ABS( ratio,s );
		if( fabs(ratio-1.0) > 0.1 )
			break;
		else
		{	r_max=ratio;		k=i-t;		nP=m-i;	}
	}
*/
	ZNAPPS( &dim,&k,&nP,shift,V,&dim,H+t*ldh+t,&ldh,resid, Q,&ldq, workL, workD );
/*	norm = BLAS_NRM2( dim,hAR_R->V,1 );
	norm = BLAS_NRM2( dim,hAR_R->V+dim*2,1 );
	BLAS_DOT_c( dim,hAR_R->V,1,hAR_R->V+dim*2,1,s);
	ABS( norm,s );*/

//	BLAS_COPY( dim,resid,1,V+dim*k,1 );
	BLAS_COPY( dim,resid,inc_1,V+dim*k,inc_1 );
	norm=BLAS_NRM2(dim,resid,inc_1 );
	SET( H[(k-1)*ldh+k],norm,0.0);
	SET(s,1.0/norm,0.0);
//	BLAS_SCAL( dim,s,V+dim*k,1 );
	BLAS_SCAL( dim,s,V+dim*k,inc_1 );
#ifdef _DEBUG
	//验证AVk=VkHk+fe*
/*	for( i = 0; i < nEV; i++ )	{
		v = hAR_R->V+i*dim;
		norm = BLAS_NRM2( dim,v,1 );
		off=MAX( off,fabs(norm-1.0) );
		n = nEV-i-1;
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,one_,v+dim,dim,v,inc_1,zero_,hAR_R->work,inc_1 );	//w=V'v
		for( j = 0; j< n; j++ ){
			ABS( norm,hAR_R->work[j] );
			ASSERT( !SCALAR_IS_NAN(norm) );
//			if( max_dot > norm )
//			{	max_dot = norm;		no=i*m+j;	}
			off=MAX( off,norm );
		}
		j=i+2;		
		ai_expand( hAR_R,hAR_R->work,i );
		BLAS_GEMV( BLAS_NOTRANS,dim,j,fuyi_,hAR_R->V,dim,hAR_R->H+i*ldh,inc_1,one_,hAR_R->work,inc_1 );		//
		norm=BLAS_NRM2(dim,hAR_R->work,1 );
		off=MAX( off,norm );
	}	
	ASSERT( off < 1.0e-12 );*/
#endif
	hAR_R->V_p=k;
}

/*
	v0.1	cys
		5/6/2007
*/
void ARE_restart( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int dim = hAR_R->dim,t=hAR_R->nCurEV,i,n=hAR_R->V_m-t,ldh=hAR_R->ldh,m_1=hAR_R->V_m;
	double norm,off;
	GRUS_FLOATING s,*rv=hAR_R->rv_A+t*dim,*v=hAR_R->V+t*dim,*ev_H=hAR_R->ev_H,
		*work=hAR_R->work,*v_0=work+dim,*w_H=hAR_R->w_H+t,*c=work+dim,*beta=c+ldh;
	clock_t start=clock( );

	DIV( s,hAR_R->w_H[hAR_R->nEV-1],hAR_R->w_H[hAR_R->nEV] );
	ABS( norm,s );
	GRUS_SE_INFO[GRUS_SE_RESTART_RATIO] = norm;
	if( hAR_R->restart_mode==GSE_RESTART_SORENSEN )	{
		ARE_Restart_Sorensen( hAR_R );
		goto END;
	}
	BLAS_COPY( dim,rv,inc_1,v,inc_1 );
	ASSIGN( s,ev_H[ldh*t+m_1-1] );
	switch( hAR_R->restart_mode  )	{
		case GSE_RESTART_RANDOM:		//bad option	7/10/2007
			GSE_random_1( dim,v );			
			break;
		case GSE_RESTART_CHEBYCHEV_0:
			GSE_CHEBYCHEV_0( hAR_R,v );
			break;
		case GSE_RESTART_JACOBI:			
			ARE_Restart_Jacobi( hAR_R,v );
			break;
		case GSE_RESTART_MOGAN:
			for( i=t; i<hAR_R->nEV; i++ )		ASSIGN( beta[i],ev_H[ldh*i+m_1-1] );
			Vandermonde_ZX_Solve( hAR_R->nEV-t,w_H,beta+t,c+t );
		break;
		default:
		break;
	}
	if( hAR_R->nEV>t+1 && hAR_R->restart_mode>0 )	{	//multiple combine
		for( i=t+1; i<hAR_R->nEV; i++ )	{
//			BLAS_GEMV( BLAS_NOTRANS,dim,m_1,one_,hAR_R->V,dim,ev_H+ldh*i,inc_1,zero_,work,inc_1 );		//
			BLAS_GEMV( BLAS_N,dim,m_1,one_,hAR_R->V,dim,ev_H+ldh*i,inc_1,zero_,work,inc_1 );		//
			BLAS_AXPY(dim,c[i],work,inc_1,v,inc_1 );
		}
	}
	if( 0 )	{
		DIV( s,ev_H[ldh*t+m_1-1],ev_H[ldh*(t+1)+m_1-1] );			
		SCALE_RECIP( s,-1.0 );
//		BLAS_GEMV( BLAS_NOTRANS,dim,m_1,one_,hAR_R->V,dim,ev_H+ldh*(t+1),inc_1,zero_,work,inc_1 );		//
		BLAS_GEMV( BLAS_N,dim,m_1,one_,hAR_R->V,dim,ev_H+ldh*(t+1),inc_1,zero_,work,inc_1 );		//
		BLAS_AXPY(dim,s,work,inc_1,v,inc_1 );
	}

	if( t > 0 )	{
		ai_update_( hAR_R,t-1 );
	}
	norm = BLAS_NRM2( dim,v,inc_1 );
	SET( s,1.0/norm,0.0 );
	BLAS_SCAL( dim,s,v,inc_1 );
	
END:
	if( hAR_R->checkOffBase )	{
		BLAS_DOT_c(dim,hAR_R->Base,inc_1,v,inc_1,s);
		ABS( off,s );
		ARE_trace( "<%.3lf>",off );
	}
	GRUS_SE_INFO[GRUS_SE_RESTART_TIME] += clock( )-start;
	GRUS_SE_INFO[GRUS_SE_RESTART_COUNT] ++;

	hAR_R->nRestart++;
}

/*
	采用normwise backward error

	v0.1	cys
		5/6/2007
	v0.2	cys
		5/11/2007

	 注意：
		1、显示的对应于A的特征值。
		2、需要区分多种情况：精度振荡，精度不变
	ISSUE-NEEDED:
		1、为啥CLEAR( hAR_R->H[t*ldh+t+1] )至关重要
*/
int ARE_is_stop( void* hSolver,double *res_norm )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,dim=hAR_R->dim,ldh=hAR_R->ldh,isConverge=0,m=hAR_R->V_m,t=hAR_R->nCurEV,n=hAR_R->V_m-t;
	double norm,res,res_y,epsi=hAR_R->tole*hAR_R->a_norm,off=0,ratio;
	GRUS_FLOATING lenda,*H,s,s_0;
	char F_norm[1]={'F'};		//Frobenius norm
	GRUS_FLOATING *rv=hAR_R->rv_A+hAR_R->nCurEV*dim,*w_H=hAR_R->w_H,*z,*work=hAR_R->work;

	ASSERT( m>=hAR_R->nEV );
//	h_norm = LAPACK_LANHS( F_norm,&m,hAR_R->H,&ldh,hAR_R->d_temp );
	ASSIGN( lenda,w_H[t] );
	switch( hAR_R->stop_mode )	{
	case GSE_STOP_DIRECT:	//lenda*e-Ae
		ARE_trans_lenda( hAR_R,&lenda );
		res = Eigen_Defect( hAR_R,lenda,rv,hAR_R->work );
		if( res < epsi )
			isConverge=1;
		break;
	case GSE_STOP_BLOCK:	
		isConverge=1;			res=0.0;
		for( i = t; i <  hAR_R->nEV; i++ )	{
			z=hAR_R->rv_A+i*dim;
//			BLAS_GEMV( BLAS_NOTRANS,dim,n,one_,hAR_R->V+t*dim,dim,hAR_R->ev_H+i*ldh,inc_1,zero_,z,inc_1 );		
			BLAS_GEMV( BLAS_N,dim,n,one_,hAR_R->V+t*dim,dim,hAR_R->ev_H+i*ldh,inc_1,zero_,z,inc_1 );		
			ASSIGN( s,w_H[i] );
			ARE_trans_lenda( hAR_R,&s );	
			res_y=Eigen_Defect( hAR_R,s,z,hAR_R->work );
			if( res_y > res )			{	res=res_y;	}		
//			if( res_y < epsi )	
//				CLEAR( hAR_R->H[i*ldh+i+1] );
		}
		isConverge = res<epsi;		
		break;	
	case GSE_STOP_KRYLOV:
		MULT( s,hAR_R->H[(m-1)*ldh+m],hAR_R->ev_H[t*ldh+m-1] );
		ABS( res,s );
//		BLAS_GEMV( BLAS_NOTRANS,hAR_R->V_m,hAR_R->V_m,one_,hAR_R->H,ldh,hAR_R->ev_H+t*ldh,inc_1,zero_,work,inc_1 );	
		BLAS_GEMV( BLAS_N,hAR_R->V_m,hAR_R->V_m,one_,hAR_R->H,ldh,hAR_R->ev_H+t*ldh,inc_1,zero_,work,inc_1 );	
		MULT( s,fuyi_,lenda );
//		BLAS_AXPY( hAR_R->V_m,s,hAR_R->ev_H+t*ldh,1,work,1 );
		BLAS_AXPY( hAR_R->V_m,s,hAR_R->ev_H+t*ldh,inc_1,work,inc_1 );
//		res_y = BLAS_NRM2( hAR_R->V_m,work,1 );	
		res_y = BLAS_NRM2( hAR_R->V_m,work,inc_1 );	
//		res=MAX( res,res_y );
		if( res < epsi/*MAX( mch_eps*h_norm,mch_eps*norm)*/ )
			isConverge=1;
		break;
	}

//	ABS( ratio,hAR_R->w_H[hAR_R->nEV-1] );		ABS( norm,hAR_R->w_H[hAR_R->nEV] );	
//	ratio/=norm;
	ASSIGN( lenda,w_H[t] );
//	ARE_trans_lenda( hAR_R,&lenda );
	if( hAR_R->V_m!=hAR_R->nGV )	//7/12/2007
		ARE_trace("[%d]",hAR_R->V_m);
	ratio = GRUS_SE_INFO[GRUS_SE_RESTART_RATIO];
	ARE_trace( "%4d:\t(%.13e,%.13e),\t%g\t%.2lf\r\n",hAR_R->nRestart+1,REAL_COMPONENT(lenda),IMAG_COMPONENT(lenda),res,ratio );
	if( isConverge )	{	//尽量不破坏V的结构
		if( hAR_R->stop_mode==GSE_STOP_BLOCK )
		{	hAR_R->nCurEV=hAR_R->nEV;	goto END;}
	/*	for( i = hAR_R->nCurEV+1; i <  hAR_R->nEV; i++ )	{
			z=hAR_R->rv_A+i*dim;
			BLAS_GEMV( BLAS_NOTRANS,dim,n,one_,hAR_R->V+t*dim,dim,hAR_R->ev_H+i*ldh,inc_1,zero_,z,inc_1 );		
			ASSIGN( s,w_H[i] );
			ARE_trans_lenda( hAR_R,&s );	
			res=Eigen_Defect( hAR_R,s,z,hAR_R->work );
			if( res <= epsi )	{
				ARE_trace( "%d Converged!!!\r\n",i+1 );
			}
		}*/
		z = hAR_R->V+(t+1)*dim;		//deflate+restart
//		BLAS_GEMV( BLAS_NOTRANS,dim,n,one_,hAR_R->V+t*dim,dim,hAR_R->ev_H+(t+1)*ldh,inc_1,zero_,rv+dim,inc_1 );		
		BLAS_GEMV( BLAS_N,dim,n,one_,hAR_R->V+t*dim,dim,hAR_R->ev_H+(t+1)*ldh,inc_1,zero_,rv+dim,inc_1 );		
//		BLAS_COPY( dim,rv,1,hAR_R->V+t*dim,1 );
		BLAS_COPY( dim,rv,inc_1,hAR_R->V+t*dim,inc_1 );
		if( t > 0 )	{			//得到schur向量
			ai_update_( hAR_R,t-1 );
		}
		ai_expand( hAR_R,z,t );	//得到H[:,t]，确认其收缩为R[:,t]。
		ai_update_( hAR_R,t );
//		ABS( norm,hAR_R->H[t*ldh+t+1] );
//		off=MAX( off,norm );
		SUB( s,hAR_R->w_H[t],hAR_R->H[t*ldh+t] );		ABS( norm,s );
		off=MAX( off,norm );
		ASSIGN( w_H[t],hAR_R->H[t*ldh+t] );		//ASSIGN( hAR_R->H[t*ldh+t],w_H[t] );		
		CLEAR( hAR_R->H[t*ldh+t+1] );
		hAR_R->nCurEV++;
	}
END:
	*res_norm =res;
	return isConverge; 
}

/*
	v0.1	cys
		6/4/2007
*/
int ARE_outputEV( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,j,dim=hAR_R->dim,nEV=hAR_R->nCurEV;
	GRUS_FLOATING *ev=hAR_R->rv_A;

	FILE *fp=fopen( "..\\GSE\\test\\eigen_pair.dat","w" );	
	fprintf( fp,"%d %d\r\n",nEV,dim );
	for( i = 0; i < nEV; i++ )	{
		fprintf( fp,"eigenvalue:\t%.13e\t%.13e\r\n",REAL_COMPONENT(hAR_R->w_H[i]),IMAG_COMPONENT(hAR_R->w_H[i]) );
		for( j = 0; j < dim; j++ )	{
			fprintf( fp,"\t%.13e\t%.13e\r\n",REAL_COMPONENT(ev[j]),IMAG_COMPONENT(ev[j]) );
		}
		ev += dim;
	}
	fclose( fp );

	return nEV;
}

/*
	v0.1	cys
		6/4/2007
*/
int ARE_loadEV( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,j,dim,nEV,ldh=hAR_R->ldh;
	GRUS_FLOATING *rv=hAR_R->rv_A,*v=hAR_R->V,lenda;
	double real,imag,off,r_lenda,i_lenda;	

	FILE *fp=fopen( ".\\GSE\\test\\eigen_pair.dat","r" );	
	fscanf( fp,"%d %d\r\n",&nEV,&dim );
	if( nEV<0 || nEV>=hAR_R->nGV-1 || dim != hAR_R->dim )
		return 0;

	ASSERT( nEV==1 );		//以下代码只适用于nEV==1
	for( i = 0; i < nEV; i++ )	{
		fscanf( fp,"eigenvalue:\t%lf %lf\r\n",&r_lenda,&i_lenda );
//		SET( hAR_R->w_H[0],real,imag );
		for( j = 0; j < dim; j++ )	{
			fscanf( fp,"%lf %lf",&real,&imag );
			SET( rv[j],real,imag );
		}
//		BLAS_COPY( dim,rv,1,v,1 );
		BLAS_COPY( dim,rv,inc_1,v,inc_1 );
		if( i > 0 )				//得到schur向量
			ai_update_( hAR_R,i-1 );
		ai_expand( hAR_R,v+dim,i );	//得到H[:,t]，确认其收缩为R[:,t]。
		ai_update_( hAR_R,i );
		ASSIGN( lenda,hAR_R->H[i*ldh+i] );
		ARE_trans_lenda( hAR_R,&lenda );
		r_lenda-=REAL_COMPONENT(lenda);		i_lenda-=IMAG_COMPONENT(lenda);
		off=sqrt( r_lenda*r_lenda+i_lenda*i_lenda );
		ASSERT( off<1.0e-10 );
		ASSIGN( hAR_R->w_H[i],hAR_R->H[i*ldh+i] );		//ASSIGN( hAR_R->H[t*ldh+t],w_H[t] );		
		ABS( off,hAR_R->H[i*ldh+i+1] );
		ASSERT( off<1.0e-10 );
		CLEAR( hAR_R->H[i*ldh+i+1] );

		hAR_R->nCurEV++;
		hAR_R->nEV ++;
		rv += dim;
	}
	fclose( fp );

	rv=hAR_R->rv_A+nEV*dim;
	GSE_random_1( dim,rv );
	if( hAR_R->checkOffBase )	{
		i=ZV_load( dim,hAR_R->Base,"..\\GSE\\test\\ev_2.dat" );
//		real=BLAS_NRM2( dim,hAR_R->Base,1 );
		real=BLAS_NRM2( dim,hAR_R->Base,inc_1 );
		ASSERT( i == 0 && fabs(real-1.0)<1.0e-14 );
	}
	ARE_restart( hAR_R );

	return nEV;
}

/*
	读入song生成的sss文件
	调用：
		GSE_ReadSSS( 477,"F:\\IEEE\\ORIGIN\\Matrix_PureA_NoSparse0.sss"	);

	v0.1	cys
		7/10/2007
*/
void 	GSE_ReadSSS( int dim,char *sPath	)	{
	int i,j,nz=0,nzMax=dim*dim,*ptr,*ind,len,ret;
	double *val,a;

	FILE *fp=fopen( sPath,"r" );

	ptr=malloc( sizeof(int)*(dim+1) );
	ind=malloc( sizeof(int)*(nzMax) );
	val=malloc( sizeof(double)*(nzMax) );
	ptr[0]=0;
	for( i = 0; i < dim; i++ )	{
		for( j = 0; j < dim; j++ )	{
			if( fscanf( fp,"%lf",&a) !=1 )
			{	ret=i+1;	goto END;	}
			if( a==0.0 && j!=i )		//保证输出对角元	7/11/2007
				continue;
			ind[nz]=j;		val[nz]=a;
			nz++;
		}
		fscanf( fp,"  ;  %d",&len);		
		ASSERT( len==dim );
		ptr[i+1]=nz;
	}
END:
	fclose( fp );

	CRS2MTX_ouput( dim,ptr, ind, val,"F:\\IEEE\\ORIGIN\\A_477_1.mtx" );	
	free(ptr);		free(ind);		free(val);

}

/*
	ISSUE-NEEDED:
		1、v1增强,参见下面的代码

	v0.1	cys
		5/2/2007
*/
int ARE_main( void* hSolver,double *w )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,j,k,dim=hAR_R->dim,ret=GSE_OK,m=hAR_R->nGV,ldh=hAR_R->ldh;
	double resi=0.0,norm,res;
	GRUS_FLOATING s,*z;
	clock_t start=clock( );

	GSE_random_1( dim,hAR_R->V );
//	i=ZV_load( dim,hAR_R->V,"H:\\GSP\\GSE\\test\\ev_1.dat" );

	if( hAR_R->sft_mode != GSE_SHIFT_NO || hAR_R->val_type!=hAR_R->A_type/*&& !IS_ZERO(hAR_R->shift)*/ )
		ARE_shift( hAR_R,0x0 );
	if( hAR_R->ivt_mode > 0 )	{
		hAR_R->hLU = GSE_invert( dim,hAR_R->ptr,hAR_R->ind,hAR_R->val,GSE_INVERT_MF_1,GSE_DATA_COMPLEX,&ret );
		if( ret!= GSE_OK )
		{	goto END;			}
	}
/*	v1增强	7/30/3007
	ARE_field( hAR_R,hAR_R->V,hAR_R->w_H );
	ARE_Restart_Jacobi( hAR_R,hAR_R->V );
	norm = BLAS_NRM2( dim,hAR_R->V,inc_1 );
	SET( s,1.0/norm,0.0 );
	BLAS_SCAL( dim,s,hAR_R->V,inc_1 );
*/
//	ARE_loadEV( hSolver );
	while( 1 )	{
		hAR_R->V_m = m;
		k = MAX( hAR_R->nCurEV,hAR_R->V_p );
		for( j = k; j < m; j++ )	{
			z = hAR_R->V+(j+1)*dim;		
			ai_expand( hAR_R,z,j );					//expand krylov space
			ai_update_( hAR_R,j );		
			ABS( norm,hAR_R->H[j*ldh+j+1] );		ASSERT(  norm > mch_eps );
/*			if( norm < epsi )		{	//find invariant space		7/4/2007
				hAR_R->V_m = j+1;
				ret=ARE_AI_CONVERGE;	goto END;		
			}*/
		}
//		ARE_verify_V( hSolver );
//		if( GSE_extraction( hSolver,GSE_MAT_HESSENBERG )!=0 )
		if( GSE_extraction( hSolver,GSE_MAT_GENERAL )!=0 )
			break;
		ARE_Verify_Resid( hSolver );
		if( ARE_is_stop( hSolver,&res )==1 )	{
			printf( "\t\t---------------------\r\n" );
		//	break;
		}
		if( hAR_R->nCurEV==hAR_R->nEV )	
			break;
		if( hAR_R->nRestart>=hAR_R->MAX_RESTART )	
			break;
		ARE_restart( hSolver );
	}
	for( i = 0; i < hAR_R->nCurEV; i++ )	{
		ARE_trans_lenda( hAR_R,&(hAR_R->w_H[i]) );
		ASSIGN( s,hAR_R->w_H[i]);
		w[2*i]=REAL_COMPONENT(s),		w[2*i+1]=IMAG_COMPONENT(s);
	}

//	ARE_outputEV( hSolver );
END:

	GRUS_SE_INFO[GRUS_SE_TOTAL_TIME] += clock( )-start;
	return ret;
}



