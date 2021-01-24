#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "..\pch\grus_def.h"
#include "..\pch\grus_version.h"
#include "..\pch\grus_version_lib.h"

#include "GSE_Krylov.h"
#include "GSE_util.h"

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static IMPLEMENT_ONE_ZERO( one_,fuyi_,zero_ );	
//DECLARE_GRUS_TRANS( )
//DECLARE_ONE_ZERO( one_,fuyi_,zero_ )

/*
	z=a*Av+b*z

	v0.1	cys
		5/4/2007

void ARE_GEMV( ARE_R_SOLVER *hAR_R,GRUS_FLOATING *a,GRUS_FLOATING *v,
			  GRUS_FLOATING *b,GRUS_FLOATING *z )	{
	GRUS_FLOATING s,t,*val;
	int i,j,dim=hAR_R->dim,*ptr=hAR_R->ptr,*ind=hAR_R->ind;
	double *d_val;

	if( IS_ZERO(*b) )	{
		for( i = 0; i < dim; i++ )	CLEAR(z[i]);
	}else if( REAL_COMPONENT(*b)==1.0 && IMAG_COMPONENT(*b)==0.0 )	{
	}else	{
		for( i = 0; i < dim; i++ )	{
			MULT( s,z[i],*b );
			ASSIGN(z[i],s);
		}
	}	
	if( hAR_R->data_type==GSE_DATA_DOUBLE )	{	//double version
		d_val = hAR_R->val;
		for( i = 0; i < dim; i++ )	{
			MULT( s,v[i],*a );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				SET( t,d_val[j],0.0 );
				MULT_ADD( z[ind[j]],t,s );	
	//			v[ind[j]] += val[j]*s;
			}
		}
	}else	{		//complex version
		val = hAR_R->val;
		for( i = 0; i < dim; i++ )	{
	//		ASSIGN( s,v[i] );
			MULT( s,v[i],*a );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				MULT_ADD( z[ind[j]],val[j],s );	
	//			v[ind[j]] += val[j]*s;
			}
		}
	}
}
*/

/*
	v0.1	cys
		5/12/2007
*/
double CCS_DLA_1( int dim,int *ptr,int *ind,double *val )	{
	int i,j;
	double a_norm = 0.0,norm;

	for( i = 0; i < dim; i++ )	{
		norm=0.0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			norm += fabs(val[j]);
		}
		if( norm > a_norm )	a_norm = norm;
	}	
	return a_norm;
};

/*
	v0.1	cys
		5/12/2007
*/
double CCS_DLA_F( int dim,int *ptr,int *ind,double *val )	{
	int i,j;
	double a_norm = 0.0;
	for( i = 0; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			a_norm += val[j]*val[j];
		}
	}
	
	a_norm=sqrt( a_norm );
	return a_norm;
};

/*
	find optimal ellipse E(c, e, a) be an ellipse containing the set w_H[m_0,m_1), 
		and having (real) center c, foci c + e, c - e, and major semi-axis a. 

	v0.1	
		cys	6/29/2007

	ref:	[The Arnoldi-Tchebycheff method for solving large complex non hermitian generalized eigenproblems.]
*/
double Optimal_Ellipse_0( int m_0,int m_1,GRUS_FLOATING *w_H,GRUS_FLOATING *e,GRUS_FLOATING *c )	{
	int i,isLarge;
	double x_1,x_0,y_1,y_0,x_d,y_d,a=0.0,b,foci,t;
	GRUS_FLOATING *lenda=w_H+m_0-1;

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

//	ASSERT( a>b );
	if( a>b )	{
		foci = sqrt( a*a-b*b );
		SET( *e,foci,0.0 );
	}else	{
		foci = sqrt( b*b-a*a );
		SET( *e,0.0,foci );
	}
	SET( *c,x_d,y_d );

	x_0=REAL_COMPONENT(*lenda);		y_0=IMAG_COMPONENT(*lenda);
	t = (x_0-x_d)*(x_0-x_d)/(a*a)+(y_0-y_d)*(y_0-y_d)/(b*b);
	return t;
}

/*
	v0.1	cys
		7/10/2007
*/
int GSE_outputEV( int nEV,int dim,GRUS_FLOATING *w,GRUS_FLOATING *ev,char* sPath )	{
	int i,j;

	FILE *fp=fopen( sPath,"w" );	
	fprintf( fp,"%d %d\r\n",nEV,dim );
	for( i = 0; i < nEV; i++ )	{
		fprintf( fp,"eigenvalue %d:\t%.15e\t%.15e\r\n",i+1,REAL_COMPONENT(w[i]),IMAG_COMPONENT(w[i]) );
		for( j = 0; j < dim; j++ )	{
			fprintf( fp,"\t%.15e\t%.15e\r\n",REAL_COMPONENT(ev[j]),IMAG_COMPONENT(ev[j]) );
		}
		ev += dim;
	}
	fclose( fp );

	return nEV;
}

/*
	v0.1	cys
		7/10/2007
*/
int ZV_load( int dim,DoubleComplex *z,char* sPath )	{
	int i,len,ret=0;
	double real,imag;

	FILE *fp=fopen( sPath,"r" );	
	fscanf( fp,"%d",&len );
	if( len!=dim )
	{	ret=-1;	goto END;	}

	for( i = 0; i < dim; i++ )	{
		if( fscanf( fp,"%lf %lf",&real,&imag ) !=2 )
		{	ret=i+1;	goto END;	}
		SET( z[i],real,imag );
	}
END:
	fclose( fp );

	return ret;
}

/*
	注意：
		暂时限定输出z的转置

	v0.1	cys
		7/12/2007
*/
int ZGE_output( int dim,DoubleComplex *z,int ldz,char* sPath,int flag )	{
	int i,j,ret=0;
	DoubleComplex *v;
	FILE *fp=fopen( sPath,"w+" );	
	char format[20]={"(%5.2lf,%5.2lf)"};

	v=z;
	if( BIT_TEST(flag,GSE_OUTPUT_SHORT) )
		strcpy( format,"(%5.2lf,%5.2lf)" );
	else
		strcpy( format,"(%.15lf,%.15lf)" );

	fprintf( fp,"%d\r\n",dim );
	for( i = 0; i < dim; i++ )	{
		for( j = 0; j < dim; j++ )	{
			fprintf( fp,format,REAL_COMPONENT(v[j*ldz]),IMAG_COMPONENT(v[j*ldz]) );
		}
		fprintf( fp,"\r\n" );
		v++;
	}
END:
	fclose( fp );

	return ret;
}

/*
	注意：
		暂时限定输出z的转置

	v0.1	cys
		7/16/2007
*/
int ZGE_load( int *n,DoubleComplex **mat,char* sPath,int flag )	{
	int i,j,ret=0,dim;
	DoubleComplex *v,*z;
	double real,imag;

	FILE *fp=fopen( sPath,"r" );			
	if( fp == 0x0 )
	{	ret=GSE_FOPEN_ERR;	goto END;	}
	fscanf( fp,"%d\r\n",&dim );

	z=(DoubleComplex*)malloc( sizeof(DoubleComplex)*dim*dim );
	v=z;
	for( i = 0; i < dim; i++ )	{
		for( j = 0; j < dim; j++ )	{
			if( fscanf( fp,"(%lf,%lf)",&real,&imag )!=2 )
			{	ret=(i+1)*dim+j+1;	goto END;	}
			SET( v[j*dim],real,imag );
		}
		fscanf( fp,"\r\n" );
		v++;
	}
END:
	fclose( fp );

	*n=dim;			*mat=z;
	return ret;
}

/*
	注意：
		暂时限定输出z的转置

	v0.1	cys
		7/12/2007
*/
int Vandermonde_ZX_Solve( int k,DoubleComplex *xita,DoubleComplex *X,DoubleComplex *c )	{
	int i,j,nrhs=1,info=0,*ipiv;
	DoubleComplex *A,*b,s;
	double norm;
	char N[1]={'N'};

	ASSERT( k >= 1 );
	if( k==1 )	{
		SET( c[0],1.0,0.0 );
		return info;
	}

	A=GRUS_alloc( sizeof(DoubleComplex)*k*(k+1) );		b=A+k*k;
	for( i = 0; i < k; i++ )	{
		CLEAR( b[i] );
		SET( A[i*k],1.0,0.0 );		ASSIGN( A[i*k+1],X[i] );
		for( j = 2; j < k; j++ )	{
			MULT( A[i*k+j],A[i*k+j-1],xita[i] );
		}
	}
	SET( b[0],1.0,0.0 );
#ifdef _DEBUG
//	ZGE_output( k,A,k,"H:\\GSP\\GSE\\test\\V.dat",GSE_OUTPUT_SHORT );
#endif

	ipiv=GRUS_alloc( sizeof(int)*k );
	zgetrf ( &k,&k,A,&k, ipiv,&info );					ASSERT( info==0 );
	zgetrs(N, &k, &nrhs, A, &k, ipiv, b, &k, &info);	ASSERT( info==0 );
	GRUS_free( ipiv );

	ABS( norm,b[0] );				ASSERT( norm != 0.0 );
	for( i = 0; i < k; i++ )
		DIV( c[i],b[i],b[0] );
//Verify
	CLEAR( s );
	for( i = 0; i < k; i++ )
		MULT_ADD( s,X[i],c[i] );
	ABS( norm,s );					ASSERT( fabs(norm) <1.0e-12 );

	GRUS_free( A );

	return info;
}

/*
	v0.1	cys
		7/31/2007

void ZV_select( int m,DoubleComplex *w,int K,int *select,int *order,int flag )	{
	int i,j,base,cur,best,isAlloc=order==0x0;
	double a,b;

	if( isAlloc )
		order=malloc( sizeof(int)*m );
	for( i = 0; i < m; i++ )	
	{	order[i]=i;		select[i]=0;		}

	for( i = 0; i < K; i++ )	{		//order the spectrum
		base = order[i];
		ABS(a,w[base]);
		best=i;
		for( j = i+1; j < m; j++ )	{
			cur = order[j];
			ABS(b,w[cur]);
			if( a < b )	{
				a = b;		best=j;
			}
		}
		if( best != i )	{
			order[i]=order[best];		order[best]=base;
		}
	}
	for( i = 0; i < K; i++ )	{		
		base = order[i];
		select[base]=1;
	}
	if( isAlloc )
		free( order );
}
*/

/*
	从m个特征值中选出K个，不一定严格排序
	与order的区别在于不变动w

	v0.1	cys
		7/31/2007
*/
int GSE_Select( int m,DoubleComplex *w,int K,int *select,int *order,DoubleComplex off,double thresh,int flag )	{
	int i,j,base,cur,best,isAlloc=order==0x0,isSame,K_new=K,old,type=(flag&GSE_ORDER_TYPE);
	double a,b;

	if( K==0 )
		goto END;
	if( isAlloc )
		order=malloc( sizeof(int)*m );
	for( i = 0; i < m; i++ )	
	{	order[i]=i;		select[i]=0;		}

	for( i = 0; i < K; i++ )	{		//order the spectrum
		base = order[i];
		OFF(a,w[base],off );		//ABS(a,w[base]);
		best=i;
		for( j = i+1; j < m; j++ )	{
			cur = order[j];
			OFF(b,w[cur],off );		//ABS(b,w[cur]);
			switch( type )	{
			case GSE_ORDER_MAX:
				if( a < b )	
				{	a = b;		best=j;		}
				break;
			case GSE_ORDER_MIN:
				if( a > b )	
				{	a = b;		best=j;		}
				break;
			default:
				break;
			}
		}
		if( best != i )	{
			order[i]=order[best];		order[best]=base;
		}
	}
	for( i = 0; i < K; i++ )	{		
		base = order[i];
		select[base]=1;
	}
//test the clustered eigenvalue
	base = order[0];		ABS(a,w[base]);
	for( j = K; j < m; j++ )	{
		cur = order[j];		ABS(b,w[cur]);
		if( fabs(a-b)<thresh*a )	{
			old=order[K_new];	order[K_new]=cur;		order[j]=old;
			select[cur]=1;
			K_new++;
		}
	}
	if( isAlloc )
		free( order );
END:
	return K_new-K;
}

/*	for( i = k; i < m; i++ )	{
		v = hAR_R->V+i*dim;
		norm = BLAS_NRM2( dim,v,1 );
		*off_OrNo=MAX( *off_OrNo,fabs(norm-1.0) );
		n = m-i-1;
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,one_,v+dim,dim,v,inc_1,zero_,work,inc_1 );	//w=V'v
		for( j = 0; j< n; j++ ){
			ABS( norm,hAR_R->work[j] );
			if( max_dot > norm )
			{	max_dot = norm;		no=i*m+j;	}
			*off_OrNo=MAX( *off_OrNo,norm );
		}
	}*/
/*	
	校验M是unitra matrix,利用U'*U=I

	v0.1	cys
		8/1/2007
*/
int ZU_verify( int m,int n,DoubleComplex* mat,int ldm,int flag )	{
	int i,j,ret=0;
	GRUS_FLOATING *I=malloc( sizeof(DoubleComplex)*n*n );
	double norm,off=0.0;
	clock_t start=clock( );

	I=malloc( sizeof(DoubleComplex)*n*n );
//	BLAS_GEMM( BLAS_CONJTRANS,BLAS_NOTRANS,n,n,m,one_,mat,ldm,mat,ldm,zero_,I,n );	
	BLAS_GEMM_0( BLAS_C,BLAS_N,n,n,m,one_,mat,ldm,mat,ldm,zero_,I,n );	
	for( i = 0; i < n; i++ )	{
		off = 0.0;
		for( j = i; j < n; j++ )	{
			ABS( norm,I[i*n+j] );
			if( i==j )	{
				off=MAX( off,fabs(norm-1.0) );
			}else	{
				off=MAX( off,norm );
			}
		}
		if( off>1.0e-5 )
		{	ret=-(i+1);			goto END;	}
	}
END:
	free( I );
	GRUS_SE_INFO[GRUS_SE_MAX_OFF_ORTHONORMAL]=MAX( GRUS_SE_INFO[GRUS_SE_MAX_OFF_ORTHONORMAL],off );
	GRUS_SE_INFO[GRUS_SE_VERIFY_TIME] += clock( )-start;

	return ret;
}

/*
	v0.1	cys
		8/3/2007
*/
int ZV_select_verify( int m,DoubleComplex *w,int K,int flag )	{
	int i,j,ret=0;
	double a,b;

	for( i = 0; i < K; i++ )	{
		ABS(a,w[i]);
		for( j = K+1; j < m; j++ )	{
			ABS( b,w[j] );
			if( a < b )	
			{	ret=-(i+1);			goto END;			}
		}
	}
END:
	return ret;
}

/*
	按与x的距离排序

	v0.1	cys
		8/7/2007
*/
void ZV_order_x( int n,DoubleComplex *w,DoubleComplex x,int flag )	{
	int i,j,best,old,type=(flag&GSE_ORDER_TYPE);
	DoubleComplex s;
	double a,b;

	for( i = 0; i < n; i++ )	{		//order the spectrum
		OFF( a,w[i],x );				//ABS(a,w[i]);
		best=i;
		for( j = i+1; j < n; j++ )	{
			OFF( b,w[j],x );			//ABS(b,w[i]);	
			switch( flag )	{
			case GSE_ORDER_MAX:
				if( a < b )	
				{	a = b;		best=j;		}
				break;
			case GSE_ORDER_MIN:
				if( a > b )	
				{	a = b;		best=j;		}
				break;
			}
		}
		if( best != i )	{
			ASSIGN(s,w[i]);			ASSIGN(w[i],w[best]);		ASSIGN(w[best],s);	
		}
	}
END:
	return;
}


/*
	v0.1	cys
		8/4/2007
	v0.2	cys		取消GSE_INVERT_MF_1
		10/24/2007
*/
void *GSE_invert( int dim,int*ptr,int *ind,GRUS_FLOATING *val,int ivt_mode,int val_type,int *ret )	{
	int info=0x0;
	clock_t start=clock( );
	void *hLU=GRUS_NULL;

	ivt_mode=GSE_INVERT_GSS;
	*ret=GSE_OK;
	switch( ivt_mode )	{
	case GSE_INVERT_GSS:
		hLU = GSS_get_solver( dim,ptr,ind,val,&info );
/*		ASSERT( val_type==GSE_DATA_DOUBLE );
		memset( GRUS_MF_CONTROL,0x0,sizeof(double)*GRUS_MF_CONTROL_ITEM );
		memset( GRUS_MF_INFO,0x0,sizeof(double)*GRUS_MF_INFO_ITEM );
		info = mf_G_InitDefault_udi( dim,dim,ptr,ind,val,GRUS_MF_UNI,GRUS_MF_AUTO );
		if( info != GRUS_OK )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}		
		hLU = mf_symbol_udi( dim,dim,ptr,ind,val );	
		if( hLU == NULL )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}		
		mf_numeric_udi( val,hLU );
		if( GRUS_MF_INFO[GRUS_MF_STATUS] != GRUS_OK )	
		{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}
		info = 0x0;*/
		break;
	case GSE_INVERT_MF_1:
/*		ASSERT( val_type==GSE_DATA_COMPLEX );
		hLU = mf_symbol_1( dim,ptr[dim],ptr,ind );			
		if( hLU == NULL )	
		{	*ret=GSE_INVERT_FAIL;		goto END;		}
		mf_numeric_1( dim,val,hLU,info );*/		
		break;
	default:		
		break;
	}	

END:
	*ret = info;

	GRUS_SE_INFO[GRUS_SE_INVERT_TIME] += clock( )-start;
	return hLU;
}

/*
	z=spectrum transform(v)

	v0.1	cys
		8/4/2007
*/
void GSE_solve( void *hLU,int dim,int *ptr,int *ind,GRUS_FLOATING *val,GRUS_FLOATING *z,int ivt_mode,int code,int *ret )	{
	int info=0;
/*	GRUS_FLOATING s,t,*val,*z_1,c;
	double *d_val,*x=d_temp,d,res;
	GRUS_FLOATING *b=work+dim;*/
	clock_t start=clock( );

	switch( ivt_mode )	{
	case GSE_INVERT_GSS:
		mf_solve( hLU,dim,ptr,ind,val,z );
/*		BLAS_COPY( dim,v,inc_1,z,inc_1 );
		if( code == GRUS_DRIVER_OPX )	{
			mf_G_solve( hLU,z );
		}else if( GRUS_DRIVER_OPTX )	{
			ASSERT( 0 );
		}*/
		break;
	case GSE_INVERT_MF_1:
/*		if( code == GRUS_DRIVER_OPX )	{
//			mf_solve_1( BLAS_NOTRANS,dim,hLU,v,z,info );
//			mf_solve_1( BLAS_N,dim,hLU,v,z,info );
		}else if( GRUS_DRIVER_OPTX )	{
//			mf_solve_1( BLAS_TRANS,dim,hLU,v,z,info );
//			mf_solve_1( BLAS_T,dim,hLU,v,z,info );
		}*/
		*ret = info;
		break;
	default:		//z=Av
		break;
	}	
	
	GRUS_SE_INFO[GRUS_SE_MV_TIME] += clock( )-start;
	GRUS_SE_INFO[GRUS_SE_MV_COUNT] ++;
}

/*
	v0.1	cys
		8/6/2007
*/
void GSE_clear( void *hLU,int ivt_mode )	{
	switch( ivt_mode )	{
	case GSE_INVERT_GSS:
		GSS_clear_solver( hLU );
		break;
	case GSE_INVERT_MF_1:
//		mf_clear_1( hLU );
		break;
	default:
		break;
	}
}

/*
	H[n,ldh]	Q[n,ldq]
	w[t:K]中包含两部分，select标记的部分排到前面，剩余的按大小排序

	与ztrsen不同，需要排序
	reorders the Schur factorization of a matrix

	v0.1	cys
		8/10/2007
*/
int GSE_MS_sort( int n,int K,int t,GRUS_FLOATING *T,int ldt,GRUS_FLOATING *Q,int ldq,GRUS_FLOATING *w,int *select,int flag )	{
	int i,j,nDo=t,info,ret=0,ifst,ilst;
	double a,b;
	char N[1]={'N'},V[1]={'V'};

	for( i = nDo; i < n; i++ )	{
		if( select[i]==0 )		continue;
		ifst=i+1;		ilst=nDo+1;	
		if( ifst!=ilst )	{
			ztrexc ( V,&n,T,&ldt,Q,&ldq,&ifst,&ilst,&info );
			ASSERT( info==0 );
		}
		nDo++;
	}
	for( i = nDo; i < K; i++ )	{
		ABS( a,T[i*ldt+i] );
		ifst=i+1;		ilst=i+1;	
		for( j = i+1; j < n; j++ )	{
			ABS( b,T[j*ldt+j] );
			if( b>a )
			{	a=b;	ifst=j+1;	}
		}			
		if( ifst!=ilst )	{
			ztrexc ( V,&n,T,&ldt,Q,&ldq,&ifst,&ilst,&info );
			ASSERT( info==0 );
		}
		nDo++;
	}
	for( i = 0; i < n; i++ )	ASSIGN( w[i],T[i*ldt+i] );	
END:
	ASSERT( nDo==K );
	ASSERT( ZV_select_verify(n,w,K,0x0)==0 );
	ret=nDo-K;

	return ret;
}

/*
	work[2*n],rwork[n]
	由于[1:t,1:t]被LOCK，注意求Q的技巧

	v0.1	cys
		8/10/2007
*/
int GSE_MG_evs( int n,int t,GRUS_FLOATING *M,int ldm,GRUS_FLOATING * w,GRUS_FLOATING *Q,int ldq,
			   GRUS_FLOATING *ev,int ldev,GRUS_FLOATING * work,double* rwork,int *flag )	{
	char N[1]={'N'},V[1]={'V'},R[1]={'R'},B[1]={'B'};
	int i,*select=GRUS_NULL,sdim,lwork,*bwork=GRUS_NULL,ldvl=n,mm,m,info;
	GRUS_FLOATING *vl=GRUS_NULL,s,*w_old;
	double norm;

	lwork=2*n;
/*	if( t > 0 )	{
		w_old = work+lwork;
		for( i = 0; i < t; i++ )	{
			ASSIGN( w_old[i],M[i*ldm+i] );
		}
	}*/
	memset( Q,0x0,sizeof(GRUS_FLOATING)*n*ldq );
//	zgees( V,N, select,&n,M,&ldm,&sdim,w,Q,&ldq,work,&lwork,rwork,bwork,&info);
//	ASSERT( info==0 && sdim==0 );
	n -= t;
	zgees( V,N, select,&n,M+t*ldm+t,&ldm,&sdim,w+t,Q+t*ldq+t,&ldq,work,&lwork,rwork,bwork,&info );
	ASSERT( info==0 && sdim==0 );
	for( i = 0; i < t; i++ )	ASSIGN(	Q[i*ldq+i],one_ );
	n+=t;
//	if( t > 0 )
//		ASSERT(	ZU_verify( n,n,Q,ldq,0x0 )==0x0 );

	mm=n;		
	ASSERT( ldq==ldev );
	memcpy( ev,Q,sizeof(GRUS_FLOATING)*n*ldq );
	ztrevc( R,B,select,&n,M,&ldm,vl,&ldvl,ev,&ldev,&mm,&m,work,rwork,&info );
	ASSERT( info==0 && m==n );

	for( i = 0; i < n; i ++ )	{
		norm = BLAS_NRM2( n,ev+ldev*i,inc_1 );			
		if( norm != 1.0 )	{
			SET(s,1.0/norm,0.0);
			BLAS_SCAL( n,s,ev+ldev*i,inc_1 );
		}
	}

	*flag=0x0;
	return *flag;
}

/*
	[REF]:	Balancing Sparse Matrices for Computing Eigenvalues
	v0.1	cys
		8/17/2007
*/
int GSE_Balance_1( int dim,int loop,double *D,double *d_temp,GRUS_FLOATING *work,int flag )	{
	int i,j,ret=GSE_OK,idist=2,iseed[4]={1,3,5,7};
	double *x=d_temp,a,d_max=0.0,d_min=1.0e100,norm;
	GRUS_FLOATING *p=work,*z=p+dim,s;

	for( j = 0; j < dim; j++ ) D[j]=1.0;
	for( i = 0; i < loop; i++ )	{
		for( j = 0; j < 3; j++ )	iseed[i] *=(i+1);
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
	v0.1	cys
		9/3/2007
*/
void GSE_driver_core( void* hSolver,void **u_0,void **u_1,int *code,int flag  )	{
	GRUS_driver_core( hSolver,u_0,u_1,code,flag  );

	return;
}

