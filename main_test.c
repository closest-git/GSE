#define _G_COMPLEX_
#include "grus_version.h"
#undef _GRUS_STATIC_
#include "GSE_Krylov.h"
#include <stdio.h>
#include <math.h>
#include "grus_debug.h"
#include "grus_def.h"

#define PI 3.1415926535897932384626433832795

static int diff_dim,*diff_flag,*diff_p,test_diff=0;
/*
	A -= shift*I

	假设	val_type==DoubleComplex
		
	v0.1	cys
		8/4/2007
*/
int SXL_D2Z_shift( int nCol,int *A_ptr,int *A_ind,double *A_val,int **AS_ptr,int **AS_ind,DoubleComplex **AS_val,double shift[2] )	{
	int i,j,nRet=GSE_OK,isFind,nX=0,nnz=A_ptr[nCol],nz,*ptr,*ind;
	double real=shift[0],imag=shift[1];
	GRUS_FLOATING *val;

	for( i = 0; i < nCol; i++ )	{
		isFind = 0;
		for( j = A_ptr[i]; j < A_ptr[i+1]; j++ )	{
			if( A_ind[j]==i )
				isFind=1;
		}
		if( isFind == 0 )
			nX++;
	}

	val=GRUS_alloc( sizeof(GRUS_FLOATING)*(nnz+nX) );		
	ptr=GRUS_alloc( sizeof(int)*(nCol+1) );		
	ind=GRUS_alloc( sizeof(int)*(nnz+nX) );		

	nz=0;
	ptr[0]=0;
	for( i = 0; i < nCol; i++ )	{
		isFind=0;
		if( test_diff==1 )	{
			if( diff_flag[i]==0 )
			{	real=shift[0],	imag=shift[1];	}
			else
			{	real=0.0;		imag=0.0;		}
		}
		for( j = A_ptr[i]; j < A_ptr[i+1]; j++ )	{
			if( A_ind[j]!=i )		{
				if( A_ind[j]>i && isFind==0 )	{
					SET( val[nz],-real,-imag );
					ind[nz]=i;		isFind=1;
					nz++;
				}
				ind[nz]=A_ind[j];
				SET( val[nz],A_val[j],0.0 );
			}else	{
				isFind=1;
				ind[nz]=A_ind[j];
				SET( val[nz],A_val[j]-real,-imag );
			}		
			nz++;
		}
		if( isFind == 0 )	{
			SET( val[nz],-real,-imag );
			ind[nz]=i;
			nz++;
		}
		ptr[i+1]=nz;
	}
//	ASSERT( nnz+nX==nz );
END:
	*AS_ptr=ptr;		*AS_ind=ind;		*AS_val=val;
	return nRet;
}

/*
	v0.1	cys
		8/5/2007
*/
void SXL_ouput( int no,int n,int dim,DoubleComplex *lenda,double shift[2],int nOPx )	{
	int i,j,ret;
	double a,b,real,imag,t;
	char sPath[80];
	FILE *fp;

	sprintf( sPath,"F:\\IEEE\\ORIGIN\\GSE_%d_out.eig",dim );
	if( no==0 )
		fp=fopen( sPath,"w" );
	else
		fp=fopen( sPath,"a" );
	ret = fseek( fp,0x0,SEEK_END );
	fprintf( fp,"\r\n*%2d位移：     %.18lf\t%.18lf\t%d\r\n",no,shift[0],shift[1],nOPx );
	for( i = 0; i < n; i++ )	{		
		real = REAL_COMPONENT(lenda[i]);		imag=IMAG_COMPONENT(lenda[i]);
	/*	a = sqrt( real*real+imag*imag );
		for( j = i+1; j < n; j++ )	{
			real = w[2*j]-shift[0];		imag=w[2*j+1]-shift[1];
			b = sqrt( real*real+imag*imag );
			if( a > b )	{
				t=w[2*i];		w[2*i]=w[2*j];			w[2*j]=t;
				t=w[2*i+1];		w[2*i+1]=w[2*j+1];		w[2*j+1]=t;
				a = b;
			}
		}*/
		fprintf( fp,"\t%.18lf\t%.18lf\t%d\r\n",real,imag,i );
	}
	fclose( fp );

}

/*
	v0.1	cys
		8/8/2007
*/
double Eigen_Defect( int nCol,int* ptr,int* ind,double *val,
			DoubleComplex lenda,DoubleComplex *ev,DoubleComplex *z )	{
	int i,j;
	double res=0.0,a;
	DoubleComplex s;

	ASSERT( test_diff==0 );
	for( i = 0; i < nCol; i++ )		CLEAR(z[i] );
	for( j = 0; j < nCol; j++ )		
		MULT_SUB( z[j],lenda,ev[j] );	
	for( i = 0; i < nCol; i++ )	{		
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			SET( s,val[j],0.0 );
			MULT_ADD( z[ind[j]],s,ev[i] );
		}
	}
	res = sqrt(res);
	for( j = 0; j < nCol; j++ )			{	
		ABS( a,z[j] );	
		res += a*a;
	}

	res = sqrt(res);

	return res;			
}

/*
	也可拓展于shift-invert 模式

	v0.1	cys
		8/9/2007
*/
double SXL_Eigen_Defect( int nCol,void *hLU,double shift[2],DoubleComplex lenda,DoubleComplex *ev,DoubleComplex *z )	{
	int i,j,cur,info;
	double res=0.0,a;
	DoubleComplex s;

	memset( z,0x0,sizeof(DoubleComplex)*nCol );
	for( j = 0; j < diff_dim; j++ )		{
		cur = test_diff==1 ? diff_p[j] : j;
		ASSIGN( z[cur],ev[j] );
	}
	GSE_solve( hLU,nCol,GSE_INVERT_MF_1,z,z,GRUS_DRIVER_OPX,&info );

	SET(s,REAL_COMPONENT(lenda)-shift[0],IMAG_COMPONENT(lenda)-shift[1] );
	RECIPROCAL( s );

	for( j = 0; j < diff_dim; j++ )		{
		cur = test_diff==1 ? diff_p[j] : j;
		MULT_SUB( z[cur],s,ev[j] );
		ABS( a,z[cur] );	
		res += a*a;
	}

	res = sqrt(res);

	return res;			
}

/*		
		v0.1	cys
		8/2/2007
*/
void GSE_test( int nCol,int *ptr,int *ind,double *val,int flag )	{
	int i,j,k,code,*AS_ptr,*AS_ind,AS_nz,info;
	int nEV,isAlwaysInit=1,nOPx=0,nOPx_total=0;
	double param[GSE_PARAM_MOST],a_norm,res,real,imag;
	double shift[15][2] = {
		{0.0,0.628000020980834960},{0.0,1.656285762786865200},{0.0,2.684571504592895500},
		{0.0,3.712857246398925800},{0.0,4.741142988204956100},{0.0,5.769428730010986300},
		{0.0,6.797714471817016600},{0.0,7.826000213623046900},{0.0,8.854285955429077100},
		{0.0,9.882571697235107400},{0.0,10.910857439041138000},{0.0,11.939143180847168000},
		{0.0,12.967428922653198000},{0.0,13.995714664459229000},{0.0,15.024000167846680000}
	};
	void *hSolver,*hLU;	
	DoubleComplex *u_1,*u_0,*AS_val,*z,*work,*lenda;
	FILE *fp_diff;
//	char gse_put[80]={"set KMP_DUPLICATE_LIB_OK=TRUE"};
//	_putenv( gse_put );

//	GSE_QR_main( nCol,ptr,ind,val,0,shift );
//	return;
	if( test_diff==1 )	{
		diff_flag = malloc( sizeof(int)*nCol*2 );		diff_p=diff_flag+nCol;
		fp_diff = fopen( "F:\\business\\SONG\\diff_flag.info","rb" );
		fread( diff_flag,sizeof(int),nCol,fp_diff );
		fclose( fp_diff );
		diff_dim=0;
		for( i = 0 ; i < nCol; i++ )	{
			if( diff_flag[i]==0 )	{
				diff_p[diff_dim]=i;
				diff_dim++;
			}
		}
	}
//必须赋值的参数
	param[GSE_METHOD] = GSE_METHOD_KSCHUR;			//ARE;RESK
	param[GSE_PARAM_NEV] = 1;
	param[GSE_PARAM_NGV] = 10;
	param[GSE_PARAM_MAX_RESTART]=10;
	param[GSE_STOP_MODE] = GSE_STOP_KRYLOV;					//BLOCK,DIRECT;
//可选的参数
	param[GSE_REORTHO_THRESH]=0.1;
	param[GSE_CONVERGE_TOL]=1.0e-15;

	param[GSE_RESTART_MODE] = GSE_RESTART_SORENSEN;				//CHEBYCHEV_0,V1;
	param[GSE_DATA_TYPE] = GSE_DATA_DOUBLE;						//GSE_DATA_COMPLEX;
	param[GSE_SHIFT_MODE] =  GSE_SHIFT_NORMAL;					//GSE_SHIFT_MA;GSE_SHIFT_NO;
	param[GSE_INVERT_MODE] = GSE_INVERT_MF_1;					//GSE_INVERT_NO;GSE_INVERT_GSS

	printf( "\tdim=%d,nnz=%d,nEV=%g,nGV=%g\r\n",
		nCol,ptr[nCol],param[GSE_PARAM_NEV],param[GSE_PARAM_NGV] );
	printf( "\tinvert mode=%d\tshift mode=%d\tstop mode=%d\r\n",
		param[GSE_INVERT_MODE],param[GSE_SHIFT_MODE],param[GSE_STOP_MODE] );

	nEV = param[GSE_PARAM_NEV];
	work = malloc( sizeof(DoubleComplex)*nCol*(nEV+1) );		z=work+nCol;
	lenda = malloc( sizeof(DoubleComplex)*nEV );

	a_norm = CCS_DLA_F( nCol,ptr,ind,val );
	for( i = 0; i < 1; i++ )	{
		printf( "\tshift=(%.12lf,%.12lf)\r\n",shift[i][0],shift[i][1] );
		SXL_D2Z_shift( nCol,ptr,ind,val,&AS_ptr,&AS_ind,&AS_val,shift[i] );
		if( test_diff==0 )	diff_dim=nCol;
		if( i == 0 || isAlwaysInit==1 )
			hSolver = GSE_Krylov_Init( diff_dim,param,shift[i],1 );
		else
			GSE_Krylov_Shift( hSolver,shift[i],0x0 );
		hLU = GSE_invert( nCol,AS_ptr,AS_ind,AS_val,0x0,GSE_DATA_COMPLEX,&info );
		if( hLU==0x0 )
			break;
		nOPx = 0;
		code=GRUS_DRIVER_START;
		while( code != GRUS_DRIVER_FINISH )	{			
//			KSchur_core( hSolver,&u_0,&u_1,&code,0x0	 );
			GSE_driver_core( hSolver,&u_0,&u_1,&code,0x0 );
			switch( code )	{
			case GRUS_DRIVER_OPX:
			case GRUS_DRIVER_OPTX:
				nOPx++;
				memset( z,0x0,sizeof(DoubleComplex)*nCol );
				if( test_diff==1 )	
				{	for( j = 0; j < diff_dim; j++ )		ASSIGN( z[diff_p[j]],u_0[j] );	}
				else
					memcpy( z,u_0,sizeof(DoubleComplex)*nCol );
				GSE_solve( hLU,nCol,GSE_INVERT_MF_1,z,z,code,&info );
				if( test_diff==1 )	
				{	for( j = 0; j < diff_dim; j++ )		ASSIGN( u_1[j],z[diff_p[j]] );	}
				else
					memcpy( u_1,z,sizeof(DoubleComplex)*nCol );
				break;				
			case GRUS_DRIVER_AX:
				for( j = 0; j < nCol; j++ )		CLEAR(u_1[j] );
				for( j = 0; j < nCol; j++ )	{		
					for( k = AS_ptr[j]; k < AS_ptr[j+1]; k++ )	{
						MULT_ADD( u_1[ AS_ind[k]],u_0[j],AS_val[k] );
					}
				}				
				break;
			case GRUS_DRIVER_ATX:
				for( j = 0; j < nCol; j++ )		CLEAR(u_1[j] );
				for( j = 0; j < nCol; j++ )	{		
					for( k = AS_ptr[j]; k < AS_ptr[j+1]; k++ )	{
						MULT_ADD( u_1[j],u_0[ AS_ind[k]],AS_val[k] );
					}
				}				
				break;
			default:
				break;
			}
		}
		nOPx_total += nOPx;
		flag = GSE_EIGEN_VALUE | GSE_RIGHT_EIGENVECTOR;
		GSE_Krylov_Post( hSolver,lenda,z,&nEV,&flag );			ASSERT( flag==0 );
		printf( "\r\n--NO--------------------LENDA-----------|Ax-lenda*x|----|lenda-shift|--\r\n" );
		for( j = 0; j < nEV; j ++ )	{
//			res = test_diff ? -1 : Eigen_Defect( nCol,ptr,ind,val,lenda[j],z+j*nCol,work );
			res = SXL_Eigen_Defect( nCol,hLU,shift[i],lenda[j],z+j*diff_dim,work );
			real=REAL_COMPONENT(lenda[j]),		imag=IMAG_COMPONENT(lenda[j]);
			printf( "%3d: (%.12lf,%.12lf)\t%g",j+1,real,imag,res );
			real-=shift[i][0];		imag-=shift[i][1];		
			printf( "\t%g\r\n",sqrt( real*real+imag*imag ) );
		}		

		GSE_Krylov_Dump( hSolver,0x0 );
		SXL_ouput( i,nEV,nCol,lenda,shift[i],nOPx );
		GSE_clear( hLU,GSE_INVERT_MF_1 );
		if( isAlwaysInit==1 )
			GSE_Krylov_Clear( hSolver );
		free(AS_ptr);		free(AS_ind);		free(AS_val);
		printf( "\r\n" );
	}
	free( work );		free(lenda);		
END:
	if( isAlwaysInit==0 )
		GSE_Krylov_Clear( hSolver );
	printf( "\r\n-nOP*x------------------------LENDA----------------------------res-------\r\n" );
	printf( "%8d\r\n",nOPx_total );
	if( test_diff==1 )
		free( diff_flag );
}
