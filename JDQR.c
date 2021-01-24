/*
	The algorithm is from "JACOBIDAVIDSON ALGORITHMS FOR VARIOUS EIGENPROBLEMS"

	v0.1	cys
		5/28/2007
*/
int JDQR_core( ARE_R_SOLVER *hJDQR,int m )	{
	int dim=hJDQR->dim,n=m,info=0,i,j;
	GRUS_FLOATING *t=hJDQR->V+m*dim,*u,*r,*M,s,es_u,es_r,lenda,sft,off,
		*w=hJDQR->work,*h=w+dim,*x_r=w,*x_u=h,epsi,*c,*val= hJDQR->val;
	double norm,res,a_1,a_2;

	if( n > 0 )	{
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,alpha_,hJDQR->V,dim,t,inc_1,gama_,h,inc_1 );	//h=V't
		BLAS_GEMV( BLAS_NOTRANS,dim,n,beta_,hJDQR->V,dim,h,inc_1,alpha_,t,inc_1 );	
	//DGKS correction]
		c = hJDQR->work;
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,alpha_,hJDQR->V,dim,t,inc_1,gama_,c,inc_1 );	//h=-V'z
		BLAS_GEMV( BLAS_NOTRANS,dim,n,beta_,hJDQR->V,dim,c,inc_1,alpha_,t,inc_1 );	
//		BLAS_AXPY( n,alpha_,c,1,h,1 );	
		norm = BLAS_NRM2( dim,t,1 );	//normalize z
		SET(s,1.0/norm,0.0);
		BLAS_SCAL( dim,s,t,1 );
	}
	ARE_AV_0( hJDQR,t,w );
//	ai_expand( hJDQR,m,w );
	M=hJDQR->H+m*hJDQR->dimH;
	n = n+1;
	BLAS_GEMV( BLAS_CONJTRANS,dim,n,alpha_,hJDQR->V,dim,w,inc_1,gama_,M,inc_1 );	//M=V'AV
	hJDQR->V_m = m+1;

	GSE_verify_V( hJDQR );
	GSE_extraction( hJDQR,GSE_MAT_HERMITIAN );
	ASSIGN( lenda,hJDQR->w_H[0] );
//	ARE_trans_lenda( hJDQR,&lenda );
	u = hJDQR->rv_A+hJDQR->nCurEV*dim;
	ARE_AV_0( hJDQR,u,x_r );						
	BLAS_SCAL( dim,beta_,x_r,1 );
	BLAS_AXPY( dim,lenda,u,1,x_r,1 );	//-r=lenda*u-Au
	res = BLAS_NRM2( dim,x_r,1 );	
	if( res < mch_eps*100*hJDQR->a_norm )
		return GSE_JDQR_CONVERGE;
//solve t from JD correction equation
	t=hJDQR->V+(m+1)*dim;
	if( 0 )	{
	ASSIGN( sft,lenda );
	SUB( off,sft,hJDQR->sft_old );
	ABS( a_1,off );			ABS( a_2,sft );
	if( a_1 > a_2*0.5 )		{
		for( i = 0; i < dim; i++ )	{
			for( j = hJDQR->ptr[i]; j < hJDQR->ptr[i+1]; j++ )	
				DECREMENT( val[j],off );	
		}
		mf_numeric_1( dim,val,hJDQR->hLU,info );		
		ASSIGN( hJDQR->sft_old,sft );
	}
		mf_solve_1( dim,hJDQR->hLU,u,x_u,info );
		mf_solve_1( dim,hJDQR->hLU,x_r,x_r,info );
		BLAS_DOT_c( dim,u,1,x_u,1,es_u );
		BLAS_DOT_c( dim,u,1,x_r,1,es_r );
		DIV( epsi,es_r,es_u );
		BLAS_SCAL( dim,beta_,x_r,1 );
		BLAS_AXPY( dim,epsi,x_u,1,x_r,1 );
		BLAS_COPY( dim,x_r,1,t,1 );		//t=epsi*M'*u-M'*r
		BLAS_DOT_c( dim,t,1,u,1,s );
	}else	{
		mf_solve_1( dim,hJDQR->hLU,u,t,info );
	}

	return GSE_JDQR_LOOP;
}

/*
	The algorithm is from "JACOBIDAVIDSON ALGORITHMS FOR VARIOUS EIGENPROBLEMS"

	v0.1	cys
		5/28/2007

int JDQR_core( ARE_R_SOLVER *hJDQR,int m )	{
	int dim=hJDQR->dim,n=m,k=hJDQR->nEV,info=0,i,j;
	GRUS_FLOATING *t=hJDQR->V+m*dim,*u,*r,*M,s,es_u,es_r,lenda,sft,off,tao,
		*w,*h=w+dim,*x_r=w,*x_u=h,epsi,*c,*val= hJDQR->val;
	double norm,res,a_1,a_2;

	w = hJDQR->W+m*dim;
	if( n > 0 )	{
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,alpha_,hJDQR->V,dim,t,inc_1,gama_,h,inc_1 );	//h=V't
		BLAS_GEMV( BLAS_NOTRANS,dim,n,beta_,hJDQR->V,dim,h,inc_1,alpha_,t,inc_1 );	
	//DGKS correction]
		c = hJDQR->work;
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,alpha_,hJDQR->V,dim,t,inc_1,gama_,c,inc_1 );	//h=-V'z
		BLAS_GEMV( BLAS_NOTRANS,dim,n,beta_,hJDQR->V,dim,c,inc_1,alpha_,t,inc_1 );	
//		BLAS_AXPY( n,alpha_,c,1,h,1 );	
		norm = BLAS_NRM2( dim,t,1 );	//normalize z
		SET(s,1.0/norm,0.0);
		BLAS_SCAL( dim,s,t,1 );
	}
	ARE_AV_0( hJDQR,t,w );
	MULT( s,tao,beta );
	BLAS_AXPY( dim,s,t,1,w,1 );		//w=A*t-tao*t
	if( k > 0 )	{	//w=w-Q(Q'w)
		BLAS_GEMV( BLAS_CONJTRANS,dim,k,alpha_,hJDQR->rv_A,dim,w,inc_1,0,h,inc_1 );	
		BLAS_GEMV( BLAS_NOTRANS,dim,n,beta_,hJDQR->rv_A,dim,h,inc_1,alpha_,w,inc_1 );	
	}
	if( n > 0 )	{	//M(*,m)=W'w;		w=w-W*M(*,m)
		MA=hJDQR->MA+m*ldm;
		BLAS_GEMV( BLAS_CONJTRANS,dim,k,alpha_,hJDQR->W,dim,w,inc_1,0,MA,inc_1 );	
		BLAS_GEMV( BLAS_NOTRANS,dim,n,beta_,hJDQR->W,dim,MA,inc_1,alpha_,w,inc_1 );	
	}
	norm=BLAS_NRM2( dim,w,1 );
	SET(hJDQR->MA+m*ldm+m,norm,0.0);
	SET(s,1.0/norm,0.0);			BLAS_SCAL( dim,s,w,1 );		//normalise w
//(3)
	if( n > 0 )	{		//M(*,m)=;		M(m,*)=
		M=hJDQR->M+m*ldm;
		BLAS_GEMV( BLAS_CONJTRANS,dim,n,alpha_,hJDQR->W,dim,v,inc_1,gama_,M,inc_1 );	
		M=hJDQR->M+m;
		BLAS_GEMV( BLAS_CONJTRANS,1,m,alpha_,w,dim,hJDQR->V,inc_1,gama_,M,ldm );	
	}
	BLAS_DOT_u( dim,w,1,v,1,hJDQR->MA+m*ldm+m );	//M(m,m)=w'*v
//(4)		QZ decomposition for( MA,M )
	JDQR_QZ( hJDQR );
//(5)
	norm=BLAS_NRM2( dim,r,1 );
	while( norm<epsi )	{
		hJDQR->nEV++;
		if( hJDQR->nEV>=hJDQR->nEV )
			return GSE_JDQR_ALL_CONVERGE;
//(6)
		
	}
//(7)
	if( )	{		//restart
	}

//(8)	solve t from JD correction equation
	t=hJDQR->V+(m+1)*dim;
	if( 0 )	{
	ASSIGN( sft,lenda );
	SUB( off,sft,hJDQR->sft_old );
	ABS( a_1,off );			ABS( a_2,sft );
	if( a_1 > a_2*0.5 )		{
		for( i = 0; i < dim; i++ )	{
			for( j = hJDQR->ptr[i]; j < hJDQR->ptr[i+1]; j++ )	
				DECREMENT( val[j],off );	
		}
		mf_numeric_1( dim,val,hJDQR->hLU,info );		
		ASSIGN( hJDQR->sft_old,sft );
	}
		mf_solve_1( dim,hJDQR->hLU,u,x_u,info );
		mf_solve_1( dim,hJDQR->hLU,x_r,x_r,info );
		BLAS_DOT_c( dim,u,1,x_u,1,es_u );
		BLAS_DOT_c( dim,u,1,x_r,1,es_r );
		DIV( epsi,es_r,es_u );
		BLAS_SCAL( dim,beta_,x_r,1 );
		BLAS_AXPY( dim,epsi,x_u,1,x_r,1 );
		BLAS_COPY( dim,x_r,1,t,1 );		//t=epsi*M'*u-M'*r
		BLAS_DOT_c( dim,t,1,u,1,s );
	}else	{
		mf_solve_1( dim,hJDQR->hLU,u,t,info );
	}

	return GSE_JDQR_LOOP;
}
*/

/*
	×¢Òâ£º
		hAR_R->ivt_modeËÆºõÃ»ÓÃ

	v0.1	cys
		5/2/2007
*/
int JDQR_engine( void* hSolver )	{
	ARE_R_SOLVER *hAR_R=hSolver;
	int i,j,dim=hAR_R->dim,p=hAR_R->nEV,ret=GSE_OK,info=0;
	double resi=0.0,norm;
	GRUS_FLOATING s,*w,*m=hAR_R->H,*v=hAR_R->V,*r=hAR_R->work,xita;

	GSE_random_1( dim,v );

	if( hAR_R->sft_mode != GSE_SHIFT_NO || hAR_R->val_type!=hAR_R->A_type/*&& !IS_ZERO(hAR_R->shift)*/ )
		ARE_shift( hAR_R,0x0 );
	ASSIGN( hAR_R->sft_old,hAR_R->shift );

	ASSERT( hAR_R->val_type==GSE_DATA_COMPLEX );
	hAR_R->hLU = mf_symbol_1( dim,hAR_R->ptr[dim],hAR_R->ptr,hAR_R->ind );			
	if( hAR_R->hLU == NULL )	
	{	info = GRUS_MF_INFO[GRUS_MF_STATUS];		goto END;		}
	mf_numeric_1( dim,hAR_R->val,hAR_R->hLU,info );		

	hAR_R->nCurEV = 0;
	for( i = 0; i < p; i++ )	{
		for( j=0; j < hAR_R->nGV;	j++ )	{
			if( JDQR_core( hSolver,j )==GSE_JDQR_CONVERGE )
				break;
			if( hAR_R->nRestart>=hAR_R->MAX_RESTART )	
				break;
		}
//		if( i+1 < p )
//			ai_deflate_( hAR_R );
		hAR_R->nCurEV++;
	}
	if( hAR_R->ivt_mode > 0 || hAR_R->sft_mode > 0 )	{	//
		for( i = 0; i < p; i++ )
			ARE_trans_lenda( hAR_R,&(hAR_R->w_H[i]) );
	}
//	ARE_verify( hSolver );
END:
	return ret;
}
