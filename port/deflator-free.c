#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

#define ds      sizeof(double)
#define zs      sizeof(doublecomplex)


void 
q(fini_df_eigcg)(struct q(DeflatorEigcg) *d_e, struct Q(State) *s)
{
    if (NULL == s || NULL == d_e)
        return;

    int vmax    = d_e->vmax;

    if (! latmat_c_is_null(&(d_e->V)))           
        q(latmat_c_free)(s, &(d_e->V));
    if (! latmat_c_is_null(&(d_e->tmp_V)))       
        q(latmat_c_free)(s, &(d_e->tmp_V));
    if (! latvec_c_is_null(&(d_e->work_c_1)))    
        q(latvec_c_free)(s, &(d_e->work_c_1));
    if (! latvec_c_is_null(&(d_e->work_c_2)))    
        q(latvec_c_free)(s, &(d_e->work_c_2));
    if (NULL != d_e->T)             
        { q(free)(s, d_e->T, vmax * vmax * zs);     d_e->T = NULL ; }
    if (NULL != d_e->hevals)        
        { q(free)(s, d_e->hevals, vmax * ds);       d_e->hevals = NULL; }
    if (NULL != d_e->hevecs2)       
        { q(free)(s, d_e->hevecs2, vmax * vmax * zs);d_e->hevecs2 = NULL; }
    if (NULL != d_e->zwork)         
        { q(free)(s, d_e->zwork, d_e->lwork * zs);  d_e->zwork = NULL ; 
          d_e->lwork = 0; }
#if defined(HAVE_LAPACK)
    if (NULL != d_e->hevecs1)       
        { q(free)(s, d_e->hevecs1, vmax * vmax * zs);d_e->hevecs1 = NULL ; }
    if (NULL != d_e->tau)           
        { q(free)(s, d_e->tau, vmax * zs);          d_e->tau = NULL ; }
    if (NULL != d_e->rwork)         
        { q(free)(s, d_e->rwork, 3 * vmax * ds);    d_e->rwork = NULL ; }
#elif defined(HAVE_GSL)
    if (NULL != d_e->gsl_T_full)    
        { gsl_matrix_complex_free(d_e->gsl_T_full); d_e->gsl_T_full = NULL ; }
    if (NULL != d_e->gsl_hevecs1)   
        { gsl_matrix_complex_free(d_e->gsl_hevecs1);d_e->gsl_hevecs1 = NULL ; }
    if (NULL != d_e->gsl_hevals1)   
        { gsl_vector_free(d_e->gsl_hevals1);        d_e->gsl_hevals1 = NULL ; }
    if (NULL != d_e->gsl_wkspace1)  
        { gsl_eigen_hermv_free(d_e->gsl_wkspace1);  d_e->gsl_wkspace1 = NULL ; }
    if (NULL != d_e->gsl_T_m1)      
        { gsl_matrix_complex_free(d_e->gsl_T_m1);   d_e->gsl_T_m1 = NULL ; }
    if (NULL != d_e->gsl_hevecs2)   
        { gsl_matrix_complex_free(d_e->gsl_hevecs2);d_e->gsl_hevecs2 = NULL ; }
    if (NULL != d_e->gsl_hevals2)   
        { gsl_vector_free(d_e->gsl_hevals2);        d_e->gsl_hevals2 = NULL ; }
    if (NULL != d_e->gsl_wkspace2)  
        { gsl_eigen_hermv_free(d_e->gsl_wkspace2);  d_e->gsl_wkspace2 = NULL ; }
    if (NULL != d_e->gsl_T_proj)    
        { gsl_matrix_complex_free(d_e->gsl_T_proj); d_e->gsl_T_proj = NULL ; }
    if (NULL != d_e->gsl_hevecs3)   
        { gsl_matrix_complex_free(d_e->gsl_hevecs3);d_e->gsl_hevecs3 = NULL ; }
    if (NULL != d_e->gsl_wkspace3)  
        { gsl_eigen_hermv_free(d_e->gsl_wkspace3);  d_e->gsl_wkspace3 = NULL ; }
    if (NULL != d_e->gsl_QR)        
        { gsl_matrix_complex_free(d_e->gsl_QR);     d_e->gsl_QR = NULL ; }
    if (NULL != d_e->gsl_Q_unpack)  
        { gsl_matrix_complex_free(d_e->gsl_Q_unpack);d_e->gsl_Q_unpack = NULL ; }
    if (NULL != d_e->gsl_tmp_MxS)   
        { gsl_matrix_complex_free(d_e->gsl_tmp_MxS);d_e->gsl_tmp_MxS = NULL ; }
    if (NULL != d_e->gsl_tau)
        { gsl_vector_complex_free(d_e->gsl_tau); d_e->gsl_tau = NULL ; }
    if (NULL != d_e->hevals_select1)
        { q(free)(s, d_e->hevals_select1, vmax * sizeof(d_e->hevals_select1[0])); d_e->hevals_select1 = NULL ; }
    if (NULL != d_e->hevals_select2)
        { q(free)(s, d_e->hevals_select2, vmax * sizeof(d_e->hevals_select2[0])); d_e->hevals_select2 = NULL ; }
#else
#  error "no linear algebra library"
#endif

}


void
q(fini_deflator)(struct Q(Deflator) *df, struct Q(State) *s)
{
    if (NULL == s || NULL == df)
        return;
    
    int umax = df->umax;

    if (! latmat_c_is_null(&(df->U)))        q(latmat_c_free)(s, &(df->U));
    if (! latvec_c_is_null(&(df->work_c_1))) q(latvec_c_free)(s, &(df->work_c_1));
    if (! latvec_c_is_null(&(df->work_c_2))) q(latvec_c_free)(s, &(df->work_c_2));
    if (NULL != df->zwork)
        { q(free)(s, df->zwork, umax * zs);         df->zwork = NULL ; }
    if (NULL != df->H)                          
        { q(free)(s, df->H, umax * umax * zs);      df->H = NULL ; }
    if (NULL != df->H_ev)                       
        { q(free)(s, df->H_ev, umax * umax * zs);   df->H_ev = NULL ; }
    if (NULL != df->hevals)                     
        { q(free)(s, df->hevals, umax * ds);        df->hevals = NULL ; }
    if (NULL != df->C)
        { q(free)(s, df->C, umax * umax * zs);      df->C = NULL ; }
#if defined(HAVE_LAPACK)                        
    if (NULL != df->eig_zwork)
        { q(free)(s, df->eig_zwork, df->lwork * zs);df->eig_zwork = NULL ; }
    if (NULL != df->eig_rwork)
        { q(free)(s, df->eig_rwork, 3 * umax * ds); df->eig_rwork = NULL ; }
#elif defined(HAVE_GSL)
    if (NULL != df->gsl_eig_wkspace)    
        { gsl_eigen_herm_free(df->gsl_eig_wkspace); df->gsl_eig_wkspace = NULL ; }
#else
#  error "no linear algebra library"
#endif
}


void
Q(free_deflator)(struct Q(Deflator) **deflator_ptr)
{
    struct Q(State) *s;

    if (deflator_ptr == 0 || *deflator_ptr == 0)
        return;

    s = (*deflator_ptr)->state;
    struct Q(Deflator) *d = *deflator_ptr;

    BEGIN_TIMING(s);

    q(fini_deflator)(d, s);
    q(free)(s, d, sizeof(struct Q(Deflator)));

    END_TIMING(s, 0, 0, 0);

    *deflator_ptr = 0;
}
