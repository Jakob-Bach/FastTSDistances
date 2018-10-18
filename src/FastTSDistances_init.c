#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _FastTSDistances_averageTimeSeries_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_averageTimeSeriesMult_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_clusterEntropy_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_conditionalEntropy_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_corDist_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_cortFactor_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_cortFactorMult_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_DTWDist_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_DTWDistMult_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_DTWDistSakoeChiba_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_DTWDistSakoeChibaMult_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_EDRDist_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_EDRDistMult_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_EDRDistSakoeChiba_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_EDRDistSakoeChibaMult_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_ERPDist(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_ERPDist_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_ERPDistMult_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_ERPDistSakoeChiba(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_ERPDistSakoeChiba_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_ERPDistSakoeChibaMult_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_fowlkesMallows_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_generalizedDB_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_generalizedDunn_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_iGeneralizedDB_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_l1Dist_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_l1DistMult_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_l2CompCorFactor_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_l2CompCorFactorMult_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_l2Dist_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_l2DistMult_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_l2Norm_fast(SEXP);
extern SEXP _FastTSDistances_lmaxDist_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_lmaxDistMult_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_PAA_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_pairCVIParameters_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_phi_fast(SEXP);
extern SEXP _FastTSDistances_PKurtAA_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_PMaxAA_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_PMedAA_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_PMinAA_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_PSDAA_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_PSkewAA_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_purity_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_randIndex_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_SAX_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_subVectorMean_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_tempCor_fast(SEXP, SEXP);
extern SEXP _FastTSDistances_vanDongen_fast(SEXP, SEXP, SEXP);
extern SEXP _FastTSDistances_vectorCrossDistMat(SEXP, SEXP);
extern SEXP _FastTSDistances_VI_fast(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_FastTSDistances_averageTimeSeries_fast",     (DL_FUNC) &_FastTSDistances_averageTimeSeries_fast,     2},
    {"_FastTSDistances_averageTimeSeriesMult_fast", (DL_FUNC) &_FastTSDistances_averageTimeSeriesMult_fast, 2},
    {"_FastTSDistances_clusterEntropy_fast",        (DL_FUNC) &_FastTSDistances_clusterEntropy_fast,        2},
    {"_FastTSDistances_conditionalEntropy_fast",    (DL_FUNC) &_FastTSDistances_conditionalEntropy_fast,    3},
    {"_FastTSDistances_corDist_fast",               (DL_FUNC) &_FastTSDistances_corDist_fast,               3},
    {"_FastTSDistances_cortFactor_fast",            (DL_FUNC) &_FastTSDistances_cortFactor_fast,            3},
    {"_FastTSDistances_cortFactorMult_fast",        (DL_FUNC) &_FastTSDistances_cortFactorMult_fast,        3},
    {"_FastTSDistances_DTWDist_fast",               (DL_FUNC) &_FastTSDistances_DTWDist_fast,               5},
    {"_FastTSDistances_DTWDistMult_fast",           (DL_FUNC) &_FastTSDistances_DTWDistMult_fast,           5},
    {"_FastTSDistances_DTWDistSakoeChiba_fast",     (DL_FUNC) &_FastTSDistances_DTWDistSakoeChiba_fast,     3},
    {"_FastTSDistances_DTWDistSakoeChibaMult_fast", (DL_FUNC) &_FastTSDistances_DTWDistSakoeChibaMult_fast, 3},
    {"_FastTSDistances_EDRDist_fast",               (DL_FUNC) &_FastTSDistances_EDRDist_fast,               4},
    {"_FastTSDistances_EDRDistMult_fast",           (DL_FUNC) &_FastTSDistances_EDRDistMult_fast,           4},
    {"_FastTSDistances_EDRDistSakoeChiba_fast",     (DL_FUNC) &_FastTSDistances_EDRDistSakoeChiba_fast,     4},
    {"_FastTSDistances_EDRDistSakoeChibaMult_fast", (DL_FUNC) &_FastTSDistances_EDRDistSakoeChibaMult_fast, 4},
    {"_FastTSDistances_ERPDist",                    (DL_FUNC) &_FastTSDistances_ERPDist,                    3},
    {"_FastTSDistances_ERPDist_fast",               (DL_FUNC) &_FastTSDistances_ERPDist_fast,               4},
    {"_FastTSDistances_ERPDistMult_fast",           (DL_FUNC) &_FastTSDistances_ERPDistMult_fast,           4},
    {"_FastTSDistances_ERPDistSakoeChiba",          (DL_FUNC) &_FastTSDistances_ERPDistSakoeChiba,          4},
    {"_FastTSDistances_ERPDistSakoeChiba_fast",     (DL_FUNC) &_FastTSDistances_ERPDistSakoeChiba_fast,     4},
    {"_FastTSDistances_ERPDistSakoeChibaMult_fast", (DL_FUNC) &_FastTSDistances_ERPDistSakoeChibaMult_fast, 4},
    {"_FastTSDistances_fowlkesMallows_fast",        (DL_FUNC) &_FastTSDistances_fowlkesMallows_fast,        2},
    {"_FastTSDistances_generalizedDB_fast",         (DL_FUNC) &_FastTSDistances_generalizedDB_fast,         2},
    {"_FastTSDistances_generalizedDunn_fast",       (DL_FUNC) &_FastTSDistances_generalizedDunn_fast,       2},
    {"_FastTSDistances_iGeneralizedDB_fast",        (DL_FUNC) &_FastTSDistances_iGeneralizedDB_fast,        2},
    {"_FastTSDistances_l1Dist_fast",                (DL_FUNC) &_FastTSDistances_l1Dist_fast,                2},
    {"_FastTSDistances_l1DistMult_fast",            (DL_FUNC) &_FastTSDistances_l1DistMult_fast,            2},
    {"_FastTSDistances_l2CompCorFactor_fast",       (DL_FUNC) &_FastTSDistances_l2CompCorFactor_fast,       2},
    {"_FastTSDistances_l2CompCorFactorMult_fast",   (DL_FUNC) &_FastTSDistances_l2CompCorFactorMult_fast,   2},
    {"_FastTSDistances_l2Dist_fast",                (DL_FUNC) &_FastTSDistances_l2Dist_fast,                4},
    {"_FastTSDistances_l2DistMult_fast",            (DL_FUNC) &_FastTSDistances_l2DistMult_fast,            4},
    {"_FastTSDistances_l2Norm_fast",                (DL_FUNC) &_FastTSDistances_l2Norm_fast,                1},
    {"_FastTSDistances_lmaxDist_fast",              (DL_FUNC) &_FastTSDistances_lmaxDist_fast,              2},
    {"_FastTSDistances_lmaxDistMult_fast",          (DL_FUNC) &_FastTSDistances_lmaxDistMult_fast,          2},
    {"_FastTSDistances_PAA_fast",                   (DL_FUNC) &_FastTSDistances_PAA_fast,                   2},
    {"_FastTSDistances_pairCVIParameters_fast",     (DL_FUNC) &_FastTSDistances_pairCVIParameters_fast,     2},
    {"_FastTSDistances_phi_fast",                   (DL_FUNC) &_FastTSDistances_phi_fast,                   1},
    {"_FastTSDistances_PKurtAA_fast",               (DL_FUNC) &_FastTSDistances_PKurtAA_fast,               4},
    {"_FastTSDistances_PMaxAA_fast",                (DL_FUNC) &_FastTSDistances_PMaxAA_fast,                2},
    {"_FastTSDistances_PMedAA_fast",                (DL_FUNC) &_FastTSDistances_PMedAA_fast,                2},
    {"_FastTSDistances_PMinAA_fast",                (DL_FUNC) &_FastTSDistances_PMinAA_fast,                2},
    {"_FastTSDistances_PSDAA_fast",                 (DL_FUNC) &_FastTSDistances_PSDAA_fast,                 3},
    {"_FastTSDistances_PSkewAA_fast",               (DL_FUNC) &_FastTSDistances_PSkewAA_fast,               3},
    {"_FastTSDistances_purity_fast",                (DL_FUNC) &_FastTSDistances_purity_fast,                2},
    {"_FastTSDistances_randIndex_fast",             (DL_FUNC) &_FastTSDistances_randIndex_fast,             2},
    {"_FastTSDistances_SAX_fast",                   (DL_FUNC) &_FastTSDistances_SAX_fast,                   2},
    {"_FastTSDistances_subVectorMean_fast",         (DL_FUNC) &_FastTSDistances_subVectorMean_fast,         3},
    {"_FastTSDistances_tempCor_fast",               (DL_FUNC) &_FastTSDistances_tempCor_fast,               2},
    {"_FastTSDistances_vanDongen_fast",             (DL_FUNC) &_FastTSDistances_vanDongen_fast,             3},
    {"_FastTSDistances_vectorCrossDistMat",         (DL_FUNC) &_FastTSDistances_vectorCrossDistMat,         2},
    {"_FastTSDistances_VI_fast",                    (DL_FUNC) &_FastTSDistances_VI_fast,                    3},
    {NULL, NULL, 0}
};

void R_init_FastTSDistances(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
