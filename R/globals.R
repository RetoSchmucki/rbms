## column names used inside data.table expressions (non-standard evaluation), declared to avoid the
## "no visible binding for global variable" notes of R CMD check

utils::globalVariables(c(
  ".", "..available_cols", "..colnames",
  "ANCHOR", "ANCHOR_PER_WEEK", "BOOTi", "COL_INDEX", "COMPLT_SEASON", "COUNT", "C_SY",
  "DATE", "DATE2", "DAY", "DAY_SINCE", "DISPERSION", "FITTED", "GLM_LEVEL", "GROUP",
  "IMPUTED_COUNT", "LOG_EFFECT", "MONTH", "M_SEASON", "M_YEAR", "M_YEAR_fac", "N", "NM",
  "NMnew", "NSITE", "NSITE_B", "NSITE_OBS", "N_SY", "SE", "SINDEX", "SINDEX_GLM", "SITE_ID",
  "SITE_SUM", "SITE_YEAR_SUM", "SPECIES", "START_END", "STRATA", "SVALUE", "S_IDX",
  "TOTAL_COUNT", "TOTAL_NM", "V1", "WDAY", "WEEK", "WEEK_COUNT", "WEEK_DAY", "WEEK_NM",
  "WEEK_SINCE", "YEAR", "Y_IDX", "ZERO", "maxCOUNT", "meanCOUNT", "occurN", "tmpV1",
  "trimDAYNO", "trimWEEKNO", "uID", "visitN", "RESP", "YEAR_F", "SITE_F", "weights",
  "HABITAT", "HAB_F", "N_HAB", "N_SITE", "WEIGHT", "S_EXP", "INDEX", "U",
  "SITE_LENGTH", "LOG_LENGTH", "N_L"
))
