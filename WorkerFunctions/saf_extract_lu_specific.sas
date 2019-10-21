/**************************************************************************************************/
/* PROGRAM: 	Extracting the donor and candidate data for the lung survival from listing model  */
/* AUTHOR(S): 	Andrew Wey                                                                        */
/* CREATED: 	08/13/2018                                                                        */
/* MODIFY:                                                                                        */
/* FUNCTION: 	This file is part of the data exporting process and is lung-specific.			  */
/* NOTES:     	Based on the file of the same name from the offer acceptance model build process  */
/*																								  */
/**************************************************************************************************/



libname allorg "\CURRENT_SAF\allorg" ACCESS=READONLY; 
libname kipa "SAF_DIRECTORY\CURRENT_SAF\kipa" ACCESS=READONLY; 
libname liin "SAF_DIRECTORY\CURRENT_SAF\liin" ACCESS=READONLY; 
libname thor "SAF_DIRECTORY\CURRENT_SAF\thor" ACCESS=READONLY; 
libname srtrexp "SAF_DIRECTORY\CURRENT_SAF\srtrexp" ACCESS=READONLY; 
libname srtrfmt 'SAF_DIRECTORY\CURRENT_SAF\srtrfmt' ACCESS=READONLY;
options fmtsearch = (srtrfmt.formats);




data work.alloc_scores;
	set srtrexp.LU_ALLOC_SCORES;
	keep CANHX_LAS_SCORE_CALC PX_ID CANHX_LAS_SCORE_BEGIN_DT;
	run;


proc datasets;
	contents data = alloc_scores;
		modify alloc_scores;
			format CANHX_LAS_SCORE_BEGIN_DT MMDDYY10.;
	contents data = alloc_scores;
	run;
quit;


PROC EXPORT DATA= WORK.alloc_scores 
            OUTFILE= "DATA_DIR\lu_alloc_scores_CURRENT_SAF.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;









data work.las_hist_alloc;
	set srtrexp.las_hist_alloc;
	keep PX_ID CANHX_GROUPING CANHX_CI CANHX_CVP_INT CANHX_VENTILATOR_USE CANHX_CALC_SERUM_CREAT CANHX_FVC_ABSOLUTE CANHX_PCO2 CANHX_SIX_MIN_WALK CANHX_SERUM_CREAT
			CANHX_BEGIN_DT;
	run;


proc datasets;
	contents data = las_hist_alloc;
		modify las_hist_alloc;
			format CANHX_BEGIN_DT MMDDYY10.;
	contents data = las_hist_alloc;
	run;
quit;


PROC EXPORT DATA= WORK.las_hist_alloc 
            OUTFILE= "DATA_DIR\lu_las_hist_alloc_CURRENT_SAF.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;







data work.ped_priority;
	set thor.PEDIATRIC_LUNG_PRIORITY;
	keep PX_ID CAN_PL_PRIORITY CAN_PL_EFFECTIVE_DT CAN_PL_PR_LIST_DT ;
	run;


proc datasets;
	contents data = ped_priority;
		modify ped_priority;
			format CAN_PL_EFFECTIVE_DT CAN_PL_PR_LIST_DT MMDDYY10.;
	contents data = ped_priority;
	run;
quit;


PROC EXPORT DATA= WORK.ped_priority 
            OUTFILE= "DATA_DIR\lu_ped_priority_CURRENT_SAF.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

