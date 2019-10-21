/**************************************************************************************************/
/* PROGRAM: 	Extracting the donor and candidate data for the liver offer acceptance model      */
/* AUTHOR(S): 	Andrew Wey                                                                        */
/* CREATED: 	07/20/2018                                                                        */
/* MODIFY:                                                                                        */
/* FUNCTION: 	This file is part of the data exporting process and is liver-specific.			  */
/*																								  */
/**************************************************************************************************/



libname allorg "SAF_DIRECTORY\CURRENT_SAF\allorg" ACCESS=READONLY; 
libname kipa "SAF_DIRECTORY\CURRENT_SAF\kipa" ACCESS=READONLY; 
libname liin "SAF_DIRECTORY\CURRENT_SAF\liin" ACCESS=READONLY; 
libname thor "SAF_DIRECTORY\CURRENT_SAF\thor" ACCESS=READONLY; 
libname srtrexp "SAF_DIRECTORY\CURRENT_SAF\srtrexp" ACCESS=READONLY; 
libname srtrfmt 'SAF_DIRECTORY\CURRENT_SAF\srtrfmt' ACCESS=READONLY;
options fmtsearch = (srtrfmt.formats);




data work.meld_hist_tmp1;
	set liin.STATHIST;
	where wl_org = "LI";
	keep PX_ID CANHX_BEGIN_DT CANHX_END_DT CAN_LISTING_DT CANHX_ALBUMIN CANHX_ASCITES CANHX_BILI CANHX_BILI_PPC CANHX_INR CANHX_STAT_CD CANHX_SRTR_LAB_MELD 
			CANHX_SRTR_LAB_MELD_TY CANHX_OPTN_LAB_MELD CANHX_OPTN_LAB_MELD_TY CANHX_EXC_DIAG_HCC1 CANHX_EXC_DIAG_HCC2 CANHX_EXC_DIAG_HCC_NOPOLICY 
			CANHX_EXC_DIAG_OTHER CANHX_EXC_FLG CANHX_SERUM_SODIUM CANHX_SERUM_CREAT;
	run;


proc datasets;
	contents data = meld_hist_tmp1;
		modify meld_hist_tmp1;
			format CANHX_BEGIN_DT CANHX_END_DT MMDDYY10.;
	contents data = meld_hist_tmp1;
	run;
quit;


PROC EXPORT DATA= WORK.meld_hist_tmp1 
            OUTFILE= "DATA_DIR\li_alloc_meld_CURRENT_SAF.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



