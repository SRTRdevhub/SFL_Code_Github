/**************************************************************************************************/
/* PROGRAM: 	Extracting the donor and candidate data for the heart survival from listing model */
/* AUTHOR(S): 	Andrew Wey                                                                        */
/* CREATED: 	08/13/2018                                                                        */
/* MODIFY:                                                                                        */
/* FUNCTION: 	This file is part of the data exporting process and is heart-specific. 			  */
/*																								  */
/**************************************************************************************************/



libname allorg "SAF_DIRECTORY\CURRENT_SAF\allorg" ACCESS=READONLY; 
libname kipa "SAF_DIRECTORY\CURRENT_SAF\kipa" ACCESS=READONLY; 
libname liin "SAF_DIRECTORY\CURRENT_SAF\liin" ACCESS=READONLY; 
libname thor "SAF_DIRECTORY\CURRENT_SAF\thor" ACCESS=READONLY; 
libname srtrexp "SAF_DIRECTORY\CURRENT_SAF\srtrexp" ACCESS=READONLY; 
libname srtrfmt 'SAF_DIRECTORY\CURRENT_SAF\srtrfmt' ACCESS=READONLY;
options fmtsearch = (srtrfmt.formats);




data work.candstat;
	set thor.STATHIST;
	where wl_org = "HR";
	keep PX_ID CANHX_BEGIN_DT CANHX_END_DT CANHX_STAT_CD;
	run;


proc datasets;
	contents data = candstat;
		modify candstat;
			format CANHX_BEGIN_DT CANHX_END_DT MMDDYY10.;
	contents data = candstat;
	run;
quit;


PROC EXPORT DATA= WORK.candstat 
            OUTFILE= "DATA_DIR\hr_stathist_CURRENT_SAF.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

