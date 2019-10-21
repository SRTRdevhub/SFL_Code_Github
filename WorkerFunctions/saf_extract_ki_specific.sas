/**************************************************************************************************/
/* PROGRAM: 	Extracting the CPRA data for the survival from listing kidney model			      */
/* AUTHOR(S): 	Andrew Wey                                                                        */
/* CREATED: 	03/15/2018                                                                        */
/* MODIFY:                                                                                        */
/* FUNCTION: 	This file is part of the data exporting process and is kidney-specific.			  */
/*																								  */
/**************************************************************************************************/



libname allorg "SAF_DIRECTORY\CURRENT_SAF\allorg" ACCESS=READONLY; 
libname kipa "SAF_DIRECTORY\CURRENT_SAF\kipa" ACCESS=READONLY; 
libname liin "SAF_DIRECTORY\CURRENT_SAF\liin" ACCESS=READONLY; 
libname thor "SAF_DIRECTORY\CURRENT_SAF\thor" ACCESS=READONLY; 
libname srtrexp "SAF_DIRECTORY\CURRENT_SAF\srtrexp" ACCESS=READONLY; 
libname srtrfmt 'SAF_DIRECTORY\CURRENT_SAF\srtrfmt' ACCESS=READONLY;
options fmtsearch = (srtrfmt.formats);





data work.stathist;
	set kipa.stathist;
	keep CAN_LISTING_DT CAN_REM_DT CANHX_BEGIN_DT CANHX_END_DT CANHX_CPRA PX_ID;
	run;


proc datasets;
	contents data = stathist;
		modify stathist;
			format CAN_LISTING_DT CAN_REM_DT CANHX_BEGIN_DT CANHX_END_DT MMDDYY10.;
	contents data = stathist;
	run;
quit;



PROC EXPORT DATA= WORK.stathist 
            OUTFILE= "DATA_DIR\ki_stathist.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




data work.ersd_dialysis_dates;
	set kipa.cand;
	keep PX_ID PERS_ESRD_FIRST_SERVICE_DT PERS_ESRD_FIRST_DIAL_DT;
	run;


proc datasets;
	contents data = ersd_dialysis_dates;
		modify ersd_dialysis_dates;
			format PERS_ESRD_FIRST_SERVICE_DT PERS_ESRD_FIRST_DIAL_DT MMDDYY10.;
	contents data = ersd_dialysis_dates;
	run;
quit;


PROC EXPORT DATA= WORK.ersd_dialysis_dates 
            OUTFILE= "DATA_DIR\ki_esrd_dial_dates.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



