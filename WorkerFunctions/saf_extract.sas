/**************************************************************************************************/
/* PROGRAM: 	Extracting candidate data for survival from listing models			              */
/* AUTHOR(S): 	Andrew Wey                                                                        */
/* CREATED: 	03/15/2018                                                                        */
/* MODIFY:                                                                                        */
/* FUNCTION: 	This is a key file called by the 'sas_saf_graf' function						  */
/*																								  */
/**************************************************************************************************/



libname allorg "SAF_DIRECTORY\CURRENT_SAF\allorg" ACCESS=READONLY; 
libname kipa "SAF_DIRECTORY\CURRENT_SAF\kipa" ACCESS=READONLY; 
libname liin "SAF_DIRECTORY\CURRENT_SAF\liin" ACCESS=READONLY; 
libname thor "SAF_DIRECTORY\CURRENT_SAF\thor" ACCESS=READONLY; 
libname srtrexp "SAF_DIRECTORY\CURRENT_SAF\srtrexp" ACCESS=READONLY; 
libname srtrfmt 'SAF_DIRECTORY\CURRENT_SAF\srtrfmt' ACCESS=READONLY;
options fmtsearch = (srtrfmt.formats);





data work.candidates;
	set ORGAN_SAF_AB.cand;
	keep RECIP_VARS;
	run;



data candidates;
	set candidates;
	can_zip = substr(CAN_PERM_ZIP, 1, 5);
	run;



PROC EXPORT DATA= WORK.candidates 
            OUTFILE= "DATA_DIR\ORGAN_AB_candidate_file.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;






data work.zipcodes;
	set SASHELP.ZIPCODE;
	keep zip x y;
	rename y = can_latitude;
	rename x = can_longitude;
	*zip = input(put(zip, $5.), 5.);
	run;


PROC EXPORT DATA= WORK.zipcodes
            OUTFILE= "DATA_DIR\zipcodes.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




data work.institute;
	set allorg.institution;
	zip = input(substr(PRIMARY_ZIP, 1, 5), 5.);
	run;


proc datasets;
	contents data = institute;
		modify institute;
			format zip z5.;
	contents data = institute;
	run;
quit;


data institute2;
	set institute;
	where CTR_ID ne .;
	run;



PROC EXPORT DATA= WORK.institute2
            OUTFILE= "DATA_DIR\institutes.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;






data work.kipa_cands;
	set kipa.cand;
	keep px_id pers_id CAN_LISTING_CTR_CD CAN_LISTING_CTR_TY CAN_LISTING_DT CAN_REM_DT WL_ORG;
	run;



PROC EXPORT DATA= WORK.kipa_cands 
            OUTFILE= "DATA_DIR\kipa_multi_candidates.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;





data work.liin_cands;
	set liin.cand;
	keep px_id pers_id CAN_LISTING_CTR_CD CAN_LISTING_CTR_TY CAN_LISTING_DT CAN_REM_DT WL_ORG;
	run;



PROC EXPORT DATA= WORK.liin_cands 
            OUTFILE= "DATA_DIR\liin_multi_candidates.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;





data work.thor_cands;
	set thor.cand;
	keep px_id pers_id CAN_LISTING_CTR_CD CAN_LISTING_CTR_TY CAN_LISTING_DT CAN_REM_DT WL_ORG;
	run;



PROC EXPORT DATA= WORK.thor_cands 
            OUTFILE= "DATA_DIR\thor_multi_candidates.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
