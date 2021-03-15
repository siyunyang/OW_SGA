/******************************************************************************
 Purpose:  This macro creates Connect-S plots used to display standardized
           differences in subgroups analysis.

 Original Author: Daniel Wojdyla
 Date: March 2021
 Version: 0.3

 In the following description a subgroup level represents a subset of the
 analysis population. For example: men or age > 55 yo.

 Data Source: the input dataset must be arranged with one row for each
 confounder by subgroup level combination which will be represented by one
 circle in the plot

 The following variables are required in the input dataset:
   a numeric variable indexing subgroup levels
   a numeric variable indexing confounders
   a numeric variable with the unweighted or weighted standardized differences

 The following variables are optional in the input dataset:
   a numeric variable with the sample size within each subgroup level
   a numeric variable with the variance inflation for each subgroup level

 The sample size and variance inflation variables are subgroup level specific
 and can be included only once (for one confounder) or multiple times (for
 every confounder)

 For example if there are 12 confounders and 15 subgroup levels, the dataset
 should contain 12 x 15 = 180 rows. Sample sizes and variance inflations can
 be included only for one confounder (any confounder works) or for all the
 confounders.

 Macro Parameters:

 Required:
 DS=         input dataset, must have structure described above

 STDVBLE=    variable containing standardized differences, unweighted or
             weighted. Either signed or absolute values are acceptable.

 CONFIDX=    numeric variable indexing the confounders (1,2,..). Confounders
             should be ordered as expected to appear in the plot

 SUBGIDX=    numeric variable indexing the subgroup levels (1,2,...) Subgroup
             levels should be ordered as expected to appear in the plot

 Optional:
 VIVBLE=     numeric variable containing variance inflation

 SAMPSIZE=   numeric variable containing sample size within each subgroup level

 GPATH=      location where the plot will be saved (default: same folder as the
             SAS program)

 IMGNAME=    name of the file containing the plots (default: connectS_plot.png)

 H=          height of the plot in inches (default: approximately 2 inches for
             every 5 subgroup levels)

 W=          width of the plot in inches (default: approximately 2 inches for
             every 5 confounders)

 TIT=        title for the plot (default: No title)

 REFS=       position for reference lines. For example, if the first subgroup
             variable has 3 levels, the first reference line should be drawn at
             3.5. If multiple lines are needed list the positions separated by
             blanks (default: no lines)

 SUBGFMT=    format for the subgroup levels. The format is used to label the
             rows (Y-axis) in the plot (default: no format, numbers are used)

 CONFFMT=    format for the confounders. The format is used to label the columns
             (X-axis) in the plot (default: no format, numbers are used)

 FS=         font size used in the plot. Certain labels are scaled from the
             FS= parameter (default: 8)

 Example call:

 %connectS(ds=selected,stdvble=w_std,confidx=vble_n,subgidx=subg_n,vivble=vi, sampsize=ss,
           gpath=/output/plots,imgname=example,h=5in,w=8in,tit=%str(Main effects Logistic Regression - IP Weights),
           refs=3.5 6.5 8.5 10.5 12.5,subgfmt=subgx.,conffmt=vblelab.,fs=8);
**************************************************************************************************************************************/
* Macro;
%macro connectS(ds=,stdvble=,confidx=,subgidx=,vivble=,sampsize=,gpath=,
                imgname=connectS_plot,h=5,w=8,tit=,refs=,subgfmt=,conffmt=,fs=8);
data ds1;
   set &ds;

   if abs(&stdvble) > 0.20 then stdcat=3;
   else if 0.15 < abs(&stdvble) <= 0.20 then stdcat=2;
   else if 0.1 < abs(&stdvble) <= 0.15 then stdcat=1;
   else if abs(&stdvble) <= 0.1 then stdcat=0;
run;

proc sort data=ds1;
   by &confidx &subgidx;
run;

data order;
   do stdcat = 0 to 3;
      output;
   end;
run;

data to_plot;
   set order
       ds1;
run;

proc sql noprint;
   select max(&subgidx) into :nsubgs from to_plot;
   select max(&confidx) into :nconf from to_plot;
quit;

proc format;
   value catlev 0="ASMD (*ESC*){unicode '2264'x} 0.10" 1="ASMD 0.10 - 0.15" 2="ASMD 0.15 - 0.20" 3="ASMD > 0.20";
run;

%let fss = %sysevalf(0.75*&fs);

%if %length(&w) = 0 %then %do;
   %let ww = %sysevalf(0.4*&nconf);
   %end;

%else %if %length(&w) > 0 %then %do;
   %let ww = &w;
   %end;

%if %length(&h) = 0 %then %do;
   %let hh = %sysevalf(0.4*&nsubgs);
   %end;

%else %if %length(&h) > 0 %then %do;
   %let hh = &h;
   %end;

ods listing image_dpi=600 gpath="&gpath";
ods graphics on / imagename="&imgname" noborder height=&hh.in width=&ww.in;
proc sgplot data=to_plot;

   %if %length(&tit) > 0 %then %do;
      title "&tit";
      %end;

   styleattrs datacolors=(white grayC0 gray70 black);
   scatter y=&subgidx x=&confidx / markerattrs=(size=12 symbol=circlefilled) filledoutlinedmarkers
                            markeroutlineattrs=(color=black thickness=0.1)
							       group=stdcat name="plot" grouporder=ascending;


   %if %length(&refs) > 0 %then %do;
      refline &refs / axis=y;
      %end;

   xaxis label="Confounder" valueattrs=(size=&fs.pt) fitpolicy=rotate type=discrete;
   yaxis reverse values=(1 to &nsubgs by 1) label=%str("Subgroups") valueattrs=(size=&fs.pt);

   %if %length(&sampsize) > 0 %then %do;
      yaxistable &sampsize. / y=&subgidx label="Sample Size" labelattrs=(size=&fss.pt) location=inside position=right valueattrs=(size=&fs.pt) valuejustify=center pad=(left=0px) stat=mean;
      %end;

   %if %length(&vivble) > 0 %then %do;
      yaxistable &vivble / y=&subgidx label="Variance Inflation" labelattrs=(size=&fss.pt) location=inside position=right valueattrs=(size=&fs.pt) valuejustify=center pad=(left=0px) stat=mean;
      %end;

   keylegend  "plot" / noborder valueattrs=(size=&fs.pt) title="Absolute Standardized Mean Difference" titleattrs=(size=&fs.pt) down=1;
   format stdcat catlev.;

   %if %length(&subgfmt) > 0 %then %do;
      format &subgidx &subgfmt;
      %end;

   %if %length(&conffmt) > 0 %then %do;
      format &confidx &conffmt;
      %end;

   %if %length(&vivble) > 0 %then %do;
      format &vivble 5.1;
      %end;
run;
%mend connectS;
