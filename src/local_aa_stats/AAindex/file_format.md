Data Format of AAindex1
************************************************************************
*                                                                      *
* Each entry has the following format.                                 *
*                                                                      *
* H Accession number                                                   *
* D Data description                                                   *
* R PMID                                                               *
* A Author(s)                                                          *
* T Title of the article                                               *
* J Journal reference                                                  *
* * Comment or missing                                                 *
* C Accession numbers of similar entries with the correlation          *
*   coefficients of 0.8 (-0.8) or more (less).                         *
*   Notice: The correlation coefficient is calculated with zeros       *
*   filled for missing values.                                         *
* I Amino acid index data in the following order                       *
*   Ala    Arg    Asn    Asp    Cys    Gln    Glu    Gly    His    Ile *
*   Leu    Lys    Met    Phe    Pro    Ser    Thr    Trp    Tyr    Val *
* //                                                                   *
************************************************************************

Data Format of AAindex2 and AAindex3
************************************************************************
*                                                                      *
* Each entry has the following format.                                 *
*                                                                      *
* H Accession number                                                   *
* D Data description                                                   *
* R PMID                                                               *
* A Author(s)                                                          *
* T Title of the article                                               *
* J Journal reference                                                  *
* * Comment or missing                                                 *
* M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV           *
*   AA                                                                 *
*   AR RR                                                              *
*   AN RN NN                                                           *
*   AD RD ND DD                                                        *
*   AC RC NC DC CC                                                     *
*   AQ RQ NQ DQ CQ QQ                                                  *
*   AE RE NE DE CE QE EE                                               *
*   AG RG NG DG CG QG EG GG                                            *
*   AH RH NH DH CH QH EH GH HH                                         *
*   AI RI NI DI CI QI EI GI HI II                                      *
*   AL RL NL DL CL QL EL GL HL IL LL                                   *
*   AK RK NK DK CK QK EK GK HK IK LK KK                                *
*   AM RM NM DM CM QM EM GM HM IM LM KM MM                             *
*   AF RF NF DF CF QF EF GF HF IF LF KF MF FF                          *
*   AP RP NP DP CP QP EP GP HP IP LP KP MP FP PP                       *
*   AS RS NS DS CS QS ES GS HS IS LS KS MS FS PS SS                    *
*   AT RT NT DT CT QT ET GT HT IT LT KT MT FT PT ST TT                 *
*   AW RW NW DW CW QW EW GW HW IW LW KW MW FW PW SW TW WW              *
*   AY RY NY DY CY QY EY GY HY IY LY KY MY FY PY SY TY WY YY           *
*   AV RV NV DV CV QV EV GV HV IV LV KV MV FV PV SV TV WV YV VV        *
* //                                                                   *
************************************************************************

