"""
Build a tabulated CSV of cometary gas production rates compiled from
Harrington Pinto, Womack, Fernandez & Bauer (2022), PSJ 3, 247
"A Survey of CO, CO2, and H2O in Comets and Centaurs"  doi:10.3847/PSJ/ac960d

Source: Table 2 (contemporaneous CO/CO2/H2O) and Table 3 (additional CO/CO2/H2O).
Production rates published in units of 1e26 molecules/s; written here in
absolute molecules/s.

Instrument / wavelength regime / methodology for Table 2 rows are taken from
HP22 Table 1 and Sections 2.1-2.3 and 4.x.  HP22 does not tabulate instruments
for Table 3, so those fields are left blank except where HP22 itself identifies
the reference (Ootsubo et al. 2012 = AKARI, HP22 Table 1).

methodology:
  direct   -> Q derived from emission of the molecule itself
              (CO2 nu3 4.26 um, CO v(1-0) 4.67 um, CO J=2-1, CO 4th Positive,
              in situ mass spectrometry)
  indirect -> Q derived from a proxy
              (CO Cameron bands, [O I] forbidden lines, CO+CO2 broadband
              photometry minus an independent CO measurement, OH as H2O proxy)
"""
import csv, os

S = 1e26  # published unit -> molecules/s

COLS = ["designation", "date", "r_hel", "ta_deg",
        "q_h2o", "q_h2o_limit", "q_co2", "q_co2_limit", "q_co", "q_co_limit",
        "q_co_q_co2", "q_co_q_co2_err", "q_co_q_h2o", "q_co_q_h2o_err",
        "q_co2_q_h2o", "q_co2_q_h2o_err", "one_over_a0",
        "reference", "instrument", "regime", "methodology",
        "source_table", "epoch_id", "notes"]

# instrument shorthands: (instrument, regime, methodology)
AKARI  = ("AKARI IRC/NC", "near-IR", "direct")
IKS    = ("Vega 1 IKS", "near-IR", "direct")
ISO    = ("ISO ISOPHOT-S", "near-IR", "direct")
DI     = ("Deep Impact/EPOXI HRI-IR", "near-IR", "direct")
SPITZ  = ("Spitzer IRAC", "near-IR", "indirect")
SMT    = ("ARO 10-m Submillimeter Telescope", "radio/mm", "direct")
SMITH  = ("Harlan J. Smith 2.7-m Telescope", "optical", "indirect")
SUBARU = ("Subaru HDS", "optical", "indirect")
ISHELL = ("IRTF iSHELL", "near-IR", "direct")
IUE_CO = ("IUE SWP", "UV", "direct")
IUE_C2 = ("IUE SWP", "UV", "indirect")
FOS_CO = ("HST FOS", "UV", "direct")
FOS_C2 = ("HST FOS", "UV", "indirect")
COS    = ("HST COS", "UV", "direct")
ACS    = ("HST ACS/SBC", "UV", "direct")
ROSINA = ("Rosetta ROSINA DFMS", "in situ", "direct")
LMI    = ("DCT LMI", "optical", "indirect")
NONE   = ("", "", "")

def row(desig, date, rh, ta, h2o, h2ol, co2, co2l, co, col,
        rcc, rcce, rch, rche, r2h, r2he, a0, ref, inst, tbl, eid, note):
    i, reg, meth = inst
    def q(v):
        return "" if v is None else "%.3e" % (v * S)
    def f(v):
        return "" if v is None else ("%g" % v)
    return dict(zip(COLS, [desig, date, f(rh), f(ta),
                           q(h2o), h2ol, q(co2), co2l, q(co), col,
                           f(rcc), f(rcce), f(rch), f(rche), f(r2h), f(r2he),
                           f(a0), ref, i, reg, meth, tbl, eid, note]))

R = []
# ----------------------------------------------------------------- Table 2
T2 = [
 # desig, date, rh, ta, h2o,lim, co2,lim, co,lim, CO/CO2,err, CO/H2O,err, CO2/H2O,err, 1/a0, ref, inst, epoch, note
 ("1P/Halley","1986/03/06",0.79,61, 4500.0,"",180.00,"",500.00,"", 2.78,0.50, 0.110,0.016, 0.04,0.005, None,"Combes et al. 1988",IKS,"T2-01",""),
 ("9P/Tempel 1","2005/07/03",1.51,359, None,"",None,"",5.00,"", None,None,None,None,None,None, None,"Feldman et al. 2006",ACS,"T2-02","CO from CO 4th Positive bands, 1 day pre-impact"),
 ("9P/Tempel 1","2005/07/04",1.51,359, 46.0,"",3.20,"",None,"", 1.88,None, 0.130,None, 0.07,None, None,"Feaga et al. 2007",DI,"T2-03","Pre-impact (July 4) water; ratios use the July 3 HST CO"),
 ("21P/Giacobini-Zinner","2018/10/03",1.07,28, None,"",None,"",None,"", None,None,None,None, 0.11,0.01, None,"Shinnaka et al. 2020",SUBARU,"T2-04","Only the CO2/H2O ratio was published (from [O I] G/R); no production rates"),
 ("21P/Giacobini-Zinner","2018/10/10",1.10,36, 203.7,"",None,"",2.55,"", 0.11,None, 0.013,0.002, None,None, None,"Roth et al. 2020",ISHELL,"T2-05","CO/CO2 formed by combining with the 2018/10/03 Subaru CO2/H2O ratio"),
 ("22P/Kopff","2009/04/22",1.61,340, 59.4,"",10.96,"",1.65,"<", 0.15,None, 0.030,None, 0.18,0.03, None,"Ootsubo et al. 2012",AKARI,"T2-06","CO upper limit; ratios involving CO are upper limits"),
 ("29P/Schwassmann-Wachmann 1","2009/11/18",6.18,135, 67.8,"",3.5,"<",291.5,"", 83.2,None, 4.630,0.708, 0.05,None, None,"Ootsubo et al. 2012",AKARI,"T2-07","CO2 upper limit; CO/CO2 is a lower limit (>83.2), CO2/H2O an upper limit (<0.05); note QCO/QH2O as published (4.630) differs from QCO/QH2O computed from the tabulated rates (4.30)"),
 ("46P/Wirtanen","2019/01/08",1.11,29, None,"",8.55,"",None,"", 0.04,None, None,None, 0.15,None, None,"McKay et al. 2021",SMITH,"T2-08","CO2 inferred from [O I] forbidden-oxygen lines; CO/CO2 is an upper limit (<0.04)"),
 ("46P/Wirtanen","2019/01/11",1.13,33, 57.0,"",None,"",0.31,"<", 0.04,None, 0.005,None, None,None, None,"McKay et al. 2021",ISHELL,"T2-09","CO upper limit; ratios involving CO are upper limits"),
 ("67P/Churyumov-Gerasimenko","2014/08/23",3.50,229, 0.5,"",0.07,"",0.21,"", 3.06,None, 0.450,None, 0.15,None, None,"Combi et al. 2020",ROSINA,"T2-10","Representative value; 67P was CO-dominant at ~3.5 au pre-perihelion"),
 ("67P/Churyumov-Gerasimenko","2015/08/13",1.24,359, 158.1,"",6.26,"",4.78,"", 0.76,None, 0.030,None, 0.04,None, None,"Combi et al. 2020",ROSINA,"T2-11","Representative value near perihelion"),
 ("67P/Churyumov-Gerasimenko","2016/01/02",2.03,89, 11.5,"",2.71,"",0.36,"", 0.13,None, 0.030,None, 0.24,None, None,"Combi et al. 2020",ROSINA,"T2-12","Representative value post-perihelion"),
 ("81P/Wild 2","2009/12/14",1.74,320, 60.0,"",8.63,"",2.24,"<", 0.26,None, 0.040,None, 0.14,0.02, None,"Ootsubo et al. 2012",AKARI,"T2-13","CO upper limit"),
 ("88P/Howell","2009/07/03",1.73,294, 33.5,"",8.65,"",2.72,"<", 0.31,None, 0.080,None, 0.26,0.04, None,"Ootsubo et al. 2012",AKARI,"T2-14","CO upper limit"),
 ("103P/Hartley 2","2010/11/02",1.06,8, 100.0,"",20.00,"",None,"", 0.01,None, 0.003,None, 0.20,None, None,"A'Hearn et al. 2011",DI,"T2-15","Ratios use the 2010/11/04 HST CO"),
 ("103P/Hartley 2","2010/11/04",1.06,8, None,"",None,"",0.26,"", None,None,None,None,None,None, None,"Weaver et al. 2011",COS,"T2-16","CO from CO 4th Positive bands"),
 ("144P/Kushida","2009/04/18",1.70,53, 37.6,"",5.57,"",1.14,"<", 0.20,None, 0.030,None, 0.15,0.02, None,"Ootsubo et al. 2012",AKARI,"T2-17","CO upper limit"),
 # --- IUE / HST FOS epochs: CO direct (4th Positive), CO2 indirect (Cameron bands) -> split by species
 ("C/1979 Y1 (Bradfield)","1980/01/10",0.71,56, 1460.0,"",None,"",51.00,"", 1.00,None, 0.030,None, None,None, 232.7,"Feldman et al. 1997; Tozzi et al. 1998",IUE_CO,"T2-18","CO from CO 4th Positive bands; CO2 for the same epoch is the companion row"),
 ("C/1979 Y1 (Bradfield)","1980/01/10",0.71,56, None,"",51.00,"",None,"", 1.00,None, None,None, 0.03,None, 232.7,"Feldman et al. 1997; Tozzi et al. 1998",IUE_C2,"T2-18","CO2 inferred from CO Cameron bands; H2O in the companion row"),
 ("C/1989 X1 (Austin)","1990/05/09",0.83,99, 1020.0,"",None,"",17.00,"", 0.81,None, 0.020,None, None,None, 0.3,"Feldman et al. 1997; Tozzi et al. 1998",IUE_CO,"T2-19","CO from CO 4th Positive bands"),
 ("C/1989 X1 (Austin)","1990/05/09",0.83,99, None,"",21.00,"",None,"", 0.81,None, None,None, 0.02,None, 0.3,"Feldman et al. 1997; Tozzi et al. 1998",IUE_C2,"T2-19","CO2 inferred from CO Cameron bands"),
 ("C/1990 K1 (Levy)","1990/08/26",1.38,291, 2480.0,"",None,"",101.00,"", 0.59,None, 0.040,None, None,None, 0.7,"Feldman et al. 1997; Tozzi et al. 1998",IUE_CO,"T2-20","CO from CO 4th Positive bands"),
 ("C/1990 K1 (Levy)","1990/08/26",1.38,291, None,"",170.00,"",None,"", 0.59,None, None,None, 0.07,None, 0.7,"Feldman et al. 1997; Tozzi et al. 1998",IUE_C2,"T2-20","CO2 inferred from CO Cameron bands"),
 ("C/1990 K1 (Levy)","1990/09/18",1.13,311, 2010.0,"",None,"",168.00,"", 0.63,None, 0.080,None, None,None, 0.7,"Feldman et al. 1997; Tozzi et al. 1998",IUE_CO,"T2-21","CO from CO 4th Positive bands"),
 ("C/1990 K1 (Levy)","1990/09/18",1.13,311, None,"",268.00,"",None,"", 0.63,None, None,None, 0.13,None, 0.7,"Feldman et al. 1997; Tozzi et al. 1998",IUE_C2,"T2-21","CO2 inferred from CO Cameron bands"),
 # --- Hale-Bopp (ISO)
 ("C/1995 O1 (Hale-Bopp)","1996/04/27",4.58,233, 131.0,"",130.00,"",700.00,"", 5.38,None, 5.340,None, 0.99,None, 38.4,"Crovisier et al. 1999a; Colom et al. 1997",ISO,"T2-22","H2O from Colom et al. 1997 OH radio observations used as a proxy"),
 ("C/1995 O1 (Hale-Bopp)","1996/09/26",2.93,248, 3300.0,"",740.00,"",2300.00,"", 3.11,None, 0.700,None, 0.22,None, 38.4,"Crovisier et al. 1999a",ISO,"T2-23",""),
 ("C/1995 O1 (Hale-Bopp)","1997/12/29",3.89,122, 280.0,"",350.00,"",1500.00,"", 4.29,None, 5.360,None, 1.25,None, 38.4,"Crovisier et al. 1999a",ISO,"T2-24",""),
 ("C/1995 O1 (Hale-Bopp)","1998/04/06",4.90,129, None,"",300.00,"",3000.00,"", 10.00,None, None,None, None,None, 38.4,"Crovisier et al. 1999a",ISO,"T2-25",""),
 # --- Hyakutake (HST FOS): CO direct, CO2 indirect
 ("C/1996 B2 (Hyakutake)","1996/06/02",0.94,170, 1530.0,"",None,"",180.00,"", 2.83,0.96, 0.111,0.007, None,None, 15.5,"McPhate 1999",FOS_CO,"T2-26","CO from CO 4th Positive bands; H2O derived from the OH band"),
 ("C/1996 B2 (Hyakutake)","1996/06/02",0.94,170, None,"",60.00,"",None,"", 2.83,0.96, None,None, 0.04,0.01, 15.5,"McPhate 1999",FOS_C2,"T2-26","CO2 inferred from CO Cameron bands"),
 # --- AKARI OCCs
 ("C/2006 OF2 (Broughton)","2008/09/16",2.43,0, 61.2,"",13.97,"",2.52,"<", 0.18,None, 0.040,None, 0.23,0.03, 0.2,"Ootsubo et al. 2012",AKARI,"T2-27","CO upper limit"),
 ("C/2006 OF2 (Broughton)","2009/03/28",3.20,59, 17.0,"",9.90,"",4.46,"<", 0.45,None, 0.260,None, 0.58,0.09, 0.2,"Ootsubo et al. 2012",AKARI,"T2-28","CO upper limit"),
 ("C/2006 Q1 (McNaught)","2008/06/03",2.78,351, 36.3,"",16.19,"",3.67,"<", 0.23,None, 0.100,None, 0.45,0.06, 0.4,"Ootsubo et al. 2012",AKARI,"T2-29","CO upper limit"),
 ("C/2006 W3 (Christensen)","2008/12/21",3.66,315, 83.0,"",84.59,"",299.30,"", 3.54,0.50, 3.610,0.519, 1.02,0.15, 3.7,"Ootsubo et al. 2012",AKARI,"T2-30",""),
 ("C/2006 W3 (Christensen)","2009/06/16",3.13,355, 201.0,"",84.66,"",197.70,"", 2.34,0.33, 0.980,0.140, 0.42,0.06, 3.7,"Ootsubo et al. 2012",AKARI,"T2-31",""),
 ("C/2007 N3 (Lulin)","2009/02/05",1.28,26, 409.1,"",48.48,"",8.67,"<", 0.18,None, 0.020,None, 0.12,0.02, 0.3,"Ootsubo et al. 2012",AKARI,"T2-32","CO upper limit"),
 ("C/2007 N3 (Lulin)","2009/03/30",1.70,64, 223.8,"",22.42,"",2.86,"<", 0.13,None, 0.010,None, 0.10,0.01, 0.3,"Ootsubo et al. 2012",AKARI,"T2-33","CO upper limit"),
 ("C/2007 Q3 (Siding Spring)","2009/03/03",3.29,292, 39.8,"",6.95,"",3.97,"<", 0.57,None, 0.100,None, 0.17,0.03, 0.3,"Ootsubo et al. 2012",AKARI,"T2-34","CO upper limit"),
 ("C/2008 Q3 (Garradd)","2009/07/05",1.81,7, 112.0,"",31.04,"",29.30,"", 0.94,0.14, 0.260,0.038, 0.28,0.04, 2.5,"Ootsubo et al. 2012",AKARI,"T2-35",""),
 ("C/2009 P1 (Garradd)","2012/03/26",2.00,56, 460.0,"",39.00,"",290.00,"", 7.44,2.45, 0.630,0.206, 0.08,0.02, 4.2,"Feaga et al. 2014",DI,"T2-36",""),
 ("C/2012 S1 (ISON)","2013/06/13",3.34,187, 100.0,"<",1.21,"",7.00,"", 5.79,None, None,None, None,None, 1.2,"Lisse et al. 2013; Meech et al. 2013",SPITZ,"T2-37","Date is for Spitzer; CO is the modelled value of Meech et al. 2013, CO2 inferred from the Spitzer CO+CO2 4.5 um excess. Water (upper limit) from May 9"),
 # --- C/2016 R2
 ("C/2016 R2 (PanSTARRS)","2018/02/12",2.76,332, None,"",90.50,"",None,"", None,None,None,None,None,None, 11.7,"McKay et al. 2019",SPITZ,"T2-38","CO2 inferred from Spitzer IRAC CO+CO2 photometry minus the independent ARO SMT CO"),
 ("C/2016 R2 (PanSTARRS)","2018/02/13",2.76,333, None,"",None,"",550.00,"", 6.08,None, None,None, None,None, 11.7,"McKay et al. 2019",SMT,"T2-39","CO J=2-1 at 230 GHz; ~11 hr from the Spitzer epoch"),
 ("C/2016 R2 (PanSTARRS)","2018/02/21",2.73,335, 3.1,"",None,"",None,"", None,None, 177.420,30.907, 29.19,None, 11.7,"McKay et al. 2019",LMI,"T2-40","H2O from DCT LMI OH narrowband imaging; closest water detection to the CO/CO2 epochs"),
 ("C/2016 R2 (PanSTARRS)","2019/06/10",4.72,84, None,"",6.21,"",None,"", 19.00,None, None,None, None,None, 11.7,"McKay, A. et al. (in prep., priv. comm.)",SPITZ,"T2-41","CO2 inferred from the Spitzer CO+CO2 proxy"),
 ("C/2016 R2 (PanSTARRS)","2019/06/12",4.73,84, None,"",None,"",118.00,"", None,None,None,None,None,None, 11.7,"Wierzchos 2019",SMT,"T2-42","CO J=2-1 at 230 GHz"),
]
for t in T2:
    (d,dt,rh,ta,h2o,h2ol,co2,co2l,co,col,rcc,rcce,rch,rche,r2h,r2he,a0,ref,inst,eid,note)=t
    R.append(row(d,dt,rh,ta,h2o,h2ol,co2,co2l,co,col,rcc,rcce,rch,rche,r2h,r2he,a0,ref,inst,"Table 2",eid,note))

# ----------------------------------------------------------------- Table 3
UNSPEC = "Instrument, regime and methodology are not specified in Harrington Pinto et al. 2022; see the original reference"
T3 = [
 # desig, date, rh, h2o,lim, co2,lim, co,lim, CO/H2O,err, CO2/H2O,err, ref, inst, epoch, note
 ("19P/Borrelly","2008/12/30",2.19, 6.4,"",1.6,"",None,"", None,None, 0.250,0.037,"Ootsubo et al. 2012",AKARI,"T3-01",""),
 ("64P/Swift-Gehrels","2009/11/23",2.27, 6.9,"",1.5,"",None,"", None,None, 0.220,0.033,"Ootsubo et al. 2012",AKARI,"T3-02",""),
 ("103P/Hartley 2","1997/12/31",1.04, None,"",12.0,"",None,"", None,None, None,None,"Crovisier et al. 1999b",NONE,"T3-03",UNSPEC),
 ("103P/Hartley 2","1998/01/01",1.04, 124.0,"",None,"",None,"", None,None, 0.097,None,"Crovisier et al. 1999b",NONE,"T3-04","CO2/H2O ratio combines this water with the 1997/12/31 CO2. "+UNSPEC),
 ("116P/Wild 4","2009/05/16",2.22, 15.8,"",0.9,"",None,"", None,None, 0.056,0.008,"Ootsubo et al. 2012",AKARI,"T3-05",""),
 ("118P/Shoemaker-Levy 4","2009/09/08",2.18, 10.0,"",3.0,"",None,"", None,None, 0.300,0.043,"Ootsubo et al. 2012",AKARI,"T3-06",""),
 ("157P/Tritton","2009/12/30",1.48, 4.5,"",0.3,"",None,"", None,None, 0.069,0.010,"Ootsubo et al. 2012",AKARI,"T3-07",""),
 ("C/2007 G1 (LINEAR)","2008/08/20",2.80, 18.4,"",4.2,"",None,"", None,None, 0.230,0.033,"Ootsubo et al. 2012",AKARI,"T3-08",""),
 ("8P/Tuttle","2008/01/27",1.03, 544.0,"",None,"",2.4,"", 0.004,0.001, None,None,"Boehnhardt et al. 2008",NONE,"T3-09",UNSPEC),
 ("8P/Tuttle","2007/12/23",1.15, 212.8,"",None,"",0.8,"<", 0.004,None, None,None,"Bonev et al. 2008",NONE,"T3-10","CO upper limit. "+UNSPEC),
 ("21P/Giacobini-Zinner","1998/10/02",1.25, 320.0,"",None,"",33.0,"", 0.100,0.058, None,None,"Mumma et al. 2000",NONE,"T3-11",UNSPEC),
 ("21P/Giacobini-Zinner","2018/07/30",1.17, 231.6,"",None,"",5.5,"", 0.024,0.006, None,None,"Faggi et al. 2019",NONE,"T3-12",UNSPEC),
 ("21P/Giacobini-Zinner","2018/10/07",1.09, 258.3,"",None,"",3.7,"", 0.014,0.008, None,None,"Faggi et al. 2019",NONE,"T3-13",UNSPEC),
 ("41P/Tuttle-Giacobini-Kresak","2017/04/26",1.06, 50.0,"",None,"",8.5,"<", 0.170,None, None,None,"Faggi et al. 2019",NONE,"T3-14","CO upper limit. "+UNSPEC),
 ("45P/Honda-Mrkos-Pajdusakova","2017/01/08",0.56, 372.0,"",None,"",1.9,"", 0.005,0.001, None,None,"DiSanti et al. 2017",NONE,"T3-15",UNSPEC),
 ("73P/Schwassmann-Wachmann 3","2006/05/27",0.95, 131.0,"",None,"",0.6,"", 0.005,0.002, None,None,"DiSanti et al. 2007",NONE,"T3-16",UNSPEC),
 ("153P/Ikeya-Zhang","2002/04/20",0.89, 2150.0,"",None,"",154.0,"", 0.072,0.004, None,None,"Lupu et al. 2007",NONE,"T3-17",UNSPEC),
 ("C/1995 O1 (Hale-Bopp)","1997/01/27",1.49, None,"",None,"",None,"", 0.267,0.029, None,None,"DiSanti et al. 2001",NONE,"T3-18","Only the CO/H2O ratio was published; no production rates. "+UNSPEC),
 ("C/1995 O1 (Hale-Bopp)","1997/03/01",1.06, None,"",None,"",None,"", 0.271,0.011, None,None,"DiSanti et al. 2001",NONE,"T3-19","Only the CO/H2O ratio was published; no production rates. "+UNSPEC),
 ("C/1995 O1 (Hale-Bopp)","1997/04/09",0.93, None,"",None,"",None,"", 0.276,0.022, None,None,"DiSanti et al. 2001",NONE,"T3-20","Only the CO/H2O ratio was published; no production rates. "+UNSPEC),
 ("C/1995 O1 (Hale-Bopp)","1997/05/01",1.06, None,"",None,"",None,"", 0.280,0.024, None,None,"DiSanti et al. 2001",NONE,"T3-21","Only the CO/H2O ratio was published; no production rates. "+UNSPEC),
 ("C/1996 B2 (Hyakutake)","1996/03/24",1.06, 2540.0,"",None,"",None,"", None,None, None,None,"Dello Russo et al. 2002",NONE,"T3-22",UNSPEC),
 ("C/1996 B2 (Hyakutake)","1996/03/24",1.06, None,"",None,"",440.0,"", 0.170,0.020, None,None,"DiSanti et al. 2003",NONE,"T3-23","CO/H2O uses the companion water of Dello Russo et al. 2002. "+UNSPEC),
 ("C/1996 B2 (Hyakutake)","1996/04/12",0.64, 3990.0,"",None,"",None,"", None,None, None,None,"Dello Russo et al. 2002",NONE,"T3-24",UNSPEC),
 ("C/1996 B2 (Hyakutake)","1996/04/12",0.64, None,"",None,"",803.0,"", 0.200,0.039, None,None,"DiSanti et al. 2003",NONE,"T3-25","CO/H2O uses the companion water of Dello Russo et al. 2002. "+UNSPEC),
 ("C/1999 H1 (Lee)","1999/08/20",1.06, 1260.0,"",None,"",23.0,"", 0.018,0.002, None,None,"Mumma et al. 2001c",NONE,"T3-26",UNSPEC),
 ("C/1999 S4 (LINEAR)","2000/07/13",0.81, 446.0,"",None,"",2.0,"", 0.004,0.003, None,None,"Mumma et al. 2001b",NONE,"T3-27",UNSPEC),
 ("C/1999 T1 (McNaught-Hartley)","2001/01/13",1.3, 820.0,"",None,"",140.0,"", 0.170,0.050, None,None,"Mumma et al. 2001a",NONE,"T3-28",UNSPEC),
 ("C/2000 WM1 (LINEAR)","2001/11/25",1.32, 177.1,"",None,"",0.9,"", 0.005,0.001, None,None,"Radeva et al. 2010",NONE,"T3-29",UNSPEC),
 ("C/2001 A2 (LINEAR)","2001/07/10",1.17, None,"",None,"",16.6,"", None,None, None,None,"Magee-Sauer et al. 2008",NONE,"T3-30",UNSPEC),
 ("C/2001 A2 (LINEAR)","2001/07/10",1.17, 430.0,"",None,"",None,"", 0.039,0.011, None,None,"Dello Russo et al. 2005",NONE,"T3-31","CO/H2O uses the companion CO of Magee-Sauer et al. 2008. "+UNSPEC),
 ("C/2001 Q4 (NEAT)","2004/04/26",1.02, 2000.0,"",None,"",176.0,"", 0.088,0.008, None,None,"Lupu et al. 2007",NONE,"T3-32",UNSPEC),
 ("C/2004 Q2 (Machholz)","2004/11/29",1.48, 1253.0,"",None,"",63.5,"", 0.051,0.005, None,None,"Bonev et al. 2009",NONE,"T3-33",UNSPEC),
 ("C/2006 M4 (SWAN)","2006/11/07",1.08, 1756.0,"",None,"",8.7,"", 0.005,0.002, None,None,"DiSanti et al. 2009",NONE,"T3-34",UNSPEC),
 ("C/2006 P1 (McNaught)","2007/01/27",0.55, 17400.0,"",None,"",341.0,"", 0.020,0.004, None,None,"Dello Russo et al. 2009",NONE,"T3-35",UNSPEC),
 ("C/2007 N3 (Lulin)","2009/02/01",1.26, 2013.0,"",None,"",43.6,"", 0.022,0.001, None,None,"Gibb et al. 2012",NONE,"T3-36",UNSPEC),
 ("C/2007 W1 (Boattini)","2008/07/10",0.90, 122.3,"",None,"",5.5,"", 0.045,0.005, None,None,"Villanueva et al. 2011",NONE,"T3-37",UNSPEC),
 ("C/2009 P1 (Garradd)","2011/09/21",2.01, 1410.0,"",None,"",67.0,"", 0.048,0.011, None,None,"McKay et al. 2015",NONE,"T3-38",UNSPEC),
 ("C/2009 P1 (Garradd)","2011/10/10",1.85, 1640.0,"",None,"",103.0,"", 0.063,0.011, None,None,"McKay et al. 2015",NONE,"T3-39",UNSPEC),
 ("C/2009 P1 (Garradd)","2012/01/25",1.62, 1050.0,"",None,"",164.0,"", 0.160,0.032, None,None,"McKay et al. 2015",NONE,"T3-40",UNSPEC),
 ("C/2009 P1 (Garradd)","2012/02/27",1.69, 1000.0,"",None,"",196.0,"", 0.196,0.050, None,None,"McKay et al. 2015",NONE,"T3-41",UNSPEC),
 ("C/2010 G2 (Hill)","2012/01/10",2.51, 61.0,"",None,"",47.0,"", 0.770,0.121, None,None,"Kawakita et al. 2014",NONE,"T3-42",UNSPEC),
 ("C/2012 F6 (Lemmon)","2013/03/31",0.75, 4590.0,"",None,"",195.0,"", 0.042,0.006, None,None,"Paganini et al. 2014b",NONE,"T3-43",UNSPEC),
 ("C/2012 S1 (ISON)","2013/11/17",0.53, 1827.0,"",None,"",22.4,"", 0.012,0.002, None,None,"DiSanti et al. 2016",NONE,"T3-44",UNSPEC),
 ("C/2013 R1 (Lovejoy)","2013/10/24",1.34, 153.2,"",None,"",15.2,"", 0.099,0.020, None,None,"Paganini et al. 2014a",NONE,"T3-45",UNSPEC),
 ("C/2013 R1 (Lovejoy)","2013/12/10",0.84, None,"",None,"",None,"", 0.123,0.017, None,None,"Dello Russo et al. 2016",NONE,"T3-46","Only the CO/H2O ratio was published; no production rates. "+UNSPEC),
 ("C/2020 F3 (NEOWISE)","2020/07/20",0.56, 4413.0,"",None,"",82.8,"", 0.019,0.001, None,None,"Faggi et al. 2021",NONE,"T3-47",UNSPEC),
 ("C/2020 F3 (NEOWISE)","2020/08/01",0.83, 1176.5,"",None,"",30.0,"", 0.026,0.001, None,None,"Faggi et al. 2021",NONE,"T3-48",UNSPEC),
 ("2I/Borisov","2019/12/19-22, 2020/01/13",2.12, None,"",None,"",None,"", 1.425,0.125, None,None,"Bodewits et al. 2020",NONE,"T3-49","Median over the listed dates; r_hel is the median. Only the CO/H2O ratio was published. "+UNSPEC),
]
for t in T3:
    (d,dt,rh,h2o,h2ol,co2,co2l,co,col,rch,rche,r2h,r2he,ref,inst,eid,note)=t
    R.append(row(d,dt,rh,None,h2o,h2ol,co2,co2l,co,col,None,None,rch,rche,r2h,r2he,None,ref,inst,"Table 3",eid,note))

out = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                   "data", "reference", "harrington_pinto2022_gas_production_rates.csv")   # notebooks/comspec/ -> project
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(R)
print("rows:", len(R))
