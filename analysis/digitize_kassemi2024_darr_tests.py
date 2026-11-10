"""Digitize Darr et al. (IJHMT 2016, Part I) 1-g vertical-upflow LN2 Tests 31 & 40 wall-T histories as replotted
in Kassemi, Hylton, Kartuzova, Ahuja, Hosangadi (NTRS 20240013323 preprint), Figs 3a,3b (Test 31) and 9a,9b (Test 40).
Embedded Excel-chart rasters (599x600 px) extracted losslessly with PyMuPDF.
Calibration: plot frame x=112.5px->0 s, x=570.5px->25 s; y=22px->300 K, y=501px->80 K (frame lines and 1 s / 5 K minor
ticks detected; tick spacing residual <0.5 px). Curve pixels: black (gray<90) or gray marker colour (110<gray<200,
low chroma); horizontal reference-line rows (dashed lines spanning plot) removed; per column the run nearest to the
previously tracked value is taken (continuity tracking), value = run midpoint.
"""
import numpy as np, pandas as pd
from PIL import Image
D='lit/raw/kassemi2024_figs/'
X0,X1,Y0,Y1=112.5,570.5,22.0,501.0
def t_of(x): return (x-X0)/(X1-X0)*25.0
def T_of(y): return 300.0-(y-Y0)/(Y1-Y0)*220.0
def digitize(f,mode):
    a=np.asarray(Image.open(D+f).convert('RGB')).astype(int); g=a.mean(2); chroma=a.max(2)-a.min(2)
    m=(g<90) if mode=='black' else ((g>110)&(g<205)&(chroma<25))
    linerows=[y for y in range(a.shape[0]) if (g[y,120:565]<235).sum()>250]
    for y in linerows: m[max(0,y-1):y+2,:]=False
    m[:int(Y0)+2,:]=False; m[int(Y1)-1:,:]=False
    prev=None; out=[]
    for x in range(int(X0)+3,int(X1)-2):
        ys=np.where(m[:,x])[0]
        if len(ys)==0: continue
        runs=np.split(ys,np.where(np.diff(ys)>4)[0]+1)
        cand=[(r.min(),r.max()) for r in runs if len(r)>=2]
        if not cand: continue
        if prev is None: c=min(cand,key=lambda r:r[0])  # topmost at start (T~290)
        else: c=min(cand,key=lambda r:abs((r[0]+r[1])/2-prev))
        if prev is not None and abs((c[0]+c[1])/2-prev)>80: continue
        prev=(c[0]+c[1])/2; out.append((t_of(x),T_of(prev),T_of(c[1]),T_of(c[0])))
    return np.array(out),linerows
PANELS={'test31_TC1-3':('p40_0.png','black'),'test31_TC4-6':('p40_2.png','gray'),
        'test40_TC1-3':('p49_1.png','black'),'test40_TC4-6':('p49_3.png','gray')}
if __name__=='__main__':
    res={k:digitize(*v) for k,v in PANELS.items()}
    for k,(o,lr) in res.items(): print(k,len(o),'t %.2f-%.2f'%(o[0,0],o[-1,0]),'T0 %.1f Tend %.1f'%(o[0,1],o[-1,1]),'refrows',lr)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,4,figsize=(20,5))
    for i,(k,(f,mode)) in enumerate(PANELS.items()):
        a=np.asarray(Image.open(D+f).convert('RGB')); ax[i].imshow(a); o=res[k][0]
        ax[i].plot(X0+o[:,0]/25*(X1-X0),Y0+(300-o[:,1])/220*(Y1-Y0),'r-',lw=0.7); ax[i].set_title(k)
    plt.savefig('data/darr2016ijhmt_1g_vertical/qc_overlay.png',dpi=90)
    import pickle; pickle.dump(res,open('/tmp/kas.pkl','wb'))

def crossing(o,T):
    i=np.where(o[:,1]<=T)[0][0]; return np.interp(T,[o[i,1],o[i-1,1]],[o[i,0],o[i-1,0]])
TABLE={'test31_TC1-3':{135:15.2,126.2:16.3,97:17.0},'test31_TC4-6':{135:16.1,126.2:17.0,97:18.0},
       'test40_TC1-3':{138:13.7,126.2:14.7,97:15.4},'test40_TC4-6':{138:14.6,126.2:15.8,97:16.7}}
BC={'31':dict(T0=289.9,Tin=86.96,Pout=314.624,Tsat_in=88.44,mdot=0.0334,Re=33940),
    '40':dict(T0=288.54,Tin=88.61,Pout=325.900,Tsat_in=88.44,mdot=0.0346,Re=35626)}
def write(res):
    A=np.pi/4*0.01168**2; errs=[]
    for test in ['31','40']:
        o1=res[f'test{test}_TC1-3'][0]; o2=res[f'test{test}_TC4-6'][0]; b=BC[test]
        lines=[]
        for k in [f'test{test}_TC1-3',f'test{test}_TC4-6']:
            for T,tt in TABLE[k].items():
                td=crossing(res[k][0],T); errs.append(td-tt); lines.append(f'{k} T={T} K: table t={tt} s, digitized t={td:.2f} s (diff {td-tt:+.2f})')
        tg=np.round(np.arange(0.15,min(o1[-1,0],o2[-1,0]),0.05),2)
        df=pd.DataFrame({'t_s':tg,'Twall_TC1-3_avg_K':np.interp(tg,o1[:,0],o1[:,1]),'Twall_TC4-6_avg_K':np.interp(tg,o2[:,0],o2[:,1])})
        hdr=f"""# SOURCE (data): Darr, S.R., Hu, H., Glikin, N.G., Hartwig, J.W., Majumdar, A.K., Leclair, A.C., Chung, J.N., 'An experimental study on terrestrial cryogenic transfer line chilldown I. Effect of mass flux, equilibrium quality, and inlet subcooling', Int. J. Heat Mass Transfer 103 (2016) 1225-1242, doi:10.1016/j.ijheatmasstransfer.2016.05.019 (Test {test}).
# SOURCE (figure digitized + BCs): Kassemi, M., Hylton, S., Kartuzova, O., Ahuja, V., Hosangadi, A., 'Predictions of Line Chilldown Boiling Regime Transitions by a Coupled CFD-Subgrid Boiling Model Validated against 1G LN2 Experiments', NASA NTRS 20240013323 (preprint, 2024), Fig. {'3a/3b' if test=='31' else '9a/9b'}; PDF lit/raw/ntrs/20240013323.pdf
# DATA AVAILABILITY: no public numeric repository found for Darr IJHMT 2016 I/II (Elsevier; no supplementary data). Digitized from the replot by Kassemi et al. (Excel chart of the measured data).
# FLUID: LN2. GRAVITY: 1 g, VERTICAL UPFLOW (Kassemi Fig.2a: g opposite to flow; 'only cases parallel to gravity considered').
# GEOMETRY: 304SS tube L=0.572 m, OD 12.70 mm, ID 11.68 mm (same UF rig as Darr 2016 npj / data/uf2016). TC1-3 at 0.149 m, TC4-6 at 0.401 m from tube inlet (3 TCs at 90 deg on OUTER wall, averaged); PT1 0.180 m, PT2 0.371 m.
#   Test section in vacuum chamber (316SS, emissivity 0.453+/-0.045, 293+/-5 K, 1.33 Pa); ~2.54 cm of tube protrudes from each chamber side. Bypass valve at inlet.
# BOUNDARY CONDITIONS (Kassemi et al. text, time-averaged): T_wall0={b['T0']} K; T_inlet={b['Tin']} K; P_out={b['Pout']} kPa; Tsat(P_in)={b['Tsat_in']} K (same value quoted for both tests; check);
#   mdot={b['mdot']} kg/s (G={b['mdot']/A:.1f} kg/m2/s); inlet liquid volume fraction=1; Re_in={b['Re']}. Mass-flow and pressure TIME HISTORIES NOT AVAILABLE (only means).
# COLUMNS: t_s [s]; Twall_TC1-3_avg_K, Twall_TC4-6_avg_K [K] outer-wall, circumferential average.
# DIGITIZATION: script analysis/digitize_kassemi2024_darr_tests.py; frame/tick calibration (0.46 K/px, 0.055 s/px); colour mask + continuity tracking; resampled to 0.05 s by linear interpolation.
# DIGITIZATION CHECK vs Kassemi Table 2/3 event times (read from the same data):
""" + ''.join(f'#   {l}\n' for l in lines) + """# UNCERTAINTY: digitization ~ +/-1 K in film boiling (marker thickness ~3 px), knee timing ~ +/-0.15 s; first ~0.3 s biased low by up to 3 K (thick markers at plot edge; true T0 given above). TC instrument ~ +/-1-2 K (Darr 2016).
"""
        with open(f'data/darr2016ijhmt_1g_vertical/test{test}_LN2_1g_upflow.csv','w') as fh:
            fh.write(hdr); df.round(3).to_csv(fh,index=False)
        print('\n'.join(lines))
    e=np.array(errs); print('event-time diff mean %.2f rms %.2f max %.2f'%(e.mean(),np.sqrt((e**2).mean()),abs(e).max()))
