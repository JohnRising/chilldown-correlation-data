"""Raster digitizer for Hartwig et al. 2022 npj Microgravity 8:33, Figs 1a, 2a-c (test F2 blue curve).
Method: figure images extracted losslessly from the article PDF (PyMuPDF, embedded 1500-px raster);
axis calibration from detected tick marks on the plot frame (linear least-squares fit pixel->value);
curve = pure-blue pixel mask (R<90,G<90,B>170) per image column, median row; legend/error-bar boxes masked.
"""
import numpy as np, sys, json
from PIL import Image
D='lit/raw/hartwig2022_figs/'
def load(f): return np.asarray(Image.open(D+f).convert('RGB')).astype(int)
def frame(dark, box):
    x0,x1,y0,y1=box  # approximate
    def refine_col(x):
        c=[(dark[y0+20:y1-20,xx].sum(),xx) for xx in range(x-4,x+5)]; return max(c)[1]
    def refine_row(y):
        c=[(dark[yy,x0+20:x1-20].sum(),yy) for yy in range(y-4,y+5)]; return max(c)[1]
    return refine_col(x0),refine_col(x1),refine_row(y0),refine_row(y1)
def ticks(dark,X0,X1,Y0,Y1):
    # left-edge ticks: inner side run of dark pixels x in [X0+2, X0+6]
    ly=[y for y in range(Y0+1,Y1) if dark[y,X0+2:X0+7].all()]
    bx=[x for x in range(X0+1,X1) if dark[Y1-6:Y1-1,x].all()]
    def group(v):
        g=[];
        for p in v:
            if g and p-g[-1][-1]<=2: g[-1].append(p)
            else: g.append([p])
        return [np.mean(q) for q in g]
    return group(ly),group(bx)
def curve(a,X0,X1,Y0,Y1,masks=()):
    R,G,B=a[...,0],a[...,1],a[...,2]
    blue=(R<90)&(G<90)&(B>170)
    for (mx0,mx1,my0,my1) in masks: blue[my0:my1,mx0:mx1]=False
    out=[]
    for x in range(X0+2,X1-1):
        ys=np.where(blue[Y0+2:Y1-1,x])[0]
        if len(ys): out.append((x,np.median(ys)+Y0+2,ys.min()+Y0+2,ys.max()+Y0+2,len(ys)))
    return np.array(out)
def fit(pix,vals):
    p=np.polyfit(pix,vals,1); res=np.array(vals)-np.polyval(p,pix); return p,np.abs(res).max()
if False:
    print('see digitize run in README')

# ---- panel definitions: (file, approx frame, y-tick (pix->value) labels, x-tick labels, masks) ----
PANELS={
 'fig1a_Twall': ('p4_img0.png',(98,693,15,434),{117.0:250,218.5:200,320.5:150},{246.0:5,395.5:10,544.0:15},[]),
 'fig2a_Twall': ('p5_img0.png',(107,706,15,436),{96.0:300,177.0:250,256.5:200,337.0:150,418.0:100},{259.0:5,409.0:10,559.0:15},[(565,680,88,175)]),
 'fig2b_Pexit': ('p5_img0.png',(884,1483,15,436),{76.0:700,134.5:600,195.0:500,255.0:400,314.0:300,374.0:200},{1029.0:5,1180.0:10,1330.0:15},[(1350,1465,60,145)]),
 'fig2c_G':     ('p5_img0.png',(107,706,523,944),{599.0:1000,675.5:800,752.0:600,828.0:400,903.0:200},{255.0:5,406.0:10,556.0:15},[(578,690,572,658)]),
}
def run():
    res={}
    for name,(f,box,yt,xt,masks) in PANELS.items():
        a=load(f); dark=a.sum(2)<250
        F=frame(dark,box)
        py,ey=fit(list(yt.keys()),list(yt.values())); px,ex=fit(list(xt.keys()),list(xt.values()))
        c=curve(a,*F,masks=masks)
        t=np.polyval(px,c[:,0]); v=np.polyval(py,c[:,1]); vlo=np.polyval(py,c[:,3]); vhi=np.polyval(py,c[:,2])
        res[name]=dict(t=t,v=v,vlo=vlo,vhi=vhi,n=c[:,4],frame=F,yfit=py.tolist(),xfit=px.tolist(),ytick_resid=ey,xtick_resid=ex,
                       dv_per_px=abs(py[0]),dt_per_px=abs(px[0]))
        print(name,'frame',F,'tick-fit max resid y=%.3g x=%.3g'%(ey,ex),'res/px: %.3g units, %.4g s'%(abs(py[0]),abs(px[0])),'npts',len(t),'t-range %.2f-%.2f'%(t.min(),t.max()))
    return res

HDR="""# SOURCE: Hartwig, J., Chung, J.N., Dong, J., Han, B., Wang, H., Darr, S., Taliaferro, M., Jain, S., Doherty, M. (author list per journal page; task brief lists Rame, E.),
#   'Nitrogen flow boiling and chilldown experiments in microgravity using pulse flow and low-thermally conductive coatings',
#   npj Microgravity 8, 33 (2022). DOI: 10.1038/s41526-022-00220-9 (CC-BY 4.0). PDF: lit/raw/hartwig2022_npj.pdf
# DATA AVAILABILITY: paper's only supplementary file is the Nature Reporting Summary (lit/raw/hartwig2022_npj_MOESM1_ESM.pdf); NO numeric time series published -> digitized from raster figures.
# TEST: F2 = Flight_Day1_2 (Table 1): bare tube, continuous flow, 0.05 g (parabolic flight, ZeroG), horizontal; inlet pressure (Table 1, time-avg at test-section inlet) 550 kPa;
#   chilldown time 14.6 s; chilldown mass 0.872 kg; steady-state Re 67324; steady-state G 533.22 kg/m2/s.
# GEOMETRY: SS304 tube, OD 12.70 mm, wall 0.51 mm (ID 11.68 mm); text says 0.914 m long test sections, Fig.5d shows 26.6 in (0.676 m) inside vacuum / '0.61 m test section' in Fig.5a.
#   Wall-TC stations (Fig.5d, distance from flange A): 4,6,8,16,24 in = 0.102,0.152,0.203,0.406,0.610 m; 3 TCs per station (bottom 1-5, side 6-10, top 11-15), soldered on OUTER wall.
#   'Averaged exit wall temperature' = mean of station-5 TCs (5,10,15) at 0.610 m (consistent with Fig.1c knee at ~13.4 s).
# EXTERNAL: test section in vacuum chamber at ~1 Pa (radiation + residual gas parasitics only). Inlet: LN2 subcooled by precooler (inlet T measured 7 cm upstream; value not reported for F2).
#   Exit pressure set by downstream needle valve NV1; flow measured as gas downstream of electric vaporizers (Alicat M3000, 2 x 3000 SLPM) -> G(t) is lagged/smoothed relative to test-section flow.
# INSTRUMENT UNCERTAINTY (paper): G 2.8%, Re 3.3%, chilldown time 2.1%, mass 2.5%; T-type TCs (error bars in Fig.1/2 ~ +/-2-3 K); 16 Hz DAQ.
"""
def write():
    import pandas as pd
    r=run()
    common="""# DIGITIZATION: images extracted losslessly from PDF (PyMuPDF, 1500-px embedded raster); axis calibration = linear LSQ fit to detected frame tick marks;
#   curve = pure-blue pixel mask (R<90,G<90,B>170), per image column the median row (legend boxes masked). Script: analysis/digitize_hartwig2022.py
"""
    for k,x in r.items():
        unit={'Twall':'K','Pexit':'kPa','G':'kg/m2/s'}[k.split('_')[1]]
        df=pd.DataFrame({'t_s':x['t'],f'value_{unit.replace("/","_")}':x['v'],'value_min_in_column':x['vlo'],'value_max_in_column':x['vhi'],'n_blue_px':x['n'].astype(int)})
        with open(f'data/hartwig2022_F2/F2_{k}_raw_columns.csv','w') as fh:
            fh.write(HDR+common+f"# PANEL: {k} (one row per image column; resolution {x['dt_per_px']:.4f} s/px, {x['dv_per_px']:.3g} {unit}/px; tick-fit max residual {x['ytick_resid']:.2g} {unit})\n")
            df.round(4).to_csv(fh,index=False)
    tg=np.round(np.arange(0.1,19.81,0.1),2)
    out=pd.DataFrame({'t_s':tg})
    out['Twall_exit_avg_K']=np.interp(tg,r['fig2a_Twall']['t'],r['fig2a_Twall']['v'])
    t1=r['fig1a_Twall']['t']; out['Twall_exit_avg_fig1a_K']=np.where(tg<=t1.max(),np.interp(tg,t1,r['fig1a_Twall']['v']),np.nan)
    out['P_exit_kPa']=np.interp(tg,r['fig2b_Pexit']['t'],r['fig2b_Pexit']['v'])
    out['G_kgm2s']=np.interp(tg,r['fig2c_G']['t'],r['fig2c_G']['v'])
    d=(out.Twall_exit_avg_K-out.Twall_exit_avg_fig1a_K).dropna(); knee=(out.t_s>=13)&(out.t_s<14)
    dd=(out.Twall_exit_avg_K-out.Twall_exit_avg_fig1a_K)
    with open('data/hartwig2022_F2/F2_bare_continuous_0p05g.csv','w') as fh:
        fh.write(HDR+common+"""# COLUMNS: t_s [s, from valve open per figure axis]; Twall_exit_avg_K [K] outer-wall T, avg of station-5 TCs, Fig.2a (primary);
#   Twall_exit_avg_fig1a_K [K] same quantity independently digitized from Fig.1a (cross-check; ends 18.8 s); P_exit_kPa [kPa abs] test-section exit pressure, Fig.2b;
#   G_kgm2s [kg/m2/s] mass flux from downstream gas flow meters, Fig.2c. Linear interpolation of column data onto 0.1 s grid.
"""+f"""# DIGITIZATION ERROR: Fig1a vs Fig2a (same data, different plots): rms {np.sqrt((d**2).mean()):.2f} K overall; rms {np.sqrt((dd[~knee].dropna()**2).mean()):.2f} K (max {dd[~knee].abs().max():.1f} K) outside quench knee;
#   rms {np.sqrt((dd[knee]**2).mean()):.1f} K (max {dd[knee].abs().max():.1f} K) for 13-14 s knee (equiv. ~0.1 s time shift on ~100 K/s slope). Resolution 0.62 K/px, 0.033 s/px.
#   P: 1.7 kPa/px, line noise/spikes t<2 s are partly unresolved (see raw column min/max). G: 2.6 kg/m2/s per px; G plotted curve is smoothed by authors.
#   Suggested total uncertainty: T +/-1.5 K (+/-3 K instrument) outside knee, +/-0.15 s on knee timing; P +/-5 kPa; G +/-5 kg/m2/s digitization (+2.8% instrument).
""")
        out.round(3).to_csv(fh,index=False)
    print(out.iloc[::20])
if __name__=='__main__' and len(sys.argv)>1 and sys.argv[1]=='write': write()
