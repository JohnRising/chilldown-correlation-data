"""Goal 2: lag between the (downstream) flow-meter G(t) and the test-section response.

tau > 0 means the flow meter lags: G_testsection(t) = G_meas(t + tau).
Estimators (all on 15 Hz data, valid window):
  A  xcorr of high-passed log h_ec (film boiling, each station) with high-passed log G(t+tau)
  B  xcorr of high-passed P1 (upstream station pressure; no transport lag) with high-passed G(t+tau)
  C  per-case regression log h = a + b log G(t+tau) + c log dT (film boiling, both stations);
     tau maximising R^2
High-pass = signal minus Savitzky-Golay trend (window 61 samples ~ 4 s, order 2).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CASES, DT, FIG, OUT, STATIONS, load, md_table, smooth
from prep import film_mask, station_data, t_valid_end

TAUS = np.arange(-2.0, 8.0 + 1e-9, DT)


def hp(y, win=61):
    return y - smooth(y, win, 2)


def shifted(t, y, tau):
    """y(t + tau) by linear interpolation; NaN outside the record."""
    return np.interp(t + tau, t, y, left=np.nan, right=np.nan)


def xcorr(t, a, b, mask, taus=TAUS):
    """corr(a(t), b(t+tau)) over samples where mask and b(t+tau) defined."""
    out = []
    for tau in taus:
        bs = shifted(t, b, tau)
        m = mask & np.isfinite(bs)
        if m.sum() < 30:
            out.append(np.nan); continue
        x = a[m] - a[m].mean(); y = bs[m] - bs[m].mean()
        out.append(np.sum(x * y) / np.sqrt(np.sum(x * x) * np.sum(y * y)))
    return np.array(out)


def peak(taus, r, lo=-1.0, hi=7.0):
    m = (taus >= lo) & (taus <= hi) & np.isfinite(r)
    k = np.where(m)[0][np.nanargmax(r[m])]
    # parabolic refinement
    if 0 < k < len(r) - 1 and np.all(np.isfinite(r[k - 1:k + 2])):
        y0, y1, y2 = r[k - 1:k + 2]
        den = y0 - 2 * y1 + y2
        dk = 0.5 * (y0 - y2) / den if den != 0 else 0
        return taus[k] + dk * DT, r[k]
    return taus[k], r[k]


def main():
    rows = []
    fig, axs = plt.subplots(2, 5, figsize=(20, 7), sharex=True)
    for ax, c in zip(axs.ravel(), CASES):
        d = load(c)
        tend = t_valid_end(d)
        t = d.t_s.values
        lG = np.log(np.clip(d.G_kgm2s.values, 5, None))
        Ghp = hp(lG)
        row = dict(case=c)
        # A: h
        rs = []
        for st in STATIONS:
            s = station_data(c, st)
            m = np.zeros(len(t), bool)
            fm = film_mask(c, st) & (s.t.values > 0.5)
            m[:len(s)] = fm
            lh = np.full(len(t), np.nan); lh[:len(s)] = np.log(np.clip(s.h_ec.values, 1, None))
            lh = np.where(np.isfinite(lh), lh, np.nanmean(lh))
            r = xcorr(t, hp(lh), Ghp, m)
            tau, rmax = peak(TAUS, r)
            row[f"tau_h_{st}"] = tau; row[f"r_h_{st}"] = rmax
            rs.append(r)
            ax.plot(TAUS, r, label=f"h_{st} vs G")
        # B: P1
        mP = (t > 0.5) & (t <= tend)
        rP = xcorr(t, hp(d.P1_kPa.values), Ghp, mP)
        tau, rmax = peak(TAUS, rP)
        row["tau_P"] = tau; row["r_P"] = rmax
        ax.plot(TAUS, rP, label="P1 vs G")
        # C: regression R2
        best = (-np.inf, np.nan)
        r2s = []
        for tauv in TAUS[(TAUS >= 0) & (TAUS <= 6)]:
            X, Y = [], []
            for st in STATIONS:
                s = station_data(c, st)
                fm = film_mask(c, st)
                Gs = shifted(t, d.G_kgm2s.values, tauv)[:len(s)]
                mm = fm & np.isfinite(Gs) & (Gs > 5)
                X.append(np.c_[np.ones(mm.sum()), np.log(Gs[mm]), np.log(s.dT.values[mm]),
                               np.full(mm.sum(), st == "dn")])
                Y.append(np.log(s.h_ec.values[mm]))
            X = np.vstack(X); Y = np.concatenate(Y)
            beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
            r2 = 1 - np.sum((Y - X @ beta) ** 2) / np.sum((Y - Y.mean()) ** 2)
            r2s.append(r2)
            if r2 > best[0]:
                best = (r2, tauv, beta[1])
        row["tau_fit"] = best[1]; row["R2_fit"] = best[0]; row["bG_fit"] = best[2]
        r2_0 = r2s[0]
        row["R2_tau0"] = r2_0
        ests = [row["tau_h_up"], row["tau_h_dn"], row["tau_P"]]
        row["tau_median"] = float(np.median(ests))
        rows.append(row)
        ax.axvline(row["tau_median"], color="k", ls="--", lw=0.8)
        ax.set_title(f"case {c}: tau_med={row['tau_median']:.2f} s")
        ax.axhline(0, color="k", lw=0.4)
        ax.grid(alpha=0.3)
        if c == 1:
            ax.legend(fontsize=7)
    for ax in axs[1]:
        ax.set_xlabel("tau [s]  (G_meas(t+tau) vs response(t))")
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig06_lag_xcorr.png", dpi=100)
    plt.close(fig)
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/lag.csv", index=False)
    # global choice
    tau_g = float(np.median(df[["tau_h_up", "tau_h_dn", "tau_P"]].values))
    with open(f"{OUT}/t02_lag.md", "w") as f:
        f.write(md_table(df.round(3), 3) + "\n\n")
        f.write(f"global median tau (all xcorr estimates) = {tau_g:.2f} s\n")
    print(df.round(2).to_string())
    print("global tau", tau_g)
    # case-10 h-peak vs G-peak
    d = load(10)
    s = station_data(10, "up")
    print("case10: t(h_up max, t<5)=", s.t[(s.t < 5)][s.h_ec[(s.t < 5)].idxmax()],
          " t(G max)=", d.t_s[d.G_kgm2s.idxmax()])


def g_ts(t, G, tau, model="delay"):
    """test-section mass flux reconstructed from the meter signal.
    delay : G_ts(t) = G_meas(t+tau) (held at the last value beyond the record)
    lag1  : meter = first-order lag of G_ts with time constant tau -> G_ts = G_m + tau dG_m/dt
            (dG/dt from a 1-s Savitzky-Golay filter), clipped >= 5."""
    if model == "delay":
        return np.interp(t + tau, t, G)
    Gs = smooth(G, 15, 2)
    return np.clip(Gs + tau * smooth(G, 15, 2, deriv=1), 5, None)


def pooled(tau_of_case, model="delay", fixed=None):
    X, Y, C = [], [], []
    for c in CASES:
        d = load(c)
        t = d.t_s.values
        G = g_ts(t, d.G_kgm2s.values, tau_of_case[c], model)
        for st in STATIONS:
            s = station_data(c, st)
            fm = film_mask(c, st)
            Gs = G[:len(s)]
            mm = fm & (Gs > 5)
            X.append(np.c_[np.ones(mm.sum()), np.log(Gs[mm]), np.log(s.dT.values[mm]),
                           np.full(mm.sum(), st == "dn")])
            Y.append(np.log(s.h_ec.values[mm]))
            C.append(np.full(mm.sum(), c))
    X = np.vstack(X); Y = np.concatenate(Y); C = np.concatenate(C)
    beta = np.linalg.lstsq(X, Y, rcond=None)[0] if fixed is None else fixed
    res = Y - X @ beta
    return beta, res, C


def startup_lags():
    rows = []
    for c in CASES:
        d = load(c)
        t = d.t_s.values
        Gs = smooth(d.G_kgm2s.values, 9, 2)
        early = t < 10
        tG50 = t[np.argmax(Gs >= 0.5 * Gs[early].max())]
        tGpk = t[early][np.argmax(Gs[early])]
        s = station_data(c, "up")
        qs = smooth(s.q_ec.values, 9, 2)
        e2 = s.t.values < 10
        tq50 = s.t.values[np.argmax(qs >= 0.5 * qs[e2].max())]
        hs = smooth(s.h_ec.values, 15, 2)
        e3 = (s.t.values > 0.3) & (s.t.values < 8) & film_mask(c, "up", t_min=0.3)
        thpk = s.t.values[e3][np.argmax(hs[e3])]
        rows.append(dict(case=c, tG50=tG50, tq50=tq50, tau_onset=tG50 - tq50, tGpk=tGpk, thpk=thpk,
                         tau_peak=tGpk - thpk))
    return pd.DataFrame(rows)


def lag_models():
    out = {}
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for k, model in enumerate(("delay", "lag1")):
        taus = np.arange(0, 6.01, 0.25)
        rms = []
        for tau in taus:
            beta, res, _ = pooled({c: tau for c in CASES}, model)
            rms.append(np.sqrt(np.mean(res**2)))
        rms = np.array(rms)
        tau_g = taus[np.argmin(rms)]
        beta_g, _, _ = pooled({c: tau_g for c in CASES}, model)
        # per-case tau with the global coefficients fixed
        per = {}
        for c in CASES:
            best = (np.inf, 0)
            for tau in taus:
                tc = {cc: tau_g for cc in CASES}; tc[c] = tau
                _, res, C = pooled(tc, model, fixed=beta_g)
                r = np.sqrt(np.mean(res[C == c] ** 2))
                if r < best[0]:
                    best = (r, tau)
            per[c] = best[1]
        out[model] = dict(taus=taus, rms=rms, tau_g=tau_g, beta=beta_g, per=per)
        ax[k].plot(taus, rms, "o-")
        ax[k].set_xlabel("tau [s]"); ax[k].set_ylabel("pooled RMS of ln h residual")
        ax[k].set_title(f"{model}: global tau={tau_g:.2f} s, b_G={beta_g[1]:.2f}, c_dT={beta_g[2]:.2f}")
        ax[k].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig07_lag_pooled.png", dpi=100)
    plt.close(fig)
    return out


if __name__ == "__main__":
    main()
    su = startup_lags()
    lm = lag_models()
    su["tau_delay_percase"] = [lm["delay"]["per"][c] for c in CASES]
    su["tau_lag1_percase"] = [lm["lag1"]["per"][c] for c in CASES]
    print(su.round(2).to_string())
    for m in lm:
        print(m, "global tau", lm[m]["tau_g"], "beta", lm[m]["beta"].round(3), "rms0", lm[m]["rms"][0].round(4), "rmsmin", lm[m]["rms"].min().round(4))
    su.to_csv(f"{OUT}/lag_startup.csv", index=False)
    with open(f"{OUT}/t02_lag.md", "a") as f:
        f.write("\n## start-up / peak lags and pooled-regression lags\n\n")
        f.write(md_table(su.round(3), 3) + "\n\n")
        for m in lm:
            f.write(f"- {m}: global tau = {lm[m]['tau_g']:.2f} s, pooled beta (a, bG, cdT, d_dn) = "
                    f"{np.round(lm[m]['beta'], 3).tolist()}, RMS(tau=0) = {lm[m]['rms'][0]:.4f}, "
                    f"RMS(min) = {lm[m]['rms'].min():.4f}\n")
