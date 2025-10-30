#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os, sys, contextlib
import argparse
import logging
import uuid
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
from array import array


from dateutil.relativedelta import relativedelta
from prettytable import PrettyTable
from tqdm import tqdm

# silence ROOT's import-time C++ spam (rootmap collisions, etc.)
def _import_root_quietly():
    stderr_fd = sys.stderr.fileno()
    saved_stderr_fd = os.dup(stderr_fd)           # duplicate real stderr
    devnull_fd = os.open(os.devnull, os.O_WRONLY) # open /dev/null
    os.dup2(devnull_fd, stderr_fd)                # redirect C-level stderr -> /dev/null
    os.close(devnull_fd)
    try:
        import ROOT as _r
    finally:
        # restore real stderr
        os.dup2(saved_stderr_fd, stderr_fd)
        os.close(saved_stderr_fd)
    return _r

r = _import_root_quietly()
r.gErrorIgnoreLevel = r.kError    # or r.kFatal
r.gROOT.SetBatch(True)


# ======================================================================
# Data classes
# ======================================================================
@dataclass
class GaussParams:
    mean: float
    sigma_percent: float
    binning: int


@dataclass
class Rates:
    total: int
    positron: int
    electron: int
    photon: int
    muon: int
    proton: int
    pion_minus: int
    pion_plus: int

    @property
    def fractions(self) -> Dict[str, float]:
        t = max(self.total, 1)
        return {
            "e+": 100.0 * self.positron / t,
            "e-": 100.0 * self.electron / t,
            "#gamma": 100.0 * self.photon / t,
            "#mu-": 100.0 * self.muon / t,
            "p": 100.0 * self.proton / t,
            "#pi-": 100.0 * self.pion_minus / t,
            "#pi+": 100.0 * self.pion_plus / t,
        }


# ======================================================================
# Helpers
# ======================================================================
def human_diff(t0: datetime, t1: datetime) -> str:
    d = relativedelta(t1, t0)
    return f"{d.hours}h {d.minutes}m {d.seconds}s"


def load_tdr_style(path: Path | None) -> None:
    if not path:
        return
    if not path.exists():
        logging.warning("TDR style file not found: %s", path)
        return
    try:
        if path.suffix.lower() == ".py":
            code = compile(path.read_text(), str(path), "exec")
            exec(code, globals(), globals())
        else:
            r.gROOT.ProcessLine(f'.x {str(path)}')
        logging.info("Loaded TDR style from %s", path)
    except Exception as exc:  # noqa: BLE001
        logging.warning("Failed to load TDR style: %s", exc)


def read_gauss_table(table_path: Path, energy_label: str) -> Optional[GaussParams]:
    """Return GaussParams from table, or None if not found."""
    if not table_path.exists():
        return None
    with table_path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                e, mean, sigma, bins = [x.strip() for x in line.split("|")]
            except ValueError:
                continue
            if e == energy_label:
                return GaussParams(float(mean), float(sigma), int(float(bins)))
    return None


def adaptive_bins(
    n_entries: int,
    x_min: float,
    x_max: float,
    mean: float,
    rms: float | None = None,
    iqr: float | None = None,
    min_bins: int = 15,
    max_bins: int = 320,
) -> int:
    span = max(x_max - x_min, 1e-9)
    N = max(int(n_entries), 0)

    # Robust sigma estimate from the narrow window
    if (rms is not None) and (rms > 0):
        sigma = rms
    elif (iqr is not None) and (iqr > 0):
        sigma = iqr / 1.349  # IQR -> sigma for normal
    else:
        sigma = span / 6.0  # assume the window roughly covers +-3sigma

    # Base: bins-per-sigma target
    bins_per_sigma = 30.0
    b = bins_per_sigma * (span / max(sigma, 1e-9))

    # Gentle sqrt(N) modulation
    if N > 0:
        scale = (N ** 0.5) / 60.0
        if scale < 0.6: scale = 0.6
        if scale > 1.8: scale = 1.8
        b *= scale

    # Absolute minimum bin width cap
    min_width_abs = max(0.02, 0.05 * sigma)
    b_cap = span / min_width_abs
    b = min(b, b_cap)

    # Final clamp & round
    b = int(round(b))
    return max(min_bins, min(max_bins, b))



def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)



def legend(x1=0.67, y1=0.80, x2=0.99, y2=0.97, text_size=0.022) -> r.TLegend:
    leg = r.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(5)
    leg.SetTextSize(text_size)
    return leg


def get_tree_any(tfile: r.TFile, names: List[str]) -> Optional[r.TTree]:
    for name in names:
        obj = tfile.Get(name)
        if isinstance(obj, r.TTree):
            return obj
    return None


# ======================================================================
# Source abstraction
# ======================================================================
class Source:
    """Wraps either a GoodParticle tree (prefixed) or a detector/VirtualDetector tree (unprefixed)."""

    def __init__(self, tree: r.TTree, det: str, mode: str):
        assert mode in ("good", "vd")
        self.tree = tree
        self.det = det
        self.mode = mode
        self.prefix = f"{det}_" if mode == "good" else ""
        # GoodParticle: GeV, VD: MeV
        self.mom_scale = 1.0 if mode == "good" else 1000.0

    def expr_ptot(self) -> str:
        p = self.prefix
        return f"(sqrt({p}Px*{p}Px + {p}Py*{p}Py + {p}Pz*{p}Pz)/{self.mom_scale})"

    def expr_x(self) -> str:
        return f"{self.prefix}x"

    def expr_y(self) -> str:
        return f"{self.prefix}y"

    def expr_initz(self) -> str:
        return f"{self.prefix}InitZ"

    def expr_pdgid(self) -> str:
        return f"{self.prefix}PDGid"

    def has_goodparticle_branches(self) -> bool:
        if self.mode != "good":
            return False
        need = [f"{self.prefix}{b}" for b in ("x", "y", "Px", "Py", "Pz", "PDGid", "InitZ")]
        blist = self.tree.GetListOfBranches()
        return all(bool(blist.FindObject(n)) for n in need)

    def draw_1d(self, expr: str, name: str, title: str, nb: int, xmin: float, xmax: float, cut: str = "") -> r.TH1:
        max_entries = int(os.environ.get("MAX_ENTRIES", "0"))
        emask = f"(Entry$<{max_entries})" if max_entries > 0 else ""
        fullcut = emask if not cut else (f"{cut} && {emask}" if emask else cut)
        h = r.TH1F(name, title, nb, xmin, xmax)
        self.tree.Draw(f"{expr}>>{name}", fullcut, "goff")
        return h

    def get_entries(self, cut: str = "") -> int:
        return int(self.tree.GetEntries(cut))


def resolve_source(tfile: r.TFile, det: str, prefer_good: bool) -> Source:
    if prefer_good:
        gp_tree = get_tree_any(tfile, ["NTuples/GoodParticle", "GoodParticle"])
        if isinstance(gp_tree, r.TTree):
            src_tmp = Source(gp_tree, det, "good")
            if src_tmp.has_goodparticle_branches():
                logging.info("Detector '%s': using GoodParticle (prefixed).", det)
                return src_tmp
            else:
                logging.info("Detector '%s': GoodParticle found but no '%s_*' branches; will try VD.", det, det)

    det_tree = get_tree_any(
        tfile,
        [
            f"VirtualDetector/{det}",
            f"NTuples/{det}",
            det,
        ],
    )
    if isinstance(det_tree, r.TTree):
        # see if this detector tree is also prefixed (rare but possible)
        blist = det_tree.GetListOfBranches()
        has_pref = all(bool(blist.FindObject(f"{det}_{b}")) for b in ("x", "y", "Px", "Py", "Pz", "PDGid", "InitZ"))
        if has_pref:
            logging.info("Detector '%s': using detector tree '%s' (prefixed).", det, det_tree.GetName())
            return Source(det_tree, det, "good")
        logging.info("Detector '%s': using detector tree '%s' (VD style).", det, det_tree.GetName())
        return Source(det_tree, det, "vd")

    raise RuntimeError(f"Detector '{det}' not found in any of the expected locations.")


# ======================================================================
# Auto Gauss (NO energy-specific branches)
# ======================================================================
def auto_gauss_from_source(
    src: Source,
    energy_hint: float,
    search_low_frac: float = 0.85,
    search_high_frac: float = 1.00,
    wide_frac: float = 1.20,
    narrow_frac: float = 0.05,
) -> GaussParams:

    E = float(energy_hint)
    if E <= 0:
        E = 100.0  # just to stay sane

    # 1) wide projection
    wide_xmin = 0.0
    wide_xmax = E * wide_frac
    nbins_wide = 240  # single constant, no energy branching
    hwide_name = f"ptot_wide_{src.det}_{uuid.uuid4().hex}"
    h_wide = src.draw_1d(src.expr_ptot(), hwide_name, "p_tot_wide", nbins_wide, wide_xmin, wide_xmax)

    if h_wide.GetEntries() == 0:
        logging.warning("(%s) Auto-Gauss: wide projection is empty. Using E=%g as mean.", src.det, E)
        return GaussParams(mean=E, sigma_percent=1.0, binning=95)

    # 2) find peak only in 0.85E -> 1.00E (or user-specified)
    search_lo = search_low_frac * E
    search_hi = search_high_frac * E
    # make sure we stay > 0
    search_lo = max(search_lo, 1.0)

    bin_lo = h_wide.FindBin(search_lo)
    bin_hi = h_wide.FindBin(search_hi)

    peak_bin = -1
    peak_val = -1.0
    for b in range(bin_lo, bin_hi + 1):
        val = h_wide.GetBinContent(b)
        if val > peak_val:
            peak_val = val
            peak_bin = b

    # fallback: if nothing in that window (e.g. someone gave E=5 GeV but file is 1 GeV),
    # expand slightly but still without E-dependent if/else trees.
    if peak_bin < 0 or peak_val <= 0.0:
        logging.warning("(%s) No entries found in [%.3f, %.3f] GeV; expanding to [0.70E, 1.05E].",
                        src.det, search_lo, search_hi)
        search_lo2 = 0.70 * E
        search_hi2 = 1.05 * E
        bin_lo2 = h_wide.FindBin(search_lo2)
        bin_hi2 = h_wide.FindBin(search_hi2)
        for b in range(bin_lo2, bin_hi2 + 1):
            val = h_wide.GetBinContent(b)
            if val > peak_val:
                peak_val = val
                peak_bin = b

    # absolute fallback -> global max
    if peak_bin < 0 or peak_val <= 0.0:
        peak_bin = h_wide.GetMaximumBin()
        peak_val = h_wide.GetBinContent(peak_bin)
        logging.warning("(%s) Using global maximum at %.3f GeV for Gauss.", src.det, h_wide.GetBinCenter(peak_bin))

    peak_x = h_wide.GetBinCenter(peak_bin)

    # 3) narrow window around the peak
    half_width = max(narrow_frac * E, 0.5)
    half_width = min(half_width, 4.0)  # don't let "narrow" be wider than +-4 GeV
    narrow_xmin = max(peak_x - half_width, 0.0)
    narrow_xmax = peak_x + half_width
    nbins_narrow = 120
    hn_name = f"ptot_narrow_{src.det}_{uuid.uuid4().hex}"
    h_narrow = src.draw_1d(src.expr_ptot(), hn_name, "p_tot_narrow", nbins_narrow, narrow_xmin, narrow_xmax)

    if h_narrow.GetEntries() == 0:
        logging.warning("(%s) Auto-Gauss: narrow projection empty. Using peak_x=%.3f.", src.det, peak_x)
        return GaussParams(mean=peak_x, sigma_percent=1.0, binning=95)

    mean = h_narrow.GetMean()
    rms = h_narrow.GetRMS()

    # robust spread via IQR (25–75%)
    probs = array('d', [0.25, 0.75])
    quans = array('d', [0.0, 0.0])
    try:
        h_narrow.GetQuantiles(2, quans, probs)
        iqr = max(quans[1] - quans[0], 0.0)
    except Exception:
        iqr = 0.0

    sigma_percent = (rms / mean * 100.0) if mean > 0 else 1.0

    # keep sigma in a sensible envelope but NOT energy-dependent
    sigma_percent = max(sigma_percent, 0.02)   # 0.02% minimum
    sigma_percent = min(sigma_percent, 5.0)    # 5% maximum

    # binning: adapt to actual stats + window size
    entries_narrow = h_narrow.GetEntries()
    logging.debug(
        "(%s auto) narrow entries=%d, window=[%.3f, %.3f], mean=%.3f",
        src.det, entries_narrow, narrow_xmin, narrow_xmax, mean
    )

    if entries_narrow <= 0:
        # fall back: estimate bins from the window even if it's empty
        binning = adaptive_bins(
            n_entries=1,  # fallback
            x_min=narrow_xmin,
            x_max=narrow_xmax,
            mean=mean,
            rms=rms,
            iqr=iqr,
            min_bins=15,
            max_bins=320,
        )
        logging.warning("(%s auto) narrow hist empty -> using fallback bins=%d", src.det, binning)
    else:
        binning = adaptive_bins(
            n_entries=entries_narrow,
            x_min=narrow_xmin,
            x_max=narrow_xmax,
            mean=mean,
            rms=rms,
            iqr=iqr,  # gentle floor
            min_bins=15,
            max_bins=320,
        )

    # if very low energy, don't go crazy
    if E < 40 and binning > 40:
        binning = 40


    logging.info("(%s auto) peak=%.3f GeV -> mean=%.3f GeV, RMS=%.4f GeV, #sigma%%=%.4f, bins=%d",
                 src.det, peak_x, mean, rms, sigma_percent, binning)

    return GaussParams(mean=mean, sigma_percent=sigma_percent, binning=binning)


# ======================================================================
# Physics / drawing
# ======================================================================
def collect_rates(src: Source, ecut_gev: float) -> Rates:
    pexpr = src.expr_ptot()
    cut_le = f"({pexpr}) <= {ecut_gev:.6f}" if ecut_gev < 9e9 else "1"
    pid = src.expr_pdgid()

    def n(sel: str) -> int:
        return src.get_entries(sel)

    total = src.get_entries()
    return Rates(
        total=total,
        positron=n(f"{pid}==-11 && ({cut_le})"),
        electron=n(f"{pid}==11  && ({cut_le})"),
        proton=n(f"{pid}==2212 && ({cut_le})"),
        photon=n(f"{pid}==22   && ({cut_le})"),
        muon=n(f"{pid}==13    && ({cut_le})"),
        pion_minus=n(f"{pid}==-211 && ({cut_le})"),
        pion_plus=n(f"{pid}==211  && ({cut_le})"),
    )


def draw_total_momentum_with_gauss(src: Source, gp: GaussParams, out_dir: Path) -> None:
    # core_half = 1sigma around the peak
    core_half = gp.mean * gp.sigma_percent / 100.0
    # show more than the full width (2*core_half)
    margin_frac = 0.75   # 50% more on each side of the core
    half_span = core_half * (1.0 + margin_frac)
    xmin = gp.mean - (2*half_span)
    xmax = gp.mean + (2*half_span)

    hname = f"h_ptot_{src.det}_{uuid.uuid4().hex}"
    hist = src.draw_1d(src.expr_ptot(), hname, f"Total Momentum @ {src.det}", gp.binning, xmin, xmax)

    c = r.TCanvas(f"c_ptot_{src.det}", f"c_ptot_{src.det}", 900, 1200)

    if hist.GetEntries() == 0:
        # widen a bit if the narrow window missed everything, but keep it centered
        span = xmax - xmin
        extra = 2.00 * span  # +30% both sides
        hist = src.draw_1d(
            src.expr_ptot(),
            hname + "_wide",
            f"Total Momentum @ {src.det}",
            gp.binning,
            xmin - (2*extra),
            xmax + extra,
        )

    hist.GetYaxis().SetTitle("Number of Events")
    hist.GetXaxis().SetTitle("P_{Total} [GeV]")
    hist.SetLineColor(r.kBlue + 1)
    hist.SetFillColor(r.kBlue + 1)
    hist.SetFillStyle(3001)

    # ------------------------------------------------------------------
    # ROBUST PEAK-LOCKED FIT
    # ------------------------------------------------------------------
    hist_sigma = hist.GetStdDev()
    nbins = hist.GetNbinsX()
    global_max = hist.GetMaximum()

    # 1) find all "very high" bins in the CURRENT window
    #    (lower this 0.90 → 0.88 if your peaks are very noisy)
    top_cut = 0.97 * global_max
    candidate_bins = [b for b in range(1, nbins + 1)
                      if hist.GetBinContent(b) >= top_cut]

    if candidate_bins:
        # take the RIGHTMOST high bin → usually the higher-mass peak
        peak_bin = max(candidate_bins)
    else:
        # fallback: ROOT's default
        peak_bin = hist.GetMaximumBin()

    peak_mu = hist.GetBinCenter(peak_bin)
    peak_height = hist.GetBinContent(peak_bin)

    # 2) get a sensible sigma scale (from auto or from data)
    peak_sigma = gp.mean * gp.sigma_percent / 100.0
    if peak_sigma <= 0:
        peak_sigma = hist_sigma if hist_sigma > 0 else 1.0
    # also guard against too-small autos
    if peak_sigma < 0.4:
        peak_sigma = 0.4

    # 3) define a *narrow* fit window around that peak
    #    (this keeps the ugly left shoulder out)
    fit_half = 1.3 * peak_sigma   # you can try 1.0 .. 1.5
    fit_min = peak_mu - fit_half
    fit_max = peak_mu + fit_half

    # 4) estimate local background AS A LINE (not const!)
    left_x  = hist.GetBinCenter(max(1, hist.FindBin(peak_mu - 1.4 * peak_sigma)))
    right_x = hist.GetBinCenter(min(nbins, hist.FindBin(peak_mu + 1.4 * peak_sigma)))
    left_y  = hist.GetBinContent(hist.FindBin(left_x))
    right_y = hist.GetBinContent(hist.FindBin(right_x))
    ped0    = min(left_y, right_y)
    span_x  = max(right_x - left_x, 1e-3)
    slope0  = (right_y - left_y) / span_x   # local slope

    # 5) Gaussian + linear background around the peak
    #    [0] amplitude, [1] mean, [2] sigma, [3] pedestal, [4] slope
    gfun = r.TF1(
        "gaussfit",
        "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*(x-[1])",
        fit_min,
        fit_max,
    )

    # initial parameters
    amp0 = max(peak_height - ped0, 1.0)
    gfun.SetParameter(0, amp0)
    gfun.SetParameter(1, peak_mu)
    gfun.SetParameter(2, peak_sigma)
    gfun.SetParameter(3, ped0)
    gfun.SetParameter(4, slope0)

    # 6) sensible limits
    # μ: don't let it run back to the shoulder
    gfun.SetParLimits(1,
                      peak_mu - 0.25 * peak_sigma,
                      peak_mu + 0.25 * peak_sigma)
    # σ: not too tiny, not too wide
    gfun.SetParLimits(2,
                      0.35 * peak_sigma,
                      2.00 * peak_sigma)
    # pedestal: non-negative, but not insane
    gfun.SetParLimits(3,
                      0.0,
                      2.0 * max(left_y, right_y, 1.0))
    # slope: allow a bit of tilt, but clamp hard so it can't fake the peak
    max_slope = 5.0 * abs(slope0) + 1e-3
    gfun.SetParLimits(4,
                      -max_slope,
                      +max_slope)

    # 7) fit
    if hist.GetEntries() > 0:
        hist.Fit(gfun, "Q")

        # 8) post-check: if μ still drifted to a non-peak place, snap back
        fit_mu = gfun.GetParameter(1)
        fit_sig = gfun.GetParameter(2)
        if abs(fit_mu - peak_mu) > 0.30 * peak_sigma:
            gfun.SetParameter(1, peak_mu)
            gfun.SetParameter(2, peak_sigma)


    ymax = hist.GetMaximum()
    if ymax > 0:
        hist.SetMaximum(ymax * 1.30)  # 30% headroom for legend

    hist.Draw("HIST")
    gfun.Draw("SAME")

    leg = legend()
    leg.AddEntry(hist, "Monte Carlo [MC]", "f")
    leg.AddEntry(hist, f"Mean [data]: {hist.GetMean():.2f}", "")
    leg.AddEntry(hist, f"Std. dev. [data]: {hist.GetStdDev():.2f}", "")
    leg.AddEntry(gfun, f"Fit #mu: {gfun.GetParameter(1):.2f}", "l")
    leg.AddEntry(gfun, f"Fit #sigma: {gfun.GetParameter(2):.4f}", "")
    leg.Draw()

    c.SaveAs(str(out_dir / f"PTotal_{src.det}.pdf"))
    c.Close()


def label_symbol(label: str) -> str:
    return {
        "Positron": "e+",
        "Electron": "e-",
        "Photon": "#gamma",
        "Proton": "p",
    }.get(label, label)


def draw_particle_stacks(src: Source, gp: GaussParams, out_dir: Path, particle_rates: Dict[str, float]) -> None:
    ths_ptot = r.THStack("PtotalStack", f"Total Momentum Spectrum at {src.det};P_{{Total}} [GeV];Number Of Events")
    ths_x = r.THStack("XStack", f"X - Beam Position at {src.det};X-Position [mm];Number Of Events")
    ths_y = r.THStack("YStack", f"Y - Beam Position at {src.det};Y-Position [mm];Number Of Events")
    ths_z = r.THStack("InitZStack", f"Beam Particle Production Areas - {src.det};Z-Position [m];Number Of Events")

    pdgs = [(-11, "Positron"), (11, "Electron"), (22, "Photon"), (2212, "Proton")]
    cols = [42, 46, 32, 36]

    leg_ptot = legend(0.42, 0.76, 0.74, 0.90)
    leg_x = legend(0.73, 0.80, 0.99, 0.97, 0.019)
    leg_y = legend(0.73, 0.80, 0.99, 0.97, 0.019)
    leg_z = legend(0.73, 0.80, 0.99, 0.97, 0.019)

    xmin = 0.0
    xmax = gp.mean + (gp.mean * gp.sigma_percent / 100.0) + 3.0
    pid = src.expr_pdgid()

    for (pdg, label), col in zip(pdgs, cols):
        h_ptot = r.TH1F(
            f"h_ptot_{src.det}_{pdg}_{uuid.uuid4().hex}",
            f"Total Momentum Spectrum of {label}",
            gp.binning,
            xmin,
            xmax,
        )
        src.tree.Draw(f"{src.expr_ptot()}>>{h_ptot.GetName()}", f"{pid}=={pdg}", "goff")
        h_ptot.SetLineColor(col); h_ptot.SetFillColor(col); h_ptot.SetFillStyle(3001); h_ptot.SetNdivisions(-502)
        ths_ptot.Add(h_ptot)
        leg_ptot.AddEntry(h_ptot, f"{label} | {particle_rates.get(label_symbol(label), 0.0):.4f}%", "f")

        h_x = r.TH1F(f"h_x_{src.det}_{pdg}_{uuid.uuid4().hex}", f"X - Beam Profile of {label}", 200, -100, 100)
        src.tree.Draw(f"{src.expr_x()}>>{h_x.GetName()}", f"{pid}=={pdg}", "goff")
        h_x.SetLineColor(col); h_x.SetFillColor(col); h_x.SetFillStyle(3001)
        ths_x.Add(h_x)
        leg_x.AddEntry(h_x, f"{label} | {particle_rates.get(label_symbol(label), 0.0):.4f}%", "f")

        h_y = r.TH1F(f"h_y_{src.det}_{pdg}_{uuid.uuid4().hex}", f"Y - Beam Profile of {label}", 200, -100, 100)
        src.tree.Draw(f"{src.expr_y()}>>{h_y.GetName()}", f"{pid}=={pdg}", "goff")
        h_y.SetLineColor(col); h_y.SetFillColor(col); h_y.SetFillStyle(3001)
        ths_y.Add(h_y)
        leg_y.AddEntry(h_y, f"{label} | {particle_rates.get(label_symbol(label), 0.0):.4f}%", "f")

        h_z = r.TH1F(f"h_z_{src.det}_{pdg}_{uuid.uuid4().hex}", f"Beam Particle Production Area of {label}", 1010, -10, 1000)
        src.tree.Draw(f"{src.expr_initz()}/1000.0>>{h_z.GetName()}", f"{pid}=={pdg}", "goff")
        h_z.SetLineColor(col); h_z.SetFillColor(col); h_z.SetFillStyle(3001)
        ths_z.Add(h_z)
        leg_z.AddEntry(h_z, f"{label} | {particle_rates.get(label_symbol(label), 0.0):.4f}%", "f")

    def save_stack(stack: r.THStack, leg: r.TLegend, fname: str, logy: bool = False, w: int = 900, h: int = 1200) -> None:
        c = r.TCanvas(f"c_{fname}", f"c_{fname}", w, h)
        if logy:
            c.SetLogy(True)
        stack.Draw("nostack")
        leg.Draw()
        c.SaveAs(str(out_dir / f"{fname}_{src.det}.pdf"))
        c.Close()

    save_stack(ths_ptot, leg_ptot, "PTotalStack", logy=True)
    save_stack(ths_x, leg_x, "XPositionStack")
    save_stack(ths_y, leg_y, "YPositionStack")
    save_stack(ths_z, leg_z, "InitZStack", w=1200, h=900)


def draw_xy_positions(src: Source, out_dir: Path, rates: Rates) -> None:
    def _draw(axis: str, fname: str, ndy: int) -> None:
        name = f"h_{axis}_{src.det}_{uuid.uuid4().hex}"
        expr = src.expr_x() if axis == "x" else src.expr_y()
        h = r.TH1F(name, f"Beam position {src.det} - All Particles", 200, -100, 100)
        src.tree.Draw(f"{expr}>>{name}", "", "goff")
        h.GetYaxis().SetTitle("Number of Events")
        h.GetXaxis().SetTitle(f"{axis.upper()} Position (mm)")
        h.GetYaxis().SetTitleOffset(1.3)
        h.GetXaxis().SetTitleOffset(0.8)
        h.SetLineColor(r.kAzure + 2)
        h.SetFillColor(r.kAzure + 2)
        h.SetFillStyle(3001)

        c = r.TCanvas(f"c_{fname}", f"c_{fname}", 900, 1200)
        h.Draw("HIST")

        leg = legend(0.68, 0.85, 0.98, 0.93, 0.022 if axis == "x" else 0.025)
        leg.AddEntry(h, f"Entries: {h.GetEntries():.0f}", "f")
        if axis == "x":
            leg.AddEntry(h, f"Mean: {h.GetMean():.3e}", "")
        else:
            leg.AddEntry(h, f"Mean: {h.GetMean():.2f}", "")
        leg.AddEntry(h, f"Std. Dev.: {h.GetStdDev():.3e}", "")
        leg.Draw()

        pave = r.TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
        pave.SetFillColor(0); pave.SetBorderSize(0); pave.SetFillStyle(0)
        pave.SetTextAlign(12); pave.SetTextSize(0.025)
        t1 = pave.AddText("Particle Rates"); t1.SetTextColor(r.kRed + 3); t1.SetTextSize(0.03)
        pave.AddText(" ")
        fr = rates.fractions
        pave.AddText(f"#gamma: {fr['#gamma']:.3f}%")
        pave.AddText(f"e^{{+}}: {fr['e+']:.3f}%")
        pave.AddText(f"e^{{-}}: {fr['e-']:.3f}%")
        pave.AddText(f"#mu^{{-}}: {fr['#mu-']:.3f}%")
        pave.AddText(f"Proton: {fr['p']:.3f}%")
        pave.Draw()

        c.SaveAs(str(out_dir / f"{fname}_{src.det}.pdf"))
        c.Close()

    _draw("x", "XPos", ndy=16)
    _draw("y", "YPos", ndy=32)


def draw_y_vs_ptot(src: Source, energy: float, gp: GaussParams, out_dir: Path, rates: Rates) -> None:
    if src.det == "ECALFront":
        bin_y = 300; ylo = -150; yhi = 150
    elif src.det == "BeamStartPoint":
        bin_y = 6; ylo = -3; yhi = 3
    else:
        bin_y = 200; ylo = -100; yhi = 100

    nb = gp.binning + 10
    name = f"h_y_vs_ptot_{src.det}_{uuid.uuid4().hex}"
    h2 = r.TH2F(name, f"Beam Total Momentum vs Beam Y Position - {src.det}", nb, 0, energy * 1.05, bin_y, ylo, yhi)
    src.tree.Draw(f"{src.expr_y()}:{src.expr_ptot()}>>{name}", "", "goff")
    h2.GetYaxis().SetTitle("Beam Y-Position [mm]")
    h2.GetXaxis().SetTitle("P_{Total} [GeV]")

    c = r.TCanvas(f"c_y_vs_ptot_{src.det}", f"c_y_vs_ptot_{src.det}", 900, 1200)
    c.SetRightMargin(0.18)
    c.SetLogz(True)
    h2.Draw("COLZ")
    leg = legend(0.65, 0.90, 0.91, 0.95)
    leg.AddEntry(h2, f"Entries: {h2.GetEntries():.0f}", "p")
    leg.Draw()

    pave = r.TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
    pave.SetFillColor(0); pave.SetBorderSize(0); pave.SetFillStyle(0)
    pave.SetTextAlign(12); pave.SetTextSize(0.025)
    t1 = pave.AddText("Particle Rates"); t1.SetTextColor(r.kRed + 3); t1.SetTextSize(0.03)
    fr = rates.fractions
    pave.AddText(" ")
    pave.AddText(f"#gamma: {fr['#gamma']:.3f}%")
    pave.AddText(f"e^{{+}}: {fr['e+']:.3f}%")
    pave.AddText(f"e^{{-}}: {fr['e-']:.3f}%")
    pave.AddText(f"#mu^{{-}}: {fr['#mu-']:.3f}%")
    pave.AddText(f"Proton: {fr['p']:.3f}%")
    pave.Draw()

    c.SaveAs(str(out_dir / f"YvsPTotal_{src.det}.pdf"))
    c.Close()


# ======================================================================
# CLI / main
# ======================================================================
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Beamline Analyzer (PyROOT)")
    p.add_argument("--initialEnergy", required=True,
                   help="Initial energy of the beam label, e.g. 300, 250, 100, 50")
    p.add_argument("--inputFile", required=True,
                   help="Input ROOT file created with G4BL")
    p.add_argument("--outputFolder", default="Results",
                   help="Output folder for PDFs")
    p.add_argument("--gaussTable", default=None,
                   help="Optional GaussTable (E|mean|sigma|binning). If missing or energy not found, auto-compute.")
    p.add_argument("--detectors", default="ECALFront",
                   help="Comma-separated detector names (ECALFront,BeamStartPoint,...)")
    p.add_argument("--tailInvestigation", action="store_true",
                   help="If set, mask events below the lower edge of the fitted Gaussian")
    p.add_argument("--tdrStyle", default="tdrStyle.py",
                   help="TDR style file name (relative to this script if not absolute)")
    p.add_argument("--goodTree", action="store_true",
                   help="Prefer NTuples/GoodParticle or top-level GoodParticle")
    p.add_argument("--verbose", action="store_true", help="Verbose logging")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    t0 = datetime.now()
    ap = build_argparser()
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(levelname)s] %(message)s",
    )

    # resolve style path relative to script
    style_path = None
    if args.tdrStyle:
        p = Path(args.tdrStyle)
        if not p.is_absolute():
            p = Path(__file__).resolve().parent / p
        style_path = p
    load_tdr_style(style_path)

    # energy hint
    energy_label = str(args.initialEnergy)
    try:
        energy_hint = float(energy_label)
    except ValueError:
        digits = "".join(c for c in energy_label if (c.isdigit() or c == "."))
        energy_hint = float(digits) if digits else 100.0

    infile = Path(args.inputFile)
    if not infile.exists():
        logging.error("Input file does not exist: %s", infile)
        return 2

    tf = r.TFile.Open(str(infile), "READ")
    if not tf or tf.IsZombie():
        logging.error("Failed to open ROOT file: %s", infile)
        return 2

    detectors = [d.strip() for d in args.detectors.split(",") if d.strip()]
    out_root = Path(args.outputFolder)

    table = PrettyTable()
    table.field_names = ["Detector", "ALL", "Positron", "Electron", "Photon", "Muon", "Proton"]

    table_path = Path(args.gaussTable) if args.gaussTable else None

    for det in tqdm(detectors, desc="Detectors"):
        det_dir = out_root / det
        ensure_dir(det_dir)

        try:
            src = resolve_source(tf, det, prefer_good=args.goodTree)
        except Exception as exc:  # noqa: BLE001
            logging.error("Detector '%s': %s - skipping", det, exc)
            continue

        # 1) try table
        gp = None
        if table_path is not None:
            gp = read_gauss_table(table_path, energy_label)

        # 2) else auto
        if gp is None:
            gp = auto_gauss_from_source(src, energy_hint)

        tail_cut = gp.mean - (gp.mean * gp.sigma_percent / 100.0) if args.tailInvestigation else 1e12

        rates = collect_rates(src, tail_cut)
        fr = rates.fractions

        table.add_row([
            det,
            rates.total,
            rates.positron,
            rates.electron,
            rates.photon,
            rates.muon,
            rates.proton,
        ])

        try:
            draw_total_momentum_with_gauss(src, gp, det_dir)
            draw_particle_stacks(src, gp, det_dir, {
                "e+": fr["e+"], "e-": fr["e-"], "#gamma": max(fr["#gamma"], 0.0001), "p": fr["p"],
            })
            draw_xy_positions(src, det_dir, rates)
            draw_y_vs_ptot(src, energy_hint, gp, det_dir, rates)
        except Exception as exc:  # noqa: BLE001
            logging.warning("Plotting failed for %s: %s", det, exc)

    print(table.get_string(title="Particles"))

    tf.Close()
    logging.info("Done in %s", human_diff(t0, datetime.now()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

