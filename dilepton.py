import argparse
import os

import ROOT

def define_m_les(input_path, tree_name):

    df = ROOT.RDataFrame(tree_name, input_path)
    df = (
        df.Define("evt", "rdfentry_")
        .Define(
            "m1sq",
            "+lepton1.energy * lepton1.energy"
            "-lepton1.px     * lepton1.px"
            "-lepton1.py     * lepton1.py"
            "-lepton1.pz     * lepton1.pz"
        )
        .Define(
            "m2sq",
            "+lepton2.energy * lepton2.energy"
            "-lepton2.px     * lepton2.px"
            "-lepton2.py     * lepton2.py"
            "-lepton2.pz     * lepton2.pz"
        )
        .Filter("m1sq >= 0 && m2sq >= 0", "reconstructable mass")
        .Define("m1", "sqrt(m1sq)")
        .Define("m2", "sqrt(m2sq)")
    )

    return df

def find_lepton_antilepton_events(df, toler):

    df = (
        df.Filter("lepton1.charge * lepton2.charge < 0", "opposite charge")
        .Filter("abs(lepton1.charge) == abs(lepton2.charge)", "equal absolute charge")
        .Filter(f"std::abs(m1 - m2) <= {toler} * std::max(m1, m2)", "same mass")
        .Define(
            "mll",
            "sqrt(std::max(0.0, "
            "+(lepton1.energy + lepton2.energy) * (lepton1.energy + lepton2.energy)"
            "-(lepton1.px     + lepton2.px)     * (lepton1.px     + lepton2.px)"
            "-(lepton1.py     + lepton2.py)     * (lepton1.py     + lepton2.py)"
            "-(lepton1.pz     + lepton2.pz)     * (lepton1.pz     + lepton2.pz)))"
        )
    )

    return df

def main():

    parser = argparse.ArgumentParser(
        description="Log lepton-antilepton events (opposite charge, matching lepton mass) "
        "from a dilepton ROOT ntuple."
    )
    parser.add_argument("input",              type=str,   default="data/dilepton.root", help="path to the ROOT file")
    parser.add_argument("output",             type=str,   default="dat/output.root",    help="output file")
    parser.add_argument("-t", "--tree",       type=str,   default="DATA",               help="name of the TTree in the file")
    parser.add_argument("-e", "--toler",      type=float, default=0.01,                 help="relative tolerance for treating two lepton masses as equal (default 1%%)")
    parser.add_argument("-b", "--bins",       type=int,   default=100,                  help="number of bins in the l_ll spectrum")
    parser.add_argument("-m", "--min",        type=float, default=0.0,                  help="lower edge of the l_ll spectrum in GeV")
    parser.add_argument("-M", "--max",        type=float, default=200.0,                help="upper edge of the l_ll spectrum in GeV")
    parser.add_argument("-g", "--l-bins",     type=int,   default=800,                  help="number of bins in the single-lepton mass spectrum")
    parser.add_argument("-j", "--l-min",      type=float, default=0.0,                  help="lower edge of the single-lepton mass spectrum in MeV")
    parser.add_argument("-J", "--l-max",      type=float, default=200.0,                help="upper edge of the single-lepton mass spectrum in MeV")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    base, _   = os.path.splitext(args.output)
    root_path = f"{base}.root"

    df_leptons = define_m_les(args.input, args.tree)
    df         = find_lepton_antilepton_events(df_leptons, args.toler)

    n_matched = df.Count()
    hist = df.Histo1D(
        ("mll_spectrum", "Dilepton invariant mass;m_{l\\bar{l}} / GeV;Events", args.bins, args.min, args.max),
        "mll",
    )

    df_leptons_mev = df_leptons.Define("m1_MeV", "m1 * 1000.0").Define("m2_MeV", "m2 * 1000.0")
    hist_m1 = df_leptons_mev.Histo1D(
        ("m_l_1", "Single-lepton mass;m_{l} / MeV;Leptons", args.l_bins, args.l_min, args.l_max),
        "m1_MeV",
    )
    hist_m2 = df_leptons_mev.Histo1D(
        ("m_l_2", "Single-lepton mass;m_{l} / MeV;Leptons", args.l_bins, args.l_min, args.l_max),
        "m2_MeV",
    )

    n_matched = n_matched.GetValue()
    n_total   = ROOT.RDataFrame(args.tree, args.input).Count().GetValue()

    bw = ROOT.TF1("bw", "[0]*TMath::BreitWigner(x,[1],[2])", args.min, args.max)
    bw.SetParameters(hist.GetMaximum(), 91.1876, 2.4952)
    bw.SetParNames("Norm", "Mass", "Width")
    hist.Fit(bw, "RQ")

    hist_lepton = hist_m1.GetValue().Clone("m_l")
    hist_lepton.Add(hist_m2.GetPtr())

    canvas = ROOT.TCanvas("c_ll", "c_ll", 800, 600)
    hist.Draw()
    bw.Draw("same")
    canvas.SaveAs(f"{base}_{hist.GetName()}.png")

    canvas_lepton = ROOT.TCanvas("c_l", "c_l", 800, 600)
    hist_lepton.Draw()
    canvas_lepton.SaveAs(f"{base}_{hist_lepton.GetName()}.png")

    out_file = ROOT.TFile(root_path, "RECREATE")
    hist.Write()
    bw.Write()
    hist_lepton.Write()
    out_file.Close()

    mass,  mass_err  = bw.GetParameter(1), bw.GetParError(1)
    width, width_err = bw.GetParameter(2), bw.GetParError(2)

    chi2    = bw.GetChisquare()
    ndf     = bw.GetNDF()
    p_value = ROOT.TMath.Prob(chi2, ndf)
    
    print()
    print(f".------------------------------------------")
    print(f"| Z boson")
    print(f"+------------------------------------------")
    print(f"| M_Z \t {mass:.2f} +- {mass_err:.2f} GeV")
    print(f"| Γ_Z \t  {width:.2f} +- {width_err:.2f} GeV")
    print()
    print(f".------------------------------------------")
    print(f"| Breit-Wigner Fit")
    print(f"+------------------------------------------")
    print(f"| chi2/ndf \t {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
    print(f"| p value  \t {p_value:.3f}")
    print()
    print(f"{n_matched} / {n_total} dilepton events")
    print()


if __name__ == "__main__":
    main()
