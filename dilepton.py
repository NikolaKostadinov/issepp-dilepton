import argparse
import os

import ROOT


def find_lepton_antilepton_events(input_path, tree_name, toler):
    
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
        .Filter("lepton1.charge * lepton2.charge < 0", "opposite charge")
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
    parser.add_argument("input",         type=str,   default="data/dilepton.root", help="path to the ROOT file")
    parser.add_argument("output",        type=str,   default="dat/output.root",    help="output file")
    parser.add_argument("-t", "--tree",  type=str,   default="DATA",               help="name of the TTree in the file")
    parser.add_argument("-e", "--toler", type=float, default=0.01,                 help="relative tolerance for treating two lepton masses as equal (default 1%%)")
    parser.add_argument("-b", "--bins",  type=int,   default=100,                  help="number of bins in the m_ll spectrum")
    parser.add_argument("-m", "--min",   type=float, default=0.0,                  help="lower edge of the m_ll spectrum, in GeV")
    parser.add_argument("-M", "--max",   type=float, default=200.0,                help="upper edge of the m_ll spectrum, in GeV")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    base, _   = os.path.splitext(args.output)
    png_path  = f"{base}.png"
    root_path = f"{base}.root"

    df = find_lepton_antilepton_events(args.input, args.tree, args.toler)
    n_matched = df.Count()
    hist = df.Histo1D(
        ("mll_spectrum", "Dilepton invariant mass;m_{ll} / GeV;Events", args.bins, args.min, args.max),
        "mll",
    )

    n_matched = n_matched.GetValue()
    n_total   = ROOT.RDataFrame(args.tree, args.input).Count().GetValue()

    bw = ROOT.TF1("bw", "[0]*TMath::BreitWigner(x,[1],[2])", args.min, args.max)
    bw.SetParameters(hist.GetMaximum(), 91.1876, 2.4952)
    bw.SetParNames("Norm", "Mass", "Width")
    hist.Fit(bw, "RQ")

    canvas = ROOT.TCanvas("c", "c", 800, 600)
    hist.Draw()
    bw.Draw("same")
    canvas.SaveAs(png_path)

    out_file = ROOT.TFile(root_path, "RECREATE")
    hist.Write()
    bw.Write()
    out_file.Close()

    mass,  mass_err  = bw.GetParameter(1), bw.GetParError(1)
    width, width_err = bw.GetParameter(2), bw.GetParError(2)

    print(f"M_Z = {mass:.3f} +- {mass_err:.3f} GeV")
    print(f"Γ_Z = {width:.3f} +- {width_err:.3f} GeV")
    print(f"{n_matched}/{n_total} events")


if __name__ == "__main__":
    main()
