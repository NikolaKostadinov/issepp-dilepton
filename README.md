# ISSEP Dilepton Project

We have a simulation data from a $p + p$ collider.
We are going to analyze dilepton events in order to show how one could rediscover the Z boson, since it has a high branching ratio to decay into a lepton-antilepton pair.

$$
Z \rightarrow l + \bar{l}
$$

We are going to filter out events where we have a dilepton pair.
This occurs when the two leptons have the same mass $m_{l} = m_{\bar{l}}$ and opposite charge $q_{l} = -q_{\bar{l}}$.
After that we are going to calculate the dilepton 4-momentum which is the same as the sum of their 4-momenta.
Since 4-momentum is conserved due to relativity we can say that if the pair was produced by an unknown particle, its 4-momentum is going to equal that of the dilepton.
In order to identify this possible unknown particle we will use the spectrum of the dilepton mass.

$$
p_{Z} = p_{l} + p_{\bar{l}}
$$

The spectrum should be a Breit-Wigner distribution from a theoretical standpoint.
This distribution can be derived from the optical theorem applied for the amplitude of the one to two decay.
For a unstable particle with mass $m$ and a decay width $\Gamma$ we can find the "probability" of the decay by squaring the Feynman amplitude which is shown in the next formula.
We are going to fit the result and extract the mass and decay width of the resonance which will be enough to prove the observation of the Z boson in the toy data.

$$
\left|\frac{i}{p^2 - m^2 + im\Gamma}\right|^2
$$

The results from the given data is a Z boson with the following mass and decay width:

$$
\begin{aligned}
m_Z      &= 91.3 \pm 0.1 \mathrm{GeV} \\
\Gamma_Z &= 2.6 \pm 0.1 \mathrm{GeV}
\end{aligned}
$$

![Dilepton Spectrum](assets/spec.png)

However the toy data is not physical.
This is because the leptons have random masses so we could not match them to the Standard model.
The signal is also very clean and one does not need to filter out the background.

## `dilepton.py`

In this implementation we have used `PyROOT`.
Use the `--help` flag to check the inline arguments.

```shell
python dilepton.py --help
```

## Misc

All `.root` data is in `data/` directory which is ignored in `git`.
