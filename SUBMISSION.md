# Space Apps 2026 — project page text (copy-paste)


---

## Project Name
Spectral-Aware Blink

## High-Level Project Summary
SPHEREx photographs the whole sky through a rainbow filter: every exposure sees a star in a different colour, so blinking SPHEREx images, the classic way to find Planet X candidates, brown dwarfs and variable stars, makes constant stars flicker. We built Spectral-Aware Blink, a pipeline and web tool that learns each source's spectrum from the other survey passes, subtracts it, and blinks only what is left: real motion and real variability. Using only public NASA SPHEREx images (258 exposures, 10 fields, ~5,100 light curves), we show that naive light curves flag 58% of constant stars as periodic (1% after correction), that the suggested detector stacks fake 7–13% variability that we can predict from metadata (r = 0.92), and that cold brown dwarfs, the objects Planet X hunters want most, are the ones colour-blind methods fake the most (26–36%). The coldest known brown dwarf, WISE 0855, appears as a moving dipole that flips sign over the year.

## Link to Project Demo
https://murphyjuan.github.io/spectral-aware-blink/demo/Spectral-Aware-Blink-30s.mp4
*(29-second narrated video with English subtitles, SI voice disclosed, plays in any browser. Alternative: 7-slide deck https://github.com/MurphyJUAN/spectral-aware-blink/blob/main/demo/Spectral-Aware-Blink-slides.pdf . Optionally upload the MP4 to YouTube and use that link.)*

## Link to Final Project
https://github.com/MurphyJUAN/spectral-aware-blink
*(Live tool: https://murphyjuan.github.io/spectral-aware-blink/web/ )*

## Detailed Project Description

**What it does.** Spectral-Aware Blink turns SPHEREx's rainbow-filtered images into a clean "blink" that shows only real change. SPHEREx observes the sky through a linear variable filter: each row of each detector sees a different wavelength, and the telescope steps along that filter between exposures to collect a 102-channel spectrum of everything. The side effect is that the same constant star is measured at a different colour in every exposure. Plotted against time it looks variable by a factor of two; plotted against wavelength it is a smooth spectrum. Our tool removes that colour effect so that citizen scientists, researchers and pipelines can trust what blinks.

**How it works.**
1. *Data access without query services.* IRSA's search services were unavailable during the hackathon, so we read cutouts directly from the public SPHEREx DR1 AWS mirror with HTTP range requests: image, variance, zodiacal-light model, quality flags and the per-pixel wavelength table. We derive the wavelength sampling of any sky position from the FITS headers alone.
2. *Preprocessing.* Zodiacal-light subtraction, a mesh background, and masking of artefacts with SPHEREx's flag bits.
3. *Spectral model with leave-one-pass-out.* For every star (forced aperture photometry) and every sky pixel, we fit a smooth spectrum per detector using all other survey passes, predict the brightness at the exact colour of each exposure, and subtract. Leaving the pass out keeps a moving or varying object from hiding in its own model. Motion then shows up as a +/− dipole and variability as a coherent residual.
4. *Fair benchmark.* We compare four methods on identical data and detectors: naive blink, same-colour pairs (used by several teams), detector-stack differencing (the challenge's suggestion), and our spectral model. Thresholds are set at equal false-alarm probability, and we inject fake moving sources and fake variability into real data.
5. *Web tool (static, no server).* Per-pass blink in four modes for short- and long-wave detectors and four fields. Click any star to see its spectrum, its light curve before and after removing the spectrum, and a metadata-only prediction of how strongly it will fake-blink in a colour-blind comparison. URLs are shareable.

**What we found** (10 independent 17′ fields, ~5,100 star–detector light curves):
- Inside a survey pass, wavelength tracks time (|ρ| = 0.97 in detectors 3 and 6, over 600 sky positions). As a result, a standard periodogram flags 50–67% of constant bright stars as periodic with fake 1–8-day periods; after removing the spectrum the rate falls to 0–5%.
- Detector-stack "variability" (6.7–12.7%) can be predicted from spectrum × wavelength sampling (r = 0.86–0.95). Our correction brings it to the ~3% photometric floor.
- Relative-colour law: after the usual zero-point match, the remaining false variability is proportional to (source colour − field colour) × the colour shift between passes. This held in all 10 fields. For the most colour-peculiar stars our model is 1.3–2.4× cleaner, and for T/Y brown dwarfs we predict 26–36% fake blinks in colour-blind comparisons.
- Same-colour pairing breaks down in detectors 1 and 4 (including the 3.3 µm methane band): at 1% tolerance, 17–23% of positions have no partner exposure at all.
- WISE 0855, the coldest known brown dwarf, appears in 3.9–5 µm residuals as a moving dipole (~8σ combined at its predicted position) that flips sign between the first and last pass.
- Honest limitation: in single-pixel blinking, SPHEREx's undersampled PSF, not colour, limits star-like movers. All colour-aware methods tie, and naive blink is 2–5× worse. Our next step is PSF forward modelling.

**Benefits.** Fewer false alarms for citizen-science blinking, trustworthy SPHEREx variability, and better sensitivity to the coldest, reddest moving objects (Planet X candidates, Y dwarfs). Every correction can be computed from public metadata before touching a pixel. We checked 13 papers. The newest SPHEREx forward-modelling work (Nov 2026) lists proper motion and variable stars as open problems, and the SPHEREx pipeline plans survey-to-survey variability flags that our relative-colour law predicts will over-select colour-peculiar sources unless they are spectrum-aware.

**What we hope to achieve.** An all-sky, spectrum-aware blink service for SPHEREx in the spirit of WiseView and Backyard Worlds, a vetted candidate list of nearby cold movers cross-matched with Gaia/CatWISE, and a contribution of the relative-colour correction to SPHEREx time-domain science.

**Tools.** Python (Astropy, NumPy, SciPy, pandas, Matplotlib, fsspec/aiohttp, Pillow); vanilla HTML/JavaScript/Canvas/SVG for the web tool; GitHub Pages hosting; a cloud GPU node used only as a fast network/CPU worker. SI tools (Claude for implementation planning, Claude Code for code, analysis and drafting) were used on top of the team's own idea (see the SI section).

## Use of Super Intelligence (SI)
Yes. The core idea of this project (removing SPHEREx's rainbow-filter flicker by subtracting a spectral model learned from the other survey passes before blinking) was conceived by Yining Juan together with a discussion partner, without SI. We then used Claude (Anthropic) to expand our idea into a detailed implementation plan, and Claude Code, an SI coding assistant, to write and debug the Python pipeline and the web tool, run and analyse the experiments, search and summarise the literature, produce the figures programmatically, and draft the documentation, slides and submission text. We reviewed the results and decided which claims to make. Every source file carries an "SI-generated" header, and the data products and web pages are labelled the same way in their metadata. Every figure is a plot made by our code from NASA data (no SI image generation). The 30-second demo video uses an SI text-to-speech voice (Microsoft Edge neural voice "Aria" via the edge-tts tool) for the narration; its visuals and its background music were generated by our SI-written code (original music, no copyrighted material). The video carries a visible "SI-assisted" watermark and says so in its file metadata.

## NASA Data
- **SPHEREx DR1 Level-2 Spectral Images** (NASA / Caltech-IPAC IRSA), the core input of every result. We used the IMAGE, VARIANCE, ZODI (zodiacal-light model), FLAGS and WCS-WAVE (per-pixel wavelength) extensions of 258 exposures covering WISE 0855−0714 (May 2025 – May 2026), read by HTTP range requests from the public AWS mirror.
  - IRSA SPHEREx mission page: https://irsa.ipac.caltech.edu/Missions/spherex.html
  - Public AWS S3 mirror: https://nasa-irsa-spherex.s3.amazonaws.com/ (prefix `dr1/level2/`)
- **IRSA SPHEREx Simple Image Access (SIA) service**, used to list the exposures covering WISE 0855 (results saved on 2025-11-14): https://irsa.ipac.caltech.edu/SIA?COLLECTION=spherex_dr1
- **IRSA SPHEREx tutorials and documentation** (FITS structure, flag bits, wavelength WCS): https://caltech-ipac.github.io/irsa-tutorials/spherex/
- **2MASS All-Sky Point Source Catalog** (NASA/IPAC), used to identify a blended star: https://irsa.ipac.caltech.edu/Missions/2mass.html (accessed through the VizieR mirror, catalog II/246)
- **NASA WISE / NEOWISE heritage**, the blink approach we extend: WiseView http://byw.tools/wiseview and Backyard Worlds: Planet 9 https://backyardworlds.org/

## Space Agency Partner Data and Other References
**Space Agency Partner data**
- **ESA Gaia DR3** astrometry, used indirectly through SIMBAD for target positions and proper motions: https://www.cosmos.esa.int/web/gaia/dr3

**Other data**
- **SIMBAD** (CDS, Strasbourg), coordinates, proper motions and parallaxes of WISE 0855 and the challenge targets, via the Sesame name resolver: https://simbad.cds.unistra.fr/ (mirror used: https://vizier.idia.ac.za/viz-bin/nph-sesame)
- **VizieR** (CDS) catalogs for the blended-star check: AAVSO VSX (B/vsx) https://www.aavso.org/vsx/ and Gaia DR3 variability classification (I/358) https://vizier.cds.unistra.fr/

**Literature** (all open access)
- Bock et al. 2025, *The SPHEREx Satellite Mission*: https://arxiv.org/abs/2511.02985
- Akeson et al. 2025, *The SPHEREx Image and Spectrophotometry Processing Pipeline*: https://arxiv.org/abs/2511.15823
- Saydjari et al. 2026, *Forward modeling linear variable filter data*: https://arxiv.org/abs/2611.12846
- Gagné et al. 2026, *A SPHEREx Pipeline and Spectral Library for Ultracool Dwarfs*: https://arxiv.org/abs/2604.22012
- Rustamkulov et al. 2026, *SPHEREx Spectra for a Sequence of Nearby Brown Dwarfs*: https://arxiv.org/abs/2607.00543
- Hui et al. 2026, *Spectral response of SPHEREx*: https://arxiv.org/abs/2602.09139
- Li et al. 2016, *Systematic Chromatic Errors*: https://arxiv.org/abs/1601.00117
- Holl et al. 2023, *Gaia DR3: scan-angle-dependent signals and spurious periods*: https://doi.org/10.1051/0004-6361/202245353
- Binnenfeld et al. 2026, *Mitigating scanning-induced systematics in Hipparcos*: https://arxiv.org/abs/2611.08282
- VanderPlas & Ivezić 2015, *Periodograms for Multiband Astronomical Time Series*: https://arxiv.org/abs/1502.01344
- Schlanger et al. 2026, *Co-addition and Subtraction of Undersampled Images*: https://arxiv.org/abs/2607.18054
- Meisner et al. 2022, *unTimely*: https://arxiv.org/abs/2209.14327
- Luhman 2014, discovery of WISE 0855−0714 (ApJL 786, L18): https://doi.org/10.1088/2041-8205/786/2/L18
- Literature search tools: arXiv API https://arxiv.org/help/api and OpenAlex https://openalex.org/

**Software and tools** (open source unless noted)
- Astropy https://www.astropy.org · NumPy https://numpy.org · SciPy https://scipy.org · pandas https://pandas.pydata.org · Matplotlib https://matplotlib.org · fsspec https://filesystem-spec.readthedocs.io · aiohttp https://docs.aiohttp.org · Pillow https://python-pillow.org · FFmpeg https://ffmpeg.org
- GitHub and GitHub Pages hosting: https://github.com/MurphyJUAN/spectral-aware-blink
- Claude / Claude Code (Anthropic), SI assistant, disclosed above: https://www.anthropic.com/claude-code
- All code in the repository was written for this project (SI-generated under the team's direction, disclosed in every file); no third-party code was copied beyond the libraries listed. No copyrighted images or music are used.
