# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Narrated, animated <30 s demo. Voice: SI text-to-speech (edge-tts, en-US-AriaNeural). Music: synthesized here (original).
Usage: python make_video2.py <ffmpeg> <voice_dir with s1..s5.mp3> <out.mp4>"""
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter
from scipy.io import wavfile

FF, VO, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "results", "figs")
W, H, FPS = 1280, 720, 25
FD = "/mnt/c/Windows/Fonts/"
fnt = lambda n, s: ImageFont.truetype(FD + n, s)
B = lambda s: fnt("segoeuib.ttf", s)
R = lambda s: fnt("segoeui.ttf", s)
M = lambda s: fnt("consolab.ttf", s)
INK, INK2, INK3, BG, GOLD = (238, 241, 247), (169, 177, 199), (111, 120, 149), (10, 15, 30), (255, 210, 138)
STOPS = [(91, 76, 240), (42, 120, 214), (31, 182, 214), (27, 175, 122), (230, 194, 41), (235, 104, 52), (201, 50, 63)]

# ---- narration timeline ----
SEG = ["Hi, I'm Yining Juan from team Jellybit.",
       "SPHEREx sees every star through a rainbow filter, so even steady stars flicker when you blink its images.",
       "Spectral-Aware Blink learns each star's spectrum, subtracts it, and blinks only what's left.",
       "Normal blinking fakes cold brown dwarfs, Planet X's best targets, by up to thirty-six percent. We cut it to three.",
       "Tested on ten fields, it catches WISE 0855 moving."]
dur = [float(subprocess.run([FF.replace("ffmpeg", "ffprobe"), "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0",
                             os.path.join(VO, f"s{i+1}.mp3")], capture_output=True, text=True).stdout) for i in range(5)]
starts, t = [], 0.15
for d in dur:
    starts.append(t); t += d + 0.1
TOTAL = min(29.7, t + 0.1)
print("segments", [round(s, 2) for s in starts], "total", TOTAL)


def lvf_color(x):
    t = min(max(x, 0), 1) * (len(STOPS) - 1)
    i = min(int(t), len(STOPS) - 2); f = t - i
    return tuple(int(STOPS[i][k] * (1 - f) + STOPS[i + 1][k] * f) for k in range(3))


def lvf(width, height):
    a = np.array([lvf_color(x / (width - 1)) for x in range(width)], np.uint8)
    return Image.fromarray(np.repeat(a[None], height, 0))


def make_base():
    im = Image.new("RGB", (W, H), BG)
    glow = Image.new("RGB", (W, H), BG); g = ImageDraw.Draw(glow)
    g.ellipse((700, -400, 1700, 300), fill=(24, 34, 74)); g.ellipse((-500, 450, 400, 1100), fill=(18, 26, 56))
    im = Image.blend(im, glow.filter(ImageFilter.GaussianBlur(160)), 0.9)
    im.paste(lvf(W, 6), (0, 0))
    d = ImageDraw.Draw(im)
    d.text((W - 36, 16), "SI-assisted video: SI voice (TTS) + code-made visuals & music · idea by the team", font=M(12), fill=INK3, anchor="ra")
    return im


BASE = make_base()
IMG = {k: Image.open(v).convert("RGB") for k, v in {"fig1": os.path.join(FIG, "fig1_problem.png"), "tool": os.path.join(FIG, "web_tool.png"),
                                                      **{f"b{k}": f"/tmp/blink_{k}.png" for k in range(4)}}.items()}


def ease(x):
    x = min(max(x, 0.0), 1.0); return x * x * (3 - 2 * x)


def ftext(im, xy, s, f, fill, t, t0, dy=18, anchor="la"):
    """fade+slide-in text starting at t0"""
    a = ease((t - t0) / 0.45)
    if a <= 0: return
    col = tuple(int(BG[k] + (fill[k] - BG[k]) * a) for k in range(3))
    ImageDraw.Draw(im).text((xy[0], xy[1] + dy * (1 - a)), s, font=f, fill=col, anchor=anchor)


def card(im, src, box, zoom=1.0, pad=10, alpha=1.0):
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0 - 2 * pad, y1 - y0 - 2 * pad
    s = min(bw / src.width, bh / src.height)
    sw, sh = int(src.width * s), int(src.height * s)
    img = src.resize((sw, sh), Image.BILINEAR)
    if zoom > 1.0:  # Ken Burns: crop centre then resize back
        cw, ch = int(sw / zoom), int(sh / zoom)
        img = img.crop(((sw - cw) // 2, (sh - ch) // 2, (sw - cw) // 2 + cw, (sh - ch) // 2 + ch)).resize((sw, sh), Image.BILINEAR)
    cx, cy = x0 + (x1 - x0 - sw - 2 * pad) // 2, y0 + (y1 - y0 - sh - 2 * pad) // 2
    layer = im.copy(); d = ImageDraw.Draw(layer)
    d.rounded_rectangle((cx, cy, cx + sw + 2 * pad, cy + sh + 2 * pad), radius=14, fill=(252, 252, 251))
    layer.paste(img, (cx + pad, cy + pad))
    return Image.blend(im, layer, alpha) if alpha < 1 else layer


def subtitle(im, t):
    for i, s0 in enumerate(starts):
        if s0 - 0.05 <= t <= s0 + dur[i] + 0.15:
            d = ImageDraw.Draw(im); f = R(22)
            words, lines, cur = SEG[i].split(), [], ""
            for w in words:
                tt = (cur + " " + w).strip()
                if d.textlength(tt, font=f) < 1080: cur = tt
                else: lines.append(cur); cur = w
            lines.append(cur)
            y = H - 22 - 30 * len(lines)
            d.rounded_rectangle((90, y - 10, W - 90, H - 16), radius=10, fill=(5, 8, 18))
            for ln in lines:
                d.text((W // 2, y), ln, font=f, fill=INK, anchor="ma"); y += 30
            return


def frame(t):
    im = BASE.copy(); d = ImageDraw.Draw(im)
    s1, s2, s3, s4, s5 = starts
    if t < s2 - 0.1:  # intro
        ftext(im, (64, 210), "NASA SPACE APPS CHALLENGE 2026 · TAICHUNG, TAIWAN", M(17), INK2, t, 0.2)
        ftext(im, (64, 250), "Hi, I'm Yining Juan.", B(70), INK, t, 0.5)
        ftext(im, (64, 345), "Team Jellybit", R(34), INK2, t, 1.1)
        ftext(im, (64, 405), "Challenge: Planet X and SPHEREx", B(30), GOLD, t, 2.6)
    elif t < s2 + 3.6:  # rainbow filter: a constant star slides across the LVF and "flickers"
        tt = t - s2
        ftext(im, (64, 120), "Planet hunters blink images to spot movers…", B(40), INK, t, s2)
        ftext(im, (64, 180), "…but SPHEREx sees every star through a rainbow filter.", R(34), INK2, t, s2 + 1.4)
        bx0, bx1, by = 64, 1216, 360
        im.paste(lvf(bx1 - bx0, 26), (bx0, by))
        d.text((bx0, by + 36), "0.75 µm", font=M(15), fill=INK3); d.text((bx1, by + 36), "5.0 µm", font=M(15), fill=INK3, anchor="ra")
        x = (tt / 3.6) % 1.0
        sed = 0.55 + 0.45 * math.sin(6.5 * x + 0.8) * math.exp(-1.2 * x) + 0.25 * math.exp(-((x - 0.62) / 0.06) ** 2)
        r = 10 + 34 * max(sed, 0.08)
        cx, cy = bx0 + x * (bx1 - bx0), 500
        glow = Image.new("RGB", (W, H), (0, 0, 0)); gd = ImageDraw.Draw(glow)
        c = lvf_color(x)
        gd.ellipse((cx - 2.2 * r, cy - 2.2 * r, cx + 2.2 * r, cy + 2.2 * r), fill=tuple(int(v * 0.5) for v in c))
        im = Image.fromarray(np.clip(np.asarray(im, np.int16) + np.asarray(glow.filter(ImageFilter.GaussianBlur(18)), np.int16), 0, 255).astype(np.uint8))
        d = ImageDraw.Draw(im)
        d.ellipse((cx - r, cy - r, cx + r, cy + r), fill=(255, 255, 255))
        d.line((cx, by + 28, cx, cy - r - 6), fill=c, width=2)
        d.text((64, 590), f"same constant star · brightness seen through this colour: {sed:4.2f}", font=M(18), fill=INK2)
    elif t < s3 - 0.1:  # real data: problem figure with slow zoom
        tt = t - (s2 + 3.6)
        ftext(im, (64, 60), "Real SPHEREx data: one steady star looks 2× variable", B(36), INK, t, s2 + 3.6)
        im = card(im, IMG["fig1"], (40, 130, 1240, 580), zoom=1.0 + 0.06 * tt, alpha=ease(tt / 0.4))
    elif t < s3 + 2.0:  # title
        ftext(im, (64, 210), "OUR PROJECT", M(17), INK2, t, s3)
        ftext(im, (64, 245), "Spectral-Aware Blink", B(86), INK, t, s3 + 0.15)
        ftext(im, (64, 370), "Remove the rainbow. Keep what really moves.", R(40), GOLD, t, s3 + 0.7)
    elif t < s4 - 0.1:  # learn spectrum, subtract, blink
        tt = t - (s3 + 2.0)
        k = min(int(tt / 0.9), 3)
        ftext(im, (64, 56), "Learn the spectrum · subtract · blink the rest", B(34), INK, t, s3 + 2.0)
        im = card(im, IMG[f"b{k}"], (40, 120, 1240, 560))
        d = ImageDraw.Draw(im)
        d.text((64, 572), f"BLINK  pass {k+1}/4", font=M(20), fill=GOLD)
    elif t < s5 - 0.05:  # KEY INSIGHT: brown dwarfs faked most
        tt = t - s4
        ftext(im, (64, 80), "THE KEY INSIGHT · 10 FIELDS · 5,100 LIGHT CURVES", M(16), INK2, t, s4)
        ftext(im, (64, 115), "Cold brown dwarfs, Planet X's best targets,", B(40), INK, t, s4 + 0.1)
        ftext(im, (64, 168), "are the ones normal blinking fakes the most.", B(40), INK, t, s4 + 0.5)
        v = 36 - 33 * ease((tt - 3.6) / 1.4)
        im.paste(lvf(6, 150), (64, 260))
        ftext(im, (90, 250), "up to", R(30), INK2, t, s4 + 1.5, dy=0)
        ftext(im, (90, 285), f"{round(v)}%", B(130), (255, 210, 138) if tt < 3.6 else INK, t, s4 + 1.5, dy=0)
        ftext(im, (92, 440), "fake blink amplitude: normal blinking  →  Spectral-Aware Blink", R(26), INK2, t, s4 + 2.0)
        ftext(im, (64, 500), "Also: 58% of constant stars look 'periodic' naively → 1% with ours", R(22), INK3, t, s4 + 4.2)
    else:  # validation + WISE 0855 + call to action
        tt = t - s5
        k = [0, 3, 0, 3, 0, 3][min(int(tt / 0.75), 5)]
        ftext(im, (64, 50), "Tested on 10 sky fields: WISE 0855 caught moving", B(34), INK, t, s5)
        im = card(im, IMG[f"b{k}"], (40, 105, 1240, 520), zoom=1.0 + 0.05 * tt)
        d = ImageDraw.Draw(im)
        d.text((64, 530), f"pass {k+1}: the dipole flips as the coldest known brown dwarf moves", font=M(19), fill=GOLD)
        ftext(im, (64, 572), "Try it online: murphyjuan.github.io/spectral-aware-blink", M(22), INK, t, s5 + 1.0, dy=0)
    subtitle(im, t)
    return im


# ---- music: soft pad Am–F–C–G + gentle arpeggio, synthesized (original) ----
SR = 44100
n = int(TOTAL * SR); tt = np.arange(n) / SR
CH = [[220.0, 261.63, 329.63], [174.61, 220.0, 261.63], [130.81, 164.81, 196.0], [196.0, 246.94, 293.66]]
mus = np.zeros(n)
seg = TOTAL / 4
for i, ch in enumerate(CH):
    i0, i1 = int(i * seg * SR), int(min((i + 1) * seg * SR + 0.6 * SR, n))
    tl = np.arange(i1 - i0) / SR
    env = np.minimum(1, tl / 1.2) * np.minimum(1, (tl[-1] - tl + 1e-3) / 1.0)
    for f in ch:
        mus[i0:i1] += env * (np.sin(2 * np.pi * f * tl) + 0.3 * np.sin(2 * np.pi * 2 * f * tl) + 0.15 * np.sin(2 * np.pi * f / 2 * tl))
    for j in range(8):  # arpeggio
        k0 = i0 + int(j * seg / 8 * SR)
        if k0 >= n: break
        L = min(int(0.6 * SR), n - k0); tn = np.arange(L) / SR
        f = ch[j % 3] * 2
        mus[k0:k0 + L] += 0.35 * np.exp(-tn * 5) * np.sin(2 * np.pi * f * tn)
mus *= np.minimum(1, tt / 1.0) * np.minimum(1, (TOTAL - tt) / 1.5)
mus = mus / np.max(np.abs(mus)) * 0.9
tmp = tempfile.mkdtemp()
wavfile.write(os.path.join(tmp, "music.wav"), SR, (mus * 32767).astype(np.int16))

# ---- frames ----
nfr = int(TOTAL * FPS)
for i in range(nfr):
    frame(i / FPS).save(os.path.join(tmp, f"f{i:05d}.png"))
print("frames", nfr)
inputs = ["-framerate", str(FPS), "-i", os.path.join(tmp, "f%05d.png"), "-i", os.path.join(tmp, "music.wav")]
fc = []
for i in range(5):
    inputs += ["-i", os.path.join(VO, f"s{i+1}.mp3")]
    ms = int(starts[i] * 1000)
    fc.append(f"[{i+2}:a]adelay={ms}|{ms},volume=1.6[v{i}]")
fc.append("[1:a]volume=0.16[m]")
fc.append("".join(f"[v{i}]" for i in range(5)) + "[m]amix=inputs=6:normalize=0:duration=longest[a]")
subprocess.run([FF, "-y", "-loglevel", "error"] + inputs + ["-filter_complex", ";".join(fc), "-map", "0:v", "-map", "[a]",
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-c:a", "aac", "-b:a", "160k", "-t", f"{TOTAL:.2f}", "-movflags", "+faststart",
                "-metadata", "title=Spectral-Aware Blink - Jellybit - NASA Space Apps 2026",
                "-metadata", "comment=SI-assisted: narration is SI text-to-speech (Microsoft Edge neural voice via edge-tts); visuals and music generated by code written with Claude Code (Anthropic); idea by the team",
                OUT], check=True)
print("wrote", OUT, TOTAL)
