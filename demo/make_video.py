# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Render the 29-second demo video frames with PIL, then encode with ffmpeg.
Usage: python make_video.py <ffmpeg> <out.mp4>"""
import os
import subprocess
import sys
import tempfile

from PIL import Image, ImageDraw, ImageFont, ImageFilter

FF, OUT = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "results", "figs")
W, H, FPS = 1280, 720, 25
FD = "/mnt/c/Windows/Fonts/"


def font(name, size):
    for n in [name, "segoeui.ttf", "DejaVuSans.ttf"]:
        for d in [FD, "/usr/share/fonts/truetype/dejavu/"]:
            if os.path.exists(d + n):
                return ImageFont.truetype(d + n, size)
    return ImageFont.load_default()


F_B = lambda s: font("segoeuib.ttf", s)
F_R = lambda s: font("segoeui.ttf", s)
F_M = lambda s: font("consolab.ttf", s)
INK, INK2, INK3, BG = (238, 241, 247), (169, 177, 199), (111, 120, 149), (10, 15, 30)
STOPS = [(91, 76, 240), (42, 120, 214), (31, 182, 214), (27, 175, 122), (230, 194, 41), (235, 104, 52), (201, 50, 63)]


def lvf(width, height):
    im = Image.new("RGB", (width, height))
    px = im.load()
    for x in range(width):
        t = x / (width - 1) * (len(STOPS) - 1)
        i = min(int(t), len(STOPS) - 2)
        f = t - i
        c = tuple(int(STOPS[i][k] * (1 - f) + STOPS[i + 1][k] * f) for k in range(3))
        for y in range(height):
            px[x, y] = c
    return im


LVF_TOP = lvf(W, 6)


def base():
    im = Image.new("RGB", (W, H), BG)
    glow = Image.new("RGB", (W, H), BG)
    g = ImageDraw.Draw(glow)
    g.ellipse((700, -400, 1700, 300), fill=(24, 34, 74))
    glow = glow.filter(ImageFilter.GaussianBlur(160))
    im = Image.blend(im, glow, 0.9)
    im.paste(LVF_TOP, (0, 0))
    d = ImageDraw.Draw(im)
    d.text((W - 40, 18), "SI-assisted video (Claude Code) · idea by the team", font=F_M(13), fill=INK3, anchor="ra")
    d.text((40, H - 30), "Jellybit · NASA Space Apps 2026 · Planet X and SPHEREx", font=F_R(15), fill=INK3)
    return im


def wrap(d, text, f, maxw):
    words, lines, cur = text.split(), [], ""
    for w in words:
        t = (cur + " " + w).strip()
        if d.textlength(t, font=f) <= maxw:
            cur = t
        else:
            lines.append(cur); cur = w
    lines.append(cur)
    return lines


def text_block(d, xy, text, f, fill, maxw, gap=1.25):
    x, y = xy
    for ln in wrap(d, text, f, maxw):
        d.text((x, y), ln, font=f, fill=fill)
        y += int(f.size * gap)
    return y


def card(im, path, box, pad=10):
    """paste image fit into box=(x0,y0,x1,y1) on a light rounded card"""
    x0, y0, x1, y1 = box
    src = Image.open(path).convert("RGB")
    bw, bh = x1 - x0 - 2 * pad, y1 - y0 - 2 * pad
    s = min(bw / src.width, bh / src.height)
    src = src.resize((int(src.width * s), int(src.height * s)), Image.LANCZOS)
    cw, ch = src.width + 2 * pad, src.height + 2 * pad
    cx, cy = x0 + (x1 - x0 - cw) // 2, y0 + (y1 - y0 - ch) // 2
    d = ImageDraw.Draw(im)
    d.rounded_rectangle((cx, cy, cx + cw, cy + ch), radius=14, fill=(252, 252, 251))
    im.paste(src, (cx + pad, cy + pad))


def kicker(d, t, y=60):
    d.text((64, y), t.upper(), font=F_M(16), fill=INK2)


# ---------------- scenes ----------------
def s_intro():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "NASA Space Apps Challenge 2026 · Taichung, Taiwan", 230)
    d.text((64, 270), "Hi, I'm Yining Juan.", font=F_B(66), fill=INK)
    d.text((64, 365), "Team Jellybit", font=F_R(30), fill=INK2)
    return im


def s_challenge():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "Our challenge: Planet X and SPHEREx", 150)
    y = text_block(d, (64, 195), "Hunters of Planet X and nearby brown dwarfs blink sky images to see what moves.", F_B(42), INK, 1120)
    y = text_block(d, (64, y + 18), "But NASA's SPHEREx sees every star through a rainbow filter.", F_R(36), INK2, 1120)
    im.paste(lvf(1152, 18), (64, y + 30))
    d.text((64, y + 58), "0.75 µm", font=F_M(15), fill=INK3)
    d.text((1216, y + 58), "5.0 µm", font=F_M(15), fill=INK3, anchor="ra")
    return im


def s_title():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "Our project", 200)
    d.text((64, 240), "Spectral-Aware Blink", font=F_B(84), fill=INK)
    d.text((64, 360), "Remove the rainbow. Keep what really moves.", font=F_R(40), fill=(255, 210, 138))
    return im


def s_problem():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "The problem", 44)
    text_block(d, (64, 78), "One constant star looks 2× variable over time, yet is a smooth spectrum in wavelength.", F_B(34), INK, 1150)
    card(im, os.path.join(FIG, "fig1_problem.png"), (40, 180, 1240, 620))
    d.text((64, 640), "Real SPHEREx data · each exposure sees the star in a different colour", font=F_R(20), fill=INK2)
    return im


def s_blink(k):
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "How it works · NASA SPHEREx DR1 images", 44)
    text_block(d, (64, 78), "We learn each star's spectrum from the other passes, subtract it, and blink only the rest.", F_B(32), INK, 1150)
    card(im, f"/tmp/blink_{k}.png", (40, 170, 1240, 610))
    d.text((64, 628), f"Blinking pass {k + 1} of 4 · circle: WISE 0855, the coldest known brown dwarf, moving 1.3 px/yr", font=F_R(20), fill=INK2)
    d.text((64, 656), "Left: naive stack · Middle: colour-blind difference · Right: our spectral residual", font=F_R(18), fill=INK3)
    return im


def s_tool():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "Free web tool", 44)
    text_block(d, (64, 78), "Click any star: its spectrum, its light curve, and how much it would fake-blink.", F_B(32), INK, 1150)
    card(im, os.path.join(FIG, "web_tool.png"), (40, 165, 1240, 640))
    d.text((64, 656), "murphyjuan.github.io/spectral-aware-blink/web", font=F_M(18), fill=INK2)
    return im


def s_stats():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "What we found · 258 exposures · 10 fields · 5,100 light curves", 90)
    rows = [("58% → 1%", "constant stars wrongly flagged as periodic"),
            ("7–13% → 3%", "fake variability in standard detector stacks"),
            ("26–36% → 3%", "fake blinks predicted for cold brown dwarfs")]
    y = 160
    for big, lab in rows:
        im.paste(lvf(6, 112), (64, y))
        d.text((90, y - 4), big, font=F_B(64), fill=INK)
        d.text((92, y + 74), lab, font=F_R(24), fill=INK2)
        y += 150
    return im


def s_future():
    im = base(); d = ImageDraw.Draw(im)
    kicker(d, "What's next", 170)
    y = text_block(d, (64, 210), "An all-sky, spectrum-aware blink for citizen scientists and SPHEREx time-domain science, to find the coldest movers near the Sun.", F_B(40), INK, 1120)
    d.text((64, y + 30), "Try it: murphyjuan.github.io/spectral-aware-blink", font=F_M(24), fill=(255, 210, 138))
    d.text((64, y + 70), "Code & findings: github.com/MurphyJUAN/spectral-aware-blink", font=F_M(20), fill=INK2)
    return im


SCENES = [(s_intro(), 2.5), (s_challenge(), 3.0), (s_title(), 3.0), (s_problem(), 4.5)] + \
         [(s_blink(k), 1.25) for k in range(4)] + [(s_tool(), 4.0), (s_stats(), 4.0), (s_future(), 3.0)]
total = sum(t for _, t in SCENES)
print("total seconds", total)
assert total <= 29.5
tmp = tempfile.mkdtemp()
n = 0
XF = 6  # crossfade frames between non-blink scenes
for i, (im, dur) in enumerate(SCENES):
    nf = int(round(dur * FPS))
    nxt = SCENES[i + 1][0] if i + 1 < len(SCENES) else None
    blink_pair = 4 <= i <= 6  # hard cuts between blink frames (that is the blink)
    for j in range(nf):
        fr = im
        if nxt is not None and not blink_pair and j >= nf - XF:
            fr = Image.blend(im, nxt, (j - (nf - XF) + 1) / (XF + 1))
        fr.save(os.path.join(tmp, f"f{n:05d}.png"))
        n += 1
print("frames", n, "=", n / FPS, "s")
subprocess.run([FF, "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", os.path.join(tmp, "f%05d.png"),
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart",
                "-metadata", "comment=SI-assisted video generated by code written with Claude Code (Anthropic); idea by the team (Jellybit)",
                "-metadata", "title=Spectral-Aware Blink - Jellybit - NASA Space Apps 2026", OUT], check=True)
SCENES[1][0].save(os.path.join(HERE, "video_poster.png"))
print("wrote", OUT)
