"""
Build the 8-10 page mini-project report as a Word document.

Everything numerical is read back out of the result files, so the report
cannot drift from the data: rerun the pipeline, rerun this, and the numbers
move with it.

Output: 07_report/GSK3B_virtual_screen_report.docx
"""
import csv, json, os, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
FIG = os.path.join(ROOT, "06_analysis", "figures")
OUT = os.path.join(ROOT, "07_report")

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ACCENT = RGBColor(0x8C, 0x2B, 0x22)
GREY = RGBColor(0x55, 0x55, 0x55)


# ----------------------------------------------------------------- data ----
def load_json(*p):
    return json.load(open(os.path.join(ROOT, *p), encoding="utf8"))


def load_csv(*p):
    return list(csv.DictReader(open(os.path.join(ROOT, *p), encoding="utf8")))



_RESIZED = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                        "06_analysis", "figures", "_embed")


def embed_copy(path, width_cm, dpi=330):
    """Downsample a figure to the resolution it is actually printed at.

    The 3Dmol captures are ~3400 px wide; at 11 cm that is over 780 dpi and
    makes the .docx ten times bigger than it needs to be.
    """
    try:
        from PIL import Image
    except ImportError:
        return path
    target = int(width_cm / 2.54 * dpi)
    try:
        im = Image.open(path)
    except Exception:
        return path
    if im.width <= target * 1.1:
        return path
    os.makedirs(_RESIZED, exist_ok=True)
    out = os.path.join(_RESIZED, os.path.basename(path))
    if not (os.path.exists(out) and os.path.getmtime(out) > os.path.getmtime(path)):
        im = im.convert("RGB")
        im.thumbnail((target, int(target * im.height / im.width)), Image.LANCZOS)
        im.save(out, "PNG", optimize=True)
    return out


# -------------------------------------------------------------- document ---
class Report:
    def __init__(self):
        d = Document()
        for name, size in (("Normal", 9.0), ):
            st = d.styles[name]
            st.font.name = "Calibri"; st.font.size = Pt(size)
        st = d.styles["Normal"].paragraph_format
        st.space_after = Pt(4); st.line_spacing = 1.0
        for s in d.sections:
            s.top_margin = Cm(1.5); s.bottom_margin = Cm(1.5)
            s.left_margin = Cm(1.8); s.right_margin = Cm(1.8)
        self.d = d
        self.fig_n = 0
        self.tab_n = 0

    # --- text -------------------------------------------------------------
    def h1(self, text):
        p = self.d.add_paragraph()
        p.paragraph_format.space_before = Pt(12); p.paragraph_format.space_after = Pt(3)
        r = p.add_run(text); r.bold = True; r.font.size = Pt(13); r.font.color.rgb = ACCENT
        return p

    def h2(self, text):
        p = self.d.add_paragraph()
        p.paragraph_format.space_before = Pt(8); p.paragraph_format.space_after = Pt(2)
        r = p.add_run(text); r.bold = True; r.font.size = Pt(10.5)
        return p

    def p(self, text, size=9.0, italic=False, align=None, space_after=3):
        par = self.d.add_paragraph()
        par.paragraph_format.space_after = Pt(space_after)
        if align:
            par.alignment = align
        self._runs(par, text, size, italic)
        return par

    def bullet(self, text, size=9.0):
        par = self.d.add_paragraph(style="List Bullet")
        par.paragraph_format.space_after = Pt(2)
        par.paragraph_format.left_indent = Cm(0.6)
        self._runs(par, text, size, False)
        return par

    @staticmethod
    def _runs(par, text, size, italic):
        """**bold** and *italic* inline markup"""
        import re
        for part in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text):
            if not part:
                continue
            if part.startswith("**") and part.endswith("**"):
                r = par.add_run(part[2:-2]); r.bold = True
            elif part.startswith("*") and part.endswith("*") and len(part) > 2:
                r = par.add_run(part[1:-1]); r.italic = True
            else:
                r = par.add_run(part)
            r.font.size = Pt(size); r.italic = r.italic or italic

    # --- objects ----------------------------------------------------------
    def figure(self, filename, caption, width_cm=16.0):
        path = os.path.join(FIG, filename)
        if not os.path.exists(path):
            print("  !! missing figure", filename); return
        self.fig_n += 1
        p = self.d.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(6); p.paragraph_format.space_after = Pt(2)
        p.add_run().add_picture(embed_copy(path, width_cm), width=Cm(width_cm))
        c = self.d.add_paragraph(); c.alignment = WD_ALIGN_PARAGRAPH.CENTER
        c.paragraph_format.space_after = Pt(8)
        r = c.add_run(f"Figure {self.fig_n}. {caption}")
        r.font.size = Pt(8); r.font.color.rgb = GREY
        return self.fig_n

    def figure_pair(self, left, left_caption, right, right_caption, width_cm=7.6):
        """two figures side by side in a borderless table, to save a page"""
        paths = [os.path.join(FIG, left), os.path.join(FIG, right)]
        if not all(os.path.exists(p) for p in paths):
            for f, c in ((left, left_caption), (right, right_caption)):
                self.figure(f, c, width_cm=width_cm + 1.0)
            return
        t = self.d.add_table(rows=1, cols=2)
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for cell, path, cap in zip(t.rows[0].cells, paths, (left_caption, right_caption)):
            cell.width = Cm(width_cm + 0.6)
            p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_after = Pt(1)
            p.add_run().add_picture(embed_copy(path, width_cm), width=Cm(width_cm))
            self.fig_n += 1
            c = cell.add_paragraph(); c.alignment = WD_ALIGN_PARAGRAPH.CENTER
            r = c.add_run(f"Figure {self.fig_n}. {cap}")
            r.font.size = Pt(8); r.font.color.rgb = GREY
        self.d.add_paragraph().paragraph_format.space_after = Pt(2)
        return t

    def table(self, headers, rows, caption=None, widths=None, size=8):
        self.tab_n += 1
        if caption:
            c = self.d.add_paragraph()
            c.paragraph_format.space_before = Pt(6); c.paragraph_format.space_after = Pt(2)
            r = c.add_run(f"Table {self.tab_n}. {caption}")
            r.font.size = Pt(8.5); r.font.color.rgb = GREY; r.bold = True
        t = self.d.add_table(rows=1, cols=len(headers))
        t.style = "Light Grid Accent 1"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        hdr = t.rows[0].cells
        for i, h in enumerate(headers):
            hdr[i].text = ""
            par = hdr[i].paragraphs[0]
            run = par.add_run(str(h)); run.bold = True; run.font.size = Pt(size)
        for row in rows:
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = ""
                par = cells[i].paragraphs[0]
                run = par.add_run("" if v is None else str(v))
                run.font.size = Pt(size)
        if widths:
            for r_ in t.rows:
                for i, w in enumerate(widths):
                    r_.cells[i].width = Cm(w)
        self.d.add_paragraph().paragraph_format.space_after = Pt(2)
        return t

    def columns(self, n=2):
        """start a continuous multi-column section (used for the references)"""
        from docx.enum.section import WD_SECTION
        sec = self.d.add_section(WD_SECTION.CONTINUOUS)
        cols = sec._sectPr.xpath("./w:cols")[0]
        cols.set(qn("w:num"), str(n))
        cols.set(qn("w:space"), "360")
        return sec

    def pagebreak(self):
        self.d.add_page_break()

    def save(self, path):
        self.d.save(path)
        print("wrote", path)
