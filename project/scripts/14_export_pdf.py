"""
Export the Word report to PDF with Word itself (keeps the layout exactly as
Word renders it) and report the page count.

Run: python 14_export_pdf.py
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DOCX = os.path.join(ROOT, "07_report", "CADD_miniproject_GSK3B_report.docx")
PDF = os.path.join(ROOT, "07_report", "GSK3B_virtual_screen_report.pdf")

WD_FORMAT_PDF = 17
WD_STAT_PAGES = 2


def main():
    try:
        import win32com.client as win32
    except ImportError:
        print("pywin32 not installed: pip install pywin32"); sys.exit(1)
    word = win32.DispatchEx("Word.Application")
    word.Visible = False
    try:
        doc = word.Documents.Open(DOCX, ReadOnly=False)
        doc.Repaginate()
        pages = doc.ComputeStatistics(WD_STAT_PAGES)
        words = doc.ComputeStatistics(0)
        doc.SaveAs(PDF, FileFormat=WD_FORMAT_PDF)
        doc.Close(False)
        print(f"{pages} pages, {words} words -> {PDF}")
    finally:
        word.Quit()


if __name__ == "__main__":
    main()
