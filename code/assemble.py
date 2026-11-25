"""Typeset the analytical R graphics in the authors' publication layout."""
import csv
import os
import sys
from itertools import groupby
from pathlib import Path

from lxml import etree as ET
import pypdfium2 as pdfium
from pypdf import PdfWriter
import vl_convert

HERE = Path(__file__).resolve().parent
SVG = "{http://www.w3.org/2000/svg}"


def layout(name):
    with (HERE / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_panel(page, source, frame, labels, family):
    graphic = ET.parse(str(source)).getroot()
    # Clip identifiers are only unique within a single svglite file.
    prefix = frame["panel"] + "_"
    for item in graphic.iter():
        if item.get("id"):
            item.set("id", prefix + item.get("id"))
        if item.get("clip-path"):
            item.set("clip-path", item.get("clip-path").replace("url(#", "url(#" + prefix))
        if "stroke-width:" in item.get("style", ""):
            style = dict(field.strip().split(": ", 1) for field in item.get("style").split(";") if field.strip())
            width = float(style["stroke-width"])
            minimum = float(frame["min_stroke_pt"])
            if style.get("stroke") != "none" and 0 < width < minimum:
                style["stroke-width"] = str(minimum)
                if "stroke-dasharray" in style:
                    style["stroke-dasharray"] = ",".join(
                        str(float(length) * minimum / width) for length in style["stroke-dasharray"].split(",")
                    )
                item.set("style", "; ".join(f"{k}: {v}" for k, v in style.items()))
    graphic.set("preserveAspectRatio", "none")
    for name in ("x", "y", "width", "height"):
        graphic.set(name, frame[name + "_pt"])
    page.append(graphic)
    raw = list(graphic.iter(SVG + "text"))
    panel_tag = frame["figure"] == "4"
    expected = labels[:-1] if panel_tag else labels
    if [e.text for e in raw] == [row["label"] for row in expected]:
        for element in raw:
            element.getparent().remove(element)
    else:
        # New datasets may generate different tick breaks. Preserve the R labels.
        print(f'{source.stem}: using R label placement for changed scales.')
        labels = labels[-1:] if panel_tag else []
    for row in labels:
        text = ET.SubElement(page, SVG + "text", {
            "transform": f'translate({row["x_pt"]} {row["baseline_y_pt"]}) rotate({row["angle_deg"]})',
            "font-family": family, "font-weight": "bold", "font-size": row["size_pt"],
            "fill": row["color"],
        })
        text.text = row["label"]


def assemble(root, formats):
    folder = root / "fig"
    output = folder / "assembled"
    output.mkdir()
    family = os.environ.get("SD_FONT_FAMILY", "Arial")
    font_path = Path(os.environ["SD_BOLD_FONT_FILE"]).read_text(encoding="utf-8-sig").strip()
    vl_convert.register_font_directory(str(Path(font_path).parent))
    labels = layout("label_layout.tsv")
    combined = PdfWriter()
    for number, rows in groupby(layout("panel_layout.tsv"), key=lambda r: r["figure"]):
        frames = list(rows)
        width, height = (frames[0][key] for key in ("page_width_pt", "page_height_pt"))
        page = ET.Element(SVG + "svg", nsmap={None: SVG[1:-1]}, width=width, height=height,
                          viewBox=f"0 0 {width} {height}")
        ET.SubElement(page, SVG + "rect", width=width, height=height, fill="white")
        for frame in frames:
            positions = [row for row in labels if row["panel"] == frame["panel"]]
            add_panel(page, folder / frame["svg"], frame, positions, family)
        dest = output / f"Fig{number}_assembled"
        # svg_to_pdf maps one user unit to one PDF point.
        dest.with_suffix(".pdf").write_bytes(vl_convert.svg_to_pdf(ET.tostring(page, encoding="unicode")))
        if "svg" in formats:
            page.set("width", width + "pt")
            page.set("height", height + "pt")
            ET.ElementTree(page).write(str(dest.with_suffix(".svg")), encoding="utf-8", xml_declaration=True)
        with pdfium.PdfDocument(dest.with_suffix(".pdf")) as document:
            image = document[0].render(scale=600/72).to_pil().convert("RGB")
            image.save(dest.with_suffix(".png"), dpi=(600, 600))
            if "tiff" in formats:
                image.save(dest.with_suffix(".tiff"), dpi=(600, 600), compression="tiff_lzw")
        combined.append(dest.with_suffix(".pdf"))
    combined.write(output / "Figures2-4.pdf")
    combined.close()
    if "svg" not in formats:
        for path in folder.glob("Fig*/*.svg"):
            path.unlink()


if __name__ == "__main__":
    assemble(Path(sys.argv[1]), sys.argv[2].split(","))
