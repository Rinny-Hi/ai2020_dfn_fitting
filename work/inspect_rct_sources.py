from pathlib import Path
import re
import zipfile
import xml.etree.ElementTree as ET

import openpyxl


downloads = Path(r"C:\Users\user\Downloads")
book = downloads / "Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx"
wb = openpyxl.load_workbook(book, read_only=True, data_only=True)
for ws in wb.worksheets:
    print("\nSHEET", ws.title)
    rows = ws.iter_rows(values_only=True)
    for index, row in enumerate(rows):
        if ws.title == "Info" or index < 3:
            print(index + 1, row)
        if ws.title != "Info" and index >= 3:
            break

needle = ("EIS-DRT", "파라미터 선정 근거", "N3/N4", "N3 + N4")
ns = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main"}
for path in downloads.glob("*.pptx"):
    try:
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                if not re.fullmatch(r"ppt/slides/slide\d+\.xml", name):
                    continue
                root = ET.fromstring(z.read(name))
                text = " ".join((t.text or "") for t in root.findall(".//a:t", ns))
                if any(part in text for part in needle):
                    print("SLIDE", path.name, name, text[:500])
    except Exception as exc:
        print("ERROR", path.name, type(exc).__name__)
