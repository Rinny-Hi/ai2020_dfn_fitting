import zipfile, xml.etree.ElementTree as ET, re
from pathlib import Path

paths = [p for p in Path('C:/Users/user/OneDrive').rglob('261021*.pptx') if 24_000_000 < p.stat().st_size < 26_000_000]
paths += list(Path('outputs').glob('261021*.pptx'))
for p in paths:
    print('\nFILE', p)
    try:
      zctx = zipfile.ZipFile(p)
    except Exception as e:
      print('UNAVAILABLE', repr(e)); continue
    with zctx as z:
        pres = ET.fromstring(z.read('ppt/presentation.xml'))
        ns = {'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
              'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
        print('SIZE', pres.find('p:sldSz', ns).attrib)
        slides = sorted((n for n in z.namelist() if re.match(r'ppt/slides/slide\d+\.xml$', n)),
                        key=lambda x:int(re.search(r'(\d+)\.xml', x).group(1)))
        print('SLIDES', len(slides))
        for s in slides:
            root = ET.fromstring(z.read(s))
            tx = [n.text for n in root.findall('.//a:t', ns) if n.text]
            pix = len(root.findall('.//p:pic', ns)); tbl = len(root.findall('.//a:tbl', ns))
            print(re.search(r'(\d+)\.xml', s).group(1), 'pics', pix, 'tables', tbl, '|', ' / '.join(tx[:12])[:220])
