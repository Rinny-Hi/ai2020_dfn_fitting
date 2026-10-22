import fs from 'node:fs/promises';
import path from 'node:path';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';

const src = path.resolve('outputs', (await fs.readdir('outputs')).find(x=>x.startsWith('261021_') && x.endsWith('.pptx')));
const out = path.resolve('work/deck_inspect');
await fs.mkdir(out,{recursive:true});
console.log('IMPORT',src);
const p = await PresentationFile.importPptx(await FileBlob.load(src));
console.log('SLIDES',p.slides.items.length);
const v = await p.inspect({select:'$.slides[*]',include:['id','index','title','elements.id','elements.kind','elements.text','elements.position','elements.style'],fileName:path.join(out,'slides.ndjson')});
await fs.writeFile(path.join(out,'slides.ndjson'),v.ndjson);
for (let i=0;i<p.slides.items.length;i++) {
  if (i<5 || i>17) continue;
  const slide=p.slides.getItem(i);
  const png=await p.export({slide,format:'png',scale:1});
  await fs.writeFile(path.join(out,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
}
console.log('DONE');
