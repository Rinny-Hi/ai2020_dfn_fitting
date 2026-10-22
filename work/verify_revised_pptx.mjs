import fs from 'node:fs/promises';
import path from 'node:path';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';
const root=path.resolve('.');
const file=path.join(root,'outputs',(await fs.readdir(path.join(root,'outputs'))).find(x=>x.endsWith('최신검증본.pptx')));
const p=await PresentationFile.importPptx(await FileBlob.load(file));
const dir=path.join(root,'work','pptx_revised_all');await fs.mkdir(dir,{recursive:true});
for(let i=0;i<p.slides.items.length;i++){
 const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});
 await fs.writeFile(path.join(dir,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await b.arrayBuffer()));
}
console.log(JSON.stringify({slides:p.slides.items.length,file}));
