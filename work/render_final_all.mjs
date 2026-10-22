import fs from 'node:fs/promises';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';
const source='C:/Users/user/Documents/Codex/2026-10-20/https-github-com-rinny-hi-ai2020-2/outputs/261021_미팅자료_DFN_검증반영_완성본.pptx';
const dir='C:/Users/user/Documents/Codex/2026-10-20/https-github-com-rinny-hi-ai2020-2/work/pptx_final_all';
await fs.mkdir(dir,{recursive:true});
const p=await PresentationFile.importPptx(await FileBlob.load(source));
for(let i=0;i<p.slides.items.length;i++){
  const blob=await p.export({slide:p.slides.items[i],format:'png',scale:1});
  await fs.writeFile(`${dir}/slide-${String(i+1).padStart(2,'0')}.png`,new Uint8Array(await blob.arrayBuffer()));
}
console.log('rendered',p.slides.items.length);
