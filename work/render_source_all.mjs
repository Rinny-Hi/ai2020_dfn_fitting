import fs from 'node:fs/promises';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';
const source='C:/Users/user/OneDrive/02. 현대자동차 미팅자료/★★현대자동차 최종 보고서 자료★★/261021 미팅자료.pptx';
const dir='C:/Users/user/Documents/Codex/2026-10-20/https-github-com-rinny-hi-ai2020-2/work/pptx_source_all';
await fs.mkdir(dir,{recursive:true});
const p=await PresentationFile.importPptx(await FileBlob.load(source));
for(let i=0;i<p.slides.items.length;i++){
  const blob=await p.export({slide:p.slides.items[i],format:'png',scale:1});
  await fs.writeFile(`${dir}/slide-${String(i+1).padStart(2,'0')}.png`,new Uint8Array(await blob.arrayBuffer()));
}
console.log('rendered',p.slides.items.length);
