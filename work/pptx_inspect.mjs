import fs from 'node:fs/promises';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';

const source = 'C:/Users/user/OneDrive/02. 현대자동차 미팅자료/★★현대자동차 최종 보고서 자료★★/261021 미팅자료.pptx';
const dir = 'C:/Users/user/Documents/Codex/2026-10-20/https-github-com-rinny-hi-ai2020-2/work/pptx_preview';
await fs.mkdir(dir, {recursive: true});
const presentation = await PresentationFile.importPptx(await FileBlob.load(source));
console.log('SLIDES', presentation.slides.items.length);
for (const number of [6,7,8,9,10,12,17,18]) {
  const slide = presentation.slides.items[number - 1];
  const blob = await presentation.export({slide, format: 'png', scale: 1});
  await fs.writeFile(`${dir}/slide-${String(number).padStart(2,'0')}.png`, new Uint8Array(await blob.arrayBuffer()));
  const layout = await slide.export({format: 'layout'});
  await fs.writeFile(`${dir}/slide-${String(number).padStart(2,'0')}.json`, await layout.text());
  console.log('RENDERED', number);
}
