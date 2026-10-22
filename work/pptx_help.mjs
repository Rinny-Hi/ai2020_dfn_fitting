import {FileBlob, PresentationFile} from '@oai/artifact-tool';
const presentation = await PresentationFile.importPptx(await FileBlob.load('C:/Users/user/OneDrive/02. 현대자동차 미팅자료/★★현대자동차 최종 보고서 자료★★/261021 미팅자료.pptx'));
for (const q of ['slide.shapes.delete|shape.delete|slide.shapes.add','slide.charts.add|slide.tables.add','speaker notes']) {
  const r = presentation.help('*',{search:q,include:['index','examples','notes'],maxChars:4000});
  console.log('QUERY',q,'\n',r.ndjson);
}
