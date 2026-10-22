import fs from 'node:fs/promises';
import fsSync from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';

const workspaceDir='C:/Users/user/Documents/Codex/2026-10-20/https-github-com-rinny-hi-ai2020-2';
const source='C:/Users/user/OneDrive/02. 현대자동차 미팅자료/★★현대자동차 최종 보고서 자료★★/261021 미팅자료.pptx';
const SKILL_DIR='C:/Users/user/.codex/plugins/cache/openai-primary-runtime/presentations/26.1015.12040/skills/presentations';
const RUNTIME_PYTHON='C:/Users/user/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe';
process.env.RUNTIME_NODE_MODULES='C:/Users/user/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const FINAL_PPTX=path.join(workspaceDir,'outputs','261021_미팅자료_DFN_검증반영_완성본.pptx');
const data=JSON.parse(await fs.readFile(path.join(workspaceDir,'work','ppt_data.json'),'utf8'));
const presentation=await PresentationFile.importPptx(await FileBlob.load(source));
const S=n=>presentation.slides.items[n-1];
const blue='#156082', blue2='#4e85bf', pale='#e8eef6', dark='#202a35', orange='#c44e32';

function erase(...ids){presentation.delete(ids)}
function text(slide, body, x,y,w,h,opts={}){
  const sh=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',
    line:{style:'solid',fill:'none',width:0}});
  sh.text=body;
  sh.text.style={fontFamily:'Pretendard',fontSize:opts.size??17,color:opts.color??dark,
    bold:opts.bold??false,alignment:opts.align??'left'};
  return sh;
}
function rule(slide,x,y,w,color=blue){slide.shapes.add({geometry:'line',position:{left:x,top:y,width:w,height:0},
  fill:'none',line:{style:'solid',fill:color,width:1}})}
function table(slide,vals,x,y,w,h,fs=12){
  const t=slide.tables.add({rows:vals.length,columns:vals[0].length,left:x,top:y,width:w,height:h,values:vals});
  t.styleOptions={headerRow:true,bandedRows:true};
  for(let i=0;i<vals.length;i++) for(let j=0;j<vals[0].length;j++){
    const c=t.getCell(i,j);
    c.fill=i===0?blue:(i%2?pale:'#f6f8fb');
    c.text.style={fontFamily:'Pretendard',fontSize:fs,color:i===0?'#ffffff':dark,bold:i===0,
      alignment:j===0?'left':'center'};
  }
  t.borders.assign({style:'solid',fill:'#ffffff',width:1});
  t.cells.block({row:0,column:0,rowCount:vals.length,columnCount:vals[0].length}).assign({
    margins:{top:1,bottom:1,left:3,right:3}
  });
  return t;
}
function img(slide,file,x,y,w,h,alt){slide.images.add({blob:fsSync.readFileSync(path.join(workspaceDir,'outputs',file)),
  contentType:'image/png',alt,fit:'contain',position:{left:x,top:y,width:w,height:h}})}
function note(slide,s){slide.speakerNotes.append('\n[2026-10-21 DFN 검증 추가] '+s)}
const metric=(ocp,rate,dir)=>data.prefit.find(r=>r.ocp===ocp&&r.rate===rate&&r.direction===dir);
const scenario=(ocp,el,phase,dir)=>data.scenario.find(r=>r.ocp===ocp&&r.electrolyte===el&&r.phase===phase&&r.direction===dir);
const fit=(ocp,el)=>data.fit.find(r=>r.ocp===ocp&&r.electrolyte===el);
const pot=(ocp,el,phase,rate)=>data.potential.find(r=>r.ocp===ocp&&r.electrolyte===el&&r.phase===phase&&r.rate===rate);
const fitted='charge_fit_discharge_validation', initial='experimental_initial';

// Slide 6: only the two blue prompts are replaced.
erase('sh/w32twb6l','sh/u1kbu1ov');
text(S(6),'• C/50 qOCV와 양극 GITT rest OCP로 window를 각각 재최적화; 음극은 Ai2020와 측정 저율 OCP 두 버전 유지.',52,135,930,34,{size:16});
text(S(6),'• 용량·전극 두께·활물질 분율에 맞춰 cₛ,max를 OCP별로 재계산. 서로 다른 window의 cₛ,max 혼용 금지.',52,172,930,34,{size:16});
text(S(6),'• 아래 수치는 동일 셀의 독립 OCP 조합. 후속 Ds·k·전해질 비교도 각각 이어서 수행.',52,209,930,34,{size:16});
const windowVals=[['음극 OCP','x₀ / x₁₀₀','y₀ / y₁₀₀','cₛ,max,n (mol m⁻³)','cₛ,max,p (mol m⁻³)']];
for(const r of data.window) windowVals.push([r.ocp,`${r.x0} / ${r.x100}`,`${r.y0} / ${r.y100}`,r.csn.toLocaleString(),r.csp.toLocaleString()]);
table(S(6),windowVals,52,265,930,130,13);
img(S(6),'../work/pptx_preview/ocp_fig_1.png',52,407,930,260,'Measured qOCV reconstruction and OCP errors');
text(S(6),'C/50 qOCV 재구성 및 오차 · 0.5C/1C/2C 동특성 그래프는 다음 장',80,669,880,20,{size:12,color:'#555'});
note(S(6),'OCP/window/csmax는 C/50 qOCV·양극 GITT 및 측정 음극 저율 OCP에서 계산. 상세 원자료·코드는 최종 간결본 notebook. 그림은 Stage 1 qOCV reconstruction 결과.');

// Slide 7: experimental Ds + k, six actual branches.
erase('sh/ydkbm5sv','sh/dcbud0ra');
text(S(7),'GITT Ds, EIS/DRT local-minimum Rct→k, 최신 입자 반경과 cₛ,max를 모두 반영',47,118,940,33,{size:16,bold:true,color:blue});
img(S(7),'Stage3_prefit_two_OCP_overlay_six_curves.png',38,155,960,455,'Experimental parameter voltage overlays');
const pre=[['OCP','0.5C 충/방 RMSE','1C 충/방 RMSE','2C 충/방 RMSE']];
for(const ocp of ['Ai2020','Measured']) pre.push([ocp,...[.5,1,2].map(rate=>`${metric(ocp,rate,'Charge').full_rmse} / ${metric(ocp,rate,'Discharge').full_rmse} mV`)]);
table(S(7),pre,80,620,880,68,12);
note(S(7),'동일 시간축 전체 구간 RMSE. 추가 MAE와 실제 full-cell SoC 10–70% 조건별 값은 부록 A2 및 최종 notebook. Dsn=4.1708e-14, Dsp=6.8882e-14 m2/s; Rct,n=27.0009 Ωcm2 (N3+N4 local-minimum), Rct,p=25.387 Ωcm2.');

// Slide 8: measured sensitivity, identifiability and restricted fit.
erase('sh/6h0fypgb');
img(S(8),'Stage2_seven_parameter_sensitivity_identifiability.png',35,155,645,500,'Seven-parameter sensitivity and Pearson identifiability');
text(S(8),'7개 후보: Dsₙ, Dsₚ, kₙ, kₚ, bₙ, bₛ, bₚ',695,164,300,34,{size:17,bold:true,color:blue});
text(S(8),'5% 교란 시 bₙ: 약 9 mV\nkₙ·kₚ: 각각 약 2 mV\nDsn·Dsp: 0.1 mV 이하',695,211,300,95,{size:17});
text(S(8),'상관: kₙ–kₚ ≈ 1.00\nbₛ–bₚ ≈ 0.99',695,320,300,72,{size:17});
text(S(8),'피팅: kₙ, bₙ만 변동\nDs·kₚ·bₛ·bₚ 고정\n충전 3율 학습 → 방전 3율 검증',695,412,300,100,{size:16,bold:true});
text(S(8),'bₙ은 porosity와 결합된 유효값. 경계 도달 시 고유 물성으로 해석하지 않음.',695,527,300,85,{size:14,color:orange});
note(S(8),'민감도/식별도 방법: COBRAPRO, J Electrochem Soc 2024 https://pangea.stanford.edu/ERE/pdf/OnoriPDF/Journals/96.pdf. Bruggeman 계수 세 개 포함. b_n은 porosity와 구조적 중복이 있으므로 유효 보정계수로만 해석.');

// Slide 9: preserve existing two tables; fill the blank post-fit area.
const rows9=[['OCP / 전해질','충전 RMSE','방전 RMSE','2C 방전 RMSE','1C 음극 최소']];
for(const ocp of ['Ai2020','Measured']) for(const el of ['Ai2020_corrected','ORegan2022']){
  const ch=scenario(ocp,el,fitted,'Charge'), dis=scenario(ocp,el,fitted,'Discharge');
  const branch2=data.branch.find(r=>r.ocp===ocp&&r.electrolyte===el&&r.phase===fitted&&r.rate===2&&r.direction==='Discharge');
  const p=pot(ocp,el,fitted,1);
  rows9.push([`${ocp} / ${el==='ORegan2022'?'ORegan':'Ai2020'}`,`${ch.full_rmse} mV`,`${dis.full_rmse} mV`,
    `${branch2.full_rmse} mV`,`${p.min_mv} mV`]);
}
table(S(9),rows9,44,572,952,109,10);
text(S(9),'상단 표는 기존 비교, 하단 표는 이번 OCP별 충전 피팅·독립 방전 검증(3율 RMS 집계)',380,543,605,25,{size:12,color:'#4b5563'});
note(S(9),'위 두 기존 표는 수정하지 않음. 아래 신규 표는 동일 실험값·두 OCP·각 전해질에 대해 충전만 피팅하고 방전은 독립 검증. 전체/SoC 10–70% RMSE·MAE 조건별 값은 notebook. ORegan2022는 전해질 bulk 수송함수 5개만 이전, 전극화학은 Enertech/Ai2020을 유지. https://www.sciencedirect.com/science/article/abs/pii/S0013468622008593');

// Slide 10: preserve existing statements and place evidence within the blank card.
erase('sh/ydcnatkn');
const oA=fit('Ai2020','ORegan2022'), oM=fit('Measured','ORegan2022');
text(S(10),'전해질 ORegan2022는 검토 후보. 충전 목적함수는 낮아졌지만 방전·음극 전위 독립 검증이 남음.',75,284,900,45,{size:18,bold:true,color:blue});
text(S(10),`ORegan 충전 피팅: Ai2020 OCP ${oA.train_rmse} mV · Measured OCP ${oM.train_rmse} mV`,75,337,900,34,{size:17});
text(S(10),'2C 방전 조기 cutoff와 1–2C 음극 전위 0 mV 하회가 확인됨. 이 조건으로 목표 달성을 선언하지 않음.',75,380,900,55,{size:17,color:orange,bold:true});
note(S(10),'현재 하단 기존 서술은 원본 유지. 이번 계산의 추가 결론은 바로 위 삽입문과 최종 notebook 기준. 충전 전압 적합만으로 도금 안전성을 검증할 수 없으며, 별도 음극 전위 측정이 필요.');

// Slide 12: OCP rationale and complete per-rate core RMSE/MAE table.
img(S(12),'../work/pptx_preview/ocp_fig_1.png',40,130,960,307,'Measured qOCV reconstruction and OCP errors');
text(S(12),'실험 파라미터 적용 후 · 전체 / 실제 풀셀 SoC 10–70% 전압오차 (RMSE/MAE, mV)',45,443,940,27,{size:15,bold:true,color:blue});
const vals12=[['조건','Ai2020 전체','Measured 전체','Ai2020 SoC 10–70','Measured SoC 10–70']];
for(const rate of [.5,1,2]) for(const dir of ['Charge','Discharge']){
 const a=metric('Ai2020',rate,dir), m=metric('Measured',rate,dir);
 vals12.push([`${rate}C ${dir==='Charge'?'충전':'방전'}`,`${a.full_rmse}/${a.full_mae}`,`${m.full_rmse}/${m.full_mae}`,
    `${a.soc_rmse}/${a.soc_mae}`,`${m.soc_rmse}/${m.soc_mae}`]);
}
table(S(12),vals12,46,477,946,167,10);
text(S(12),'* 2C 충전은 cutoff로 SoC 70% 미도달: Ai2020 ~59%, Measured ~63%까지만 비교',55,675,940,20,{size:11,color:orange});
note(S(12),'qOCV 그림은 원본 Stage 1 notebook 계산 결과. 표는 수정된 실제 full-cell SoC 10–70% 재계산(실험 capacity coulomb counting, branch 초기 rest 전압 SoC). 2C 충전 비교 끝 SoC는 시나리오별로 달라짐; 상세 범위 CSV 참고.');

// Slide 17: electrolyte provenance and the limits of transfer.
erase('sh/rm5czq50','sh/6lwvql4v');
text(S(17),'ORegan2022에서 이전한 항목 (전해질 bulk 함수만)',48,128,940,32,{size:19,bold:true,color:blue});
table(S(17),[['항목','Ai2020 기준','ORegan2022 반영'],
  ['Dₑ(cₑ,T)','Ai2020 함수','ORegan2022 함수'],
  ['κₑ(cₑ,T)','Ai2020 함수','ORegan2022 함수'],
  ['t₊, χ(cₑ,T)','Ai2020 함수','ORegan2022 함수'],
  ['초기 cₑ','1000 mol m⁻³','ORegan2022 값 (1000)']],70,188,900,205,14);
text(S(17),'전극 OCP, 기하, cₛ,max, GITT Ds, EIS/DRT k 초기값은 Enertech 셀 값 유지.',68,420,900,38,{size:17,bold:true});
text(S(17),'ORegan2022 원논문은 LG M50 계열 NMC/Si-graphite 셀이다. 따라서 전체 parameter set을 가져오면 셀 화학이 달라진다.',68,472,900,70,{size:16});
text(S(17),'전압 개선과 음극 전위의 0 mV 통과를 별도로 평가하며, 전해질 변경만으로 plating 위험 예측을 검증했다고 보지 않는다.',68,561,900,70,{size:16,color:orange});
text(S(17),'변경 이유: Ai2020 전해질에서 1C ≈ −22 mV, 2C ≈ −78 mV (두 OCP 공통). ORegan2022도 재검증 필요.',68,637,900,45,{size:15,bold:true,color:blue});
note(S(17),'ORegan et al. 2022 https://www.sciencedirect.com/science/article/abs/pii/S0013468622008593. PyBaMM ORegan2022 함수 중 Electrolyte diffusivity, conductivity, Cation transference number, Thermodynamic factor, Initial concentration만 전달. 기타 전극 parameter는 현재 셀 측정값 유지.');

// Slide 18: preserve the existing discussion, replace the blue prompt with evidence.
erase('sh/je5knut0');
img(S(18),'Stage3_ORegan_two_OCP_anode_potential_PPT.png',48,395,945,264,'Two-OCP anode potential at separator interface after charge fit');
text(S(18),'모델 예측: 1–2C에서 0 mV 하회. 실측 reference electrode 검증 필요.',260,671,630,22,{size:12,color:orange});
note(S(18),'그림은 DFN negative electrode surface potential difference at separator interface, Ai2020 anode OCP + ORegan electrolyte, charge-fitted parameters. 0 mV는 열역학적 참고선이며 도금 여부의 직접 실험 증거는 아님.');

// Render changed slides for manual QA, then run the required finalizer.
const preview=path.join(workspaceDir,'work','pptx_final_preview');
await fs.mkdir(preview,{recursive:true});
for(const n of [6,7,8,9,10,12,17,18]){
  const blob=await presentation.export({slide:S(n),format:'png',scale:1});
  await fs.writeFile(path.join(preview,`slide-${String(n).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL_DIR,'container_tools','artifact_tool_utils.mjs')).href);
const stagingDir=path.join(workspaceDir,'.codex-finalizer');
await fs.mkdir(stagingDir,{recursive:true});
const candidatePath=path.join(stagingDir,'261021_meeting_candidate.pptx');
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);
const sha=crypto.createHash('sha256').update(await fs.readFile(source)).digest('hex');
const result=await finalizePresentation({
  explicitTotalSlideCount:21,
  workspaceDir,candidatePath,finalPath:FINAL_PPTX,pythonExecutable:RUNTIME_PYTHON,
  integrityValidatorPath:path.join(SKILL_DIR,'container_tools','inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(SKILL_DIR,'container_tools','inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','9906000,6858000','--validate-heading-fit'],
  verifyArtifactToolImport:true,
  receiptPath:path.join(stagingDir,'261021_meeting_validation.json'),
});
console.log(JSON.stringify({result,finalPath:FINAL_PPTX,preview},null,2));
