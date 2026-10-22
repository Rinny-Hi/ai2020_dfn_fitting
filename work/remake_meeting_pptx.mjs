import fs from 'node:fs/promises';
import fsSync from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';

const root=path.resolve('.');
const outDir=path.join(root,'outputs');
const files=await fs.readdir(outDir);
const src=path.join(outDir,files.find(x=>x.startsWith('261021_')&&x.endsWith('완성본.pptx')));
const final=path.join(outDir,'261021_미팅자료_Ai2020_ORegan_최신검증본.pptx');
const skill='C:/Users/user/.codex/plugins/cache/openai-primary-runtime/presentations/26.1015.12040/skills/presentations';
const py='C:/Users/user/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe';
process.env.RUNTIME_NODE_MODULES='C:/Users/user/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const p=await PresentationFile.importPptx(await FileBlob.load(src));
const S=n=>p.slides.items[n-1];
const blue='#156082', medium='#4E85BF', dark='#262626', pale='#E8EEF6', light='#F6F8FB', orange='#B84426', green='#117A64';
const font='Pretendard';

function tx(s,v,x,y,w,h,z={}){
 const q=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{style:'solid',fill:'none',width:0}});
 q.text=v; q.text.insets=0;
 q.text.style={fontFamily:font,fontSize:z.size??16,color:z.color??dark,bold:z.bold??false,alignment:z.align??'left'};
 return q;
}
function line(s,x,y,w,c=blue,t=2){s.shapes.add({geometry:'line',position:{left:x,top:y,width:w,height:0},fill:'none',line:{style:'solid',fill:c,width:t}})}
function hd(s,sub,app=false){
 s.background.fill='#FFFFFF';
 tx(s,app?'Appendix. 파라미터 선정 근거':'Ai2020 모델 Fitting',44,22,945,48,{size:30,bold:true});
 line(s,38,77,966,'#263A85',2.2);
 tx(s,sub,47,86,946,43,{size:22,bold:true,color:blue});
}
function table(s,rows,x,y,w,h,z={}){
 const t=s.tables.add({rows:rows.length,columns:rows[0].length,left:x,top:y,width:w,height:h,values:rows});
 t.styleOptions={headerRow:true,bandedRows:false};
 for(let i=0;i<rows.length;i++)for(let j=0;j<rows[0].length;j++){
  const c=t.getCell(i,j); c.fill=i===0?blue:(i%2?pale:light);
  c.text.style={fontFamily:font,fontSize:z.size??13,color:i===0?'#FFFFFF':dark,bold:i===0,alignment:j===0?'left':'center'};
 }
 t.borders.assign({style:'solid',fill:'#FFFFFF',width:1});
 t.cells.block({row:0,column:0,rowCount:rows.length,columnCount:rows[0].length}).assign({margins:{top:1,bottom:1,left:4,right:4}});
 return t;
}
function img(s,name,x,y,w,h,alt){
 const file=path.join(outDir,name); s.images.add({blob:fsSync.readFileSync(file),contentType:'image/png',alt,fit:'contain',position:{left:x,top:y,width:w,height:h}});
}
function notes(s,v){s.speakerNotes.text=v;}
const map=await p.inspect({select:'$.slides[*]',include:['index','elements.id']});
for(const n of [6,7,8,9,10,12,13,15,16,17,18]){
 const ids=map.records[n-1].elements.map(e=>e.id).filter(Boolean);
 if(ids.length)p.delete(ids);
}

// 6 — OCP alternatives and stoichiometry-specific concentrations.
{
 const s=S(6);hd(s,'[1단계] OCP 2가지 버전과 stoichiometry window');
 tx(s,'같은 셀의 C/50 qOCV를 각 OCP로 재구성하고, window에 맞춰 cₛ,max를 각각 계산',50,135,940,34,{size:17,bold:true,color:blue});
 table(s,[['음극 OCP','x₀ / x₁₀₀','y₀ / y₁₀₀','cₛ,max,n (mol m⁻³)','cₛ,max,p (mol m⁻³)'],
 ['Ai2020','0.003803 / 0.954167','0.970716 / 0.431646','25,184.17','49,143.14'],
 ['Measured','0.002150 / 0.739784','0.976291 / 0.448293','32,447.1','50,174']],49,190,944,152,{size:13});
 const pic=path.join(root,'work','pptx_preview','ocp_fig_1.png');
 if(fsSync.existsSync(pic))s.images.add({blob:fsSync.readFileSync(pic),contentType:'image/png',alt:'C/50 qOCV reconstruction and OCP comparison',fit:'contain',position:{left:70,top:364,width:900,height:270}});
 tx(s,'현재 선정 분기: Ai2020 음극 + 측정 양극 OCP. 이후 kₙ·kₚ 계산에는 이 분기의 cₛ,max 사용.',52,659,930,28,{size:14,color:orange,bold:true});
 notes(s,'C/50 qOCV와 GITT rest OCP 기반. 두 음극 OCP는 비교 후보. 현 선택은 Ai2020 graphite OCP. c_s,max는 각 window와 Q_cell=2.356528 Ah로 계산. 출처: 최종 간결본 및 Stage 1 노트북.');
}
// 7 — stage A/B/C comparison with all six branch curves.
{
 const s=S(7);hd(s,'[2–3단계] 실험 초기값 → Bruggeman 고정 → 충전 kₙ 피팅');
 tx(s,'A: EIS kₙ, bₙ=2.914   →   B: bₙ=2.5 고정   →   C: 충전 3율만 kₙ 피팅·방전 3율 독립 검증',49,128,960,34,{size:15,bold:true,color:blue});
 img(s,'ORegan2022_Ai2020_피팅전후_충방전6곡선.png',33,162,973,468,'Six voltage curves comparing experiment and stages A, B, C');
 table(s,[['6곡선 전체 RMSE','A · 실험 초기','B · bₙ=2.5','C · 최종'],['mV','66.69','39.70','37.14']],155,632,730,59,{size:14});
 notes(s,'Cell 6-6 초기 RPT, 0.5/1/2C. A k_n=6.307602e-7 b_n=2.914; B same k_n b_n=2.5; C k_n=9.702793e-7 b_n=2.5. All stages use latest k_p=4.768216e-7. D_s from GITT, electrolyte bulk from ORegan2022. Exact data: ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv.');
}
// 8 — sensitivity, identifiability and method.
{
 const s=S(8);hd(s,'민감도·식별도 분석과 피팅 변수 선정');
 img(s,'Stage2_seven_parameter_sensitivity_identifiability.png',35,148,620,485,'Sensitivity and identifiability for seven DFN parameters');
 tx(s,'후보 7개',675,160,310,30,{size:20,bold:true,color:blue});
 tx(s,'Dₛ,n · Dₛ,p · kₙ · kₚ\nbₙ · bₛ · bₚ',675,196,300,68,{size:17});
 tx(s,'5% 교란: bₙ ≈ 9 mV; kₙ·kₚ ≈ 2 mV\nDₛ,n·Dₛ,p < 0.1 mV',675,282,300,72,{size:16});
 tx(s,'kₙ–kₚ 상관 ≈ 1.00\nbₛ–bₚ 상관 ≈ 0.99',675,373,300,70,{size:16,bold:true});
 tx(s,'자유 kₙ+bₙ 피팅은 kₙ 상한에 도달하고 2C 방전 RMSE 122.06 mV. 따라서 bₙ=2.5 고정 후 kₙ만 피팅.',675,467,305,125,{size:16,color:orange,bold:true});
 tx(s,'충전 3곡선 학습 / 방전 3곡선 검증. PSO 대신 1변수 bounded scalar 탐색.',54,655,925,30,{size:15,bold:true,color:blue});
 notes(s,'COBRAPRO 자료와 Stage2_seven_parameter_sensitivity_identifiability.png 기반. 7개 후보와 3개 Bruggeman 계수 포함. 현재 1변수 objective: 140 fixed time points/rate to 98% experimental end; separate 20 early cutoff penalty terms/rate of 300 mV*shortfall fraction. This penalty is heuristic, not uncertainty-calibrated. SciPy bounded scalar: https://docs.scipy.org/doc/scipy/reference/optimize.minimize_scalar-bounded.html ; PyBOP: https://github.com/pybop-team/PyBOP.');
}
// 9 — final curves.
{
 const s=S(9);hd(s,'최종 선정 모델: 충전 피팅 3곡선 · 방전 독립 검증 3곡선');
 tx(s,'Ai2020 흑연 OCP + 측정 양극 OCP + ORegan2022 전해질 벌크 물성',51,129,940,31,{size:17,bold:true,color:blue});
 img(s,'ORegan2022_Ai2020_최종_충방전6곡선.png',36,160,969,476,'Final DFN voltage compared with six experiment curves');
 tx(s,'전체 RMSE: 충전 41.85 mV · 방전 31.74 mV     |     kₙ=9.7028×10⁻⁷, bₙ=2.5, kₚ=4.7682×10⁻⁷',55,651,930,39,{size:16,bold:true,color:blue});
 notes(s,'Black=experiment; orange=selected DFN. Charge 0.5/1/2C fitted, discharge 0.5/1/2C held out. Rate aggregation is RMS across three branch RMSE values, not pooled point RMSE. Final per-rate metrics on slide 12.');
}
// 10 — selected model and potential criterion.
{
 const s=S(10);hd(s,'최종 모델 선정: 1C까지의 SoC별 음극 전위 확인');
 img(s,'ORegan2022_Ai2020_최종_SoC별음극전위.png',30,150,690,460,'Selected DFN charging anode potential versus full-cell SOC');
 tx(s,'전위 최저값',735,162,270,35,{size:21,bold:true,color:blue});
 table(s,[['충전율','최저 (mV)'],['0.5C','+36.72'],['1C','−3.33'],['2C','−49.55']],728,207,260,160,{size:16});
 tx(s,'1C: SoC 70% +19.54 mV; 75% +8.94 mV. 0 mV 미만은 약 79–80% SoC의 끝부분에 한정.',735,385,266,125,{size:16,bold:true,color:dark});
 tx(s,'2C는 약 34.5–61.1% SoC에서 음수. 현재 모델 선정의 허용 기준은 1C까지.',735,520,265,95,{size:16,color:orange,bold:true});
 tx(s,'음극 전위는 DFN 분리막 계면 예측값이며 Li plating 실측 판정은 아님.',48,650,950,38,{size:15,color:orange,bold:true});
 notes(s,'Negative electrode surface potential difference at separator interface. 1C min -3.3335 mV at ~80.1% SoC; below zero ~79.0-80.1%. Mesh 15 to 25 points: -3.334 to -3.332 mV. 2C reaches only ~61.1% SoC. No reference electrode or plating measurement.');
}
// 12 — exact per-rate metrics and tradeoffs.
{
 const s=S(12);hd(s,'A2. 최종 충·방전 전압오차 및 중앙 SoC 검증',true);
 tx(s,'전압오차 mV; 중앙 = 실제 풀셀 SoC 10–70% 중 실험·모델 공통 관측 구간',49,135,945,34,{size:16,bold:true,color:blue});
 table(s,[['조건','역할','전체 RMSE / MAE','중앙 RMSE / MAE','종료용량 오차'],
 ['0.5C 충전','피팅','35.94 / 28.97','37.60 / 35.16','−1.84%'],
 ['0.5C 방전','검증','35.60 / 27.26','30.01 / 26.28','+0.46%'],
 ['1C 충전','피팅','42.18 / 28.57','26.21 / 24.00','−1.78%'],
 ['1C 방전','검증','36.01 / 27.68','29.52 / 25.07','+0.91%'],
 ['2C 충전','피팅','46.74 / 23.13','10.69 / 9.67','−3.40%'],
 ['2C 방전','검증','21.41 / 11.30','8.31 / 4.95','+0.03%']],48,191,945,292,{size:12.5});
 tx(s,'B→C 전체 RMSE 변화 (mV)',52,511,930,32,{size:18,bold:true,color:blue});
 table(s,[['0.5C 충 / 방','1C 충 / 방','2C 충 / 방'],['32.53→35.94 / 29.72→35.60','44.09→42.18 / 32.73→36.01','57.10→46.74 / 35.21→21.41']],52,555,936,76,{size:12.5});
 tx(s,'C는 2C와 종료용량을 개선하지만 0.5C 두 방향은 B보다 악화. 2C 충전 중앙값은 SoC 10.8–59.2%만 평가.',51,652,940,39,{size:14,color:orange,bold:true});
 notes(s,'Numbers from ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv. Metrics on common elapsed-time segment. Central actual full-cell SoC by experimental capacity. 2C charge stops before 70%; central 10.84–59.22% only. Early cutoff separately reflected in capacity error.');
}
// 13 — current OCP and capacity basis.
{
 const s=S(13);hd(s,'A2. OCP window와 cₛ,max 계산 근거',true);
 tx(s,'실험 평균용량 Qcell = (2.355758 + 2.357298) / 2 = 2.356528 Ah',50,139,940,31,{size:18,bold:true,color:blue});
 tx(s,'Qₙ=Qcell/(x₁₀₀−x₀),  Qₚ=Qcell/(y₀−y₁₀₀),  cₛ,max,k=3600 Qₖ/(F A Lₖ εₛ,k)',52,185,938,41,{size:17});
 table(s,[['항목','음극 n','양극 p'],['공통 overlap 면적','787.236 cm²','787.236 cm²'],['전극 두께','76.5 µm','68.0 µm'],['활물질 부피분율','0.610','0.620'],['Ai2020 OCP window','0.003803–0.954167','0.970716–0.431646'],['cₛ,max','25,184.17 mol m⁻³','49,143.14 mol m⁻³']],92,262,855,292,{size:15});
 tx(s,'측정 음극 OCP 분기는 x₀/x₁₀₀=0.002150/0.739784, cₛ,max,n=32,447.1 mol m⁻³로 별도 계산.',61,591,920,46,{size:16,bold:true,color:blue});
 tx(s,'주의: 과거 슬라이드의 Qcell 2.35628 Ah는 자릿수 누락; kₙ·kₚ 환산에 다른 OCP 분기의 cₛ,max를 혼용하지 않음.',61,648,920,44,{size:14,color:orange,bold:true});
 notes(s,'Original A2 slide had typographic omission 2.35628 for 2.356528 Ah. Current window obtained from latest C/50 qOCV fit. F=96485.33212 C/mol. Source: final compact notebook, Anode_DRT_Rct_local_minimum_계산요약.md.');
}
// 15 — latest cathode kinetics.
{
 const s=S(15);hd(s,'A4. 양극 EIS–DRT 기반 kₚ 최신값',true);
 tx(s,'양극 charge-transfer 피크 P4의 반복 재현성을 근거로 Rct=25.387 Ω cm² 적용',50,139,935,39,{size:18,bold:true,color:blue});
 table(s,[['계산 입력','적용값','근거'],['Rct,p','25.387 Ω cm²','양극 EIS–DRT P4'],['Lₚ · εₛ,p · rₚ','68.0 µm · 0.620 · 4.239 µm','최신 전극 기하'],['cₛ,max,p','49,143.14 mol m⁻³','현재 OCP window'],['cₑ','1000 mol m⁻³','초기 전해질 농도'],['j₀,p','0.339 A m⁻²','Rct 환산'],['kₚ','4.768216×10⁻⁷','현재 모든 A/B/C 계산에 고정']],55,212,930,305,{size:14});
 tx(s,'과거 PPT의 4.37×10⁻⁷ 및 10/20 공동 피팅의 5.89952×10⁻⁷은 이번 모델의 입력값이 아님.',59,550,920,56,{size:17,color:orange,bold:true});
 tx(s,'k 단위: A m⁻² (m³ mol⁻¹)¹·⁵   |   비교 시 OCP에 맞는 cₛ,max와 전극 기하를 함께 적용',59,631,920,49,{size:15,color:blue});
 notes(s,'Latest k_p=4.768215585284518e-7 in all A/B/C calculations. Derived from cathode Rct 25.387 Ω cm2 and current stoichiometry / c_s,max. Old 10/20 joint-fit k_p=5.899520e-7 only for historical restoration. See ORegan2022_Ai2020_최종결과_방법검토.md.');
}
// 16 — anode DRT local-minimum computation and fit value.
{
 const s=S(16);hd(s,'A4. 음극 N3+N4 local-minimum 적분과 kₙ',true);
 img(s,'Anode_DRT_local_minimum_50SOC.png',42,154,568,358,'Anode DRT 50% SOC local-minimum bounds');
 table(s,[['50% SoC 항목','값'],['N4 기여','9.9608 Ω cm²'],['N3 기여','17.0401 Ω cm²'],['N3+N4 Rct','27.0009 Ω cm²'],['EIS 초기 kₙ','6.307602×10⁻⁷'],['충전 피팅 kₙ','9.702793×10⁻⁷']],620,165,372,298,{size:14});
 tx(s,'옛 N2/N3 ↔ 현재 N3/N4: 고주파 접촉 피크를 현재 N1로 따로 표기한 결과',54,533,932,45,{size:17,bold:true,color:blue});
 tx(s,'골 위치는 smoothing으로 탐색하고, 원본 γ를 ln f에 대해 적분. 음극 half-cell Rct는 초기값이며 확정 물성으로 해석하지 않음.',54,595,930,80,{size:16,color:orange,bold:true});
 notes(s,'At 50% SOC N4/LF=0.2861 Hz; N3/N4=5.0447 Hz; N2/N3=86.4315 Hz. Rct=27.0009 Ωcm2, j0=0.2509 A/m2. N3+N4 from original DRT gamma over d ln f; smoothing only for local minima. See outputs/Anode_DRT_Rct_local_minimum_계산요약.md.');
}
// 17 — transfer and objective.
{
 const s=S(17);hd(s,'A5. ORegan2022 전해질 이식과 최적화 방식',true);
 table(s,[['ORegan2022에서 이식','대상 셀 값으로 유지'],['Dₑ(cₑ,T), κₑ(cₑ,T)','전극 OCP · stoichiometry window'],['t₊(cₑ,T), χ(cₑ,T), 초기 cₑ','기하 · cₛ,max · GITT Dₛ · EIS k 초기값']],55,160,930,170,{size:14});
 tx(s,'ORegan2022 전체 파라미터셋이 아니라 전해질 벌크 함수 5개만 이식',60,350,925,42,{size:18,bold:true,color:blue});
 table(s,[['피팅 구분','10/20 과거','현재 선정'],['학습 데이터','충전+방전 6곡선','충전 3곡선'],['검증 데이터','독립 방전 없음','방전 3곡선 보류'],['변수 / 탐색','kₙ+bₙ / PSO 7×7','bₙ=2.5, kₙ / bounded scalar'],['kₚ','5.89952×10⁻⁷','4.768216×10⁻⁷']],58,416,924,207,{size:13});
 tx(s,'현재 조기 종료 벌점(300 mV×시간 부족률)은 경험적 가중치. 중앙 SoC 악화 원인 검토 및 반복성 기반 비용함수 개선이 후속 과제.',58,640,926,52,{size:14,color:orange,bold:true});
 notes(s,'ORegan et al. 2022 original parameter set targets LG M50 NMC/Si-graphite; only electrolyte bulk transport functions were transferred. Past-vs-current methodology: ORegan2022_과거현재_방식_파라미터_차이.md. Objective: fixed time voltage samples (140 per rate, to 98% experimental end) and 20 separate early-cutoff terms/rate. SciPy bounded scalar official docs https://docs.scipy.org/doc/scipy/reference/optimize.minimize_scalar-bounded.html. PyBOP architecture https://github.com/pybop-team/PyBOP.');
}
// 18 — past and current potential.
{
 const s=S(18);hd(s,'A6. 과거 공동 피팅 대비 음극 전위 변화',true);
 img(s,'ORegan2022_과거공동피팅_vs_현재선정_음극전위.png',42,150,955,380,'Historical and current model charging anode potential comparison');
 table(s,[['충전율','과거 최저','현재 최저','변화'],['0.5C','+31.43 mV','+36.72 mV','+5.29 mV'],['1C','−11.20 mV','−3.33 mV','+7.87 mV'],['2C','−59.74 mV','−49.55 mV','+10.19 mV']],109,537,820,131,{size:13.5});
 tx(s,'같은 SoC 전위 상승은 주로 kₙ 증가·bₙ 감소에 기인; kₚ 변화는 주로 컷오프 시점에 작용.',51,678,941,31,{size:13.5,color:orange,bold:true});
 notes(s,'Old 10/20 model: joint six-curve fit with k_n=8.726541e-7, b_n=2.533583, k_p=5.899520e-7, r_n=3.7542 µm. Current: charge-only fit k_n=9.702793e-7, b_n=2.5, k_p=4.768216e-7, r_n=3.692 µm. Same OCP family and electrolyte. Changes cannot be attributed to one parameter alone. See ORegan2022_과거현재_방식_파라미터_차이.md and contribution CSV.');
}

const preview=path.join(root,'work','pptx_revised_preview');await fs.mkdir(preview,{recursive:true});
for(const n of [6,7,8,9,10,12,13,15,16,17,18]){
 const png=await p.export({slide:S(n),format:'png',scale:1});
 await fs.writeFile(path.join(preview,`slide-${String(n).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
}
const staging=path.join(root,'work','pptx_revised_build');await fs.mkdir(staging,{recursive:true});
const draft=path.join(staging,'candidate.pptx');await(await PresentationFile.exportPptx(p)).save(draft);
const {finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools','artifact_tool_utils.mjs')).href);
const result=await finalizePresentation({explicitTotalSlideCount:21,workspaceDir:root,candidatePath:draft,finalPath:final,pythonExecutable:py,
 integrityValidatorPath:path.join(skill,'container_tools','inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools','inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','9906000,6858000','--validate-heading-fit',...[6,7,10,12,13,15,16,17,18].flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:[6,7,10,12,13,15,16,17,18],verifyArtifactToolImport:true,
 receiptPath:path.join(staging,'validation.json')});
console.log(JSON.stringify({final,preview,result},null,2));
