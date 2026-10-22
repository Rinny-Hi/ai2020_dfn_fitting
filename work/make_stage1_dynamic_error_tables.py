"""Slide-ready error tables for the Stage-1 three-way dynamic comparison."""
from pathlib import Path
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'
data=pd.read_csv(glob.glob(str(OUT/'Stage1*3단계*동적전압_RMSE_MAE.csv'))[0])
font_manager.fontManager.addfont('C:/Windows/Fonts/malgun.ttf')
font_manager.fontManager.addfont('C:/Windows/Fonts/malgunbd.ttf')
plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False})
stages=list(data['stage'].drop_duplicates())
labels=['A  Ai2020 공개 window','B  Ai2020 OCP\nwindow 재최적화','C  측정 양극 OCP\nwindow 재최적화']
conditions=[(.5,'Charge','0.5C 충전'),(.5,'Discharge','0.5C 방전'),(1.0,'Charge','1C 충전'),(1.0,'Discharge','1C 방전'),(2.0,'Charge','2C 충전'),(2.0,'Discharge','2C 방전')]

def draw(kind):
    if kind=='full':
        title='1단계 동적 전압 오차 · 전체 공통 시간 구간'
        rmse='full_RMSE_mV';mae='full_MAE_mV'
        foot='각 조건 400개 공통 경과시간점. 마지막 행은 6조건 RMSE의 제곱평균(RMS).'
        out=OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_전체오차표.png'
    else:
        title='1단계 동적 전압 오차 · 중앙 전달용량 10–70%'
        rmse='center_10_70_RMSE_mV';mae='center_10_70_MAE_mV'
        foot='중앙 구간은 간결본의 전달용량 기준 10–70% 시간 구간. 마지막 행은 6조건 RMSE의 제곱평균(RMS).'
        out=OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_중앙오차표.png'
    rows=[]
    for rate,direction,label in conditions:
        vals=[]
        for stage in stages:
            r=data[(data.stage==stage)&(data.C_rate==rate)&(data.direction==direction)].iloc[0]
            vals.append(f'{r[rmse]:.2f} / {r[mae]:.2f}')
        rows.append([label,*vals])
    overall=[float(np.sqrt(np.mean(data.loc[data.stage==stage,rmse]**2))) for stage in stages]
    rows.append(['6곡선 종합 RMSE',*[f'{v:.2f}' for v in overall]])

    fig,ax=plt.subplots(figsize=(14.8,6.5),dpi=180)
    fig.patch.set_facecolor('white');ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    ax.text(.035,.925,title,fontsize=22,fontweight='bold',color='#156082',va='center')
    ax.text(.035,.872,'단위: mV  |  셀 안 숫자: RMSE / MAE',fontsize=15,color='#344251',va='center')
    x=[.035,.235,.49,.745, .965]
    top=.805;rowh=.091
    header=['조건',*labels]
    for j in range(4):
        ax.add_patch(Rectangle((x[j],top-rowh),x[j+1]-x[j],rowh,facecolor='#156082',edgecolor='white',linewidth=1))
        ax.text((x[j]+x[j+1])/2,top-rowh/2,header[j],ha='center',va='center',color='white',fontsize=11.5,fontweight='bold',linespacing=1.15)
    for i,row in enumerate(rows):
        y=top-(i+2)*rowh
        fill='#DCE7F2' if i%2==0 else '#F3F6FA'
        if i==len(rows)-1:fill='#C4DAEB'
        for j,val in enumerate(row):
            ax.add_patch(Rectangle((x[j],y),x[j+1]-x[j],rowh,facecolor=fill,edgecolor='white',linewidth=1))
            ax.text((x[j]+x[j+1])/2,y+rowh/2,val,ha='center',va='center',fontsize=15.5 if j else 14.5,
                    color='#202A35',fontweight='bold' if i==len(rows)-1 else 'normal')
    ax.text(.035,.047,foot,fontsize=12.5,color='#4A5560',va='center')
    fig.savefig(out,dpi=200,facecolor='white',bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    print(out)

draw('full')
draw('center')
