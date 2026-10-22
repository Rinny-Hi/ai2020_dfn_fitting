from pathlib import Path
import numpy as np
import openpyxl
from scipy.signal import find_peaks, savgol_filter

p=Path(r"C:\Users\user\Downloads\Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx")
w=openpyxl.load_workbook(p,read_only=True,data_only=True)
rows=list(w["Origin_Paste"].iter_rows(values_only=True))
headers=rows[0]
for col in range(0,18,2):
    soc=int(str(headers[col]).split('%')[0])
    arr=np.array([(r[col],r[col+1]) for r in rows[3:] if len(r)>col+1 and isinstance(r[col],(int,float)) and isinstance(r[col+1],(int,float))],float)
    arr=arr[np.argsort(arr[:,0])]
    f,y=arr[:,0],arr[:,1]
    sm=savgol_filter(y,21,3)
    peaks,_=find_peaks(sm,prominence=max(sm)*0.04,distance=10)
    mins,_=find_peaks(-sm,prominence=max(sm)*0.02,distance=10)
    print(f'SOC {soc}: n={len(f)}, f={f.min():.3g}-{f.max():.3g}, gamma={y.min():.4g}-{y.max():.4g}')
    print(' peaks',[(round(f[i],4),round(sm[i],3)) for i in peaks])
    print(' minima',[(round(f[i],4),round(sm[i],3)) for i in mins])
