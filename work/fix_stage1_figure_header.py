"""Finalize scientific plot header with a font supporting Korean labels."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import font_manager

root=Path(__file__).resolve().parents[1]
out=root/'outputs'
png=out/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_qOCV_비교.png'
svg=out/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_qOCV_비교.svg'
im=plt.imread(png)
h,w=im.shape[:2]
fp=font_manager.FontProperties(fname='C:/Windows/Fonts/malgunbd.ttf',size=15)
fig=plt.figure(figsize=(w/220,h/220),dpi=220)
ax=fig.add_axes([0,0,1,1])
ax.imshow(im)
ax.add_patch(patches.Rectangle((0,0),w,112,facecolor='white',edgecolor='none',zorder=3))
ax.text(145,64,'공통 입력: Ai2020 De 함수 × 1e-4   |   rn = 3.692 µm   |   rp = 4.239 µm',
        fontproperties=fp,color='#156082',va='center',zorder=4)
ax.set_xlim(0,w);ax.set_ylim(h,0);ax.axis('off')
fig.savefig(png,dpi=220,facecolor='white',pad_inches=0)
fig.savefig(svg,facecolor='white',pad_inches=0)
plt.close(fig)
print(png)
