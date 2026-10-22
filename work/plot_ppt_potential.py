from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

OUT=Path(__file__).resolve().parents[1]/'outputs'
curves=pd.read_csv(OUT/'Stage3_charge_fit_curves.csv')
fig,axes=plt.subplots(1,2,figsize=(15,4.3),sharey=True,layout='constrained')
for ax,(ocp,label) in zip(axes,[('Ai2020_anode__measured_cathode','Ai2020 anode OCP'),
                                ('measured_anode__measured_cathode','Measured anode OCP')]):
    for rate,color in [(.5,'#156082'),(1.,'#c44e32'),(2.,'#6a4c93')]:
        d=curves.loc[curves.OCP_case.eq(ocp)&curves.Electrolyte.eq('ORegan2022')
            &curves.Phase.eq('charge_fit_discharge_validation')&curves.Direction.eq('Charge')
            &curves.C_rate.eq(rate)&curves.Series.eq('DFN')]
        ax.plot(d.SOC_pct,d.Anode_potential_mV,color=color,lw=2.4,label=f'{rate:g}C')
    ax.axhline(0,color='black',lw=1.2)
    ax.axhline(20,color='#777',ls='--',lw=.9)
    ax.set(xlim=(10,80),ylim=(-125,210),xlabel='Full-cell SoC (%)',title=label)
    ax.grid(alpha=.2)
    ax.legend(frameon=False,ncol=3,loc='upper right')
axes[0].set_ylabel('Anode potential vs electrolyte (mV)')
fig.suptitle('ORegan2022 electrolyte | charge-fitted DFN | separator interface',fontsize=14)
fig.savefig(OUT/'Stage3_ORegan_two_OCP_anode_potential_PPT.png',dpi=170)
print(OUT/'Stage3_ORegan_two_OCP_anode_potential_PPT.png')
