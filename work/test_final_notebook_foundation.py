import json, os
from pathlib import Path
root=Path(__file__).resolve().parents[1]
os.chdir(root/'outputs')
nb=json.loads((root/'outputs'/'261021_Ai2020_DFN_최종_간결본.ipynb').read_text(encoding='utf-8'))
scope={'__name__':'__main__'}
for i in (3,4,5,6,7,8,9,10,11,13):
    source=''.join(nb['cells'][i]['source'])
    if i==3: source='\n'.join(x for x in source.splitlines() if not x.startswith('%pip'))
    print('EXEC',i,flush=True)
    exec(compile(source,f'notebook_cell_{i}','exec'),scope)
    print('DONE',i,flush=True)
case=scope['cases']()['Ai2020_anode__measured_cathode']
print('READY',scope['CELL_CAPACITY_AH'],scope['kinetic_constants'](case),flush=True)
